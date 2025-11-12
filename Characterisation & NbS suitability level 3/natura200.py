import os
import geopandas as gpd
import matplotlib.pyplot as plt
import fiona
from shapely.geometry import shape, mapping
import rasterio
from rasterio.features import rasterize
import pandas as pd

# === File Paths ===

current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
input_dir = os.path.join(parent_dir, "data", "input", "natura2000", "natura2000", "Natura2000_end2023_epsg4326.shp")
print(input_dir)
out_dir = os.path.join(parent_dir, "data", "output")
boundaries_dir = os.path.join(parent_dir, "Demo Boundaries")
# Create output folders
output_dir = os.path.join(out_dir, "natura2000")
os.makedirs(output_dir, exist_ok=True) #make sure output directory exists

# === Define Regions ===
regions_info = [
    {'region': 'Split_Dalmatia', 'country': 'Croatia', 'shapefile': os.path.join(boundaries_dir, "Croatia", "Split_Dalmatia.shp")},
    {'region': 'Macedonia_Thrace', 'country': 'Greece', 'shapefile': os.path.join(boundaries_dir, "Greece", "Macedonia_Thrace.shp")},
    {'region': 'Potenza', 'country': 'Italy', 'shapefile': os.path.join(boundaries_dir, "Italy", "Potenza.shp")},
    {'region': 'Corse_du_Sud', 'country': 'France', 'shapefile': os.path.join(boundaries_dir, "France", "Corse_du_Sud.shp")},
    {'region': 'Sardegna', 'country': 'Italy', 'shapefile': os.path.join(boundaries_dir, "Italy", "Sardegna.shp")},
    {'region': 'Beiras_Centro', 'country': 'Portugal', 'shapefile': os.path.join(boundaries_dir, "Portugal", "Beiras_Centro.shp")},
    {'region': 'Nicosia', 'country': 'Cyprus', 'shapefile': os.path.join(boundaries_dir, "Cyprus", "Nicosia.shp")},
    {'region': 'Valencia', 'country': 'Spain', 'shapefile': os.path.join(boundaries_dir, "Spain", "Valencia.shp")}
]

# === Process Each Region ===
for item in regions_info:
    region = item['region']
    shapefile_path = item['shapefile']
    print(f" Displaying Natura 2000 map for {region}...")

    try:
        region_gdf = gpd.read_file(shapefile_path).to_crs("EPSG:4326")
    except UnicodeDecodeError:
        region_gdf = gpd.read_file(shapefile_path, encoding="ISO-8859-1").to_crs("EPSG:4326")

    # Get bounding box for region
    bbox = tuple(region_gdf.total_bounds)  # (minx, miny, maxx, maxy)

    # Filter features inside bbox using fiona
    with fiona.open(input_dir, 'r') as src:
        filtered_features = [
            {**feature, 'geometry': feature['geometry']}
            for feature in src.filter(bbox=bbox)
        ]


    # Convert to GeoDataFrame and clip to actual region geometry
    natura_gdf = gpd.GeoDataFrame.from_features(filtered_features, crs="EPSG:4326")
    natura_clipped = gpd.clip(natura_gdf, region_gdf)

    # # Save to Shapefile
    shp_path = os.path.join(output_dir, f"{region}_natura2000.shp")
    natura_clipped.to_file(shp_path)

    # Save to GeoTIFF (rasterize)
    try:
        bounds = natura_clipped.total_bounds  # minx, miny, maxx, maxy
        resolution = 0.01  # degrees per pixel (adjust as needed)
        width = int((bounds[2] - bounds[0]) / resolution)
        height = int((bounds[3] - bounds[1]) / resolution)

        transform = rasterio.transform.from_origin(bounds[0], bounds[3], resolution, resolution)

        shapes = [(mapping(geom), 1) for geom in natura_clipped.geometry if geom is not None]

        raster = rasterize(
            shapes,
            out_shape=(height, width),
            transform=transform,
            fill=0,
            dtype='uint8'
        )

        tiff_path = os.path.join(output_dir, f"{region}_natura.tiff")
        with rasterio.open(
            tiff_path,
            'w',
            driver='GTiff',
            height=height,
            width=width,
            count=1,
            dtype='uint8',
            crs='EPSG:4326',
            transform=transform
        ) as dst:
            dst.write(raster, 1)

    except Exception as e:
        print(f"Could not save TIFF for {region}: {e}")


    # Plot
    fig, ax = plt.subplots(figsize=(10, 10))
    region_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=1)
    natura_clipped.plot(ax=ax, color='green', alpha=0.5, label='Natura 2000 Areas')

    ax.set_title(f"Natura 2000 Areas in {region.replace('_', ' ')}", fontsize=14)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.legend()
    ax.tick_params(labelsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()



#%%
# Create list to store area stats
area_stats = []

# Replace your plotting loop with this (summary only shown)
for item in regions_info:
    region = item['region']
    shapefile_path = item['shapefile']

    print(f" Processing area for {region}...")

    try:
        region_gdf = gpd.read_file(shapefile_path).to_crs("EPSG:4326")
    except UnicodeDecodeError:
        region_gdf = gpd.read_file(shapefile_path, encoding="ISO-8859-1").to_crs("EPSG:4326")

    bbox = tuple(region_gdf.total_bounds)

    with fiona.open(natura_path, 'r') as src:
        filtered_features = [
            {**feature, 'geometry': feature['geometry']}
            for feature in src.filter(bbox=bbox)
        ]

    if not filtered_features:
        print(f" No Natura 2000 areas found in {region}")
        area_stats.append({'Region': region, 'Area_km2': 0})
        continue

    natura_gdf = gpd.GeoDataFrame.from_features(filtered_features, crs="EPSG:4326")
    natura_clipped = gpd.clip(natura_gdf, region_gdf)

    # Reproject to a metric CRS for accurate area (e.g., EPSG:3035 – Europe LAEA)
    natura_clipped = natura_clipped.to_crs("EPSG:3035")
    natura_clipped['area_km2'] = natura_clipped.geometry.area / 10**6

    total_area = natura_clipped['area_km2'].sum()
    area_stats.append({'Region': region, 'Area_km2': total_area})



# Create DataFrame
area_df = pd.DataFrame(area_stats).sort_values(by='Area_km2', ascending=False)
print(area_df)

# Save to CSV
csv_path = os.path.join(output_dir, f"natura2000_km2.csv")
area_df.to_csv(csv_path, index=False)


# Plot
plt.figure(figsize=(10, 6))
plt.barh(area_df['Region'], area_df['Area_km2'], color='green')
plt.xlabel("Natura 2000 Area (km²)")
plt.title("Natura 2000 Coverage by Region")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
