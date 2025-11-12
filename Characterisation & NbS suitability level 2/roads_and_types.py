import fiona
import geopandas as gpd
from shapely.geometry import shape
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
import os

# -------------------------------
# Road type and availability labels
# -------------------------------
road_type_labels = {
    1: "Highways", 2: "Primary", 3: "Secondary", 4: "Tertiary",
    5: "Local", 0: "Unspecified"
}
road_availability_labels = {
    1: "Seasonal", 2: "All year", 0: "Unspecified"
}
road_colors = {
    "Highways": "red", "Primary": "orange", "Secondary": "yellow",
    "Tertiary": "green", "Local": "blue", "Unspecified": "gray"
}

# -------------------------------
# Region info + appropriate projected CRS
# -------------------------------
# Get the directory of the current script


current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
out_dir = os.path.join(parent_dir, "data", "output")
boundaries_dir = os.path.join(parent_dir, "Demo Boundaries")
roads_path = os.path.join(parent_dir, "data", "input", "GRIP4_Region4_vector_shp", "GRIP4_region4.shp")
# Create output folders
output_dir = os.path.join(out_dir, "road_outputs")
os.makedirs(output_dir, exist_ok=True) #make sure output directory exists


# Load shapefiles
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

# -------------------------------
# Try encoding helper
# -------------------------------
def safe_read_shapefile(filepath):
    encodings_to_try = ['utf-8', 'ISO-8859-1', 'latin1']
    for enc in encodings_to_try:
        try:
            return gpd.read_file(filepath, encoding=enc)
        except UnicodeDecodeError:
            print(f" Encoding failed: {enc}")
    raise RuntimeError(" Failed to decode shapefile with known encodings.")

# -------------------------------
# Process each region
# -------------------------------
for info in regions_info:
    region_name = info['region']
    region_shapefile = info['shapefile']
    print(f"\n Processing region: {region_name}")

    try:
        # Load region shapefile
        region_gdf = safe_read_shapefile(region_shapefile)

        # Read road CRS from Fiona
        with fiona.open(roads_path) as road_src:
            roads_crs = road_src.crs_wkt  # safer full definition
            src_crs = road_src.crs  # for fallback

        # Reproject region to match roads CRS for intersection
        region_gdf_match = region_gdf.to_crs(roads_crs)
        region_geom_match = region_gdf_match.geometry.union_all()

    except Exception as e:
        print(f" Failed to load region or reproject: {e}")
        continue

    # -------------------------------
    # Select intersecting roads
    # -------------------------------
    selected_features = []
    with fiona.open(roads_path) as src:
        for feat in src:
            if feat['geometry'] is None:
                continue
            try:
                geom = shape(feat['geometry'])
                if geom.is_empty:
                    continue
                if geom.intersects(region_geom_match):
                    selected_features.append({
                        'geometry': geom,
                        'properties': feat['properties']
                    })
            except Exception:
                continue

    if not selected_features:
        print(f" No valid roads found in {region_name}")
        continue

    # -------------------------------
    # Convert to GeoDataFrame
    # -------------------------------
    roads_gdf = gpd.GeoDataFrame.from_features(selected_features, crs=src_crs)

    # Reproject both region and roads to EPSG:3035 for accurate measurement
    target_crs = "EPSG:3035"
    region_gdf = region_gdf.to_crs(target_crs)
    roads_gdf = roads_gdf.to_crs(target_crs)

    # -------------------------------
    # Map attributes and compute length
    # -------------------------------
    roads_gdf["road_type"] = roads_gdf["GP_RTP"].map(road_type_labels)
    roads_gdf["availability"] = roads_gdf["GP_RAV"].map(road_availability_labels)
    roads_gdf["length_km"] = roads_gdf.geometry.length / 1000  # In km

    # -------------------------------
    # Summary statistics
    # -------------------------------
    type_summary = roads_gdf.groupby("road_type")["length_km"].sum().sort_values(ascending=False)
    print(f" Road Type Summary (km):\n{type_summary.round(2)}")

    # -------------------------------
    # Plot Map
    # -------------------------------
    fig, ax = plt.subplots(figsize=(10, 10))
    region_gdf.boundary.plot(ax=ax, color='black', linewidth=1)

    for rtype, subset in roads_gdf.groupby("road_type"):
        color = road_colors.get(rtype, "gray")
        subset.plot(ax=ax, color=color, linewidth=0.6, label=rtype)

    ax.set_title(f"{region_name} - GRIP4 Roads by Type", fontsize=14)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.legend(title="Road Type", fontsize=9)
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{region_name}_roads_map.png"), dpi=300)
    plt.show()

    # -------------------------------
    # Bar Chart
    # -------------------------------
    type_summary.plot(kind='bar', color=[road_colors.get(rt, "gray") for rt in type_summary.index])
    plt.title(f"{region_name} - Total Road Length by Type")
    plt.ylabel("Length (km)")
    plt.xlabel("Road Type")
    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{region_name}_roads_chart.png"), dpi=300)
    plt.show()