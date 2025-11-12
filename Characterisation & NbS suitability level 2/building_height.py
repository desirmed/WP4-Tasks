import os
import geopandas as gpd
import ee
import matplotlib.pyplot as plt
import rioxarray
from urllib.request import urlretrieve
import numpy as np
import pandas as pd

# Initialize Earth Engine
try:
    ee.Initialize(project='xxxxx')  # Replace with your GEE project ID
except Exception:
    ee.Authenticate()
    ee.Initialize(project='xxxxx')  # Replace with your GEE project ID

current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
out_dir = os.path.join(parent_dir, "data", "output")
boundaries_dir = os.path.join(parent_dir, "Demo Boundaries")
# Create output folders
output_dir = os.path.join(out_dir, "ghs_height_rasters")
os.makedirs(output_dir, exist_ok=True) #make sure output directory exists


ghs_obat_iso_assets = {
    'Italy':     'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_ITA_E2020_R2024A_V1_0',
    'France':    'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_FRA_E2020_R2024A_V1_0',
    'Spain':     'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_ESP_E2020_R2024A_V1_0',
    'Portugal':  'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_PRT_E2020_R2024A_V1_0',
    'Croatia':   'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_HRV_E2020_R2024A_V1_0',
    'Cyprus':    'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_CYP_E2020_R2024A_V1_0',
    'Greece':    'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_GRC_E2020_R2024A_V1_0'
}

# Load shapefiles
shapefiles = []
for root, dirs, files in os.walk(boundaries_dir):
    for file in files:
        if file.endswith('.shp'):
            full_path = os.path.join(root, file)
            region = os.path.splitext(file)[0]
            country = os.path.basename(os.path.dirname(full_path))
            shapefiles.append((country, region, full_path))

def gdf_to_ee(gdf):
    return ee.Geometry(gdf.geometry.union_all().__geo_interface__)

results = []
# Loop through all regions
for country, region, path in shapefiles:
    asset_path = ghs_obat_iso_assets.get(country)

    if not asset_path:
        print(f" No GHS-OBAT asset for {country}, skipping {region}")
        continue

    try:
        gdf = gpd.read_file(path, encoding='ISO-8859-1', engine='fiona')
        if gdf.empty:
            raise ValueError("Empty shapefile")

        ee_geom = gdf_to_ee(gdf)
        ee_geom_bbox = ee_geom.bounds()

        # Filter buildings with height data -> extract all buildings within the region that have a positive height attribute
        fc = ee.FeatureCollection(asset_path)\
               .filterBounds(ee_geom)\
               .filter(ee.Filter.neq('height', None))\
               .filter(ee.Filter.gt('height', 0))

        count = fc.size().getInfo()
        if count == 0:
            print(f" No height data for {region}")
            continue

        # Create and download building height raster for the region
        image = fc.reduceToImage(['height'], ee.Reducer.first()).clip(ee_geom_bbox)
        filename = f"{region}_height.tif".replace(" ", "_")
        tif_path = os.path.join(output_dir, filename)

        if os.path.exists(tif_path):
            os.remove(tif_path)  # delete old 500m file

        url = image.getDownloadURL({
            'scale': 150,
            'region': ee_geom_bbox,
            'format': 'GEO_TIFF'
        })
        urlretrieve(url, tif_path)

        # if not os.path.exists(tif_path):
        #     url = image.getDownloadURL({
        #         'scale': 150,
        #         'region': ee_geom_bbox,
        #         'format': 'GEO_TIFF'
        #     })
        #     urlretrieve(url, tif_path)

        # Open and mask raster
        da = rioxarray.open_rasterio(tif_path).squeeze()
        da.rio.write_crs("EPSG:4326", inplace=True)
        masked = da.rio.clip(gdf.geometry.values, gdf.crs, drop=False)
        
        # Compute building height statistics (ignore NaNs)
        data = masked.values.flatten()
        data = data[~np.isnan(data)]

        if len(data) == 0:
            print(f" No valid height pixels for {region}")
            continue

        mean_height = float(np.mean(data))
        median_height = float(np.median(data))
        max_height = float(np.max(data))
        min_height = float(np.min(data))
        perc95_height = float(np.percentile(data, 95))
        pixel_count = int(len(data))

        # Store results
        results.append({
            'Country': country,
            'Region': region,
            'Pixel_Count': pixel_count,
            'Mean_Height': mean_height,
            'Median_Height': median_height,
            'Min_Height': min_height,
            'Max_Height': max_height,
            'P95_Height': perc95_height
        })

        # Plot
        # Fix visualization range
        # vmax = float(np.nanpercentile(masked.values, 95))
        # vmin = 1
        vmin, vmax = 0, 10

        fig, ax = plt.subplots(figsize=(7, 6))
        im = masked.plot(
            ax=ax,
            cmap='RdYlGn_r',
            vmin=vmin,
            vmax=vmax,
            add_colorbar=False,
            alpha=1.0
        )

        gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.7)
        ax.set_title(f"{region} Building Heights", fontsize=13, fontweight='bold')
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_facecolor('white')

        cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.03)
        cbar.set_label("Building Height (m)", fontsize=10)

        plt.axis('equal')
        plt.tight_layout()
        out_path = f"{region}_building_heights.png".replace(" ", "_") 
        plt.savefig(os.path.join(output_dir, out_path), dpi=300)
       # plt.show()

        print(f" Saved: {out_path}")

    except Exception as e:
        print(f" Error for {region}: {e}")

# Save statistics to CSV
if results:
    df = pd.DataFrame(results)
    csv_path = os.path.join(output_dir, "building_height_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n📄 Saved summary CSV: {csv_path}")