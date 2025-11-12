import os
import geopandas as gpd
import ee
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import rioxarray
import numpy as np
from urllib.request import urlretrieve

# Initialize Earth Engine
try:
    ee.Initialize(project='xxxxx')  # Replace with your GEE project ID
except Exception:
    ee.Authenticate()
    ee.Initialize(project='xxxxx')  # Replace with your GEE project ID

# GHS-OBAT assets per country
ghs_obat_iso_assets = {
    'Italy':     'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_ITA_E2020_R2024A_V1_0',
    'France':    'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_FRA_E2020_R2024A_V1_0',
    'Spain':     'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_ESP_E2020_R2024A_V1_0',
    'Portugal':  'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_PRT_E2020_R2024A_V1_0',
    'Croatia':   'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_HRV_E2020_R2024A_V1_0',
    'Cyprus':    'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_CYP_E2020_R2024A_V1_0',
    'Greece':    'projects/sat-io/open-datasets/JRC/GHS-OBAT/GHS_OBAT_GPKG_GRC_E2020_R2024A_V1_0'
}


current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
out_dir = os.path.join(parent_dir, "data", "output")
boundaries_dir = os.path.join(parent_dir, "Demo Boundaries")
# Create output folders
output_dir = os.path.join(out_dir, "ghs_epoch_rasters")
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


def gdf_to_ee(gdf):
    return ee.Geometry(gdf.geometry.union_all().__geo_interface__)

# Epoch color map
epoch_colors = {
    1: '#40a6ff',  # <1980
    2: '#2aff80',  # 1980–1990
    3: '#ffee40',  # 1990–2000
    4: '#ff7f0e',  # 2000–2010
    5: '#b2182b'   # 2010–2020
}
cmap = plt.matplotlib.colors.ListedColormap([epoch_colors[i] for i in epoch_colors])
bounds = list(epoch_colors.keys()) + [6]
norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap.N)

#  Loop through all regions
for item in regions_info:
    region = item['region']
    country = item['country']
    shapefile_path = item['shapefile']

    print(f" Processing {region} in {country}...")

    try:
        try:
            gdf = gpd.read_file(shapefile_path)
        except UnicodeDecodeError:
            gdf = gpd.read_file(shapefile_path, encoding='ISO-8859-1')  # Fallback if UTF-8 fails

        ee_geom = gdf_to_ee(gdf)
        bbox = ee_geom.bounds()

        # Get FeatureCollection and reduce to image
        asset_path = ghs_obat_iso_assets[country]
        fc = ee.FeatureCollection(asset_path).filterBounds(ee_geom).filter(ee.Filter.neq('epoch', None))
        image = fc.reduceToImage(['epoch'], ee.Reducer.first()).clip(bbox)

        # Download GeoTIFF
        filename = f"{region}_epoch.tif".replace(" ", "_")
        tif_path = os.path.join(output_dir, filename)

        if not os.path.exists(tif_path):
            url = image.getDownloadURL({
                'scale': 150,
                'region': bbox,
                'format': 'GEO_TIFF',
                'maxPixels': 1e9
            })
            urlretrieve(url, tif_path)

        # Load raster and reproject
        da = rioxarray.open_rasterio(tif_path, masked=True).squeeze().astype("float32")
        da = da.rio.reproject("EPSG:4326")
        masked_da = da.where(da > 0)

        # Reproject shapefile
        gdf = gdf.to_crs("EPSG:4326")

        # Plot
        fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')
        masked_da.plot.imshow(ax=ax, cmap=cmap, norm=norm, add_colorbar=False)
        gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.8)

        ax.set_title(f"{region} - Building Epochs", fontsize=14, fontweight='bold')
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

        # Legend
        labels = ['1 (<1980)', '2 (1980–1990)', '3 (1990–2000)', '4 (2000–2010)', '5 (2010–2020)']
        patches = [mpatches.Patch(color=epoch_colors[i], label=labels[i - 1]) for i in epoch_colors]
        ax.legend(handles=patches, title='GHS-OBAT Epoch', loc='lower right', frameon=True)

        # Save + Show
        output_png = os.path.join(output_dir, f"{region}_Epoch_LonLat.png")
        plt.tight_layout()
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        plt.show()

    except Exception as e:
        print(f" Error processing {region}: {e}")