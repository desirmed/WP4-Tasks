import os
import geopandas as gpd
import ee
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import rioxarray
import numpy as np
from urllib.request import urlretrieve
import pandas as pd

# === Initialize Earth Engine ===
try:
    ee.Initialize(project='xxxxx')  # Replace with your GEE project ID
except Exception:
    ee.Authenticate()
    ee.Initialize(project='xxxxx') # Replace with your GEE project ID

# === Define SARL Asset ===
sarl = ee.Image("projects/sat-io/open-datasets/SARL")

# === Class Info (excluding 0, 5, 6) ===
sarl_classes = {
    1: 'Permanent River',
    2: 'Permanent Lake',
    3: 'Seasonal River',
    4: 'Seasonal Lake'
}
sarl_colors = {
    1: '#FFD700',  # Yellow - Permanent River
    2: '#00FFFF',  # Cyan   - Permanent Lake
    3: '#0000FF',  # Blue   - Seasonal River
    4: '#6A0DAD'   # Purple - Seasonal Lake
}


current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
out_dir = os.path.join(parent_dir, "data", "output")
boundaries_dir = os.path.join(parent_dir, "Demo Boundaries")
input_dir= os.path.join(parent_dir, "data", "input", "HydroRIVERS_v10_eu_shp", "HydroRIVERS_v10_eu_shp")
# Create output folders
output_dir = os.path.join(out_dir, "SARL")
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



# === Function to Convert GeoDataFrame to EE Geometry ===
def gdf_to_ee(gdf):
    return ee.Geometry(gdf.geometry.union_all().__geo_interface__)

# === Visualization Setup ===
palette = [sarl_colors[k] for k in sorted(sarl_colors)]
cmap = plt.matplotlib.colors.ListedColormap(palette)
bounds = [0.5, 1.5, 2.5, 3.5, 4.5]
norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap.N)

# #%% HydroRIVERS dataset are overlayed on top of SARL data for visualization purposes
# # === Load HydroRIVERS dataset ===
# rivers_path = os.path.join(input_dir, "HydroRIVERS_v10_eu.shp")
# rivers_gdf = gpd.read_file(rivers_path)
# rivers_gdf = rivers_gdf.to_crs("EPSG:4326")

# # === Process Each Region ===
for item in regions_info:
    region = item['region']
    shapefile_path = item['shapefile']
    print(f" Processing {region}...")

    try:
        gdf = gpd.read_file(shapefile_path)
    except UnicodeDecodeError:
        print(f" Encoding issue detected in {region}, retrying with ISO-8859-1...")
        gdf = gpd.read_file(shapefile_path, encoding='ISO-8859-1')

    gdf = gdf.to_crs("EPSG:4326")
    ee_geom = gdf_to_ee(gdf)
    bbox = ee_geom.bounds()

    # Select SARL year band
    year_band = 'Y2021'
    band_img = sarl.select(year_band).clip(bbox)

    # Mask out background (0), no-data (5, 6)
    mask = band_img.gt(0).And(band_img.lt(5))
    clean_img = band_img.updateMask(mask)

    # Download as GeoTIFF
    filename = f"{region}_SARL_{year_band}.tif".replace(" ", "_")
    tif_path = os.path.join(output_dir, filename)

    if not os.path.exists(tif_path):
        url = clean_img.getDownloadURL({
            'scale': 100,
            'region': bbox,
            'format': 'GEO_TIFF'
        })
        urlretrieve(url, tif_path)

    # Clip river layer to current region
    rivers_clipped = gpd.clip(rivers_gdf, gdf)

    # Read and Plot
    da = rioxarray.open_rasterio(tif_path, masked=True).squeeze()
    fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')
    da.plot.imshow(ax=ax, cmap=cmap, norm=norm, add_colorbar=False)
    gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    # Plot rivers
    if not rivers_clipped.empty:
        rivers_clipped.plot(ax=ax, color='steelblue', linewidth=0.8, alpha=0.8, label='Rivers')

    # Title & Labels
    ax.set_title(f"{region} - SARL {year_band}", fontsize=14, fontweight='bold')
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # Legend
    labels = [sarl_classes[i] for i in sorted(sarl_classes)]
    patches = [mpatches.Patch(color=sarl_colors[i], label=labels[i - 1]) for i in sorted(sarl_classes)]
    if not rivers_clipped.empty:
        patches.append(mpatches.Patch(color='steelblue', label='Rivers'))
    ax.legend(handles=patches, title='Water Class', loc='lower right', frameon=True)

    # Save
    out_png = os.path.join(output_dir, f"{region}_SARL_{year_band}.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.show()

#%% Save SARL data to .tif and .png without overlaying rivers data on top of it

# Create colormap & normalization
cmap = mcolors.ListedColormap([sarl_colors[i] for i in sarl_colors])
norm = mcolors.BoundaryNorm(list(sarl_classes.keys()) + [5], cmap.N)

# =========================
#   MAIN LOOP PER REGION
# =========================
for item in regions_info:
    region = item['region']
    shapefile_path = item['shapefile']
    print(f"Processing {region}...")

    # Load region shapefile
    #gdf = gpd.read_file(shapefile_path).to_crs("EPSG:4326")
    # Load region shapefile safely (handles encoding issues)
    try:
        gdf = gpd.read_file(shapefile_path)
    except UnicodeDecodeError:
        print(f"Encoding issue in {region} — retrying with ISO-8859-1...")
        gdf = gpd.read_file(shapefile_path, encoding='ISO-8859-1')

    gdf = gdf.to_crs("EPSG:4326")

    # Convert to Earth Engine geometry
    ee_geom = gdf_to_ee(gdf)
    bbox = ee_geom.bounds()

    # Select SARL year
    year_band = "Y2021"
    band_img = sarl.select(year_band).clip(bbox)

    # Mask out classes 0 (background), 5 (No Data Lake), 6 (No Data River)
    mask = band_img.gt(0).And(band_img.lt(5))  # Keep only 1–4
    clean_img = band_img.updateMask(mask)

    # Download GeoTIFF if not exists
    #tif_name = f"{region}_SARL_{year_band}.tif".replace(" ", "_")
    #tif_path = os.path.join(output_dir, tif_name)
    tif_path = os.path.join(output_dir, f"{region}_SARL_{year_band}.tif")

    if not os.path.exists(tif_path):
        url = clean_img.getDownloadURL({
            'scale': 100,
            'region': bbox,
            'format': 'GEO_TIFF'
        })
        urlretrieve(url, tif_path)

    # Read raster with rioxarray
    da = rioxarray.open_rasterio(tif_path, masked=True).squeeze()

    # === PLOT SARL WATER ONLY ===
    fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')
    da.plot.imshow(ax=ax, cmap=cmap, norm=norm, add_colorbar=False)

    # Plot region boundary
    gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    # Title and labels
    ax.set_title(f"{region} - SARL {year_band}", fontsize=14, fontweight='bold')
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # === Legend Outside ===
    patches = [mpatches.Patch(color=sarl_colors[i], label=sarl_classes[i]) for i in sarl_classes]

    ax.legend(
        handles=patches,
        title='SARL Water Class',
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        frameon=True
    )

    out_png = os.path.join(output_dir, f"{region}_SARL_{year_band}.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")  # ensures legend is included
    plt.show()
#%%

# Output CSV
csv_path = os.path.join(output_dir, "SARL_summary_table.csv")

# Collect results
results = []

for item in regions_info:
    region = item['region']
    shapefile_path = item['shapefile']
    print(f"Processing {region}...")

    # Load region shapefile safely
    try:
        gdf = gpd.read_file(shapefile_path)
    except UnicodeDecodeError:
        gdf = gpd.read_file(shapefile_path, encoding='ISO-8859-1')
    gdf = gdf.to_crs("EPSG:4326")

    # Convert to EE geometry
    ee_geom = gdf_to_ee(gdf)

    # Select SARL year
    year_band = "Y2021"
    band_img = sarl.select(year_band)

    # Mask out 0, 5, 6
    mask = band_img.gt(0).And(band_img.lt(5))
    clean_img = band_img.updateMask(mask)

    # Count pixels per class
    hist = clean_img.reduceRegion(
        reducer=ee.Reducer.frequencyHistogram(),
        geometry=ee_geom,
        scale=30,
        maxPixels=1e13
    ).get(year_band)

    hist_dict = ee.Dictionary(hist).getInfo()  # Convert to Python dict

    # Total masked pixels
    total_pixels = sum([hist_dict.get(str(cls), 0) for cls in water_classes])

    # Store results
    row = {"Region": region}
    for cls in water_classes:
        count = hist_dict.get(str(cls), 0)
        percent = (count / total_pixels * 100) if total_pixels > 0 else 0
        row[f"{water_classes[cls]}_Pixels"] = count
        row[f"{water_classes[cls]}_Percent"] = percent
    results.append(row)

# Convert to DataFrame
df_summary = pd.DataFrame(results)

# Optional: reorder columns nicely
cols_order = ["Region"]
for cls in water_classes.values():
    cols_order += [f"{cls}_Pixels", f"{cls}_Percent"]
df_summary = df_summary[cols_order]

# Save CSV
df_summary.to_csv(csv_path, index=False)
print(f"Summary table saved to {csv_path}")
# %%
df = pd.read_csv(csv_path)
print(df.columns)

# Get all classes from column names automatically
percent_cols = [col for col in df.columns if col.endswith("_Percent")]
classes = [col.replace("_Percent", "") for col in percent_cols]

# Now use 'classes' and 'percent_cols' in plotting
for region in df["Region"].unique():
    site_df = df[df["Region"] == region]
    percentages = [site_df[col].values[0] for col in percent_cols]

# Directory to save plots
plot_dir = os.path.join(output_dir, "SARL_percentage_plots")
os.makedirs(plot_dir, exist_ok=True)

# Loop through each region/site
for region in df["Region"].unique():
    site_df = df[df["Region"] == region]

    # Extract percentages
    percentages = [site_df[f"{cls}_Percent"].values[0] for cls in classes]

    # Create bar plot
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(classes, percentages, color=["#08306b", "#4292c6", "#a6bddb", "#c6dbef"])
    
    # Add percentage labels above bars
    for bar, pct in zip(bars, percentages):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, f"{pct:.1f}%", 
                ha='center', va='bottom', fontsize=10)

    ax.set_title(f"{region} - SARL Water Classes Percentage", fontsize=14, fontweight='bold')
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 100)
    ax.set_xlabel("Water Class")
    
    plt.tight_layout()
    
    # Save figure
    plot_path = os.path.join(output_dir, f"{region}_SARL_percentage.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    print(f"Plot saved for {region}: {plot_path}")
