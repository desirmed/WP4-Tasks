#%%


# # Define Earth Engine ERA5 dataset and filter by date
# era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') \
#     .filter(ee.Filter.date('2000-01-01', '2025-01-01')) \
#     .select('total_precipitation_sum')

# # Set spatial extent to cover all of Croatia with some buffer
# # Automatically zoom to each country — uncomment as needed

# # Portugal
# # plt.xlim(-10.0, -6.0)   # Longitude range
# # plt.ylim(36.5, 42.5)    # Latitude range

# # Italy
# plt.xlim(6.0, 19.0)
# plt.ylim(36.0, 47.5)

# # France
# # plt.xlim(-5.0, 9.5)
# # plt.ylim(41.0, 51.5)

# # Croatia
# # plt.xlim(13.0, 20.5)
# # plt.ylim(42.0, 47.0)

# # Cyprus
# # plt.xlim(32.0, 34.0)
# # plt.ylim(34.5, 35.7)

# # Spain
# # plt.xlim(-9.5, 4.5)
# # plt.ylim(36.0, 44.5)

# # plt.suptitle("Standardized Precipitation Index", fontsize=16)
# plt.show()

#%%
#======  New script cell  ======    
# The previous script would not work since python relies on xarray and would require a custom ee-xarray engine.
# A workaround is to export the data from GEE as GeoTIFFs and then process them locally with xarray.
# To do that, you can use the following code ERA5_Italy.js to export the Geotiffs, then you need to process them using those lines of code
# This is not the most elegant way, consider to downlooad ERA5 netcdf directly

import xarray as xr
import os
import glob
import rioxarray
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
from climate_indices import indices, compute

# ---------------------- 1. Load GeoTIFF ----------------------
# load .tif from your repository - to get the data use the ERA5_Italy.js code

ds = rioxarray.open_rasterio(r"C:\Users\ERA5_Italy_2000_2024.tif")

# Remove extra dimension if exists
da = ds.squeeze()

# Rename spatial dims from x,y → lon,lat (important!)
if 'x' in da.dims and 'y' in da.dims:
    da = da.rename({'x': 'lon', 'y': 'lat'})

# Rename band → time
da = da.rename({'band': 'time'})

# Create proper monthly timestamps
da = da.assign_coords(time=pd.date_range('2000-01-01', periods=da.sizes['time'], freq='M'))

# ---------------------- 2. SPI Function ----------------------
scale = 3
distribution = indices.Distribution.gamma
data_start_year = 2000
calibration_year_initial = 2000
calibration_year_final = 2024
periodicity = compute.Periodicity.monthly

def spi_func(precip_series):
    precip_series = np.array(precip_series, dtype=float)
    return indices.spi(
        precip_series,
        scale,
        distribution,
        data_start_year,
        calibration_year_initial,
        calibration_year_final,
        periodicity
    )

# ---------------------- 3. Apply to every pixel ----------------------
start = time.time()

spi_da = xr.apply_ufunc(
    spi_func,
    da,
    input_core_dims=[["time"]],
    output_core_dims=[["time"]],
    vectorize=True,               # Apply to each pixel independently
    dask="parallelized",          # Optional — speeds up
    output_dtypes=[float],
)

spi_da = spi_da.assign_coords(time=da.time)
spi_da.name = "SPI"

end = time.time()
print(f"SPI computation took {end - start:.2f} seconds")

# ---------------------- 4. Plot 2024 maps ----------------------
spi_2024 = spi_da.sel(time='2024')

spi_2024.plot(
    cmap='RdBu',
    col='time',
    col_wrap=4,
    levels=[-3, -2, -1, 0, 1, 2, 3],
    aspect=1.2
)

plt.xlim(6.0, 19.0)
plt.ylim(36.0, 47.5)

# Now you can save it as .csv, .tiff, .nc (preferably)
current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
out_dir = os.path.join(parent_dir, "data", "output")
spi_dir = os.path.join(out_dir, "spi") #create output directory
os.makedirs(spi_dir, exist_ok=True) #make sure output directory exists

df = spi_da.to_dataframe(name='SPI').reset_index()
csv_path = os.path.join(spi_dir, "SPI_Italy_2000_2024.csv")
df.to_csv(csv_path, index=False)

output_nc = "SPI_Italy_2000_2024.nc"
# Add metadata
spi_da.attrs["title"] = "Standardized Precipitation Index (SPI)"
spi_da.attrs["summary"] = "SPI calculated from ERA5 monthly precipitation for Italy (2000–2024)"
spi_da.attrs["source"] = "ERA5-Land Monthly Aggregated Precipitation"
spi_da.attrs["methodology"] = "SPI computed using gamma distribution and 3-month scale via climate_indices"
spi_da.attrs["projection"] = "EPSG:4326 (WGS84)"
spi_da.attrs["history"] = "Created with Python and xarray"
nc_path = os.path.join(spi_dir, output_nc)
spi_da.to_netcdf(nc_path)


# save all years as single GeoTIFF
#output_tif = "SPI_Italy_2000_2024.tif"
# output_path = "SPI_Italy_2000_2024.tif"
# spi_da.rio.to_raster(output_path)
# save only one year as GeoTIFF


# Select SPI for 2024
spi_2024 = spi_da.sel(time='2024')

# Rename spatial dimensions and assign CRS
spi_2024 = spi_2024.rename({'lat': 'y', 'lon': 'x'}).rio.write_crs("EPSG:4326")

# Reorder dimensions: ('time', 'y', 'x')
spi_2024 = spi_2024.transpose('time', 'y', 'x')

# Loop over each month
for i, t in enumerate(spi_2024.time.values):
    month_da = spi_2024.isel(time=i)
    month_str = pd.to_datetime(str(t)).strftime("%Y_%m")
    output_tif = f"SPI_Italy_{month_str}.tif"
    tif_path = os.path.join(spi_dir, output_tif)
    
    # Remove existing file if present
    if os.path.exists(tif_path):
        try:
            os.remove(tif_path)
        except PermissionError:
            raise PermissionError(f"Close the file {tif_path} if it's open.")
    
    # Save GeoTIFF
    month_da.rio.to_raster(tif_path)
    print(f"✅ Saved: {tif_path}")

# %%
