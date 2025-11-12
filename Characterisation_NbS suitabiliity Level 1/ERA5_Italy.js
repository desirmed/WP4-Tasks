// This script allows to retrieve the ERA5_Land monthly precipitation sum, used as input to compute the SPI
// ----------------------------------
// 1. Load Italy boundaries
// ----------------------------------
var adm0 = ee.FeatureCollection("projects/sat-io/open-datasets/geoboundaries/CGAZ_ADM0");
var italy  = ee.FeatureCollection("projects/sat-io/open-datasets/geoboundaries/CGAZ_ADM0").filter('shapeISO == "ITA"')

// ----------------------------------
// 2. Load ERA5 monthly precipitation
// ----------------------------------
var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
            .filterDate('2000-01-01', '2024-12-31') // Calibration period
            .select('total_precipitation_sum')
            .map(function(img) {
              return img.clip(italy);
            });

// ----------------------------------
// 3. Stack all months into one image
// ----------------------------------
var era5Stacked = era5.toBands();

// Rename bands to include date for clarity
var dateStrings = era5.aggregate_array('system:index');
era5Stacked = era5Stacked.rename(dateStrings);

// ----------------------------------
// 4. Export to Google Drive
// ----------------------------------
Export.image.toDrive({
  image: era5Stacked,
  description: 'ERA5_Italy_2000_2024',
  folder: 'ERA5_Italy',  // change to your desired folder in Drive
  fileNamePrefix: 'ERA5_Italy_2000_2024',
  scale: 10000,          // 10 km resolution
  region: italy.geometry(),
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

print('Stacked ERA5 image ready for export');
