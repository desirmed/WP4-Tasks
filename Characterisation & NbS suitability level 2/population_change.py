import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import rioxarray as rxr
import xrspatial
import zipfile

data_folder = 'data'
output_folder = 'output'

if not os.path.exists(data_folder):
    os.mkdir(data_folder)
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

def download(url):
    filename = os.path.join(data_folder, os.path.basename(url))
    if not os.path.exists(filename):
        from urllib.request import urlretrieve
        local, _ = urlretrieve(url, filename)
        print('Downloaded ' + local)

data_url = 'https://github.com/spatialthoughts/geopython-tutorials/releases/download/data/'

gpw_pop_2010 = 'gpw-v4-population-count-rev11_2010_2pt5_min_tif.zip'
gpw_pop_2020 = 'gpw-v4-population-count-rev11_2020_2pt5_min_tif.zip'

download(data_url + gpw_pop_2010)
download(data_url + gpw_pop_2020)






#Zoom to meditaranean
#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import rioxarray as rxr
import os
import xrspatial

# File paths
current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
inp_dir = os.path.join(parent_dir, "data", "input")
pop2010_path = os.path.join(inp_dir, 'gpw-v4-population-count-rev11_2010_2pt5_min_tif', 'gpw_v4_population_count_rev11_2010_2pt5_min.tif')
pop2020_path = os.path.join(inp_dir, 'gpw-v4-population-count-rev11_2020_2pt5_min_tif', 'gpw_v4_population_count_rev11_2020_2pt5_min.tif')
out_dir = os.path.join(parent_dir, "data", "output")
output_dir = os.path.join(out_dir, "population_change")
os.makedirs(output_dir, exist_ok=True) #make sure output directory exists


# Load rasters
pop2010 = rxr.open_rasterio(pop2010_path, mask_and_scale=True).squeeze()
pop2020 = rxr.open_rasterio(pop2020_path, mask_and_scale=True).squeeze()
change = pop2020 - pop2010

# Clip to Mediterranean bounding box using rio.clip_box
change_med = change.rio.clip_box(minx=-10, maxx=35, miny=30, maxy=50)

# Reclassify the clipped data
class_bins = [-100, 100, 1000, np.inf]
class_values = [1, 2, 3, 4]
change_class = xrspatial.classify.reclassify(change_med, bins=class_bins, new_values=class_values)

# Define colors and labels
colors = ['#3288bd', '#e0e0e0', '#fdae61', '#d7191c']
labels = ['Decline', 'Neutral', 'Growth', 'High Growth']

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

change_class.plot.imshow(
    ax=ax,
    add_colorbar=False,
    levels=[1, 2, 3, 4, 5],
    colors=colors
)

patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(4)]
ax.legend(handles=patches, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.1))

ax.set_title(" Mediterranean Population Change Classes (2010–2020)", y=1.02)
ax.axis('off')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "population_change_mediterranean.png"), dpi=300)
plt.show()