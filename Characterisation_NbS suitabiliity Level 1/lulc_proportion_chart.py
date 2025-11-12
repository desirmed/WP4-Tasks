import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# To run this script, you first need to download the CORINE LULC zonal statistics CSV using the lulc_proportion.js script
# https://github.com/desirmed/WP4-Tasks/blob/main/Characterisation_NbS%20suitabiliity%20Level%201/lulc_proportions.js

# Load your data
file_path = r"C:\Users\Gebruiker\OneDrive\DesirMED info\Paper\Cleaned_CORINE_Zonal_Stats.csv" # change it with your directory

df = pd.read_csv(file_path)

# Clean column names if needed
df.columns = [col.strip() for col in df.columns]

# Melt to long format
df_melted = df.melt(id_vars='Region', 
                    var_name='LULC Class', 
                    value_name='Proportion')

# Plot settings
sns.set(style="whitegrid")
regions = df['Region'].tolist()

n_cols = 4
n_rows = int(len(regions) / n_cols) + (len(regions) % n_cols > 0)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 12), sharex=False)

axes = axes.flatten()

for i, region in enumerate(regions):
    ax = axes[i]
    region_data = df_melted[df_melted['Region'] == region]
    
    barplot = sns.barplot(
        data=region_data, 
        y='LULC Class', 
        x='Proportion', 
        ax=ax, 
        palette='pastel'
    )
    
    ax.set_title(region, fontsize=18)
    ax.set_xlabel("Proportion (%)", fontsize=16)
    ax.set_ylabel("LULC Class", fontsize=16)
    ax.set_xlim(0, region_data["Proportion"].max() * 1.25)



    for bar in ax.patches:
        width = bar.get_width()
        ax.text(
            width + 0.5,  # Always offset to the right
            bar.get_y() + bar.get_height() / 2,
            f'{width:.2f}',
            va='center', ha='left', fontsize=8
        )

# Remove extra axes
for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

import matplotlib.pyplot as plt



fig.suptitle("LULC Proportions per Region (CORINE Level 1)", fontsize=18)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

