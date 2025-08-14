import numpy as np
import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt


input_file = "normalised_SRAs.xlsx"
sheet_name = "Sheet4"


# Dataframe generation

df= pd.read_excel(input_file, sheet_name, index_col=0)

print(df.head())


# Heatmap 
plt.figure(figsize=(28, 40))
ax = sns.heatmap(df, annot=False, cmap='coolwarm', linewidths=0.5, cbar_kws={'label': 'Hit Count'},
            xticklabels=True, yticklabels=True, square=False)

# Increase colorbar label and tick font size
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)  
cbar.set_label('Hit Count', fontsize=20, fontweight='bold')  

# Labels and titles

plt.title('Probe Hits in AD SRA samples', fontsize=22, fontweight='bold')
plt.xlabel('Probes', fontsize=22, fontweight='bold')
plt.ylabel('SRA Samples', fontsize=22, fontweight='bold')
plt.xticks(rotation=45, ha="right", fontsize=20, fontweight='bold')  # Rotate x-axis labels 45 degrees
plt.yticks(rotation=0, fontsize=20, fontweight='bold')

# Plot
plt.tight_layout()
plt.savefig('AD_heatmap_v5.png', dpi=300, bbox_inches='tight')
plt.show()



