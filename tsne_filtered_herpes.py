import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import TSNE
from collections import Counter

# File paths and parameters
n_clusters = 6
distance_matrix_file = "protein_distance_matrix_filtered.csv"
output_plot = "tsne_filtered.png"
output_csv = "tsne_filtered.csv"

# Load distance matrix 
dist_matrix = pd.read_csv(distance_matrix_file, index_col=0)

# Perform Agglomerative Clustering 
clusterer = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
labels = clusterer.fit_predict(dist_matrix.values)

# Run t-SNE embedding 
tsne = TSNE(n_components=2, metric='precomputed', init='random', random_state=42)
embedding = tsne.fit_transform(dist_matrix.values)

# Prepare consistent colors for clusters 
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in range(n_clusters)]

# Plot the t-SNE results 
plt.figure(figsize=(12, 8))
plt.scatter(embedding[:, 0], embedding[:, 1], c=[colors[label] for label in labels], s=50, alpha=0.7)

# Add legend based on color mapping
handles = [plt.Line2D([0], [0], marker='o', color='w',
                      label=f'Cluster {i}',
                      markerfacecolor=colors[i],
                      markersize=10) for i in range(n_clusters)]
plt.legend(handles=handles, title='Clusters', loc='best')

plt.title('t-SNE Visualisation of Filtered HHV Protein Sequences')
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')

# Save the plot 
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved t-SNE plot as {output_plot}")

# Save cluster labels and coordinates to CSV 
tsne_df = pd.DataFrame({
    'Sequence_ID': dist_matrix.index,
    'Cluster_Label': labels,
    'TSNE_1': embedding[:, 0],
    'TSNE_2': embedding[:, 1]
})
tsne_df.to_csv(output_csv, index=False)
print(f"Saved t-SNE coordinates and cluster labels to {output_csv}")

# Print top 3 virus/gene labels per cluster - adds context to t-SNE
print("\nCluster composition summary (top virus/gene per cluster):")
for cluster_id in range(n_clusters):
    cluster_seqs = dist_matrix.index[labels == cluster_id]
    simplified_labels = []
    for seq in cluster_seqs:
        parts = seq.split("_")
        if len(parts) >= 3:
            simplified_labels.append(f"{parts[0]}_{parts[2]}")
        else:
            simplified_labels.append(seq)
    most_common = Counter(simplified_labels).most_common(3)
    print(f"Cluster {cluster_id}:")
    for label, count in most_common:
        print(f"  {label}: {count} sequences")
    print()
