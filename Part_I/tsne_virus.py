import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import TSNE
from collections import Counter

# Input files
distance_files = [
    "hadv_distance_matrix.csv",
    "hcov_distance_matrix.csv",
    "hpv_distance_matrix.csv"
]

default_n_clusters = 6  # Best output after manual testing

for dist_file in distance_files:
    print(f"\nProcessing {dist_file}...")
    
    prefix = dist_file.replace("_distance_matrix.csv", "")
    output_plot = f"{prefix}_tsne_clusters.png"
    output_csv = f"{prefix}_tsne_clusters.csv"

    # Load distance matrix
    dist_matrix = pd.read_csv(dist_file, index_col=0)

    # Agglomerative clustering
    clusterer = AgglomerativeClustering(
        n_clusters=default_n_clusters,
        metric='precomputed',
        linkage='average'
    )
    labels = clusterer.fit_predict(dist_matrix.values)

    # t-SNE embedding
    tsne = TSNE(n_components=2, metric='precomputed', init='random', random_state=42)
    embedding = tsne.fit_transform(dist_matrix.values)

    # Consistent colors for clusters
    cmap = plt.get_cmap('tab10')
    colors = [cmap(i) for i in range(default_n_clusters)]

    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(embedding[:, 0], embedding[:, 1], c=[colors[label] for label in labels], s=50, alpha=0.7)

    # Legend handles
    handles = [plt.Line2D([0], [0], marker='o', color='w',
                          label=f'Cluster {i}',
                          markerfacecolor=colors[i],
                          markersize=10) for i in range(default_n_clusters)]
    plt.legend(handles=handles, title='Clusters', loc='best')

    plt.title(f't-SNE Visualisation of {prefix.upper()} Sequences Coloured by Cluster')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')

    # Removed point annotations for cluster labels

    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved t-SNE plot to: {output_plot}")

    # Save t-SNE + cluster data
    tsne_df = pd.DataFrame({
        'Sequence_ID': dist_matrix.index,
        'Cluster_Label': labels,
        'TSNE_1': embedding[:, 0],
        'TSNE_2': embedding[:, 1]
    })
    tsne_df.to_csv(output_csv, index=False)
    print(f"Saved t-SNE data to: {output_csv}")

    # Cluster composition summary
    print("\nCluster composition summary (top virus/gene per cluster):")
    for cluster_id in range(default_n_clusters):
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

print("\nAll distance matrices processed.")
