import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

# File paths
blast_file = "blastp_results_mc.tsv"
output_similarity = "glycoprotein_similarity_matrix_mc.csv"
output_distance = "glycoprotein_distance_matrix_mc.csv"
output_tree = "herpesvirus_phylogenetic_tree_mc.nwk"

# BLAST column names specified
cols = [
    "query_id", "subject_id", "%identity", "align_length", "mismatches",
    "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"
]

# BLAST dataframe filtering
blast_df = pd.read_csv(blast_file, sep="\t", names=cols)
blast_df = blast_df[blast_df["query_id"] != blast_df["subject_id"]]

# Relaxed thresholds for divergent sequences
blast_df = blast_df[
    (blast_df["%identity"] >= 20) &
    (blast_df["evalue"] <= 1e-3)
]

# Create symmetric similarity matrix (bit scores)
pivot_df = blast_df.pivot_table(index="query_id", columns="subject_id", values="bit_score", fill_value=0)

# Ensure all IDs are included in both axes
all_ids = sorted(set(pivot_df.index).union(set(pivot_df.columns)))
pivot_df = pivot_df.reindex(index=all_ids, columns=all_ids, fill_value=0)

# Save the similarity matrix
similarity_matrix = pivot_df.combine(pivot_df.T, func=np.maximum)
similarity_matrix.to_csv(output_similarity)
print(f"Similarity matrix saved to: {output_similarity}")

# Convert to distance matrix
bit_score_matrix = similarity_matrix.values
valid_scores = bit_score_matrix[bit_score_matrix > 0]
max_bit_score = np.max(valid_scores)

# Replace 0s (no hit) with max distance to avoid errors
distance_matrix = np.where(bit_score_matrix == 0, max_bit_score, max_bit_score - bit_score_matrix)
np.fill_diagonal(distance_matrix, 0)

# Save distance matrix 
distance_df = pd.DataFrame(distance_matrix, index=all_ids, columns=all_ids)
distance_df.to_csv(output_distance)
print(f"Distance matrix saved to: {output_distance}")

# Build tree (Neighbor-Joining) - this step might require a separate conda environment
lower_triangular = []
for i in range(len(all_ids)):
    row = [float(distance_matrix[i][j]) for j in range(i + 1)]
    lower_triangular.append(row)

dm = DistanceMatrix(all_ids, lower_triangular)
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Save tree in Newick format 
Phylo.write(tree, output_tree, "newick")
print(f"Phylogenetic tree saved to: {output_tree}")

# References
# https://pandas.pydata.org/docs/reference/api/pandas.pivot_table.html

