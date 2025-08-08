import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

# List of virus-specific BLAST result files 
blast_files = {
    "hadv": "hadv_blastp_results.tsv",
    "hcov": "hcov_blastp_results.tsv",
    "hpv": "hpv_blastp_results.tsv"
}

# BLAST column names specified
cols = [
    "query_id", "subject_id", "%identity", "align_length", "mismatches",
    "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"
]

# Process each virus-specific BLAST file 
for virus, blast_file in blast_files.items():
    print(f"Processing {virus} ({blast_file})...")

    # BLAST dataframe filtering
    blast_df = pd.read_csv(blast_file, sep="\t", names=cols)
    blast_df = blast_df[blast_df["query_id"] != blast_df["subject_id"]]

    # Relaxed filters for divergent sequences
    blast_df = blast_df[
        (blast_df["%identity"] >= 20) & 
        (blast_df["evalue"] <= 1e-3)
    ]

    # Pivot into bit score matrix
    pivot_df = blast_df.pivot_table(index="query_id", columns="subject_id", values="bit_score", fill_value=0)

    # Ensure all protein IDs are present in both axes
    all_ids = sorted(set(pivot_df.index).union(set(pivot_df.columns)))
    pivot_df = pivot_df.reindex(index=all_ids, columns=all_ids, fill_value=0)

    # Save the similarity matrices (use maximum value between ij and ji)
    similarity_matrix = pivot_df.combine(pivot_df.T, func=np.maximum)
    sim_file = f"{virus}_similarity_matrix.csv"
    similarity_matrix.to_csv(sim_file)
    print(f"Saved similarity matrix to: {sim_file}")

    # Convert to distance matrix: max(bit_score) - bit_score
    bit_score_matrix = similarity_matrix.values
    valid_scores = bit_score_matrix[bit_score_matrix > 0]
    max_bit_score = np.max(valid_scores) if len(valid_scores) > 0 else 1  
    distance_matrix = np.where(bit_score_matrix == 0, max_bit_score, max_bit_score - bit_score_matrix)
    np.fill_diagonal(distance_matrix, 0)

    # Save distance matrices
    dist_file = f"{virus}_distance_matrix.csv"
    distance_df = pd.DataFrame(distance_matrix, index=all_ids, columns=all_ids)
    distance_df.to_csv(dist_file)
    print(f"Saved distance matrix to: {dist_file}")

    # Build tree using Neighbor-Joining - might need to be run in a conda environment
    lower_triangular = []
    for i in range(len(all_ids)):
        row = [float(distance_matrix[i][j]) for j in range(i + 1)]
        lower_triangular.append(row)

    dm = DistanceMatrix(all_ids, lower_triangular)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Save tree to Newick file
    tree_file = f"{virus}_tree.nwk"
    Phylo.write(tree, tree_file, "newick")
    print(f"Saved phylogenetic tree to: {tree_file}\n")

print("All virus matrices and trees have been generated.")
