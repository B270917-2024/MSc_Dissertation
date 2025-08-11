import pandas as pd
import numpy as np
from Bio import SeqIO

# Load t-SNE cluster assignments 
clusters_df = pd.read_csv("tsne_clusters_mc_final_v1.csv")

# Load similarity matrix (assuming bitscores, will normalise) 
similarity_matrix = pd.read_csv("glycoprotein_similarity_matrix_mc.csv", index_col=0)

# Load major capsid protein IDs from FASTA 
mcp_fasta = "major_capsid_proteins.fasta"
mcp_ids = {record.id for record in SeqIO.parse(mcp_fasta, "fasta")}

# Identify protein type (glycoprotein or MCP) with stricter matching 
def get_protein_type(seq_id):
    seq_id_lower = seq_id.lower()
    protein_map = {"gb": "gB", "gd": "gD", "gh": "gH", "gl": "gL", "gm": "gM", "gn": "gN"}
    for tag, full_tag in protein_map.items():
        if f"_{tag}_" in seq_id_lower or f"_{tag}$" in seq_id_lower:
            return full_tag
    if seq_id in mcp_ids:
        return "MCP"
    return "Unknown"

clusters_df["protein_type"] = clusters_df["Sequence_ID"].apply(get_protein_type)

# Group by protein type 
protein_groups = clusters_df.groupby("protein_type")["Sequence_ID"].apply(list).to_dict()

# Average pairwise similarity within group (normalize bitscores to 0-100) 
def avg_pairwise_similarity(protein_list, matrix):
    if len(protein_list) < 1:
        return 0.0
    try:
        submatrix = matrix.loc[protein_list, protein_list]
        # Normalize bitscore to percentage (assuming max bitscore ~1000, adjust if needed)
        max_bitscore = submatrix.max().max()  # Estimate max bitscore
        if max_bitscore > 0:
            submatrix = (submatrix / max_bitscore) * 100
        upper_triangle = submatrix.where(np.triu(np.ones(submatrix.shape), k=1).astype(bool))
        mean_sim = upper_triangle.stack().mean()
        return round(mean_sim, 2) if not np.isnan(mean_sim) else 0.0
    except Exception as e:
        print(f"Error computing similarity for group: {e}")
        return 0.0

# Count unique clusters for group 
def num_clusters(protein_list, df):
    return df[df["Sequence_ID"].isin(protein_list)]["Cluster_Label"].nunique()

# Score each protein type 
score_data = []
for ptype, proteins in protein_groups.items():
    if ptype == "Unknown":
        continue
    avg_sim = avg_pairwise_similarity(proteins, similarity_matrix)
    cluster_count = num_clusters(proteins, clusters_df)
    score = avg_sim / cluster_count if cluster_count > 0 else 0
    score_data.append({
        "protein_type": ptype,
        "num_proteins": len(proteins),
        "avg_similarity": avg_sim,
        "num_clusters": cluster_count,
        "score": round(score, 2)
    })

# Create and sort the scoring DataFrame 
score_df = pd.DataFrame(score_data)
score_df = score_df.sort_values(by="score", ascending=False)

# Display results 
print("\nRanked Glycoproteins and MCP Based on Conservation and Diversity:\n")
print(score_df.to_string(index=False))

# References 

#Upper triangle method
#https://numpy.org/doc/stable/reference/generated/numpy.triu.html
#Stacking
#https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.stack.html 
