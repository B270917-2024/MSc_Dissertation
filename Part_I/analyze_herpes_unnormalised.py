import pandas as pd
import numpy as np
from Bio import SeqIO

# Load t-SNE cluster assignments 
clusters_df = pd.read_csv("tsne_clusters_mc_final_v1.csv")  

# Load similarity matrix 
similarity_matrix = pd.read_csv("glycoprotein_similarity_matrix_mc.csv", index_col=0)

# Load major capsid protein IDs from FASTA 
mcp_fasta = "major_capsid_proteins.fasta"
mcp_ids = {record.id for record in SeqIO.parse(mcp_fasta, "fasta")}  

# Identify protein type (glycoprotein or MCP) 
def get_protein_type(seq_id):
    for tag in ["gB", "gD", "gH", "gL", "gM", "gN"]:
        if tag.lower() in seq_id.lower():
            return tag
    if seq_id in mcp_ids:
        return "MCP"
    return "Unknown"

clusters_df["protein_type"] = clusters_df["Sequence_ID"].apply(get_protein_type)

# Group by protein type 
protein_groups = clusters_df.groupby("protein_type")["Sequence_ID"].apply(list).to_dict()

# Average pairwise similarity within group 
def avg_pairwise_similarity(protein_list, matrix):
    try:
        if len(protein_list) < 2:
            return 0.0
        submatrix = matrix.loc[protein_list, protein_list]
        upper_triangle = submatrix.where(np.triu(np.ones(submatrix.shape), k=1).astype(bool))
        return upper_triangle.stack().mean()
    except Exception as e:
        print(f"Error computing similarity for group: {e}")
        return 0.0

# Count unique clusters for group 
def num_clusters(protein_list, df):
    return df[df["Sequence_ID"].isin(protein_list)]["Cluster_Label"].nunique()

# Score each protein type 
score_data = []
for ptype, proteins in protein_groups.items():
    if ptype == "Unknown" or len(proteins) < 2:
        continue
    avg_sim = avg_pairwise_similarity(proteins, similarity_matrix)
    cluster_count = num_clusters(proteins, clusters_df)
    score = avg_sim / cluster_count if cluster_count > 0 else 0
    score_data.append({
        "protein_type": ptype,
        "num_proteins": len(proteins),
        "avg_similarity": round(avg_sim, 2),
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
