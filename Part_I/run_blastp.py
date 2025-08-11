import subprocess
import time

# Input and Output Files 
input_fasta = "filtered_proteins.fasta"  # Use the HHV shortlisted structural protein file 
blast_db = "filtered_db"  
blast_output = "blastp_filtered.tsv"


# Run makeblastdb
print(f"Creating BLAST database from {input_fasta}...")
subprocess.run([
    "makeblastdb",
    "-in", input_fasta,
    "-dbtype", "prot",
    "-out", blast_db
], check=True)
print(f"BLAST database created: {blast_db}")

# Time and run blastp 
start_time = time.time()
print(f"Running BLASTP on {input_fasta} with {blast_db}...")

subprocess.run([
    "blastp",
    "-query", input_fasta,
    "-db", blast_db,
    "-out", blast_output,
    "-evalue", "1e-5",
    "-outfmt", "6",
    "-num_threads", "8"
], check=True)

elapsed = time.time() - start_time
print(f"BLASTP completed in {elapsed:.2f} seconds")
print(f"BLASTP output saved to: {blast_output}")

