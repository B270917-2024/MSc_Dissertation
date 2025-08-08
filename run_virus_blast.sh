#!/bin/bash

LOGFILE="viruswise_blastp.log"
echo "Viruswise BLASTP pipeline started at $(date)" > "$LOGFILE"

start_all=$(date +%s)

declare -A viruses=(
    ["hadv"]="hadv_structural_proteins.fasta"
    ["hpv"]="hpv_structural_proteins_exp.fasta"
    ["hcov"]="hcov_structural_proteins_exp.fasta"
)

for virus in "${!viruses[@]}"; do
    input_fasta="${viruses[$virus]}"
    db_name="${virus}_blast_db"
    output_tsv="${virus}_blastp_results.tsv"

    echo "[$(date)] Processing $virus..." | tee -a "$LOGFILE"

    start_virus=$(date +%s)

    # Create a local database
    makeblastdb -in "$input_fasta" \
                -dbtype prot \
                -out "$db_name" >> "$LOGFILE"

    echo "[$(date)] Created BLAST DB: $db_name" | tee -a "$LOGFILE"

    # Set thresholds and metrics 
    blastp -query "$input_fasta" \
           -db "$db_name" \
           -out "$output_tsv" \
           -evalue 1e-5 \
           -outfmt 6 \
           -num_threads 8 >> "$LOGFILE"

    end_virus=$(date +%s)
    elapsed_virus=$((end_virus - start_virus))

    echo "[$(date)] Finished BLASTP for $virus in ${elapsed_virus}s. Output: $output_tsv" | tee -a "$LOGFILE"
    echo | tee -a "$LOGFILE"
done

end_all=$(date +%s)
elapsed_all=$((end_all - start_all))

echo "Pipeline finished at $(date)" | tee -a "$LOGFILE"
echo "Total elapsed time: ${elapsed_all}s" | tee -a "$LOGFILE"
