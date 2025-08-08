#!/usr/bin/env python3
from Bio import SeqIO

print("Starting HPV probe design...")

fasta_path = "hpv_structural_exp_cds_cleaned.fasta"  # Input HPV FASTA
output_path = "probes.fasta"

probe_length = 120
step_size = 30
gc_min = 0.40
gc_max = 0.60

def simpson_score(seq):
    length = len(seq)
    counts = {nt: seq.count(nt) for nt in "ACGT"}
    freqs = [count / length for count in counts.values()]
    return sum(f ** 2 for f in freqs)

picked = set()  # track (virus, protein) pairs for which probe was chosen
probe_count = 0

with open(output_path, "w") as outfile:
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        rec_id = record.id

        # Parse virus and protein from record id:
        # Example ID: Human_papillomavirus_16_L1_2972221061
        parts = rec_id.split('_')
        if len(parts) < 4:
            print(f"Skipping unexpected ID format: {rec_id}")
            continue
        virus = "_".join(parts[:3])  # e.g. Human_papillomavirus_16
        protein = parts[3]           # e.g. L1 or L2

        if (virus, protein) in picked:
            continue  # Already selected probe for this virus+protein

        best_probe = None
        best_score = 1.0  # Simpson score max is 1

        for i in range(0, len(seq) - probe_length + 1, step_size):
            probe = seq[i:i+probe_length]
            if 'N' in probe:
                continue
            gc_content = (probe.count('G') + probe.count('C')) / probe_length
            if gc_min <= gc_content <= gc_max:
                score = simpson_score(probe)
                if score < best_score:
                    best_score = score
                    best_probe = probe

        if best_probe:
            probe_count += 1
            header = f">{virus}_{protein}_PROBE_{probe_count}"
            outfile.write(f"{header}\n{best_probe}\n")
            print(f"Selected probe for {virus} {protein} [Simpson={best_score:.3f}]")

            picked.add((virus, protein))

print(f"HPV probe file generation completed. Total probes: {probe_count}")
