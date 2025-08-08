#!/usr/bin/env python3
from Bio import SeqIO

print("Starting Adenovirus probe design...")

fasta_path = "pVII_rep.fasta"  # change this per protein file (hexon, penton base, pVII, fiber)
output_path = "pVII_probes.fasta" # will also change depending on the protein

probe_length = 120
step_size = 30
gc_min = 0.40
gc_max = 0.60

def simpson_score(seq):
    length = len(seq)
    counts = {nt: seq.count(nt) for nt in "ACGT"}
    freqs = [count / length for count in counts.values()]
    return sum(f ** 2 for f in freqs)

picked = set()  # track virus IDs for which probe was chosen
probe_count = 0

with open(output_path, "w") as outfile:
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        rec_id = record.id

        # For adenovirus IDs, e.g., "Human_adenovirus_B_hexon_1846449183"
        parts = rec_id.split('_')
        if len(parts) < 4:
            print(f"Skipping unexpected ID format: {rec_id}")
            continue
        virus = "_".join(parts[:3])  # e.g. Human_adenovirus_B
        protein = parts[3]           # e.g. hexon

        if virus in picked:
            continue  # Already selected probe for this virus (one probe per virus here)

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

            picked.add(virus)

print(f"Adenovirus probe file generation completed. Total probes: {probe_count}")
