#!/usr/bin/env python3
from Bio import SeqIO

print("Starting coronavirus spike probe design...")

fasta_path = "hcov_structural_exp_cds.fasta"
output_path = "corona_probes.fasta"

probe_length = 120
step_size = 30
gc_min = 0.40
gc_max = 0.60
max_probes_per_virus = 3

# Simpson Diversity Index
def simpson_score(seq):
    length = len(seq)
    counts = {nt: seq.count(nt) for nt in "ACGT"}
    freqs = [count / length for count in counts.values()]
    return sum(f ** 2 for f in freqs)

probes_by_virus = {}

for record in SeqIO.parse(fasta_path, "fasta"):
    seq = str(record.seq).upper()
    rec_id = record.id
    parts = rec_id.split('_')
    if len(parts) < 3:
        print(f"Skipping unexpected ID format: {rec_id}")
        continue
    virus = "_".join(parts[:3])
    # GC content and step size of 30 nt 
    candidates = []
    for i in range(0, len(seq) - probe_length + 1, step_size):
        probe = seq[i:i+probe_length]
        if 'N' in probe:
            continue
        gc_content = (probe.count('G') + probe.count('C')) / probe_length
        if gc_min <= gc_content <= gc_max:
            score = simpson_score(probe)
            candidates.append((score, probe))

    # Sort candidates by Simpson score ascending (best first)
    candidates.sort(key=lambda x: x[0])

    if virus not in probes_by_virus:
        probes_by_virus[virus] = []

    # Add up to max_probes_per_virus best candidates (if not already added)
    for score, probe in candidates:
        if len(probes_by_virus[virus]) >= max_probes_per_virus:
            break
        if probe not in probes_by_virus[virus]:
            probes_by_virus[virus].append(probe)

print(f"Selected up to {max_probes_per_virus} probes per virus, writing output...")

probe_count = 0
with open(output_path, "w") as outfile:
    for virus, probes in probes_by_virus.items():
        for idx, probe in enumerate(probes, 1):
            probe_count += 1
            header = f">{virus}_S_PROBE_{idx}"
            outfile.write(f"{header}\n{probe}\n")

print(f"Coronavirus probe file generation completed. Total probes: {probe_count}")
