#!/usr/bin/env python3
from Bio import SeqIO

print("Starting Herpesvirus probe design with top 5 probes ranked by Simpson score...")

viruses = {
    "HSV-1": "NC_001806",
    "HSV-2": "NC_001798",
    "VZV": "NC_001348",
    "EBV": "NC_007605",
    "CMV": "NC_006273",
    "HHV-6A": "NC_001664",
    "HHV-6B": "NC_000898",
    "KSHV": "NC_009333",
    "HHV-7": "U43400.1"
}

probe_length = 120
step_size = 30
gc_min = 0.40
gc_max = 0.60

# Simpson Diversity Index 
def simpson_score(seq):
    length = len(seq)
    counts = {nt: seq.count(nt) for nt in "ACGT"}
    freqs = [count / length for count in counts.values()]
    return sum(f ** 2 for f in freqs)

output_path = "herpesvirus_output/probes_herpes_top5.fasta"
probe_count = 0

# Uses plus strands generated in the previous step
with open(output_path, "w") as outfile:
    for protein in ["gB", "gH", "MCP"]:
        print(f"Processing protein: {protein} ...")
        fasta_path = f"herpesvirus_output/{protein}_plus.fasta"  
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            rec_id = record.id

            virus = None
            for v, acc in viruses.items():
                if acc in rec_id or v in rec_id:
                    virus = v
                    break
            if not virus:
                print(f"Warning: Virus not identified in record {rec_id}, skipping.")
                continue
	    # Use G and C count & step size of 30nt here 
            candidate_probes = []
            for i in range(0, len(seq) - probe_length + 1, step_size):
                probe = seq[i:i+probe_length]
                if 'N' in probe:
                    continue
                gc_content = (probe.count('G') + probe.count('C')) / probe_length
                if gc_min <= gc_content <= gc_max:
                    score = simpson_score(probe)
                    candidate_probes.append((score, probe))

            # Sort probes by Simpson score (ascending)
            candidate_probes.sort(key=lambda x: x[0])

            # Take top 5 probes
            top_probes = candidate_probes[:5]

            for idx, (score, probe) in enumerate(top_probes, 1):
                probe_count += 1
                header = f">{virus}_{protein}_PROBE_{idx}"
                outfile.write(f"{header}\n{probe}\n")
                print(f"Selected probe {idx} for {virus} {protein} [Simpson={score:.3f}]")

print(f"Probe file generation completed. Total probes: {probe_count}")
