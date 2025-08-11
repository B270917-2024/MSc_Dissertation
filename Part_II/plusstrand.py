#!/usr/bin/env python3
import subprocess
from Bio import SeqIO

# Convert minus strands to plus strands
FASTA_FILES = {
    "gB": "gB_cds.fasta",
    "gH": "gH_cds.fasta",
    "MCP": "major_capsid_cds.fasta"
}

for protein, fasta in FASTA_FILES.items():
    blast_out = f"herpesvirus_output/{protein}_blast.txt"
    plus_fasta = f"herpesvirus_output/{protein}_plus.fasta"
    cmd = f"blastn -query {fasta} -db herpesvirus_output/genbank_refs/herpesviruses_db -out {blast_out} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand'"
    subprocess.run(cmd, shell=True, check=True)
    seen_ids = set()
    plus_sequences = []
    with open(blast_out) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 13: continue
            qseqid, strand = fields[0], fields[12]
            if qseqid in seen_ids: continue
            seen_ids.add(qseqid)
            for record in SeqIO.parse(fasta, "fasta"):
                if record.id == qseqid:
                    if strand == "minus":
                        record.seq = record.seq.reverse_complement()
                        print(f"{protein}: {qseqid} converted to plus-strand")
                    else:
                        print(f"{protein}: {qseqid} is plus-strand")
                    plus_sequences.append(record)
                    break
    SeqIO.write(plus_sequences, plus_fasta, "fasta")
