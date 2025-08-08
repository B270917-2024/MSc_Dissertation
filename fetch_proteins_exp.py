import os
import time
import logging
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Set up logging  
logging.basicConfig(filename='fetch_viral_cds.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Set Entrez credentials. Insert your NCBI email and key intoo corresponding brackets
Entrez.email = os.getenv('NCBI_EMAIL', '')
Entrez.api_key = os.getenv('NCBI_API_KEY', '8')
host_term = "Homo sapiens"

# Start timer
start_time = time.time()
logging.info("Entrez email set to: %s", Entrez.email)

# Function for Adenovirus
def fetch_adenovirus_cds():
    hadv_proteins = ["hexon", "penton base", "fiber", "pVII"]
    hadv_taxids = {
        "Human_adenovirus_A": "129952",
        "Human_adenovirus_B": "108098",
        "Human_adenovirus_C": "129951",
        "Human_adenovirus_D": "130293",
        "Human_adenovirus_E": "130294",
        "Human_adenovirus_F": "130295",
        "Human_adenovirus_G": "28285"
    }
    all_records = []

    for virus, taxid in hadv_taxids.items():
        for protein in hadv_proteins:
            search_term = f"txid{taxid}[Organism] AND {host_term}[Host] AND {protein}[Gene] AND cds[Feature Key]"
            logging.info("Searching adenovirus %s: %s", virus, search_term)
            try:
                handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1000)
                record = Entrez.read(handle)
                handle.close()
                ids = record["IdList"]
                logging.info("Found %d entries for %s %s", len(ids), virus, protein)
                for seq_id in ids:
                    try:
                        fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                        seq_record = SeqIO.read(fetch_handle, "genbank")
                        fetch_handle.close()
                        for feature in seq_record.features:
                            if feature.type == "CDS" and "gene" in feature.qualifiers:
                                if feature.qualifiers["gene"][0].lower() == protein.lower():
                                    cds_seq = feature.extract(seq_record.seq)
                                    record_id = f"{virus}_{protein}_{seq_id}"
                                    all_records.append(SeqRecord(cds_seq, id=record_id, description=""))
                                    logging.info("Fetched: %s", record_id)
                    except Exception as e:
                        logging.warning("Error fetching adenovirus %s %s: %s", virus, seq_id, str(e))
                    time.sleep(0.2)
            except Exception as e:
                logging.error("Error searching adenovirus %s %s: %s", virus, protein, str(e))
    return all_records

# Function for Papillomavirus
def fetch_papillomavirus_cds():
    hpv_proteins = ["L1", "L2"]
    hpv_taxids = {
        "Human_papillomavirus_6": "31542",
        "Human_papillomavirus_11": "31544",
        "Human_papillomavirus_16": "333760",
        "Human_papillomavirus_18": "333761",
        "Human_papillomavirus_31": "10585",
        "Human_papillomavirus_33": "31639",
        "Human_papillomavirus_45": "337043"
    }
    all_records = []

    for virus, taxid in hpv_taxids.items():
        for protein in hpv_proteins:
            search_term = f"txid{taxid}[Organism] AND {host_term}[Host] AND {protein}[Gene] AND cds[Feature Key]"
            logging.info("Searching papillomavirus %s: %s", virus, search_term)
            try:
                handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1000)
                record = Entrez.read(handle)
                handle.close()
                ids = record["IdList"]
                logging.info("Found %d entries for %s %s", len(ids), virus, protein)
                for seq_id in ids:
                    try:
                        fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                        seq_record = SeqIO.read(fetch_handle, "genbank")
                        fetch_handle.close()
                        for feature in seq_record.features:
                            if feature.type == "CDS" and "gene" in feature.qualifiers:
                                if feature.qualifiers["gene"][0].lower() == protein.lower():
                                    cds_seq = feature.extract(seq_record.seq)
                                    record_id = f"{virus}_{protein}_{seq_id}"
                                    all_records.append(SeqRecord(cds_seq, id=record_id, description=""))
                                    logging.info("Fetched: %s", record_id)
                    except Exception as e:
                        logging.warning("Error fetching papillomavirus %s %s: %s", virus, seq_id, str(e))
                    time.sleep(0.2)
            except Exception as e:
                logging.error("Error searching papillomavirus %s %s: %s", virus, protein, str(e))
    return all_records

# Function for Human Coronaviruses 
def fetch_coronavirus_cds():
    hcov_proteins = ["S"]
    hcov_taxids = {
        "SARS-CoV-2": "2697049",
        "SARS-CoV": "694009",
        "HCoV-229E": "11137",
        "HCoV-NL63": "277944",
        "HCoV-OC43": "31631",
        "HCoV-HKU1": "290028"
    }
    all_records = []

    for virus, taxid in hcov_taxids.items():
        for protein in hcov_proteins:
            search_term = f"txid{taxid}[Organism] AND {host_term}[Host] AND {protein}[Gene] AND cds[Feature Key]"
            logging.info("Searching coronavirus %s: %s", virus, search_term)
            try:
                handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1000)
                record = Entrez.read(handle)
                handle.close()
                ids = record["IdList"]
                logging.info("Found %d entries for %s %s", len(ids), virus, protein)
                for seq_id in ids:
                    try:
                        fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                        seq_record = SeqIO.read(fetch_handle, "genbank")
                        fetch_handle.close()
                        for feature in seq_record.features:
                            if feature.type == "CDS" and "gene" in feature.qualifiers:
                                if feature.qualifiers["gene"][0].lower() == protein.lower():
                                    cds_seq = feature.extract(seq_record.seq)
                                    record_id = f"Human_coronavirus_{virus}_{protein}_{seq_id}"
                                    all_records.append(SeqRecord(cds_seq, id=record_id, description=""))
                                    logging.info("Fetched: %s", record_id)
                    except Exception as e:
                        logging.warning("Error fetching coronavirus %s %s: %s", virus, seq_id, str(e))
                    time.sleep(0.2)
            except Exception as e:
                logging.error("Error searching coronavirus %s %s: %s", virus, protein, str(e))
    return all_records

#  Execute fetching
hadv_records = fetch_adenovirus_cds()
hpv_records = fetch_papillomavirus_cds()
hcov_records = fetch_coronavirus_cds()

# Save to separate FASTA files
output_files = {
    "hadv": "hadv_structural_exp_cds.fasta",
    "hpv": "hpv_structural_exp_cds.fasta",
    "hcov": "hcov_structural_exp_cds.fasta"
}

with open(output_files["hadv"], "w") as handle:
    SeqIO.write(hadv_records, handle, "fasta")
with open(output_files["hpv"], "w") as handle:
    SeqIO.write(hpv_records, handle, "fasta")
with open(output_files["hcov"], "w") as handle:
    SeqIO.write(hcov_records, handle, "fasta")

# Summary 
end_time = time.time()
logging.info("Saved %d hAdV sequences to %s", len(hadv_records), output_files["hadv"])
logging.info("Saved %d HPV sequences to %s", len(hpv_records), output_files["hpv"])
logging.info("Saved %d hCoV sequences to %s", len(hcov_records), output_files["hcov"])
logging.info("Total execution time: %.2f minutes", (end_time - start_time) / 60)
print(f"Saved {len(hadv_records)} hAdV sequences to {output_files['hadv']}")
print(f"Saved {len(hpv_records)} HPV sequences to {output_files['hpv']}")
print(f"Saved {len(hcov_records)} hCoV sequences to {output_files['hcov']}")
