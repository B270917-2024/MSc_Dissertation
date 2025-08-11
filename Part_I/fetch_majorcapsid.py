import os
import time
import logging
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Set up logging
logging.basicConfig(filename='$HOME/dissertation/trial/fetch_major_capsid.log',
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Start timer
start_time = time.time()

# Set Entrez credentials. Insert your NCBI email and key into corresponding brackets
Entrez.email = os.getenv('NCBI_EMAIL', '')
Entrez.api_key = os.getenv('NCBI_API_KEY', '')
logging.info("Entrez email set to: %s", Entrez.email)

# Herpesviruses with major capsid gene names and synonyms
herpesviruses = {
    'HSV-1': {'taxid': '10298', 'genes': {'major_capsid': ['UL19', 'major capsid protein']}},
    'HSV-2': {'taxid': '10310', 'genes': {'major_capsid': ['UL19', 'major capsid protein']}},
    'VZV': {'taxid': '10335', 'genes': {'major_capsid': ['ORF40', 'major capsid protein']}},
    'EBV': {'taxid': '10376', 'genes': {'major_capsid': ['BLLF1', 'major capsid protein']}},
    'CMV': {'taxid': '10359', 'genes': {'major_capsid': ['UL86', 'major capsid protein']}},
    'HHV-6A': {'taxid': '32603', 'genes': {'major_capsid': ['U57', 'major capsid protein']}},
    'HHV-6B': {'taxid': '32604', 'genes': {'major_capsid': ['U57', 'major capsid protein']}},
    'HHV-7': {'taxid': '32606', 'genes': {'major_capsid': ['U57', 'major capsid protein']}},
    'KSHV': {'taxid': '37296', 'genes': {'major_capsid': ['ORF25', 'major capsid protein']}},
}

output_file = '$HOME/dissertation/trial/herpesvirus_alignments2/majorcapsid/major_capsid.fasta'

def fetch_cds(taxid, gene, synonyms, virus, protein_type):
    sequences = []
    try:
        for term in [gene] + synonyms:
            search_term = f"txid{taxid}[Organism] ({term}[Gene Name] OR {term}[Protein Name])"
            handle = Entrez.esearch(db='nucleotide', term=search_term, retmax=100)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(0.1)  # NCBI rate limit
            
            for nuc_id in record['IdList']:
                try:
                    handle = Entrez.efetch(db='nucleotide', id=nuc_id, rettype='gb', retmode='text')
                    genbank = SeqIO.read(handle, 'genbank')
                    handle.close()
                    
                    for feature in genbank.features:
                        if feature.type == 'CDS':
                            gene_match = feature.qualifiers.get('gene', [''])[0].lower() == gene.lower()
                            product_match = any(term.lower() in feature.qualifiers.get('product', [''])[0].lower() for term in [gene]+synonyms)
                            if gene_match or product_match:
                                seq = feature.extract(genbank.seq)
                                strain = genbank.annotations.get('strain', genbank.id)
                                seq_id = f"{virus}_{gene}_{protein_type}_{strain}"
                                sequences.append(SeqRecord(seq=seq, id=seq_id, description=''))
                                logging.info("Fetched %s", seq_id)
                except Exception as e:
                    logging.error("Error processing nucleotide %s: %s", nuc_id, e)
                time.sleep(0.1)
    except Exception as e:
        logging.error("Error fetching %s %s: %s", virus, gene, e)
    return sequences

def main():
    unique_sequences = []
    for virus, info in herpesviruses.items():
        taxid = info['taxid']
        gene_terms = info['genes']['major_capsid']
        gene = gene_terms[0]
        synonyms = gene_terms[1:]
        logging.info(f"Fetching major capsid {virus} ({gene}, {synonyms})...")
        records = fetch_cds(taxid, gene, synonyms, virus, 'major_capsid')
        for record in records:
            seq_key = (str(record.seq), record.id)
            if seq_key not in [(str(r.seq), r.id) for r in unique_sequences]:
                unique_sequences.append(record)
    
    # Make sure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w') as out_handle:
        SeqIO.write(unique_sequences, out_handle, 'fasta')

    logging.info("Saved %d sequences to %s", len(unique_sequences), output_file)
    end_time = time.time()
    logging.info("Execution time: %.2f seconds", end_time - start_time)
    print(f"Saved {len(unique_sequences)} major capsid sequences to {output_file}")
    print(f"Execution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
