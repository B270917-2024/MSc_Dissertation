import os
import time
import logging
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Set up logging
logging.basicConfig(filename='$HOME/dissertation/trial/fetch_glycoproteins.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Start timer
start_time = time.time()

# Set Entrez credentials. Insert your NCBI email and key into corresponding brackets 
Entrez.email = os.getenv('NCBI_EMAIL','')
Entrez.api_key = os.getenv('NCBI_API_KEY','')
logging.info("Entrez email set to: %s", Entrez.email)

# Define herpesviruses and glycoproteins with gene names and synonyms
herpesviruses = {
    'HSV-1': {'taxid': '10298', 'genes': {
        'gB': ['UL27', 'glycoprotein B'], 
        'gH': ['UL22', 'glycoprotein H'], 
        'gD': ['US6', 'glycoprotein D'], 
        'gL': ['UL1', 'glycoprotein L'], 
        'gM': ['UL10', 'glycoprotein M'], 
        'gN': ['UL49A', 'glycoprotein N']
    }},
    'HSV-2': {'taxid': '10310', 'genes': {
        'gB': ['UL27', 'glycoprotein B'], 
        'gH': ['UL22', 'glycoprotein H'], 
        'gD': ['US6', 'glycoprotein D'], 
        'gL': ['UL1', 'glycoprotein L'], 
        'gM': ['UL10', 'glycoprotein M'], 
        'gN': ['UL49A', 'glycoprotein N']
    }},
    'VZV': {'taxid': '10335', 'genes': {
        'gB': ['ORF31', 'glycoprotein B'], 
        'gH': ['ORF60', 'glycoprotein H'], 
        'gD': ['ORF68', 'glycoprotein D'], 
        'gL': ['ORF37', 'glycoprotein L'], 
        'gM': ['ORF50', 'glycoprotein M'], 
        'gN': ['ORF9A', 'glycoprotein N']
    }},
    'EBV': {'taxid': '10376', 'genes': {
        'gB': ['BALF4', 'glycoprotein B'], 
        'gH': ['BXLF2', 'glycoprotein H'], 
        'gL': ['BDLF2', 'glycoprotein L'], 
        'gM': ['BBRF2', 'glycoprotein M'], 
        'gN': ['BMLF2', 'glycoprotein N'], 
        'gp42': ['BZLF2', 'glycoprotein 42']
    }},
    'CMV': {'taxid': '10359', 'genes': {
        'gB': ['UL55', 'glycoprotein B'], 
        'gH': ['UL75', 'glycoprotein H'], 
        'gL': ['UL115', 'glycoprotein L'], 
        'gM': ['UL100', 'glycoprotein M'], 
        'gN': ['UL73', 'glycoprotein N'], 
        'gO': ['UL74', 'glycoprotein O']
    }},
    'HHV-6A': {'taxid': '32603', 'genes': {
        'gB': ['U39', 'glycoprotein B'], 
        'gH': ['U48', 'glycoprotein H'], 
        'gL': ['U50', 'glycoprotein L'], 
        'gM': ['U72', 'glycoprotein M'], 
        'gN': ['U100', 'glycoprotein N'], 
        'gO': ['U47', 'glycoprotein O']
    }},
    'HHV-6B': {'taxid': '32604', 'genes': {
        'gB': ['U39', 'glycoprotein B'], 
        'gH': ['U48', 'glycoprotein H'], 
        'gL': ['U50', 'glycoprotein L'], 
        'gM': ['U72', 'glycoprotein M'], 
        'gN': ['U100', 'glycoprotein N'], 
        'gO': ['U47', 'glycoprotein O']
    }},
    'HHV-7': {'taxid': '32606', 'genes': {
        'gB': ['U39', 'glycoprotein B'], 
        'gH': ['U48', 'glycoprotein H'], 
        'gL': ['U50', 'glycoprotein L'], 
        'gM': ['U72', 'glycoprotein M'], 
        'gN': ['U100', 'glycoprotein N'], 
        'gO': ['U47', 'glycoprotein O']
    }},
    'KSHV': {'taxid': '37296', 'genes': {
        'gB': ['ORF8', 'glycoprotein B'], 
        'gH': ['ORF22', 'glycoprotein H'], 
        'gL': ['ORF47', 'glycoprotein L'], 
        'gM': ['ORF39', 'glycoprotein M'], 
        'gN': ['ORF53', 'glycoprotein N']
    }}
}

output_file = '$HOME/dissertation/trial/herpesvirus_alignments2/majorcapsid/new_glycoproteins.fasta'

# Store unique sequences
unique_sequences = []

def fetch_cds(taxid, gene, synonyms, virus, glycoprotein):
    sequences = []
    try:
        # Search for CDS using gene name and synonyms
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
                    
                    # Extract CDS
                    for feature in genbank.features:
                        if feature.type == 'CDS':
                            gene_match = feature.qualifiers.get('gene', [''])[0].lower() == gene.lower()
                            product_match = any(term.lower() in feature.qualifiers.get('product', [''])[0].lower() for term in [gene] + synonyms)
                            if gene_match or product_match:
                                seq = feature.extract(genbank.seq)
                                strain = genbank.annotations.get('strain', genbank.id)
                                seq_id = f"{virus}_{gene}_{glycoprotein}_{strain}"
                                sequences.append(SeqRecord(seq=seq, id=seq_id, description=''))
                                logging.info("Fetched %s", seq_id)
                except Exception as e:
                    logging.error("Error processing nucleotide %s: %s", nuc_id, e)
                time.sleep(0.1)
    except Exception as e:
        logging.error("Error fetching %s %s: %s", virus, gene, e)
    return sequences

# Fetch sequences
for virus, info in herpesviruses.items():
    taxid = info['taxid']
    for glycoprotein, terms in info['genes'].items():
        gene = terms[0]
        synonyms = terms[1:]
        logging.info("Fetching %s %s (%s, %s)...", virus, glycoprotein, gene, synonyms)
        records = fetch_cds(taxid, gene, synonyms, virus, glycoprotein)
        for record in records:
            seq_key = (str(record.seq), record.id)
            if seq_key not in [(str(r.seq), r.id) for r in unique_sequences]:
                unique_sequences.append(record)

# Write to FASTA
with open(output_file, 'w') as handle:
    SeqIO.write(unique_sequences, handle, 'fasta')

# Calculate and log execution time
end_time = time.time()
execution_time = end_time - start_time
logging.info("Saved %d sequences to %s", len(unique_sequences), output_file)
logging.info("Total execution time: %.2f seconds (%.2f minutes)", execution_time, execution_time / 60)

print(f"Saved {len(unique_sequences)} sequences to {output_file}")
print(f"Total execution time: {execution_time:.2f} seconds ({execution_time / 60:.2f} minutes)")
