# eToL-V: development of a rapid method for detecting viruses  
# Dissertation

## Project Overview

The project developed a Python-based workflow to detect four human viral classes in RNA-seq data by targeting their structural proteins. It was validated using human brain samples obtained from the Edinburgh Brain Bank (EBB), including controls and samples diagnosed with neurodegenerative disorders such as **Alzheimer's disease, Lewy body disease, and vascular dementia**.

## Project Structure  
The project is divided into two main parts:  
- **Part I:** Structural protein homology analysis  
- **Part II:** Probe design and testing

## Viral Classes Covered  
- Human herpesviruses — HSV-1, HSV-2, VZV, EBV, CMV, HHV-6A, HHV-6B, HHV-7, KSHV  
- Human adenoviruses — HAdV -A, -B, -C, -D, -E, -F, -G  
- Human papillomaviruses — HPV-6, -11, -16, -18, -31, -33, -45  
- Human coronaviruses — SARS-CoV, SARS-CoV-2, HCoV-229E, HCoV-NL63, HCoV-OC43, HCoV-HKU1

# Scripts

## Part I: Structural Protein Homology Analysis Scripts

- `fetch_majorcapsid.py` and `fetch_proteins_exp.py`: Use to acquire coding region sequences of herpesviral major capsid proteins and structural proteins of adenoviruses, papillomaviruses, and coronaviruses, respectively.  
- `run_virus_blast.sh`: Runs iterative BLASTp analysis through adenoviruses, papillomaviruses, and coronaviruses (manual run needed for herpesviruses).  
- `glycoprotein_similaritydis_matrix.py` and `similaritydis_matrix.py`: Generate similarity and distance matrices for viral structural proteins and glycoproteins, assembling neighbour-joining phylogenetic trees in Newick format.  
- `tsne_herpes.py` and `tsne_virus.py`: Perform agglomerative clustering and non-linear dimensionality reduction analysis on herpesvirus and other viral class sequences.  
- `analyze_herpes.py` and `analyze_herpes_unnormalised.py`: Use to shortlist best scoring herpesvirus structural proteins from normalised and unnormalised datasets, respectively.  
- `run_blastp.py`, `filtered_similarity_matrix.py`, and `tsne_filtered_herpes.py`: Scripts for additional BLASTp analysis, filtered similarity matrix generation, and clustering/dimensionality reduction on shortlisted herpesvirus proteins.  

## Part II: Probe Design and Testing Scripts

Use the following scripts for probe design targeting the four viral classes:  
- `probes_final.py`: Generates probes for herpesviruses.  
- `probeshadv.py`: Generates probes for human adenoviruses.  
- `probeshpv2.py`: Generates probes for human papillomaviruses.  
- `probescorona.py`: Generates probes for human coronaviruses.  

Additional script:  
- `plusstrand.py`: Parses BLAST results to identify strand orientations and converts sequences to plus strand via reverse complementation.  

Heatmap scripts for results visualisation:

- `heatmap_SRAs.py`: Visualises the initial probe hit results, displaying unvalidated probe signals across all EBB samples.  
- `heatmap.py`: Generates the validated hit heatmap, showing only confirmed probe signals after artefact removal.  


`Final_probes.fasta` is the final probe list used for performance testing against EBB samples. `Flagged_probes.fasta` includes probe sequences that generated positive signals, many of which were validated as artefacts.

