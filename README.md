# eToL-V: development of a rapid method for detecting viruses  
# Dissertation

## 1. Structural protein homology analysis  
The scripts `fetch_glycoproteins.py` and `fetch_majorcapsid.py` are used to acquire coding region sequences of the herpesviral structural proteins. `fetch_proteins_exp.py` is used to fetch structural protein sequences of the other three viral classes analysed in this project.

`run_virus_blast.sh` runs BLASTp analysis iteratively through adenoviruses, papillomaviruses, and coronaviruses. Analogous manual run can be done for herpesviruses. `glycoprotein_similaritydis_matrix.py` and `similaritydis_matrix.py` are used to generate similarity and distance matrices for the viral classes, alongside assembling neighbour-joining based phylogenetic trees in Newick (.nwk) format. `tsne_herpes.py` & `tsne_virus.py` perform agglomerative clustering of the sequences, alongside non-linear dimensional reduction analysis. `analyze_herpes.py` and `analyze_herpes_unnormalised.py` can be utilised to shortlist the best scoring structural proteins to be used in Part II of the project.

`run_blastp.py` and `filtered_similarity_matrix.py` and `tsne_filtered_herpes.py` are scripts for additional analyses of the shortlisted herpesvirus structural proteins.

## 2. Probe design and testing  
The second part of the project included a pipeline to design and test the performance of the probe sequences targeting viral genomes in the RNA-seq samples of the human brain. `plusstrand.py` script parses BLAST results to find the strand orientations and converts them to plus strands, using reverse complementation. The following scripts can be used to generate probes for the four chosen viral classes: `probes_final.py` for herpesviruses, `probeshadv.py` for human adenoviruses, `probeshpv2.py` for human papillomaviruses, and `probescorona.py` for human coronaviruses.

`final_probes.fasta` is the final probe list that was used for the performance testing against EBB samples. `flagged_probes.fasta` includes the probe sequences that generated the positive signals of the viral genome presence in the samples. These sequences were then used to validate their results, many of which turned out to be artefacts.

**Heatmap scripts:**  
* `heatmap_SRAs.py` visualises the **initial probe hit results**, displaying unvalidated probe signals across all RNA-seq samples.  
* `heatmap.py` was used to generate the **validated hit heatmap**, showing only confirmed probe signals after artefact removal.

The objective of the project was composing a Python-based workflow for accurate and efficient detection of the viral genome in human samples. The method included detection of four human viral classes in RNA-seq data by targeting viral structural proteins. Its effectiveness was then tested in human brain samples obtained from Edinburgh Brain Bank (EBB), some of which were controls, and others had positive diagnoses of neurodegenerative disorders like **Alzheimer's disease, Lewy body disease, and vascular dementia**.

**Viral classes covered in this project**:  
- Human herpesviruses - HSV-1, HSV-2, VZV, EBV, CMV, HHV-6A, HHV-6B, HHV-7, KSHV  
- Human adenoviruses - HAdV -A, -B, -C, -D, -E, -F, -G  
- Human papillomaviruses - HPV-6, -11, -16, -18, -31, -33, -45  
- Human coronaviruses - SARS-CoV, SARS-CoV-2, HCOV-229E, HCoV-NL63, HCoV-OC43, HCoV-HKU1  

# Scripts

The overall project was divided into two parts: Part I focused on the analysis of structural protein homologies, while Part II designed and tested probes and their performance in viral genome detection.

## 1. Structural protein homology analysis  

The scripts `fetch_glycoproteins.py` and `fetch_majorcapsid.py` are used to acquire coding region sequences of the herpesviral structural proteins. `fetch_proteins_exp.py` is used to fetch structural protein sequences of the other three viral classes analysed in this project.  

`run_virus_blast.sh` runs BLASTp analysis iteratively through adenoviruses, papillomaviruses, and coronaviruses. Analogous manual run can be done for herpesviruses. `glycoprotein_similaritydis_matrix.py` and `similaritydis_matrix.py` are used to generate similarity and distance matrices for the viral classes, alongside assembling neighbour-joining based phylogenetic trees in Newick (.nwk) format. `tsne_herpes.py` & `tsne_virus.py` perform agglomerative clustering of the sequences, alongside non-linear dimensional reduction analysis. `analyze_herpes.py` and `analyze_herpes_unnormalised.py` can be utilised to shortlist the best scoring structural proteins to be used in Part II of the project.  

`run_blastp.py` and `filtered_similarity_matrix.py` and `tsne_filtered_herpes.py` are scripts for additional analyses of the shortlisted herpesvirus structural proteins.  

## 2. Probe design and testing  

The second part of the project included a pipeline to design and test the performance of the probe sequences targeting viral genomes in the RNA-seq samples of the human brain. `plusstrand.py` script parses BLAST results to find the strand orientations and converts them to plus strands, using reverse complementation. The following scripts can be used to generate probes for the four chosen viral classes: `probes_final.py` for herpesviruses, `probeshadv.py` for human adenoviruses, `probeshpv2.py` for human papillomaviruses, and `probescorona.py` for human coronaviruses.  

`final_probes.fasta` is the final probe list that was used for the performance testing against EBB samples. `flagged_probes.fasta` includes the probe sequences that generated the positive signals of the viral genome presence in the samples. These sequences were then used to validate their results, many of which turned out to be artefacts.

**Heatmap scripts:**  
* `heatmap_SRAs.py` visualises the **initial probe hit results**, displaying unvalidated probe signals across all RNA-seq samples.  
* `heatmap.py` was used to generate the **validated hit heatmap**, showing only confirmed probe signals after artefact removal.
