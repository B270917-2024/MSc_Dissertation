# eToL-V: development of a rapid method for detecting viruses
# Dissertation

## 1. Structural protein homology analysis
The scripts fetch_glycoproteins.py and fetch_majorcapsid.py are used to acquire coding region sequences of the herpesviral structural proteins. Fetch_proteins_exp.py is used to fetch structural protein sequences of the other three viral classes analysed in this project.

Run_virus_blast.sh runs BLASTp analysis iteratevely through adenoviruses, papillomaviruses, and coronaviruses. Analogous manual run can be done for herpesviruses. Glycoprotein_similaritydis_matrix.py and similaritydis_matrix.py are used to generate similarity and distance matrices for the the viral classes, alongside assembling neighbour-joining based phylogenetic trees in Newick (.nwk) format. Tsne_herpes.py & tsne_virus.py perform agglomerative clustering of the sequences, alongside non-linear dimensional reduction analysis. Analyze_herpes.py and analyze_herpes_unnormalised.py can be utilised to shortlist the best scoring structural proteins to be used in Part II of the project.

Run_blastp.py and filtered_similarity_matrix.py and tsne_filtered_herpes.py are scripts fro adcdtioanl analyses of the shortlisted herpesvirus structural proteins.

## 2. Probe design and testing
The second part of the project included a pipeline to design and test the performance of the probe sequences targeting viral genomes in the RNA-seq samples of the human brain. Plusstrand.py script parses BLAST results to find the strand orientations and converts them to plus strands, using reverse complementation. The following scripts can be used to generate probes for the four chosen viral classes: probes_final.py for herpesviruses, probeshadv.py for human adenoviruses, probeshpv2.py for human papillomaviruses, and probescorona.py for human coronaviruses.

Final_probes.fasta is the final probe list that was used for the performance testing against EBB samples. Flagged_probes.fasta includes the probe sequences that generated the positive signals of the viral genome presence in the samples. These sequences were then used to validate their results, many of which turned to be artifact.

=======
The objective of the project was composing a Python-based workflow for accurate and efficient detection of the viral genome in the human samples.


=======
The objective of the project was composing a Python-based workflow for accurate and efficient detection of the viral genome in the human samples. The method included detection of four human viral classes in RNA-seq data by targeting viral structural proteins. Its effectiveness was thhen tested in the human brain samples obtained from Edinburgh Brain Bank (EBB), some of which were controls, and others had positive diagnosis in neurodegenerative disorders like **Alzheimer's disease, Lewis body disease, and Vascular dementia**.

**Viral classes covered in this project**:
- Human herpesviruses - HSV-1, HSV-2, VZV, EBV, CMV, HHV-6A, HHV-6B, HHV-7, KSHV
- Human adenoviruses - HAdV -A, -B, -C, -D, -E, -F, -G
- Human papillomaviruses - HPV-6, -11, -16, -18, -31, -33, -45
- Human coronaviruses - SARS-CoV, SARS-CoV-2, HCOV-229E, HCoV-NL63, HCoV-OC43, HCoV-HKU1


# Scripts

The overall project was divided into two parts: Part I focused on the analysis of structural protein homologies, while Part II designed and tested probes and their performance in viral genome detection.


## 1. Structural protein homology analysis 

The scripts fetch_glycoproteins.py and fetch_majorcapsid.py are used to acquire coding region sequences of the herpesviral structural proteins. Fetch_proteins_exp.py is used to fetch structural protein sequences of the other three viral classes analysed in this project. 

Run_virus_blast.sh runs BLASTp analysis iteratevely through adenoviruses, papillomaviruses, and coronaviruses. Analogous manual run can be done for herpesviruses. Glycoprotein_similaritydis_matrix.py and similaritydis_matrix.py are used to generate similarity and distance matrices for the the viral classes, alongside assembling neighbour-joining based phylogenetic trees in Newick (.nwk) format. Tsne_herpes.py & ... perform agglomerative clustering of the sequences, alongside non-linear dimensional reduction analysis. Analyze_herpes.py and analyze_herpes_unnormalised.py can be utilised to shortlist the best scoring structural proteins to be used in Part II of the project. 

Run_blastp.py and filtered_similarity_matrix.py and tsne_filtered_herpes.py are scripts fro adcdtioanl analyses of the shortlisted herpesvirus structural proteins. 

## 2. Probe design and testing 

The second part of the project included a pipeline to design and test the performance of the probe sequences targeting viral genomes in the RNA-seq samples of the human brain. Plusstrand.py script parses BLAST results to find the strand orientations and converts them to plus strands, using reverse complementation. The following scripts can be used to generate probes for the four chosen viral classes: probes_final.py for herpesviruses, probeshadv.py for human adenoviruses, probeshpv2.py for human papillomaviruses, and probescorona.py for human coronaviruses. 

Final_probes.fasta is the final probe list that was used for the performance testing against EBB samples. Flagged_probeds.fasta includes the probe sequences that generated the positive signals of the viral genome presence in the samples. These sequences were then used to validate their results, many of which turned to be artifacts. 
