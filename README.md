# Tissue-expression-Nt-proteoforms
This repository contains code and output of the analysis of N-terminal proteoform expression in human tissues. Input dataset and large output files can be downloaded from PRIDE repository project PXD039339. This repository is part of a publication:

**[N-terminal proteoforms may engage in different protein complexes](https://pubmed.ncbi.nlm.nih.gov/37316325/)**

Annelies Bogaert 1,2 , Daria Fijalkowska 1,2, An Staes 1,2, Tessa Van de Steene 1,2, Marnik Vuylsteke 3, Charlotte Stadler 4 , Sven Eyckerman 1,2, Kerstin Spirohn 5,6,7, Tong Hao 5,6,7,  Michael A. Calderwoo d5,6,7 and Kris Gevaert 1,2,*

1 VIB Center for Medical Biotechnology, VIB, Ghent, B9052, Belgium

2 Department of Biomolecular Medicine, Ghent University, Ghent, B9052, Belgium

3 Gnomixx, Melle, Belgium

4 Department of Protein Science, KTH Royal Institute of Technology and Science for Life Laboratories, Stockholm, Sweden

5 Center for Cancer Systems Biology (CCSB), Dana-Farber Cancer Institute, Boston, USA

6 Department of Genetics, Blavatnik Institute, Harvard Medical School, Boston, USA

7 Department of Cancer Biology, Dana-Farber Cancer Institute, Boston, USA

\* To whom correspondence should be addressed. Tel: +32 (0) 9 224 98 35; Email: kris.gevaert@vib-ugent.be



## Project description
To evaluate the expression of N-terminal proteoforms (protein isoforms with different N-terminus expressed from the same gene) in healthy human tissues, we re-analyzed public proteomics data of the draft map of the human proteome developed by the Pandey group (PRIDE project PXD000561). Using ionbot search engine and a custom-build protein sequence database (composed of UniProt and Ribo-seq derived proteoforms) we confirmed that N-terminal proteoforms are not only expressed in histologically healthy human tissues, but also display tissue specificity. Our results were deposited to PRIDE with a project identifier PXD039339.

## Sample processing
Mass spectrometry data of the draft human proteome map developed by the Pandey group (1), composed of 30 histologically normal human samples including 17 adult tissues, 7 fetal tissues and 6 purified primary hematopoietic cells, were downloaded from PRIDE project PXD000561. Of the 30 samples, each was processed by several sample preparation methods and MS acquisition pipelines to generate 84 technical replicates.
Data processing
RAW files were converted to MGF format using ThermoRawFileParser (2) with default settings and searched with ionbot version 0.8.0 (3). We first generated target and decoy databases from our custom-build database containing UniProt canonical and isoform entries appended with Ribo-seq derived protein sequences (4). Next, we searched the mass spectrometry data with semi-tryptic specificity, DeepLC retention time predictions (5) and protein inference enabled, precursor mass tolerance set to 10 ppm and q-value filter of 0.01. Carbamidomethylation of cysteines was set as a fixed modification, oxidation of methionines and N-terminal acetylation were set as variable modifications and open-modification search was disabled. Downstream analysis was performed in R version 4.1.0 using dplyr (1.0.9), Biostring (2.26.0), GenomicRanges (1.46.1) and biomaRt (2.50.3, release 98). To constrict the results to first-ranked PSM per spectrum, we used ionbot.first.csv output and filtered out decoy hits, common contaminants that do not overlap with target fasta and used psm q-value <= 0.01. Due to the complexity of our custom protein database, most PSMs were associated with several protein accessions. We sorted accessions to prioritize UniProt canonical followed by UniProt isoforms, followed by Ribo-seq, higher peptide count (in the whole sample), start (smallest first), accession (alphabetically). These steps yielded a filtered PSM table. Subsequently, we sorted PSMs by N-terminal modification (to prioritize N-terminally acetylated peptidoforms) and highest PSM score. Sorted PSMs were grouped by matched peptide sequence yielding a unique peptide table. Peptides were grouped by sorted accession to generate a protein table, complemented with sample and protein metadata (such as gene and protein names, descriptions). Per sample and replicate, we obtained a unique peptide count, spectral count and NSAF (normalized spectral abundance factor) quantification.

## Statistical analysis and visualization
Differential expression analysis across all tissues was performed using limma (3.50.3) based on log2NSAF values only for proteoforms found in all replicates of at least one tissue (9,644/26,159 proteoforms). Importantly, for non-canonical proteoforms, only peptides that do not map to canonical UniProt proteins via dbtoolkit (version 4.2.5) were considered for the quantification. We extracted pairwise contrasts adult vs. adult; fetal vs. fetal and adult vs. fetal of the same tissue, considering only significant differences in expression with a Benjamini-Hochberg adjusted p-value of 0.05. Boxplots were created using ggplot2 (3.3.6), whereas heatmaps of Nt-proteoform expression were generated using pheatmap (1.0.12). To determine the row clustering, we used log2NSAF values converted to binary data as input for MONothetic Analysis (cluster version 2.1.3).

1.	M. S. Kim et al., A draft map of the human proteome. Nature 509, 575-581 (2014).
2.	N. Hulstaert et al., ThermoRawFileParser: Modular, Scalable, and Cross-Platform RAW File Conversion. J Proteome Res 19, 537-542 (2020).
3.	S. Degroeve et al., ionbot: a novel, innovative and sensitive machine learning approach to LC-MS/MS peptide identification. bioRxiv, 2021.2007.2002.450686 (2022).
4.	A. Bogaert et al., Limited Evidence for Protein Products of Noncoding Transcripts in the HEK293T Cellular Cytosol. Mol Cell Proteomics 21, 100264 (2022).
5.	R. Bouwmeester, R. Gabriels, N. Hulstaert, L. Martens, S. Degroeve, DeepLC can predict retention times for peptides that carry as-yet unseen modifications. Nat Methods 18, 1363-1369 (2021).


