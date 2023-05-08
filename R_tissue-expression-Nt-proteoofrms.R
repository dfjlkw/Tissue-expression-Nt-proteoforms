library(tidyr)
library(stringr)
library(ggplot2)
library(Biostrings)
library(dplyr)
library(biomaRt)
library(GenomicRanges)
library(naniar)
library(limma)
library(corrplot)
library(pheatmap)
library(dendextend)
library(scales)
library(RColorBrewer)
library(reshape2)
library(imputeLCMD)
library(M3C)
library(mclust)
library(cluster)


#Analysis of pandey data search with custom db Ribo-seq + UniPRot+Iso
path.input="D:/MySharedFolder/Annelies/02082022_ionbot_pandey"

#intermediate output
dir.create(path="intermediate_results_02082022/")
path.out="intermediate_results_02082022/"

filenames=list.files(path = path.input, pattern = "ionbot.first.csv", full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
foldernames=str_split_fixed(filenames, "/", 6)[,5]
meta=as.data.frame(str_split_fixed(foldernames, "_", 5))
colnames(meta) <- c("tissue_development","tissue_type","sample_type", "ms_instrument", "mgf_file")
meta$file_number=as.numeric(meta$mgf_file)
meta$ms_sample=foldernames
meta$PSM_file=filenames

meta %>%
  arrange(tissue_development,tissue_type,sample_type, file_number ) %>%
  group_by(tissue_development,tissue_type) %>%
  mutate(rep=row_number(), count_rep=n()) %>%
  dplyr::select (-mgf_file)-> meta

#target protein fasta
fasta=readAAStringSet("D:/MySharedFolder/Annelies/Protein_fasta/combination LeeGao with UniProt canANDiso/proteoformer_LeeGao_uniprot_plusisoform_2019_04.fasta", format="fasta")
fasta.db=data.frame(len=nchar(fasta),name=names(fasta))
fasta.db$accession=str_split_fixed(fasta.db$name, "\\|", 3)[,2]
fasta.db$identifier=str_split_fixed(str_split_fixed(fasta.db$name, "\\|", 3)[,3], " ", 2)[,1]

#contaminants fasta (parsing taken from https://github.ugent.be/compomics/ionbot/blob/dev/ionbot/parse_fasta.py)
contaminants=readAAStringSet("ionbot_contaminants.fasta", format="fasta")
names_contaminants=names(contaminants)
names_contaminants=str_split_fixed(names_contaminants, " ", n=3)[,1]
contaminant_proteins=vector()
for (i in 1:length(names_contaminants)) {
  if (length(unlist(str_split(names_contaminants[i], "\\|"))) > 2) { 
    contaminant_proteins[i]= unlist(str_split(names_contaminants[i], "\\|"))[3]
  } else {
    contaminant_proteins[i]=names_contaminants[i]} 
}
contaminant_proteins_overlapping_target_fasta= contaminant_proteins[contaminants %in% fasta]
contaminant_proteins=contaminant_proteins[!contaminant_proteins %in% contaminant_proteins_overlapping_target_fasta]
rm(contaminants, names_contaminants)
#do not select column names, take all
#ionbot_match_id	spectrum_title	scan	spectrum_file	precursor_mass	peptide_mass	observed_retention_time	charge	database_peptide	matched_peptide	modifications	modifications_delta	corrected_retention_time	unexpected_modification	database	psm_score	q-value	PEP	proteins
#QUESTIONS: is unique ionbot_match_id? YES
#semi_tryptic instead of ragging
#why modifications_delta is always 0 for acetyl, but oxidation/carbramido has a value, but not always?
#REMARKS: output PSM file contains decoys AND contaminants
#REMARKS:contaminats name parsing doesn't work
#REMARKS:there is no column to easily filter contaminants

#REMARKS: protein group info in separate table ionbot.first.proteins.csv; protein_length & protein_description DO NOT MATCH uniprot_id !!! 
#REMARKS: 27214 PSMs and only 19166 have a protein group
#REAMRKS:not all peptides map to all protein in a protein group, 
#QUESTION: what are the protein inference rules?
#QUESTION: is it intentional that protein groups are named by something_HUMAN? It is the same value for uniprot isoforms, so naturally they will get the same protein group name; maybe protein group ID should be added just in case these two isoforms should be in separate protein groups?


#sel_cols=c("title","identified_first","database_peptide","matched_peptide","modifications","charge","precursor_mass","peptide_mass","retention_time_corrected","retention_time_predicted","intensity_explained","unexpected_modification","psm_q.value_first","psm_q.value","proteins","protein_group_first","protein_group_q.value_first","protein_group","protein_group_q.value")


PSM.in=list()
PSM.filt=list()
prot.group=list()
PSM.without.prot.group=list()
for (j in meta$file_number){
#for (j in c(1)){
  
  #get protein groups
  prot.group[[j]] =read.csv(file=paste0(path.input,"/",meta$ms_sample[meta$file_number == j],"/", pattern = "ionbot.first.proteins.csv"))
  prot.group[[j]]  %>%
    group_by(uniprot_id) %>% #group by protein name
    mutate(psm_count_protein=n()) %>% #get psm count for that protein
    ungroup %>%
    group_by(protein_group) %>% #group by protein group
    arrange(protein_group, desc(psm_count_protein), uniprot_id) %>% #get names of all proteins in protein group, order by PSM count and then alphabetically
    mutate(concat=paste0(protein, "((", uniprot_id, "))"),
           protein_group_members=paste(unique(concat), collapse="|"),
           protein_group_total_psm_count=length(unique(ionbot_match_id))) %>%
    ungroup() %>%
    group_by(ionbot_match_id, protein_group) %>%
    summarise_all(first) %>%
    ungroup() %>%
    group_by(ionbot_match_id) %>%
    arrange(ionbot_match_id, desc(protein_group_total_psm_count)) %>% #for each PSM with multiple protein groups, sort protein groups based on total PSM count
    summarise_all(paste, collapse = ";") %>%
    mutate(is_shared_peptide = ifelse(is_shared_peptide=="False", FALSE, TRUE))-> prot.group[[j]] 
  
  #get PSM
  PSM.in[[j]]=read.csv(file=meta$PSM_file[meta$file_number==j]) 
  PSM.in[[j]]$ms_sample=str_split_fixed(PSM.in[[j]]$spectrum_file, ".mgf", 2)[,1]
  PSM.in[[j]]$PSM_ID=1:nrow(PSM.in[[j]])
  
  #add protein group columns to PSM table
  PSM.in[[j]]$is_shared_peptide=prot.group[[j]]$is_shared_peptide[match( PSM.in[[j]]$ionbot_match_id, prot.group[[j]]$ionbot_match_id)]
  PSM.in[[j]]$protein_group=prot.group[[j]]$protein_group[match( PSM.in[[j]]$ionbot_match_id, prot.group[[j]]$ionbot_match_id)]
  PSM.in[[j]]$protein_group_members=prot.group[[j]]$protein_group_members[match( PSM.in[[j]]$ionbot_match_id, prot.group[[j]]$ionbot_match_id)]
  PSM.in[[j]]$protein_group_total_psm_count=prot.group[[j]]$protein_group_total_psm_count[match( PSM.in[[j]]$ionbot_match_id, prot.group[[j]]$ionbot_match_id)]
  PSM.in[[j]]$protein_group_q.value=prot.group[[j]]$protein_group_q.value[match( PSM.in[[j]]$ionbot_match_id, prot.group[[j]]$ionbot_match_id)]
  
  #mark contaminants that do not overlap with target fasta
  PSM.in[[j]]$contaminants=apply(sapply(contaminant_proteins, function(x) grepl(x, PSM.in[[j]]$proteins)), 1, any)
  
  #filter only valid PSMs
  #if you want to consider non-first ranked chimeric PSMS, use ionbot.lower.csv
  #If you want to constrict the results to first-ranked PSMs only, use ionbot.first.csv and filter out decoy hits (database == T) and contaminants==FALSE and use psm q-value <= 0,01
  PSM.filt[[j]]=PSM.in[[j]][PSM.in[[j]]$database=="T" & PSM.in[[j]]$contaminants==FALSE & PSM.in[[j]]$q.value<=0.01,]
  PSM.in[[j]]=nrow(PSM.in[[j]])
  row.count=nrow(PSM.filt[[j]])
  
  #in the modifications column, "" doesn't get grouped together, change to "unmodified"
  PSM.filt[[j]]$modifications=ifelse(nchar(PSM.filt[[j]]$modifications)==0,"unmodified", PSM.filt[[j]]$modifications)
  
  #mark contaminants that overlap with target fasta, do not remove them yet
  PSM.filt[[j]]$contaminants=apply(sapply(contaminant_proteins_overlapping_target_fasta, function(x) grepl(x, PSM.filt[[j]]$proteins)), 1, any)
  
  #split protein group into rows
  PSM.filt[[j]] %>% tidyr::separate_rows('proteins', sep="\\|\\|") -> PSM.filt[[j]]
  #get accession start and end
  PSM.filt[[j]]$accession=gsub("\\)\\)","",str_split_fixed(PSM.filt[[j]]$proteins, "\\(\\(", 3)[,3])
  PSM.filt[[j]]$start.end=gsub("\\)\\)","",str_split_fixed(PSM.filt[[j]]$proteins, "\\(\\(", 3)[,2])
  PSM.filt[[j]]$start=as.numeric(str_split_fixed(PSM.filt[[j]]$start.end ,"-",2)[,1])
  PSM.filt[[j]]$end=as.numeric(str_split_fixed(PSM.filt[[j]]$start.end ,"-",2)[,2])
  #correct start position for N-terminal ragging
  PSM.filt[[j]]$start_after_ragging=ifelse(grepl("semi_tryptic",PSM.filt[[j]]$modifications),PSM.filt[[j]]$start+nchar(PSM.filt[[j]]$database_peptide) - nchar(PSM.filt[[j]]$matched_peptide),PSM.filt[[j]]$start)
  
  #bug ionbot: peptide end position is shifted +1 @@@@ new on 15092022
  PSM.filt[[j]]$end=PSM.filt[[j]]$end-1
  PSM.filt[[j]]$protein_identifier=str_split_fixed(PSM.filt[[j]]$proteins, "\\(\\(", 3)[,1]
  PSM.filt[[j]]$proteins=paste0(PSM.filt[[j]]$protein_identifier, "((",PSM.filt[[j]]$start, "-", PSM.filt[[j]]$end, "))((",PSM.filt[[j]]$accession,"))" )
  
  #give a categorical order of importance to accessions
  PSM.filt[[j]]$cat_db=NA
  PSM.filt[[j]]$cat_db[grep("_aTIS_",PSM.filt[[j]]$accession)]<-3
  PSM.filt[[j]]$cat_db[grep("_CDS_",PSM.filt[[j]]$accession)]<-4
  PSM.filt[[j]]$cat_db[grep("_5UTR_",PSM.filt[[j]]$accession)]<-5
  PSM.filt[[j]]$cat_db[grep("_3UTR_",PSM.filt[[j]]$accession)]<-6
  PSM.filt[[j]]$cat_db[grep("_ntr_",PSM.filt[[j]]$accession)]<-7
  PSM.filt[[j]]$cat_db[grep("-",PSM.filt[[j]]$accession)]<-2
  PSM.filt[[j]]$cat_db[is.na(PSM.filt[[j]]$cat_db)]<-1
  PSM.filt[[j]]$cat_db[apply(sapply(contaminant_proteins_overlapping_target_fasta, function(x) grepl(x, PSM.filt[[j]]$proteins)), 1, any)]<-8 #potential contaminants with identical sequence in target fasta, get lowest priority
  #mark if a peptide is a TIS in any accession
  PSM.filt[[j]]$TIS_in_any_DB=ifelse(PSM.filt[[j]]$start_after_ragging %in% c(1,2), PSM.filt[[j]]$accession, NA)
  #calculate total peptide count per accession
  PSM.filt[[j]] %>%
    group_by(accession) %>%
    mutate(peptide_count=length(unique(matched_peptide))) %>%
    ungroup() %>%
  #sort accessions by PSM_ID, cat_db (uniprot first followed by uniprot isoforms, followed by ribo-seq), peptide_count(in the whole ms_sample), start (smallest first), accession (alphabetically) 
  #group by PSM_ID to get back to the start file
    arrange(PSM_ID, cat_db, desc(peptide_count), start, accession) %>%
    group_by(PSM_ID) %>%
    mutate(proteins=paste(proteins,collapse="||"),
           first_accession=first(accession),
           first_accession_peptide_count=first(peptide_count),
           peptide_counts=paste(peptide_count, collapse=";"),
           accessions=paste(accession, collapse=";"),
           TIS_in_any_DB=paste(unique(TIS_in_any_DB)[!is.na(unique(TIS_in_any_DB))],collapse=";")
           ) %>%
    dplyr::select(- start.end, -accession, -peptide_count, -protein_identifier) %>%
    summarise_all(first) -> PSM.filt[[j]]
  
  #bug in ionbot 0.8.0 database_peptide sequence has Leu substituted by Ile
  peptide_ranges <- GRanges(
    seqnames = fasta.db$name[match(PSM.filt[[j]]$first_accession, fasta.db$accession)],
    ranges = IRanges(start=PSM.filt[[j]]$start, end=PSM.filt[[j]]$end),
    strand = rep("*", nrow(PSM.filt[[j]])))
  PSM.filt[[j]]$database_peptide=as.character(fasta[peptide_ranges])
  #which(end(peptide_ranges) >width(fasta[as.character(seqnames(peptide_ranges))]))
  #peptide_ranges[94]
  #fasta[as.character(seqnames( peptide_ranges[94]))]
  #PSM.filt[[j]][94,]
  #FIXED!peptide end position outside protein sequence (protein Q9UQE7 length 1217, peptide DFVEDDTTHG end 1218 D:\MySharedFolder\Annelies\02082022_ionbot_pandey\Adult_Adrenalgland_bRP_Velos_1\ionbot.first.csv)
  #peptide end is always +1 more than should be based on peptide sequence reported
  
  #saveRDS(PSM.filt[[j]], file=paste0(path.out, "PSMs.filtered.accessions.sorted-", j, ".rds"))
  write.table(PSM.filt[[j]], file=paste0(path.out, "PSMs.filtered.accessions.sorted-", j, ".tsv"), sep="\t", row.names = FALSE)
  PSM.without.prot.group[[j]]=sum(is.na(PSM.filt[[j]]$protein_group))/nrow(PSM.filt[[j]])
  PSM.filt[[j]]=ifelse(nrow(PSM.filt[[j]])==row.count, nrow(PSM.filt[[j]]), "error")
  prot.group[[j]]="done"
  
}

PSM.excluded=unlist(PSM.in) - unlist(PSM.filt)

################ save image
save.image(file = "backup_01092022.RData")

rm(fasta)
gc()

pep=list()
for (j in meta$file_number){
#for (j in c(1)){
  PSM.filt[[j]]=read.table(file=paste0(path.out, "PSMs.filtered.accessions.sorted-", j, ".tsv"), sep="\t", header=TRUE)
  #PSM.filt[[j]]=readRDS(file=paste0(path.out, "PSMs.filtered.accessions.sorted-", j, ".rds"))
  PSM.filt[[j]] %>%
    #sort by seq, modif,  psm_score (highest) 
    arrange(matched_peptide, modifications, desc(psm_score)) %>%
    #group by    seq, modif,  get spectral_count per N-term mod
    group_by(matched_peptide, modifications) %>%
    mutate(modif_peptide_spectral_counts=n()) %>%
    summarise_all(first) %>%
    #N-term-mod to numerical category
    ungroup() %>%
    mutate(Nterm.Ace=ifelse(grepl("0\\|\\[1\\]acetyl", modifications ), 1,0)) %>%
    #sort by N-term order: 1:Ace 2:(anything else) 	, then sort by highest scoring
    arrange(matched_peptide, desc(Nterm.Ace), desc(psm_score)) %>%
    #group by sequence => UNIQUE  SEQUENCES
    group_by(matched_peptide) %>%
    #make a separate column for Ace peptide counts
    #calc % Ace  per uniq seq; get spectral counts
    mutate(
      Ace_spectral_counts=ifelse(Nterm.Ace==1, sum(modif_peptide_spectral_counts), 0),
      total_peptide_spectral_counts=sum(modif_peptide_spectral_counts),
      perc.Ace.peptide=round(100*Ace_spectral_counts/total_peptide_spectral_counts, 3),
      all_peptidoforms=paste0(matched_peptide,": ", paste(paste0("modif=", modifications, " spectral_count=",modif_peptide_spectral_counts), collapse = ";")),
      PSM_ID=paste(PSM_ID, collapse=";")) %>%
    summarise_all(first) -> pep[[j]]
  pep[[j]]$peptide_ID=1:nrow(pep[[j]])
  #do not group by start + accession, not compatible with protein level quantification later
  write.table(pep[[j]], file=paste0(path.out, "peptides-", j, ".tsv"), sep="\t", row.names = FALSE)
  pep[[j]]=nrow(pep[[j]])
  PSM.filt[[j]]=nrow(PSM.filt[[j]])
}

################ save image
save.image(file = "backup_01092022.RData")

pep.db=list()
for (j in meta$file_number){
pep.db[[j]]=read.table(file=paste0(path.out, "peptides-", j, ".tsv"), sep="\t", header=TRUE)
}


  
pep.db=as.data.frame(do.call(rbind,pep.db))

#################################################
#bug in ionbot 0.8.0 database_peptide sequence has Leu substituted by Ile
fasta=readAAStringSet("D:/MySharedFolder/Annelies/Protein_fasta/combination LeeGao with UniProt canANDiso/proteoformer_LeeGao_uniprot_plusisoform_2019_04.fasta", format="fasta")
matched_peptide_ranges <- GRanges(
  seqnames = fasta.db$name[match(pep.db$first_accession, fasta.db$accession)],
  ranges = IRanges(start=pep.db$start_after_ragging , end=pep.db$end),
  strand = rep("*", nrow(pep.db)))
pep.db$matched_peptide_cor_seq=as.character(fasta[matched_peptide_ranges])
rm(fasta)
######################################## map via dbtoolkit, use the corrected matched peptide! 12102022
#write.table(unique(pep.db$database_peptide), file="dbtoolkit_input.txt", row.names = FALSE,  quote=FALSE, sep="\t", col.names = FALSE)
write.table(unique(pep.db$matched_peptide_cor_seq), file="dbtoolkit_input.txt", row.names = FALSE,  quote=FALSE, sep="\t", col.names = FALSE)
map_uni<-read.csv(file="pandey_map_to_Unicanon.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
map_uni=map_uni[map_uni$Sequence!="",]
map_uni$Isoforms=sapply(map_uni$Isoforms, function(x) gsub("\\^A",";", x))
map_uni$all=paste0(map_uni$Accession, " (", map_uni$Start,"-",map_uni$Stop,");",map_uni$Isoforms)
map_uni$all=ifelse(map_uni$Isoforms=="",sapply(map_uni$all, function(x) gsub(";","", x)),map_uni$all)
#seems to be selected on start and alphabetical

pep.db$peptide_maps_to_UniProtCanon=NA
pep.db$peptide_maps_to_UniProtCanon=map_uni$all[match(pep.db$matched_peptide_cor_seq,map_uni$Sequence)]
pep.db$peptide_maps_to_UniProtCanon=gsub(" ", "", pep.db$peptide_maps_to_UniProtCanon)

#read R_pandey workspace
load("../backup_01092022.RData")
#at this stage we want to keep spectral counts/NSAF scores for non-canonical proteoforms only from their peptides not-matching UniProt via dbtoolkit
pep.db %>%
  mutate(protein_group=ifelse(is.na(protein_group),"undefined",protein_group)) %>%
  group_by(ms_sample, first_accession, matched_peptide) %>%
  mutate(matched_peptide_cor_seq_info = paste0(matched_peptide_cor_seq, "(", start_after_ragging, "-", end, ")"),
         matched_modifications = paste0(matched_peptide_cor_seq, "(", paste(unique(modifications), collapse = "_"), ")"),
         peptide_maps_to_UniProtCanon = ifelse( is.na(peptide_maps_to_UniProtCanon),
                                                NA,
                                                paste0(matched_peptide_cor_seq, "(", peptide_maps_to_UniProtCanon, ")")),
         peptide_not_in_UniProtCanon = ifelse( is.na(peptide_maps_to_UniProtCanon),
                                               paste0(matched_peptide_cor_seq, "(","none", ")"),
                                               NA)) %>%
  ungroup() %>%
  group_by(first_accession) %>%
  mutate(keep = case_when(sum(!is.na(peptide_not_in_UniProtCanon))>0 & !is.na(peptide_not_in_UniProtCanon) ~ TRUE, #if protein has non canonical peptides (any sample) AND this row is non canonical, keep it
                          sum(!is.na(peptide_not_in_UniProtCanon))==0 ~ TRUE, #if protein has only canonical peptides, keep all rows
                          sum(!is.na(peptide_not_in_UniProtCanon))>0 & is.na(peptide_not_in_UniProtCanon) ~ FALSE)) %>% #if protein has non canonical peptides BUT this row is canonical, remove it
  filter(keep == TRUE) %>%
  ungroup() %>%
  arrange(first_accession, start_after_ragging) %>%
  group_by(ms_sample, first_accession) %>%
  summarise(spectral_count=sum(total_peptide_spectral_counts),#spectra of unique peptides per first accession
            unique_peptide_count=n(),
            razor_and_unique_peptide_count=first(first_accession_peptide_count),
            cat_db=first(cat_db),
            TIS_in_any_DB=any((sum(TIS_in_any_DB!="")>0)),
            peptides=paste(unique(matched_peptide_cor_seq_info), collapse = " "),
            database_peptides = paste(unique(database_peptide), collapse = " "),
            matched_peptides = paste(unique(matched_peptide), collapse = " "),
            matched_modifications=paste(unique(matched_modifications), collapse = " "),
            PSM_IDs=paste(PSM_ID, collapse = ";"),
            peptide_maps_to_UniProtCanon=paste(unique(peptide_maps_to_UniProtCanon), collapse = " "),
            peptide_not_in_UniProtCanon=paste(unique(peptide_not_in_UniProtCanon), collapse = " "),
            accessions=paste(unique(unlist(str_split(accessions,";"))), collapse = " "),
            first_protein_group=first(protein_group)
            )-> prot.db
################ save image
save.image(file = "backup_07122022.RData")

#add protein length
prot.db$protein_length=fasta.db$len[match(prot.db$first_accession, fasta.db$accession)]
#SAF score is calculated as = spectral count of the protein/the length of the protein
#This is then normalised by = individual SAF/sum of all the SAF values in the sample
prot.db %>%
  group_by(ms_sample, first_accession) %>%
  mutate(SAF=spectral_count/protein_length) %>%
  ungroup() %>%
  group_by(ms_sample) %>%
  mutate(sum_SAF=sum(SAF))%>% #na.rm=TRUE means that protein length was not recovered from fasta, could not be
  ungroup() %>%
  mutate(NSAF=SAF/sum_SAF) -> prot.db

############## Merge with sample metadata
prot.db=merge(prot.db,meta,by="ms_sample")
prot.db$ms_sample=factor(prot.db$ms_sample, levels=meta$ms_sample)

############## add protein metadata
A1=data.frame(Accession = unique(prot.db$first_accession))

biomart_full<-read.csv(file="C:/Users/DariaF/Documents/Annelies/4. ANALYSIS_MERGED_DATAC13/Annotation biomart UniProt/mart_export_fullENS.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
biomart_uni<-read.csv(file="C:/Users/DariaF/Documents/Annelies/4. ANALYSIS_MERGED_DATAC13/Annotation biomart UniProt/mart_export_withUniProt.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
uni_annot<-read.csv(file="C:/Users/DariaF/Documents/Annelies/4. ANALYSIS_MERGED_DATAC13/Annotation biomart UniProt/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review-- (1).tab",header=TRUE, stringsAsFactors = FALSE, sep='\t')

A1$source=ifelse(grepl("ENST", A1$Accession), "ENST", "UniProt")

A1$UniProt.stable.ID=ifelse(A1$source!="ENST",sapply(A1$Accession , function(x) sapply(x, function(x) as.matrix(strsplit(x, "-")[[1]][1]))),NA)
A1$ENST=ifelse(A1$source=="ENST",sapply(A1$Accession , function(x) sapply(x, function(x) as.matrix(strsplit(x, "_")[[1]][1]))),NA)

A1$ENSG=ifelse(A1$source!="ENST",biomart_uni$Gene.stable.ID[match(A1$UniProt.stable.ID , biomart_uni$UniProtKB.Swiss.Prot.ID)],biomart_full$Gene.stable.ID[match(A1$ENST , biomart_full$Transcript.stable.ID)])
A1$Gene.name=ifelse(A1$source!="ENST",biomart_uni$Gene.name[match(A1$UniProt.stable.ID , biomart_uni$UniProtKB.Swiss.Prot.ID)],biomart_full$Gene.name[match(A1$ENST , biomart_full$Transcript.stable.ID)])
A1$ENSG_biotype=biomart_full$Gene.type[match(A1$ENSG , biomart_full$Gene.stable.ID)]

A1$ENST_biotype=ifelse(A1$source=="ENST",biomart_full$Transcript.type[match(A1$ENST , biomart_full$Transcript.stable.ID)],NA)
A1$ENST_TSL=ifelse(A1$source=="ENST",biomart_full$Transcript.support.level..TSL.[match(A1$ENST , biomart_full$Transcript.stable.ID)],NA)

mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",  host = "http://sep2019.archive.ensembl.org")
listAttributes(mart, what = c("name","description","page"))
ENS_additional=getBM(attributes = c("description", "ensembl_gene_id"),
                     filters="ensembl_gene_id",
                     values=A1$ENSG,
                     mart=mart)
A1$Description=ifelse(A1$source=="ENST",ENS_additional$description[match(A1$ENSG,ENS_additional$ensembl_gene_id)],
                      uni_annot$Protein.names[match(A1$UniProt.stable.ID,uni_annot$Entry)])

#add gene names/descriptions from UniProt annotation

A2=A1[is.na(A1$Gene.name),]
A3=A2[,colnames(A2) %in% c( "Accession", "UniProt.stable.ID")]
A3$Gene.name=uni_annot$Gene.names[match(A3$UniProt.stable.ID,uni_annot$Entry)]
A3$Description=uni_annot$Protein.names[match(A3$UniProt.stable.ID,uni_annot$Entry)]

#some proteins without genes fix manually
A3$Description[A3$Accession=="P59089"] <- "Putative uncharacterized protein encoded by LINC00205" 
A3$Gene.name[A3$Accession=="P59089"] <- "LINC00205"
A3$Gene.name[A3$Gene.name == ""] <- NA
A3$Description[A3$Description == ""] <- NA


#integrate the ENSG and Gene.name to A1 table
A1$Gene.name=ifelse(is.na(A1$Gene.name),A3$Gene.name[match(A1$Accession, A3$Accession)],A1$Gene.name)
A1$Description=ifelse(is.na(A1$Gene.name),A3$Description[match(A1$Accession, A3$Accession)],A1$Description)
rm (A2, A3)


############## Save output table
write.table(merge(A1,
                  prot.db[,!colnames(prot.db) %in% c("SAF","sum_SAF", "sample_type", "ms_instrument","count_rep", "rep")],
                  by.y="first_accession", by.x="Accession", all.y=TRUE),
            file=paste0("proteins.tsv"), sep="\t", row.names = FALSE)

write.table(merge(pep.db,
                   A1,by.x="first_accession", by.y="Accession", all.x=TRUE),
             file=paste0(path.out, "peptides.tsv"), sep="\t", row.names = FALSE)
 
############# Tar output files that are too big
 
tar(tarfile = "peptides_perMGF.tar.gz", files = list.files(path = path.out, pattern = "peptides-", full.names = TRUE), compression = "gzip")
tar(tarfile = "PSMs_perMGF.tar.gz", files = list.files(path = path.out, pattern = "PSMs.filtered.accessions.sorted-", full.names = TRUE), compression = "gzip")
tar(tarfile = "peptides.tar.gz", files = paste0(path.out, "peptides.tsv"), compression = "gzip")
 
########## Remove files if they were compressed
if (all(file.exists(list.files(path = path.out, pattern = "peptides-", full.names = TRUE))) & file.exists(paste0(path.out, "peptides.tsv")) & file.exists("peptides.tar.gz") & file.exists("peptides_perMGF.tar.gz")) {
   #Delete file if it exists
   file.remove(paste0(path.out, "peptides.tsv"))
   file.remove(list.files(path = path.out, pattern = "peptides-", full.names = TRUE))
 }
 
 if (all(file.exists(list.files(path = path.out, pattern = "PSMs.filtered.accessions.sorted-", full.names = TRUE))) & file.exists("PSMs_perMGF.tar.gz")) {
   #Delete file if it exists
   file.remove(list.files(path = path.out, pattern = "PSMs.filtered.accessions.sorted-", full.names = TRUE))
 }
 

############## plot before normalization
ggplot(data=prot.db, aes(ms_sample, log2(NSAF)))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=90, size=4))
ggplot(data=prot.db, aes(log2(NSAF)))+
  geom_histogram()+
  geom_text(data=data.frame() ,aes(x=-7, y=c(40000, 37000), label= c("distribution is not normal","Shapiro test p-value 1.639e-10")))+
  ylab("count")
# r sample multiple times without replacement
shapiro.test(sample (log2(prot.db$NSAF),  size=5000, replace =F))
#p-value = 1.639e-10

############ go to wide format
prot.db %>% 
  arrange(ms_sample) %>%
  pivot_wider(id_cols = first_accession , names_from= ms_sample, values_from= unique_peptide_count) -> prot.pep.count
prot.pep.count=as.data.frame(prot.pep.count)
row.names(prot.pep.count)<-prot.pep.count[,1]
prot.pep.count=as.matrix(prot.pep.count[-1])

prot.db %>% 
  arrange(ms_sample) %>%
  pivot_wider(id_cols = first_accession , names_from= ms_sample, values_from= spectral_count) -> prot.spectral.count
prot.spectral.count=as.data.frame(prot.spectral.count)
row.names(prot.spectral.count)<-prot.spectral.count[,1]
prot.spectral.count=as.matrix(prot.spectral.count[-1])


prot.db %>% 
  arrange(ms_sample) %>%
  pivot_wider(id_cols = first_accession , names_from= ms_sample, values_from= NSAF) -> prot.nsaf
prot.nsaf=as.data.frame(prot.nsaf)
row.names(prot.nsaf)<-prot.nsaf[,1]
prot.nsaf=as.matrix(prot.nsaf[-1])


colnames(prot.spectral.count)==meta$ms_sample

########
prot.db%>%
  mutate(razor_and_unique_peptide_count=paste0(ms_sample,":", razor_and_unique_peptide_count),
         PSM_IDs=paste0(ms_sample,":", PSM_IDs)) %>%
  arrange(ms_sample) %>%
  dplyr::select(-c(SAF, sum_SAF, file_number, rep, ms_instrument, sample_type, count_rep)) %>%
  pivot_wider(id_cols = c(first_accession,cat_db) , unused_fn = function(x) paste(unique(unlist(str_split(x,"[ ]+"))), collapse = ";"), names_from= ms_sample, values_from= c(NSAF, unique_peptide_count, spectral_count)) %>%#id_cols Defaults to all columns in data except for the columns specified in names_from and values_from
  group_by(first_accession) %>%
  mutate(sample_count=length(unlist(str_split(PSM_file, ';'))), 
         TIS_in_any_DB=grepl("TRUE",TIS_in_any_DB )) -> prot.wide 

#pretty peptide_maps_to_UniProtCanon and peptide_not_in_UniProtCanon
#prot.db$peptide_maps_to_UniProtCanon[prot.db$first_accession =="ENST00000397561_9_14922268_ntr_100db1"]
#prot.wide$peptide_maps_to_UniProtCanon[prot.wide$first_accession =="ENST00000397561_9_14922268_ntr_100db1"]
prot.wide$peptide_maps_to_UniProtCanon=gsub( "^NA;", "",
      gsub(";NA;", ";",
           gsub(";NA$", "", prot.wide$peptide_maps_to_UniProtCanon)))
prot.wide$peptide_not_in_UniProtCanon=gsub( "^NA;", "",
            gsub(";NA;", ";",
                 gsub(";NA$", "", prot.wide$peptide_not_in_UniProtCanon)))

################ save image
save.image(file = "backup_07122022.RData")

############## log2NSAF

############## #normalize the medians in log2NSAF?
pdf(file="1.plots_normalization.pdf", width = 15, height = 4)
par(cex=0.4,mar = c(15,2,5,2) + 0.1)
boxplot(prot.nsaf,las=2, main="NSAF distribution BEFORE log2")
prot.nsaf=log2(prot.nsaf)
boxplot(prot.nsaf,las=2, main="log2(NSAF) distribution BEFORE median normalization")
M=apply(prot.nsaf, 2, median, na.rm=TRUE)
M.adj=M-M[which.min(abs(median(M)-M))]
prot.nsaf.norm=prot.nsaf
for(i in 1:ncol(prot.nsaf.norm)){prot.nsaf.norm[,i]=prot.nsaf.norm[,i]-M.adj[i]}
boxplot(prot.nsaf.norm,las=2, main="log2(NSAF) distribution AFTER median normalization")
dev.off()
#decided median normalization is not needed
rm(prot.nsaf.norm)

############count missing values
#install.packages("naniar")

pdf(file="2.plots_missing_values.pdf", width = 8, height = 4)
v<-vis_miss(as.data.frame(prot.nsaf), 'warn_large_data' = FALSE, cluster = TRUE) 
v + 
  theme(axis.text.x = element_text(angle = 90, size = 6))

hist(rowSums(is.na(prot.nsaf)), main="number of missing NSAF values per protein", xlab ="# missing values/protein", breaks=84)
ggplot(data=prot.db, aes(log2(unique_peptide_count),fill=cat_db==1))+
  geom_density(position="identity", alpha=0.5)+
  theme_classic()+
  labs(fill="UniProt canonical")+
  ggtitle("no major bias in peptide count for UniProt canonical vs. other proteins\nwe can use one threshold to filter low expressed")

#how many proteins below 50% missing data
sum(rowSums(is.na(prot.nsaf))<50) #4437
sum(rowSums(is.na(prot.nsaf))<50)/nrow(prot.nsaf) #17%

#### keep proteins with at least two peptide IDs in at least one ms_sample
keep <- rowSums(prot.pep.count >= 2, na.rm=TRUE) >= 1 #corrected
prot.nsaf[keep,]
  
hist(rowSums(is.na(prot.nsaf[keep,])), main="number of missing NSAF values per protein\nkeep proteins with: 2 peptides in 1 sample",
     xlab ="# missing values/protein", breaks=84)
#how many proteins below 50% missing data
sum(rowSums(is.na(prot.nsaf[keep,]))<50) #4391
sum(rowSums(is.na(prot.nsaf[keep,]))<50)/nrow(prot.nsaf[keep,]) #34%


#### keep proteins found in all replicates of at least one tissue_developement+tissue_type group
prot.db %>% 
  group_by(tissue_development,tissue_type,first_accession) %>%
  summarise(coverage.per.sample.group=n()/first(count_rep),
           count.per.sample.group=n(),
           count_rep=first(count_rep)) -> per.sample.group
length(unique(per.sample.group$first_accession[per.sample.group$coverage.per.sample.group==1])) #9727 proteins covered in 100% replicates of at least one tissue
keep2=row.names(prot.nsaf) %in% unique(per.sample.group$first_accession[per.sample.group$coverage.per.sample.group==1])

hist(rowSums(is.na(prot.nsaf[keep2,])), main="number of missing NSAF values per protein\nkeep proteins present in all replicates of 1 tissue",
     xlab ="# missing values/protein", breaks=84)
#how many proteins below 50% missing data
sum(rowSums(is.na(prot.nsaf[keep2,]))<50) #4436
sum(rowSums(is.na(prot.nsaf[keep2,]))<50)/nrow(prot.nsaf[keep2,]) #46%


prot.nsaf[keep & keep2,]

hist(rowSums(is.na(prot.nsaf[keep & keep2,])), main="number of missing NSAF values per protein\nkeep proteins with: 2 peptides in 1 sample\n& present in all replicates of 1 tissue",
     xlab ="# missing values/protein", breaks=84)
#how many proteins below 50% missing data
sum(rowSums(is.na(prot.nsaf[keep & keep2,]))<50) #4390
sum(rowSums(is.na(prot.nsaf[keep & keep2,]))<50)/nrow(prot.nsaf[keep & keep2,]) #49%
dev.off()


############## Save output table filtered
prot.wide$non_canonical_proteoforms=! (prot.wide$cat_db %in% c(1,8))
prot.wide$filter_min2peptides_1sample=keep
prot.wide$filter_all_replicates_1tissue=keep2


prot.wide.export=merge(prot.wide,
                       A1,by.x="first_accession", by.y="Accession", all.x=TRUE)
expression_columns = grep("NSAF_|unique_peptide_count_|spectral_count", colnames(prot.wide.export))
other_columns=(1:ncol(prot.wide.export))[!(1:ncol(prot.wide.export) %in% expression_columns)]
prot.wide.export=prot.wide.export[,c(other_columns, expression_columns)]
colnames(prot.wide.export)


write.table(prot.wide.export, file=paste0("proteins_wide.tsv"), sep="\t", row.names = FALSE, quote=TRUE)

################ save image
save.image(file = "backup_07122022.RData")
########

############### END HERE
#limma without imputation
#impute data from truncated distribution
#

#entire pandey, run limma,  select non-canonical  and save plots
if(TRUE){
  #differentially expressed in tissues? consider filter_all_replicates_1tissue
  sum(keep2) #9644
  sum(keep2 & prot.wide$peptide_not_in_UniProtCanon != "NA") #511
  sum(keep2 & keep & prot.wide$peptide_not_in_UniProtCanon != "NA") #291
  #how much missing data per protein?
  median(rowSums(is.na(prot.nsaf))/ncol(prot.nsaf)) #median 94% missing in whole pandey
  median(rowSums(is.na(prot.nsaf[keep2,]))/ncol(prot.nsaf)) #63%
  median(rowSums(is.na(prot.nsaf[keep2 & keep,]))/ncol(prot.nsaf)) #61% second filter doesn't help much
  #decide to do limma in keep2=filter_all_replicates_1tissue
  #prot.nsaf is already log2 transformed
  #library(limma)
  f = as.factor(apply(str_split_fixed(colnames(prot.nsaf), pattern = "_", 3)[,1:2],1, paste, collapse="_"))
  design <- model.matrix(~0+f)
  colnames(design) <- levels(f)
  # make only contrasts adult vs. adult; fetal vs. fetal; adult the same tissue vs. fetal the same tissue
  contrasts =as.data.frame( t(combn(levels(f),2)))
  contrasts$both=grepl("Adult_",contrasts[,1] ) == grepl("Adult_",contrasts[,2] )
  contrasts$same_tissue=str_split_fixed(contrasts[,1], pattern = "_", 2)[,2] == str_split_fixed(contrasts[,2], pattern = "_", 2)[,2]
  contrasts=contrasts[contrasts$both | contrasts$same_tissue,1:2]
  contrasts= apply(contrasts, 1, paste, collapse = "-")
  #fit model
  fit <- lmFit(prot.nsaf[keep2,], design)
  contMatrix <- makeContrasts(contrasts=contrasts, levels=design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  summ <- summary(decideTests(fit2, adjust.method="BH", p.value=0.05))
  
  ##How much filtering to do? The mean-variance plots will guide you. For limma and log-counts, look at plotSA(fit) where fit is output from lmFit and eBayes. The trend line should be smooth. If it increases sharply at the low count end, then you haven't filtered enough.
  results.anova = topTable(fit2,sort="none",n=Inf,adjust="BH")
  pdf(file="3.QualityControlLimmaPlots.pdf", width = 8, height = 8)
  limma::plotMA(fit2, main="limma ANOVA-like", status=results.anova$adj.P.Val<0.05)
  hist(results.anova$P.Value, main="p-values histogram\nlimma ANOVA-like")
  plotSA(fit2,main="limma ANOVA-like")
  plotMDS(prot.nsaf[keep2,],main="limma ANOVA-like", top=500, gene.selection = "common", labels = str_split_fixed(colnames(prot.nsaf), pattern = "_", 3)[,2] )
  summ= as.data.frame(summ)
  colnames(summ)<-c("test_result","contrasts","count")
  summ %>%
    filter(test_result!="NotSig")%>%
    group_by(contrasts) %>%
    mutate(n=sum(count)) %>%
    arrange(desc(n)) %>%
    ungroup() %>%
    mutate(pair1=str_split_fixed(contrasts, "-",2)[,1],
           pair2=str_split_fixed(contrasts, "-",2)[,2]) %>%
    select(c(pair1, pair2, n)) %>%
    distinct()%>%
    pivot_longer(-n, names_to = "pair", values_to="sample")%>%
    group_by(sample)%>%
    summarise(n=sum(n)) %>%
    arrange(desc(n)) %>%
    mutate(sample = factor(sample, levels = unique(sample))) -> summ
  plot(ggplot(data = summ, aes(sample, n))+
    geom_bar(position="stack", stat="identity")+ 
    theme(axis.text.x = element_text(angle = 90, size = 8))+
    ggtitle("Samples ~ significant test result count"))
  
  
  dev.off()
  rm(results.anova)
  
  #get pvalues from pairwise contrasts only for significant proteins
  decide = as.data.frame(decideTests(fit2, adjust.method="BH", p.value=0.05) )
  decide$protein=rownames(decide)
  decide %>%
    pivot_longer(-protein, names_to = "contrast", values_to = "test_p005", values_drop_na = TRUE ) %>%
    filter(test_p005!=0) -> decide
  decide$adj.P.Val = NA
  decide$logFC = NA
  for (i in 1:nrow(decide)){
    res=topTable(fit2,sort="none",n=Inf,adjust="BH", coef = decide$contrast[i])
    decide$adj.P.Val[i] = signif( res$adj.P.Val[rownames(res) == decide$protein[i] ], 3)
    decide$logFC[i] = signif (res$logFC[rownames(res) == decide$protein[i] ], 3)
    rm(res)
  }
  decide$tissues_logFC_adjPval=paste0(decide$contrast, " ", decide$logFC, " ", decide$adj.P.Val)
  decide$tissues_logFC_adjPval=sub("-", "/", decide$tissues_logFC_adjPval)
  
  decide %>%
    group_by(protein) %>%
    summarize (tissuesCompared_logFC_adjPval=paste(tissues_logFC_adjPval, collapse = ";"), 
               count_significant_comparisons= n()) -> decide.short
  #export table pandey after limma
  limma.export= prot.nsaf
  colnames(limma.export)=paste0("log2NSAF_", colnames(limma.export))
  select_cols=c("first_accession","Gene.name",	"Description" ,"ENSG" ,"peptides", "peptide_maps_to_UniProtCanon",	"peptide_not_in_UniProtCanon","matched_modifications","cat_db", 	"accessions",	"filter_min2peptides_1sample",	"filter_all_replicates_1tissue","sample_count")		
  
  limma.export=merge(prot.wide.export[,select_cols],
                     limma.export, by.x="first_accession", by.y="row.names", all.x=TRUE)
  limma.export=merge(decide.short,
                     limma.export,
                     by.x="protein", by.y="first_accession", all=TRUE)
  
  
  write.table(limma.export, file=paste0("limma.export.tsv"), sep="\t", row.names = FALSE, quote=TRUE)
  
  
  
  
  #NSAF values of significant proteins
  limma.sign = prot.nsaf[rownames(prot.nsaf) %in% decide.short$protein,]
  #correlation plot for samples
  cor=cor(prot.nsaf, method = "pearson", use = "complete.obs")
  cor.sign=cor(limma.sign, method = "pearson", use = "complete.obs")
  col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                "#4393C3", "#2166AC", "#053061")))
  col2 <- colorRampPalette(rev(c("#B2182B", "#FDDBC7", "#FFFFFF", "#D1E5F0",  "#2166AC")))
  
  pdf(file="4.pearson_correlation.pdf", width = 8, height = 8)
  corrplot(cor, title = "based on all data", type = "upper", order = "hclust", 
           tl.col = "black", col=col2(256), tl.cex = 0.4)
  corrplot(cor.sign,title = "based on differential proteins",
           type = "upper", order = "hclust", 
           tl.col = "black", col=col2(256), tl.cex = 0.4)
  dev.off()
  
  
  #heatmap for significant
  {
    cal_z_score <- function(x){
      (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
    }
    table(rowSums(is.na(limma.sign)))#max 81/84 cols NA
    sort(apply(limma.sign, 2, function(x) sum(is.na(x)))/nrow(limma.sign))#max 80% missing per sample
    #pheatmap(limma.sign) #fails
    #are there rows or cols with 0 variance?
    rownames(limma.sign)[apply(limma.sign, 1, function(x) var(x, na.rm=TRUE)) == 0] #no
    colnames(limma.sign)[apply(limma.sign, 2, function(x) var(x, na.rm=TRUE)) == 0] #no
    rem=which(rowSums(is.na(limma.sign)) >=50)
    #50 missing is a problem
    #have to filter to get any histogram
    nrow(limma.sign)#5409 significant
    length(rem) #filter out 1808 with >= 50/84 values missing
    pheatmap(limma.sign[-rem,])
    
    
    #impute per row(per protein) from truncated distribution
    #limma.sign_imp.backup=limma.sign_imp
    limma.sign_imp=t(impute.QRILC(t(limma.sign))[[1]])
    #normality test
    norm.test=sapply(as.data.frame(t(limma.sign)), shapiro.test)
    norm.test.imp=sapply(as.data.frame(t(limma.sign_imp)), shapiro.test)
    table(norm.test["p.value",]>0.05) #before imputation 3387/3397+2022 are normally distributed
    table(norm.test.imp["p.value",]>0.05) #after imputation 3711/ are normally distributed
    #imputed values vs real values
    only.imp=limma.sign_imp
    only.imp[!is.na(limma.sign)] <- NA
    #per protein how often imputed value max is higher than real value min or median
    table(apply(only.imp, 1, max, na.rm=TRUE) >   apply(limma.sign, 1, min, na.rm=TRUE)) # often! 4679/5409
    table(apply(only.imp, 1, max, na.rm=TRUE) >   apply(limma.sign, 1, median, na.rm=TRUE)) #never
    #histogram after imputation per protein
    pheatmap(limma.sign_imp)
    #histogram after imputation and zscoring
    limma.sign_norm <- t(apply(limma.sign_imp, 1, cal_z_score))
    pheatmap(limma.sign_norm)
    
    
    #limma.sign_norm <- t(apply(limma.sign, 1, cal_z_score))
    #pheatmap(limma.sign_norm[-rem,])
    
    
    #get tissue order in columns
    col_cluster=pheatmap(limma.sign_norm)
    col_cluster=col_cluster$tree_col
    col_cluster= str_split_fixed(col_cluster$labels[col_cluster$order], "_", 3)[,2]
    #split similar tissues into systems
    systems=data.frame(systems=c( "white blood cells","white blood cells","white blood cells","white blood cells","white blood cells",
                                  "heart",
                                  "digestive system", "digestive system", "digestive system","digestive system", "digestive system","digestive system", 
                                  "liver", "lung",
                                  "nervous system","nervous system","nervous system","nervous system",
                                  "reproductive system","reproductive system","reproductive system",
                                  "adrenalgland",   "platelets",
                                  "urinary system","urinary system",
                                  "placenta"),
                       tissue=c("NKcells"   ,     "Monocytes"    ,  "CD8Tcells"   ,   "Bcells", "CD4Tcells",
                                "Heart",
                                "Pancreas","Colon" ,"Gut","Gallbladder" ,   "Rectum", "Esophagus" ,
                                "Liver","Lung",
                                "Spinalcord" ,"Frontalcortex", "Retina" ,"Brain",
                                "Prostate" ,"Ovary"     ,     "Testis",
                                "Adrenalgland",   "Platelets",
                                "Urinarybladder", "Kidney",
                                "Placenta"))
    
    
    #hierarchical clustering
    my_hclust_gene <- hclust(dist(limma.sign_norm), method = "complete")
    
    #dendogram
    as.dendrogram(my_hclust_gene) %>%
      plot(horiz = TRUE)
    # select the number of clusters
    CLUSTERS=10
    my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = CLUSTERS)
    my_gene_col = data.frame(Cluster=paste0("Cluster ",my_gene_col), row.names = names(my_gene_col))
    #annotate rows - choose cluster colors
    hex_codes <- hue_pal()(CLUSTERS)
    names(hex_codes)=unique(my_gene_col$Cluster)
    # annotate columns
    annotation_col = as.data.frame(str_split_fixed(colnames(limma.sign_norm), "_",3)[,1:2]) #!!!!!! velos vs elite effect (first low abundant samples arevelos)
    rownames(annotation_col)=colnames(limma.sign_norm)
    colnames(annotation_col)=c("Development", "Tissue")
    annotation_col$Tissue_type=systems$systems[match(annotation_col$Tissue, systems$tissue)]
    annotation_col=annotation_col[,c("Development","Tissue_type" )]
    #colors
    ann_colors = list(
      Cluster = hex_codes,
      Development = c("white", "pink"),
      Tissue_type = RColorBrewer::brewer.pal(length(unique(annotation_col$Tissue_type )), "Paired")#hue_pal(l = 90)(length(unique(annotation_col$Tissue_type )))
    )
    names(ann_colors[["Development"]]) <- sort(unique(annotation_col$Development ))
    names(ann_colors[["Tissue_type"]]) <- unique(systems$systems)
    #enhance expression colors
    breaksList = seq(-2, 2, by = 0.1)
    pal=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
    
     
    ph5=pheatmap(limma.sign_norm, 
             annotation_row = my_gene_col,
             annotation_colors = ann_colors,
             show_colnames=FALSE, show_rownames=FALSE,annotation_col=annotation_col,
             cutree_rows = CLUSTERS,
             color = pal,
             breaks = breaksList)
    dev.off()
    pdf(paste("5.heatmap.significant.NSAF.pdf", sep=""),width=5,height=8) 
    print(ph5)
    dev.off()
    
    #save row order based on clusters
    row.order=data.frame(name=my_hclust_gene$labels[c(my_hclust_gene$order)])
    row.order$cluster = my_gene_col$Cluster[match(row.order$name, rownames(my_gene_col))]
    row.order$number=as.numeric(str_split_fixed(row.order$cluster, " ", 2)[,2])
    #save column order
    col.order=ph5$tree_col$labels[ph5$tree_col$order ]
    
    #histogram only for non-canonical significant
    non.can=prot.wide.export$first_accession[prot.wide.export$peptide_not_in_UniProtCanon!="NA"]
    non.can=which(rownames(limma.sign_norm) %in% non.can) #101
    limma.sign_norm.noncan=limma.sign_norm[non.can,]
    
    pdf(paste("6.heatmap.noncanonical.significant.NSAF.pdf", sep=""),width=5,height=8)  
    order=row.order$name[row.order$name %in% rownames(limma.sign_norm.noncan)]
    pheatmap(limma.sign_norm.noncan[order,col.order], 
             annotation_row = my_gene_col,
             annotation_colors = ann_colors,
             show_colnames=FALSE, show_rownames=FALSE,annotation_col=annotation_col,
             #cutree_rows = CLUSTERS,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = pal,
             border_color = FALSE,
             breaks = breaksList)
    dev.off()
    
    #histogram without imputation, after zscoring
    limma.sign_norm_noimp <- t(apply(limma.sign, 1, cal_z_score))
    
    #plot heatmap with missing values using row clustering of the imputed data
    #dev.off()
    pdf(paste("7.heatmap.noImputation.significant.NSAF.pdf", sep=""),width=5,height=8)  
    pheatmap(limma.sign_norm_noimp[row.order$name,col.order], 
             annotation_row = my_gene_col,
             annotation_colors = ann_colors,
             show_colnames=FALSE, show_rownames=FALSE,annotation_col=annotation_col,
             #cutree_rows = CLUSTERS,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = pal,
             border_color = FALSE,
             breaks = breaksList)
    dev.off()
    
    #plot heatmap noncanonical with missing values using row clustering of the imputed data
    limma.sign_norm_noimp.noncan=limma.sign_norm_noimp[non.can,]
    #dev.off()
    pdf(paste("8.heatmap.noncanonical.noImputation.significant.NSAF.pdf", sep=""),width=5,height=8)  
    
    pheatmap(limma.sign_norm_noimp.noncan[order,col.order], 
             annotation_row = my_gene_col,
             annotation_colors = ann_colors,
             show_colnames=FALSE, show_rownames=FALSE,annotation_col=annotation_col,
             #cutree_rows = CLUSTERS,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = pal,
             border_color = FALSE,
             breaks = breaksList)
    dev.off()

  }

}


#N-terminal COFRADIC data overlap pandey, save plots
if(TRUE){
  #on 02092022 load R output of KNIME pipeline Annelies
  "C:/Users/DariaF/Documents/R/R_Annelies_C13_1extrapeptide/12082022-C13-FullMergedResults-R_output.txt"
  Nter =read.csv(file="C:/Users/DariaF/Documents/R/R_Annelies_C13_1extrapeptide/12082022-C13-FullMergedResults-R_output.txt", sep = "\t", header = TRUE)
  
  #overlap if pandey data has either:
  #A) the same peptide as COFRADIC
  Nter.A=sapply(pep.db$matched_peptide_cor_seq, function(x) x %in%  Nter$Sequence)
  #B) the same accession + start combination
  a = paste0(pep.db$first_accession,"_", pep.db$start_after_ragging)
  b = paste0(Nter$Accession,"_",Nter$Start)
  Nter.B=sapply(a, function(x) x %in% b)
  rm(a,b)
  #final accession selection
  Nter.acc.sel = pep.db$first_accession[Nter.A | Nter.B]
  Nter.prot.wide.all=prot.wide.export[prot.wide.export$first_accession %in% Nter.acc.sel,]
  #ALL N-term COFRADIC overlap with pandey
  nrow(Nter.prot.wide.all) #897 proteins
  #N-term COFRADIC overlap with pandey not mapping to UniProtCanon
  sum(Nter.prot.wide.all$peptide_not_in_UniProtCanon != "NA") #24 proteins
  Nter.noncan.acc.sel=Nter.prot.wide.all$first_accession[Nter.prot.wide.all$peptide_not_in_UniProtCanon != "NA"]
  #plot N-term COFRADIC overlap with pandey not mapping to UniProtCanon #NSAF, no limma, no imputation
  Nter.nsaf=prot.nsaf[Nter.noncan.acc.sel,]
  #cluster as binary data (to avoid imputation)
  Nter.noncan.order=cluster::mona(!is.na(Nter.nsaf[,colSums(!is.na(Nter.nsaf))>0]))
  pdf(paste("9A.heatmap.Nter.noncanon.noImputation.NSAF.pdf", sep=""),width=10,height=5)  
  
  pheatmap(Nter.nsaf[Nter.noncan.order$order, rownames(annotation_col)[order(annotation_col$Tissue_type)]], 
           annotation_colors = ann_colors,
           show_colnames=FALSE, show_rownames=TRUE,annotation_col=annotation_col,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = FALSE)
  dev.off()
  #try spectral counts instead
  Nter.counts=prot.spectral.count[Nter.noncan.acc.sel,]
  pdf(paste("9B.heatmap.Nter.noncanon.noImputation.spectral.counts.pdf", sep=""),width=10,height=5)  
  pheatmap(Nter.counts[Nter.noncan.order$order, rownames(annotation_col)[order(annotation_col$Tissue_type)]], 
           annotation_colors = ann_colors,
           show_colnames=FALSE, show_rownames=TRUE,annotation_col=annotation_col,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = FALSE )
  dev.off()
  
  #plot heatmap Nter with missing values using row clustering of the imputed data
  Nter.sel=which(rownames(limma.sign_norm) %in% Nter.acc.sel) #length(Nter.sel) #582
  Nter.limma.sign_norm_noimp=limma.sign_norm_noimp[Nter.sel,]
  Nter.row.order=row.order$name[row.order$name %in% rownames(Nter.limma.sign_norm_noimp)]
  #dev.off()
  pdf(paste("9.heatmap.Nter.noImputation.significant.NSAF.pdf", sep=""),width=5,height=8)  
  
  pheatmap(Nter.limma.sign_norm_noimp[Nter.row.order,col.order], 
           annotation_row = my_gene_col,
           annotation_colors = ann_colors,
           show_colnames=FALSE, show_rownames=FALSE,annotation_col=annotation_col,
           #cutree_rows = CLUSTERS,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = pal,
           border_color = FALSE,
           breaks = breaksList)
  dev.off()
  
  #plot heatmap Nter with missing values using row clustering of the imputed data
  Nter.noncan.sel=which(rownames(limma.sign_norm) %in% Nter.noncan.acc.sel) #length(Nter.noncan.sel) #3
  Nter.noncan.limma.sign_norm_noimp=limma.sign_norm_noimp[Nter.noncan.sel,]
  Nter.noncan.row.order=row.order$name[row.order$name %in% rownames(Nter.noncan.limma.sign_norm_noimp)]
  #dev.off()
  pdf(paste("9C.heatmap.Nter.noncanon.noImputation.significant.NSAF.pdf", sep=""),width=5,height=8)  
  
  pheatmap(Nter.noncan.limma.sign_norm_noimp[Nter.noncan.row.order,col.order], 
           annotation_row = my_gene_col,
           annotation_colors = ann_colors,
           show_colnames=FALSE, show_rownames=FALSE,annotation_col=annotation_col,
           #cutree_rows = CLUSTERS,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = pal,
           border_color = FALSE,
           breaks = breaksList)
  dev.off()
  
 
  
  
}


################ save image
save.image(file = "backup_07122022.RData")

#load a list of proteoforms selected for virotrap
#"C:/Users/DariaF/Documents/Annelies/4. ANALYSIS_MERGED_DATAC13/virotrap/virotrap_selection.txt"
#viro =read.csv(file="C:/Users/DariaF/Documents/Annelies/4. ANALYSIS_MERGED_DATAC13/virotrap/virotrap_selection.txt", sep = "\t", header = TRUE)
viro =read.csv(file="C:/Users/DariaF/Documents/Annelies/4. ANALYSIS_MERGED_DATAC13/virotrap/85SelectedGenesFull.txt", sep = "\t", header = TRUE)
nrow(viro) #174 proteoforms
length(unique(viro$Gene.name)) #85

#Virotrap overlap if pandey data has either:
#A) the same peptide as COFRADIC
viro.A=sapply(pep.db$matched_peptide_cor_seq, function(x) x %in%  viro$Sequence)
#B) the same accession + start combination
#C) the same accession as the virotrap data
viro.C = pep.db$first_accession %in%  viro$Accession
viro.acc.sel = unique(pep.db$first_accession[viro.A | viro.C])
sum(prot.wide.export$first_accession %in% viro.acc.sel) #100 proteoforms
length(unique(prot.wide.export$Gene.name[prot.wide.export$first_accession %in% viro.acc.sel])) #83 genes
length(prot.wide.export$first_accession[(prot.wide.export$first_accession %in% viro.acc.sel) 
                                        & prot.wide.export$peptide_not_in_UniProtCanon != "NA"]) #12 proteoforms with peptides outside UniProt canonical

 #UMAP and tSNE
if(TRUE){
  #Plot non-linear dimensionality reduction
  set.seed(42)
  U1=umap(limma.sign_norm, labels = str_split_fixed(colnames(limma.sign_norm), "_",5)[,2])
  U2=umap(limma.sign_norm, labels = str_split_fixed(colnames(limma.sign_norm), "_",5)[,1])
  U3=umap(limma.sign_norm, labels = str_split_fixed(colnames(limma.sign_norm), "_",5)[,4])
  U4=umap(limma.sign_norm, labels = str_split_fixed(colnames(limma.sign_norm), "_",5)[,3])
  U5=umap(limma.sign_norm, labels = annotation_col$Tissue_type)
  T5=tsne(limma.sign_norm, labels = annotation_col$Tissue_type , perplex = 20)
  umap.sample.label=my_gene_col$Cluster
  umap.sample.label=factor(umap.sample.label, levels = unique(my_gene_col$Cluster))
  U6=umap(t(limma.sign_norm), labels = umap.sample.label , dotsize = 1)
  T6=tsne(t(limma.sign_norm), labels = umap.sample.label , dotsize = 1, perplex = 15)
  
  pdf(file="10.UMAP_samples.pdf", width=8, height=6)
  plot(U5)
  dev.off()
  pdf(file="11.UMAP_tSNE_proteins.pdf", width=8, height=6)
  plot(U6)
  plot(T6)
  dev.off()
  #GMM
  #UMAP GMM
  BIC.proteins <- mclustBIC(U6$data)
  plot(BIC.proteins)
  summary(BIC.proteins)
  gmm.proteins = Mclust(U6$data)
  gmm.proteins.fix = Mclust(U6$data,10)
  
  pdf(file=paste0("12.GMM_clustering_proteins_UMAP.pdf"), width = 8, height = 8)
  plot(BIC.proteins)
  plot(gmm.proteins, what = "density", addEllipses=FALSE)
  plot(gmm.proteins, what=c("classification"),addEllipses=TRUE)
  mtext(paste0("GMM clustering proteins after UMAP - optimal clusters: ",gmm.proteins$G), side = 1, line=4)
  dev.off()
  
  plot(gmm.proteins.fix, what=c("classification"),addEllipses=TRUE)
  mtext(paste0("GMM clustering proteins after UMAP - fixed clusters: ",gmm.proteins.fix$G), side = 1, line=4)
  
  #tSNE GMM
  BIC.proteins.tSNE <- mclustBIC(T6$data)
  gmm.proteins.tsne = Mclust(T6$data)
  gmm.proteins.tSNE.fix = Mclust(T6$data,7)
  
  pdf(file=paste0("13.GMM_clustering_proteins_tSNE.pdf"), width = 8, height = 8)
  plot(BIC.proteins.tSNE)
  plot(gmm.proteins.tsne, what = "density", addEllipses=FALSE)
  plot(gmm.proteins.tsne, what=c("classification"),addEllipses=TRUE)
  mtext(paste0("GMM clustering proteins after tSNE - optimal clusters: ",gmm.proteins.tsne$G), side = 1, line=4)
  dev.off()
  
  plot(gmm.proteins.tSNE.fix, what=c("classification"),addEllipses=TRUE)
  mtext(paste0("GMM clustering proteins after tSNE - fixed clusters: ",gmm.proteins.tSNE.fix$G), side = 1, line=4)
  
  #summary: best is UMAP 9 clusters, export them to result table of limma
  
}

#add cluster number, cluster presumed content and NterCOFRADIC intersection infor to limma export table
cluster.enriched.tissue=c("heart", "digestive and urinary system","white blood cells",  "nervous system", "liver","reproductive system and blood","platelets","undetermined","platelets and digestive system", "undetermined") 
names(cluster.enriched.tissue) =unique(my_gene_col$Cluster)
decide.short$hclust = my_gene_col$Cluster[match(decide.short$protein, rownames(my_gene_col))]
decide.short$hclust.enriched.tissue = cluster.enriched.tissue[match(decide.short$hclust, names(cluster.enriched.tissue))]
decide.short$umap.gmm.cluster = gmm.proteins$classification[match(decide.short$protein, names(gmm.proteins$classification))]

################ save image
save.image(file = "backup_07122022.RData")

########################################## enrichment analysis on clusters 
library(clusterProfiler)
library(enrichplot)
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#http://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html
#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html
#reminder: A1=data.frame(Accession = unique(prot.db$first_accession))
# Convert ENSG to Entrez
id1<-bitr(A1$ENSG, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
id1 = id1[!duplicated(id1[c("ENSEMBL")]),]
# Convert UniProt to Entrez
id2<-bitr(A1$UniProt.stable.ID, fromType = "UNIPROT", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
id2 = id2[!duplicated(id2[c("UNIPROT")]),]
A1$Entrez=id1$ENTREZID[match(A1$ENSG, id1$ENSEMBL)]
A1$Entrez=ifelse(is.na(A1$Entrez), id2$ENTREZID[match(A1$UniProt.stable.ID, id2$UNIPROT)],A1$Entrez)
sum(is.na(A1$Entrez))/nrow(A1) #4%

#select input data for GO enrichment
go.input=merge(decide.short, A1, by.x="protein", by.y="Accession", all.x=TRUE)
go.input=go.input[!grepl( "-|_",go.input$protein), ] #only uniprot canonical proteins
sum(is.na(go.input$Entrez))/nrow(go.input) #0.7% missing Entrez
go.input=go.input[!is.na(go.input$Entrez),] #only known Entrez
go.input$hclust=factor(go.input$hclust, levels = unique(my_gene_col$Cluster))

######### GO enrichment for hclust
if (TRUE) {
  compare_GO <- compareCluster(Entrez~hclust, data=go.input, fun="enrichGO", 
                               OrgDb= "org.Hs.eg.db",
                               ont           = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.1,
                               readable      = TRUE)
  compare_GO@compareClusterResult$hclust = factor(compare_GO@compareClusterResult$hclust, levels=unique(my_gene_col$Cluster))
  
  hclust.list=split.data.frame(go.input,go.input$hclust)
  ego=list()
  for (i in 1:length(hclust.list)){
    ego[[i]] <- enrichGO(gene          = hclust.list[[i]]$Entrez,
                         universe      = unique(A1$Entrez),
                         OrgDb         = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.1,
                         readable      = TRUE)
  }
  
  pdf("14.enriched_GOterms_in_hclust.pdf",width=5,height=4)
  for (i in 1:length(ego)){
    print(dotplot(ego[[i]], showCategory=10, title = paste("Enriched GO terms\n", names(hclust.list)[i], sep=""), font.size = 6))
  }
  dev.off()
  pdf("15.compare_enriched_GOterms_in_hclust.pdf",width=8,height=10)
  dotplot(compare_GO, x = "hclust",  font.size = 6)
  dev.off()
  
}
######### GO enrichment for UMAP
if (TRUE) {
  compare_GO.umap <- compareCluster(Entrez~umap.gmm.cluster, data=go.input, fun="enrichGO", 
                                    OrgDb= "org.Hs.eg.db",
                                    ont           = "ALL",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.1,
                                    readable      = TRUE)
  
  umap.list=split.data.frame(go.input,as.factor(go.input$umap.gmm.cluster))
  ego.umap=list()
  for (i in 1:length(umap.list)){
    ego.umap[[i]] <- enrichGO(gene          = umap.list[[i]]$Entrez,
                              universe      = unique(A1$Entrez),
                              OrgDb         = org.Hs.eg.db,
                              keyType = "ENTREZID",
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.1,
                              readable      = TRUE)
  }
  
  pdf("16.enriched_GOterms_in_UMAP.pdf",width=5,height=4)
  for (i in 1:length(ego.umap)){
    print(dotplot(ego.umap[[i]], showCategory=10, title = paste("Enriched GO terms\n", names(umap.list)[i], sep=""), font.size = 6))
  }
  dev.off()
  pdf("17.compare_enriched_GOterms_in_UMAP.pdf",width=8,height=10)
  dotplot(compare_GO.umap, x = "umap.gmm.cluster",  font.size = 6)
  dev.off()
}

################ save image
save.image(file = "backup_07122022.RData")

######export tables
limma.export= prot.nsaf
colnames(limma.export)=paste0("log2NSAF_", colnames(limma.export))
select_cols=c("first_accession","Gene.name",	"Description" ,"ENSG" ,"peptides", "peptide_maps_to_UniProtCanon",	"peptide_not_in_UniProtCanon","matched_modifications","cat_db", 	"accessions",	"filter_min2peptides_1sample",	"filter_all_replicates_1tissue","sample_count")		

limma.export=merge(prot.wide.export[,select_cols],
                   limma.export, by.x="first_accession", by.y="row.names", all.x=TRUE)
limma.export=merge(decide.short,
                   limma.export,
                   by.x="protein", by.y="first_accession", all=TRUE)
limma.export$NterCOFRADIC_intersect=limma.export$protein %in% Nter.acc.sel
limma.export$virotrap_intersect=limma.export$protein %in% viro.acc.sel
limma.export$hclust=factor(limma.export$hclust, levels = unique(my_gene_col$Cluster))
limma.export %>%
  arrange(hclust ) -> limma.export

write.table(limma.export, file=paste0("limma.export.tsv"), sep="\t", row.names = FALSE, quote=TRUE)
write.table(limma.export[limma.export$NterCOFRADIC_intersect,],
            file=paste0("Nter.limma.export.tsv"), sep="\t", row.names = FALSE, quote=TRUE)
write.table(limma.export[limma.export$peptide_not_in_UniProtCanon!="NA"& !is.na(limma.export$tissuesCompared_logFC_adjPval),],
            file=paste0("noncanonical_significant.limma.export.tsv"), sep="\t", row.names = FALSE, quote=TRUE)

limma.export %>%
  arrange(Gene.name) -> limma.export
write.table(limma.export[limma.export$virotrap_intersect,],
            file=paste0("virotrap.limma.export.tsv"), sep="\t", row.names = FALSE, quote=TRUE)

#find significant canonical with the same gene name as significant noncanonical
limma.export[!is.na(limma.export$tissuesCompared_logFC_adjPval),] %>%
  filter(!is.na(Gene.name)) %>%
  group_by(Gene.name) %>%
  mutate(count.gene=n()) %>%
  filter(count.gene>1) %>%
  arrange(Gene.name, peptide_not_in_UniProtCanon !="NA")%>%
  mutate(distinct_clusters= (n_distinct(hclust)>1 & n_distinct(umap.gmm.cluster)>1),#and proteoforms in different clusters
         gene_with_noncanonical_proteoform = sum(peptide_not_in_UniProtCanon !="NA")) %>%
  filter(distinct_clusters ==TRUE, gene_with_noncanonical_proteoform > 0 ) %>%
  dplyr::select(-c("count.gene",	"distinct_clusters",	"gene_with_noncanonical_proteoform"))->significant.pairs.limma

write.table(significant.pairs.limma,
            file=paste0("significant.distinct.proteoform.pairs.limma.export.tsv"), sep="\t", row.names = FALSE, quote=TRUE)

length(unique(significant.pairs.limma$Gene.name)) #26 genes with proteoforms that have distinct expression patters accross tissues (significant limma and different cluster)

##create a simple table deciding if protein in tissue is up or down
decide %>%
  mutate(upregulated.tissue = ifelse(test_p005 == 1, str_split_fixed(contrast, "-", 2)[,1],
                              str_split_fixed(contrast, "-", 2)[,2]),
         downregulated.tissue = ifelse(test_p005 == 1, str_split_fixed(contrast, "-", 2)[,2],
                                     str_split_fixed(contrast, "-", 2)[,1]))  %>%
  select (protein, upregulated.tissue, downregulated.tissue)   %>%
  pivot_longer(!protein, names_to = "differential_expression", values_to = "tissue") %>%
  group_by(protein,differential_expression, tissue ) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  group_by(protein, tissue ) %>%
  arrange(desc(n) ) %>%
  summarize_all(first) %>%
  mutate(differential_expression=ifelse(differential_expression=="upregulated.tissue", "up", "down")) %>% select(!n)->decide.plot

##plot per gene
significant.pairs.limma %>%
  pivot_longer(cols =grep( "log2NSAF", colnames(significant.pairs.limma) ),names_to = "sample", values_to = "log2NSAF")-> significant.pairs.limma.long
significant.pairs.limma.long$tissue=  apply(str_split_fixed(significant.pairs.limma.long$sample, pattern = "_", 4)[,2:3],1, paste, collapse="_")

significant.pairs.limma.long$differential.expr=mapply(function(x,y) grepl(x,y), significant.pairs.limma.long$tissue, significant.pairs.limma.long$tissuesCompared_logFC_adjPval)
significant.pairs.limma.long=merge(significant.pairs.limma.long, decide.plot, by.x=c("protein","tissue" ), by.y=c("protein", "tissue"), all.x=TRUE)
significant.pairs.limma.long$differential_expression=factor(ifelse(is.na(significant.pairs.limma.long$differential_expression), "not signif.", significant.pairs.limma.long$differential_expression), levels=c("up", "down", "not signif." ))
significant.pairs.limma.long$gene.desc=A1$Description[match(significant.pairs.limma.long$protein, A1$Accession)]

pdf("18.differentially_expressed_proteoforms_per_gene_plots.pdf",width=8,height=6)
for (i in unique(significant.pairs.limma.long$Gene.name)) {
 plot(ggplot(data=significant.pairs.limma.long[significant.pairs.limma.long$Gene.name==i,], aes(x=tissue, y=log2NSAF))+
        geom_boxplot(aes(fill = differential_expression), alpha=0.6,  outlier.shape = NA)+
        geom_jitter(aes(color = differential_expression))+
        ggtitle(paste0(i,"\n" , unique(significant.pairs.limma.long$gene.desc[significant.pairs.limma.long$Gene.name==i])))+
        facet_grid(paste0("Hierarchical ", hclust) + paste0("UMAP Cluster ",umap.gmm.cluster) + protein~.)+
        theme_minimal()+
        scale_color_manual(values=c("orange", "darkblue","grey60"))+
        scale_fill_manual(values=c("orange", "darkblue", "grey60"))+
        theme(text=element_text(size=8), axis.text.x = element_text(angle=90, size=6, vjust = 0.4 ))
      )
  
}
dev.off()

################ save image
save.image(file = "backup_07122022.RData")

############ @@@@@@ boxplot example gene

limma.export %>%
  pivot_longer(cols =grep( "log2NSAF", colnames(limma.export) ),names_to = "sample", values_to = "log2NSAF")-> boxplot.long
boxplot.long$tissue=  apply(str_split_fixed(boxplot.long$sample, pattern = "_", 4)[,2:3],1, paste, collapse="_")
boxplot.long.copy=boxplot.long

boxplot.long$differential.expr=mapply(function(x,y) grepl(x,y), boxplot.long$tissue, boxplot.long$tissuesCompared_logFC_adjPval)
boxplot.long=merge(boxplot.long, decide.plot, by.x=c("protein","tissue" ), by.y=c("protein", "tissue"), all.x=TRUE)
boxplot.long$differential_expression=factor(ifelse(is.na(boxplot.long$differential_expression), "not signif.", boxplot.long$differential_expression), levels=c("up", "down", "not signif." ))
boxplot.long$gene.desc=A1$Description[match(boxplot.long$protein, A1$Accession)]


pdf("19.Nt-COFRADIC_proteoforms_tissue_expression_boxplots.pdf",width=8,height=6)
for (i in c("TPM3", "EPB41L3")) {
  plot(ggplot(data=boxplot.long[boxplot.long$Gene.name==i & boxplot.long$NterCOFRADIC_intersect,], aes(x=tissue, y=log2NSAF))+
         geom_boxplot(aes(fill = differential_expression), alpha=0.6, outlier.shape = NA)+
         geom_jitter(aes(color = differential_expression))+
         ggtitle(paste0(i,"\n" , unique(boxplot.long$gene.desc[boxplot.long$Gene.name==i & !is.na(boxplot.long$Gene.name)])))+
         facet_grid(paste0("Hierarchical ", hclust) + paste0("UMAP Cluster ",umap.gmm.cluster) + protein~.)+
         theme_minimal()+
         scale_color_manual(values=c("orange","darkblue", "grey60"))+
         scale_fill_manual(values=c("orange","darkblue", "grey60"))+
         theme(text=element_text(size=8), axis.text.x = element_text(angle=90, size=6, vjust = 0.4 ))
  )
  
}
dev.off()

#based on plot 19, refine selection
boxplot.selection=list()
boxplot.selection[[1]]=list(Gene.name=c("TPM3"),
                             protein = c("P06753", "P06753−2", "P06753−3"),
                             tissue = c("Adult_Colon", "Adult_Esophagus", "Adult_Pancreas", "Adult_Rectum", "Adult_Frontalcortex", "Fetal_Brain"))
boxplot.selection[[2]]=list(Gene.name=c("EPB41L3"),
                            protein = c("Q9Y2J2", "Q9Y2J2−2"),
                            tissue = c("Adult_Adrenalgland", "Adult_Frontalcortex", "Adult_Spinalcord", "Fetal_Liver"))



for (i in 1:length(boxplot.selection)) {
  pdf(paste0("19.", i, ".Nt-COFRADIC_proteoforms_tissue_expression_boxplots_selected.pdf"),
      width=1.2*(length(boxplot.selection[[i]]$protein) + length(boxplot.selection[[i]]$tissue)),height=5)
  plot(ggplot(data=boxplot.long[boxplot.long$Gene.name==boxplot.selection[[i]]$Gene.name & 
                                  boxplot.long$protein %in% boxplot.selection[[i]]$protein & 
                                  boxplot.long$tissue %in% boxplot.selection[[i]]$tissue,],
              aes(x=tissue, y=log2NSAF))+
         geom_boxplot(aes(fill = protein), alpha=0.6, outlier.shape = NA, position = position_dodge2(preserve = "single"))+
         geom_point(aes(color = protein), position=position_jitterdodge(jitter.width=0.15))+
         ggtitle(paste0(boxplot.selection[[i]]$Gene.name,"\n" , 
                        str_split(unique(boxplot.long$gene.desc[boxplot.long$Gene.name==boxplot.selection[[i]]$Gene.name & !is.na(boxplot.long$Gene.name)]), " \\[")[[1]])
                        )+
         facet_grid(.~tissue, scale="free_x")+
         theme_minimal()+
         scale_color_manual(values=c("orange","darkblue", "grey60"))+
         scale_fill_manual(values=c("orange","darkblue", "grey60"))+
         theme(strip.text.x = element_blank())
  )
  dev.off()
}

################ save image
save.image(file = "backup_07122022.RData")

# cluster.enriched.tissue "heart", "digestive and urinary system","white blood cells",  "nervous system", "liver","reproductive system and blood","platelets","undetermined","platelets and digestive system", "undetermined"
# NKcells"   ,     "Monocytes"    ,  "CD8Tcells"   ,   "Bcells", "CD4Tcells",
#                                 "Heart",
#                                 "Pancreas","Colon" ,"Gut","Gallbladder" ,   "Rectum", "Esophagus" ,
#                                 "Liver","Lung",
#                                 "Spinalcord" ,"Frontalcortex", "Retina" ,"Brain",
#                                 "Prostate" ,"Ovary"     ,     "Testis",
#                                 "Adrenalgland",   "Platelets",
#                                 "Urinarybladder", "Kidney",
#                                 "Placenta")
################### 08032023 tissue meeting make nicer plots for candidate genes
unique(boxplot.long$tissue)
#based on plot 18, refine selection
boxplot.selection=list()
boxplot.selection[[1]]=list(Gene.name=c("HNRNPC"),
                            protein = c("P07910", "P07910−2"),
                            tissue = c("Adult_Monocytes","Adult_Bcells","Adult_CD4Tcells","Adult_CD8Tcells","Adult_NKcells"))
boxplot.selection[[2]]=list(Gene.name=c("CDC42"),
                            protein = c("P60953", "P60953−1"),
                            tissue = c("Adult_CD8Tcells", "Adult_Frontalcortex", "Adult_Spinalcord", "Adult_Liver","Fetal_Liver"))
boxplot.selection[[3]]=list(Gene.name=c("ADARB1"),
                            protein = c("P78563", "P78563−4"),
                            tissue = c("Adult_Retina", "Adult_CD4Tcells","Adult_CD8Tcells", "Adult_NKcells","Adult_Liver"))
boxplot.selection[[4]]=list(Gene.name=c("SLC25A3"),
                            protein = c("Q00325", "Q00325−2"),
                            tissue = c("Adult_Heart","Adult_Kidney","Adult_Liver"))
boxplot.selection[[5]]=list(Gene.name=c("SPTAN1"),
                            protein = c("Q13813", "Q13813-2", "Q13813-3","ENST00000630866_9_128566741_aTIS_111db1"),
                            tissue = c("Adult_Heart","Fetal_Heart","Adult_Frontalcortex","Adult_Spinalcord","Fetal_Brain","Adult_NKcells", "Adult_Platelets"))



for (i in 1:length(boxplot.selection)) {
  pdf(paste0("20.", i, ".select_proteoforms_tissue_expression_boxplots.pdf"),
      width=1.2*(length(boxplot.selection[[i]]$protein) + length(boxplot.selection[[i]]$tissue)),height=5)
  temp.df=boxplot.long[boxplot.long$Gene.name==boxplot.selection[[i]]$Gene.name & 
                         boxplot.long$protein %in% boxplot.selection[[i]]$protein & 
                         boxplot.long$tissue %in% boxplot.selection[[i]]$tissue,]
  temp.df$tissue=factor(temp.df$tissue, levels=boxplot.selection[[i]]$tissue)
  plot(ggplot(data=temp.df,
              aes(x=tissue, y=log2NSAF))+
         geom_boxplot(aes(fill = protein), alpha=0.6, outlier.shape = NA, position = position_dodge2(preserve = "single"))+
         geom_point(aes(color = protein), position=position_jitterdodge(jitter.width=0.15))+
         ggtitle(paste0(boxplot.selection[[i]]$Gene.name,"\n" , 
                        str_split(unique(boxplot.long$gene.desc[boxplot.long$Gene.name==boxplot.selection[[i]]$Gene.name & !is.na(boxplot.long$Gene.name)]), " \\[")[[1]])
         )+
         facet_grid(.~tissue, scale="free_x")+
         theme_minimal()+
         scale_color_manual(values=c("orange","darkblue", "grey60", "brown1"))+
         scale_fill_manual(values=c("orange","darkblue", "grey60", "brown1"))+
         theme(strip.text.x = element_blank())
  )
  dev.off()
}


################ save image
save.image(file = "backup_08032023.RData")


########### PhD thesis responses

#make figure with expression ratio non-canon vs. canon across tissues

Nter.noncan.acc.sel
#plot N-term COFRADIC overlap with pandey not mapping to UniProtCanon #NSAF, no limma, no imputation
#Nter.nsaf=prot.nsaf[Nter.noncan.acc.sel,]
Nter.noncan.gene.sel=vector()
Nter.canon.nsaf=list()
for (i in 1:length(Nter.noncan.acc.sel)){
  Nter.noncan.gene.sel[i]=unique(prot.wide.export$Gene.name[prot.wide.export$first_accession == Nter.noncan.acc.sel[i]])
  if(Nter.noncan.gene.sel[i] == "ACTBP11") {Nter.noncan.gene.sel[i] <-"ACTB" }
  Nter.canon.nsaf[[i]] = prot.wide.export[prot.wide.export$Gene.name == Nter.noncan.gene.sel[i] & !is.na(prot.wide.export$Gene.name) & prot.wide.export$cat_db==1 ,]
  if (nrow(Nter.canon.nsaf[[i]])==0) {
    Nter.canon.nsaf[[i]] =as.data.frame(matrix(nrow=1, ncol=ncol(prot.wide.export)))
    colnames(Nter.canon.nsaf[[i]]) <- colnames(prot.wide.export)
  }
}
Nter.canon.nsaf=do.call(rbind, Nter.canon.nsaf)
Nter.canon.nsaf.rows=make.unique( Nter.canon.nsaf$first_accession )
Nter.canon.nsaf.rows[is.na(Nter.canon.nsaf.rows)] <- "none"
rownames(Nter.canon.nsaf)= Nter.canon.nsaf.rows
Nter.noncan.gene.sel == Nter.canon.nsaf$Gene.name
Nter.canon.nsaf=as.matrix(Nter.canon.nsaf[,grepl( "NSAF", colnames(Nter.canon.nsaf))])
Nter.canon.nsaf=log2(Nter.canon.nsaf)
colnames(Nter.canon.nsaf)=gsub("NSAF_","",colnames(Nter.canon.nsaf))
str(Nter.nsaf)
str(Nter.canon.nsaf)
Nter.compare=Nter.nsaf-Nter.canon.nsaf
Nter.noncan.unique=!is.na(Nter.nsaf) & is.na(Nter.canon.nsaf)
colSums(Nter.noncan.unique)
for (i in 1: ncol(Nter.compare)){
  Nter.compare[,i][Nter.noncan.unique[,i]]<- 6
}
rownames(Nter.compare)<-paste0(Nter.noncan.gene.sel," ", rownames(Nter.nsaf))

#cluster as binary data (to avoid imputation)
#Nter.noncan.order=cluster::mona(!is.na(Nter.nsaf[,colSums(!is.na(Nter.nsaf))>0]))

pdf(paste("21.heatmap.Nter.compare.noncan.vs.canon.log2FC.pdf", sep=""),width=11,height=5.5)  

pheatmap(Nter.compare[Nter.noncan.order$order, rownames(annotation_col)[order(annotation_col$Tissue_type)]], 
         annotation_colors = ann_colors,
         show_colnames=FALSE, show_rownames=TRUE,annotation_col=annotation_col,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = FALSE,
         main = "log2FC non-canonical vs. canonical proteoform")
dev.off()

pdf(paste("22.heatmap.Nter.only.canon.log2NSAF.pdf", sep=""),width=11,height=5.5)  

pheatmap(Nter.canon.nsaf[Nter.noncan.order$order, rownames(annotation_col)[order(annotation_col$Tissue_type)]], 
         annotation_colors = ann_colors,
         show_colnames=FALSE, show_rownames=TRUE,annotation_col=annotation_col,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = FALSE,
         main = "log2NSAF canonical counteraprt of 24 non-canon.")
dev.off()


################ save image
save.image(file = "backup_12042023.RData")