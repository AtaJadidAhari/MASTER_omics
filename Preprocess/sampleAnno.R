library(dplyr)
library(stringr)

analyzed_bams <- basename(scan("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/strandedness/bam_paths.txt", character()))
analyzed_pids <- unlist(strsplit(analyzed_bams, "_"))[3 * 1:length(analyzed_bams) - 1] %>%
  gsub(".*-","",.)

unique_pids <- unique(analyzed_pids)


sample_annot_RNA <- read.csv("//omics/odcf/analysis/hipo/hipo_021/outlier_analysis//sample_data/Basic_Clinical_Info_DKTK-P115_ClinicalInformation.csv",header = T, sep = ";") %>%
  filter(RNA.data=="yes")


# proteomics sa


proteomics_new_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/proteomics_TUM/MASTER_preprocessed_fp_with_annot_scores_298.csv")
proteomics_new_sa[, first_fasta_header := tstrsplit(`Fasta headers`, ";", keep = 1)]
proteomics_new_sa[, first_protein_name := tstrsplit(`Protein IDs`, ";", keep = 1)]


proteomics_new_sa[, c("protein_database", "uniprot_id") := tstrsplit(`first_fasta_header`, "\\|", keep = 1:2)]
proteomics_new_sa[, gene_name := sub(".*GN=([^ ]+).*", "\\1", `first_fasta_header`)]


table(proteomics_new_sa$first_protein_name == proteomics_new_sa$gene_name) # first `PROTEIN IDs` is the same as first part of fasta header, and then chekcing the gene nanme

proteomics_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/proteomics_TUM/MASTER_subset_pancancer_Batch298.csv")


proteomics_new_sa$gene_name == proteomics_sa$`Gene names`
proteomics_new_sa$gene_name == proteomics_new_sa$`Protein IDs`


aditional_info <- proteomics_new_sa[, .SD, .SDcols=c("Fasta headers", "Protein IDs", "Majority protein IDs", "Score", "Peptide counts (unique)", "Q-value")]
fwrite(aditional_info, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/additional_protein_data.tsv", sep="\t")
aditional_info