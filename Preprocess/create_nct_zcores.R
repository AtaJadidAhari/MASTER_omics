library(data.table)
library(dplyr)
library(stringr)

or_res_all <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv")
nct_z <- fread("/omics/odcf/analysis/hipo/hipo_021/cohort_analysis/ReeA_exp_profile/250526215122_RNAsampleInfo_TPM.cohort_info.Q4.thresholds-1_2.all.expression_profile.v5.tsv")
head(nct_z)

or_res_genes <- unique(or_res_all$geneID)

nct_z <- nct_z[gene_id %in% or_res_genes]

drop_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")


tumor_keys <- c("T", paste0("T", 1:15))
tumor_values <- c("tumor", sprintf("tumor%02d", 1:15))

# Create metastasis entries
met_keys <- c("M", paste0("M", 1:15))
met_values <- c("metastasis", sprintf("metastasis%02d", 1:15))


# Combine into a named vector
lookup_map <- c(setNames(tumor_keys, tumor_values),
                setNames(met_keys, met_values))

lookup_map["T1"] <- "tumor"
lookup_map["M1"] <- "metastasis"



nct_z <- nct_z %>%
  mutate(pid = sapply(str_split(PID, "\\."), `[`, 2)) %>%
  mutate(sample_type = sapply(str_split(PID, "\\."), `[`, 3))

nct_res <- left_join(nct_z, drop_sa[, c("pid", "sample_type", "Diag", "Oncotree Code")], by=c("pid", "sample_type"))


 
nct_res_non_na <- nct_res[!is.na(Diag)]

nct_res_non_na[, zScore := FINAL_ZSCORE]
nct_res_non_na[,aberrant := FALSE]
nct_res_non_na[abs(zScore) >= 5, aberrant := TRUE]
table(nct_res_non_na$aberrant)

nct_res_non_na[, abs_zScore := abs(zScore)]

nct_res_non_na <- nct_res_non_na[order(-abs_zScore)]
setnames(nct_res_non_na, "gene_id", "geneID")
nct_res_non_na[, sampleID := pid]
nct_res_non_na <- nct_res_non_na[!is.na(FINAL_ZSCORE)]

fwrite(nct_res_non_na[!is.na(FINAL_ZSCORE)], "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/nct_zScores_or_gene_subset.tsv", sep="\t")


hist(nct_res_non_na, breaks=100)

nct_res_subset <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/nct_zScores_or_gene_subset.tsv")




