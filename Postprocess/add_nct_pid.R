library(OUTRIDER)
library(data.table)
library(reshape2)

drop_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")


nct_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/nct_zScores.tsv")


#### export OR normalized counts ####

or_res_all <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv")
or_res_norm_counts <- or_res_all[, c("geneID", "sampleID", "normcounts")]
setnames(or_res_norm_counts, "sampleID", "pid")

or_res_norm_counts <- left_join(or_res_norm_counts, drop_sa[, c("pid", "RNA_BAM_FILE", "sample_type")], )
or_res_norm_counts[,PID := paste0(pid, ".", sample_type)]

normalized_counts <- dcast(or_res_norm_counts, geneID ~ PID, value.var = "normcounts")

fwrite(normalized_counts, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/aberrant_expression/v19/outrider/ae_normalized_counts_all.tsv", sep="\t")



or_res_dir <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider"
or_res_files <- list.dirs(or_res_dir, full.names = TRUE, recursive = FALSE)


#### create splicing outlier results ####

fr_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_results.tsv")

fr_res <- left_join(fr_res, drop_sa[, c("pid", "RNA_BAM_FILE", "sample_type")], by=c("sampleID" = "pid"))
fr_res[,PID := paste0(sampleID, ".", sample_type)]

gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))


fr_res <- left_join(fr_res, gene_annot_dt[, c("gene_id", "gene_name")], by=c("hgncSymbol" = "gene_name"))
setnames(fr_res, "gene_id", "geneID")
fr_res <- fr_res[deltaPsi >= 0.3]
fr_res <- fr_res[order(pValue),]

fwrite(fr_res, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_outliers_dpsi_gt_0_3.tsv", sep="\t")

