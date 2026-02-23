library(data.table)
library(dplyr)
library(ggplot2)
library(UpSetR)
library(grid)

gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))

cgc <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/Cosmic_CancerGeneCensus_v102_GRCh37.tsv")
cgc_tsg <- cgc[grepl("TSG", ROLE_IN_CANCER), ]
cgc_tsg <- inner_join(gene_annot_dt, cgc_tsg, by=c("gene_name"="GENE_SYMBOL"))
# loosing 3 genes

cgc_oncogene <- cgc[grepl("oncogene", ROLE_IN_CANCER), ]
cgc_oncogene <- inner_join(gene_annot_dt, cgc_oncogene, by=c("gene_name"="GENE_SYMBOL"))


pr_output_name <- "zScore_gt5"
pr_res <- fread(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_", pr_output_name, "/protrider_summary.csv"))

protein_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protrider_sa.tsv")

pr_res <- left_join(pr_res, protein_sa, by=c("sampleID"="full_name"))

pr_res_aberrant <- pr_res[PROTEIN_outlier == T]





or_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_outliers.tsv")

fr_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_outliers.tsv")


rna_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")
  



stranded_rna_samples <- rna_sa[Diag != "Unstranded_data", pid]

pr_samples <- protein_sa$pid


both_pr_rna_samples <- intersect(stranded_rna_samples, pr_samples)


fr_res_common <- fr_res[sampleID %in% both_pr_rna_samples & abs(deltaPsi) >= 0.3]
fr_res_common[, gene_sample := paste0(sampleID, "_", hgncSymbol)]

or_res_common <- or_res[sampleID %in% both_pr_rna_samples]
or_res_common[, gene_sample := paste0(sampleID, "_", hgncSymbol)]

pr_res_common <- pr_res_aberrant[sampleID %in% both_pr_rna_samples]
pr_res_common[, gene_sample := paste0(pid, "_", proteinID)]



listInput <- list(`protein(|zScore|>5)` = pr_res_common$gene_sample, expression = or_res_common$gene_sample, 
                  splicing = fr_res_common$gene_sample)

p <- upset(fromList(listInput), order.by = "freq", text.scale = 1.5)
p
p <- p + grid.text("All outliers",x = 0.65, y=0.95, gp=gpar(fontsize=20))



listInput <- list(`protein(|zScore|>5)` = pr_res_common[PROTEIN_ZSCORE < 0]$gene_sample, expression = or_res_common[zScore < 0]$gene_sample, 
                  splicing = fr_res_common$gene_sample)

p <- upset(fromList(listInput), order.by = "freq", text.scale = 1.5)
p
p + grid.text("Underexpression outliers",x = 0.65, y=0.95, gp=gpar(fontsize=20))



listInput <- list(`protein(|zScore|>5)` = pr_res_common[PROTEIN_ZSCORE > 0]$gene_sample, expression = or_res_common[zScore > 0]$gene_sample, 
                  splicing = fr_res_common$gene_sample)

p <- upset(fromList(listInput), order.by = "freq", text.scale = 1.5)
p
p + grid.text("Overexpression outliers",x = 0.65, y=0.95, gp=gpar(fontsize=20))


