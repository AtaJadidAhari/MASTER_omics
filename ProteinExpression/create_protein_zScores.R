library(data.table)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(scales)

zScore_cutoff <- 5

res_dir <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/results/"
results <- list.files(res_dir, full.names = TRUE, recursive = FALSE)
results <- c("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/results//protrider_all_gaussian(gs).csv")

# Between sample normalisation
size_factor_normlize <- function(data) {

  sf= DESeq2::estimateSizeFactorsForMatrix(data)
  data.norm <- as.data.frame(t(t(data)/sf))
  return(data.norm)
}


robust_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)  # constant=1 = unscaled MAD
  if (mad_val == 0) return(rep(0.0005, length(x)))   # avoid division by zero
  return(0.6745 *(x - med) / mad_val)
}

calculate_zScores <- function(intensities){
  drop_group_name <- "protein_all"
  #fwrite(zScore_res_sorted, paste0(drop_group, "/zScores.tsv"))
  rownames(intensities) <- intensities$protein_ID
  raw_intensities <- intensities[, -c("protein_ID")]
  raw_intensities <- size_factor_normlize(raw_intensities)
  log2_intensities <- log2(raw_intensities)
  log2_intensities_zcores <- as.data.table(t(apply(log2_intensities, 1, scale)))
  
  log2_intensities_robust_zcores <- as.data.table(t(apply(log2_intensities, 1, robust_z))) 
  setnames(log2_intensities_robust_zcores, colnames(log2_intensities))
  
  rownames(log2_intensities_zcores) <- rownames(intensities)
  setnames(log2_intensities_zcores, colnames(log2_intensities))
  

  log2_intensities_zcores[, proteinID := rownames(intensities)]
  zScore_res <- melt(log2_intensities_zcores, value.name = "zScore", variable.name = "sampleID")
  
  # log2_intensities_robust_zcores <- as.data.table(log2_intensities_robust_zcores)
  # robust_res <- melt(log2_intensities_robust_zcores, id.vars = "geneID",
  #                    value.name = "robust_zScore", variable.name = "sampleID")
  # 
  # # Merge both by geneID and sampleID
  # zScore_res <- merge(zScore_res, robust_res, by = c("geneID", "sampleID"))
  # 
  zScore_res[, abs_zScore := abs(zScore)]
  zScore_res_sorted <- zScore_res[order(zScore_res$abs_zScore, decreasing = TRUE), ]
  
  zScore_res_sorted[, DROP_GROUP:= drop_group_name]
  
  p_values <- 2 * pnorm(-abs(zScore_res_sorted$abs_zScore))
  zScore_res_sorted$pvalue <- p_values
  
  p_adj <- p.adjust(p_values, method = "BY")
  zScore_res_sorted$PROTEIN_PADJ <- p_adj
  
  zScore_res_sorted[padjust < 0.05, PROTEIN_outlier := TRUE]
  

}


intensity_subset <- fread( "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protrider_intensities.tsv")

protein_na_threshold <- 0.3
total_cols <- ncol(intensity_subset) - 1

na_sum_dt <- data.table(na_count =rowSums(is.na(intensity_subset)))
na_sum_dt[, na_0 := na_count <= 0]
na_sum_dt[, na_30 := na_count <= 0.3 * total_cols]
na_sum_dt[, na_50 := na_count <= 0.5 * total_cols]





intensity_subset_nonna <- intensity_subset[which(na_sum_dt$na_30 == TRUE)]
nrow(intensity_subset_nonna)


zScore_res_sorted <- calculate_zScores(intensity_subset_nonna)


gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))


zScore_res_sorted <- zScore_res_sorted %>%
  mutate(full_name = sampleID) %>% 
  mutate(sampleID = sapply(str_split(sampleID, "-"), `[`, 2))



zScore_res_sorted <- left_join(zScore_res_sorted, gene_annot_dt, by=c("proteinID"= "gene_name"))
setnames(zScore_res_sorted, "gene_id", "geneID")

