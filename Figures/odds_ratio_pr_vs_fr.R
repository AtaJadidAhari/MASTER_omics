library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(arrow)

source("~/Scripts/ProteinExpression/load_pr_data.R")


prepare_rna_protein_dt <- function(pr_output_name, or_res, method_name){
  pr_res <- load_pr_data(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_", pr_output_name, "/protrider_summary.csv"))
 
  
  protrider_proteins <- unique(pr_res$geneID_short)
  protrider_samples <- unique(pr_res$sampleID)

  or_res <- or_res[sampleID %in% protrider_samples]
  or_res[, geneID:= sub("\\..*", "", geneID)]
  
  or_genes <- unique(or_res[, geneID])
  
  all_genes <- intersect(or_genes, protrider_proteins)
  
  merged_df <- full_join(or_res[geneID %in% all_genes], pr_res[geneID_short %in% all_genes], by=c("sampleID", "geneID"="geneID_short"))
  
  setnames(merged_df, c("aberrant.x", "aberrant.y"), c("RNA_aberrant", "protein_aberrant"))
  if (!"PROTEIN_ZSCORE" %in% names(merged_df)) {
    setnames(merged_df, c("zScore.y"), c("PROTEIN_ZSCORE"))
  }
  
  
  merged_df[, annotation := "Non outlier"]
  merged_df[RNA_aberrant == TRUE & protein_aberrant == TRUE, annotation := "RNA + protein outlier"]
  merged_df[(is.na(RNA_aberrant) | RNA_aberrant == FALSE) & protein_aberrant == TRUE, annotation := "Protein outlier"]
  merged_df[RNA_aberrant == TRUE & (is.na(protein_aberrant) | protein_aberrant == FALSE), annotation := "RNA outlier"]
  
  merged_df <- merged_df %>%
    mutate(
      protein_direction = case_when(
        PROTEIN_ZSCORE > 0 ~ "Protein overexpression",
        PROTEIN_ZSCORE < 0 ~ "Protein underexpression",
        TRUE ~ NA_character_  # Exclude 0s (neutral)
      )
    ) %>%
    filter(!is.na(protein_direction))  # Remove neutral rows
  
  merged_df[, method := method_name]
  return (merged_df)
}

# Function to compute fisher test stats per group
get_fisher_stats <- function(data) {
  print(head(data))
  cont_table <- table(data$RNA_aberrant, data$protein_aberrant)
  print(cont_table)
  fisher <- fisher.test(cont_table)
  print(fisher$p.value)
  data.frame(
    odds_ratio = fisher$estimate,
    ci_low = fisher$conf.int[1],
    ci_high = fisher$conf.int[2],
    p_value = ifelse(fisher$p.value == 0, "P < 1e-308", sprintf("P = %.2e", fisher$p.value))
  )
}

# or_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv")

or_res <- data.table(read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust.parquet", package="arrow"))
print(head(or_res))


pr_or_dt <- list()
pr_output_name <- "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"
pr_or_dt[[pr_output_name]] <- prepare_rna_protein_dt(pr_output_name, or_res, "PROTRIDER")

#pr_output_name <- "t_gs_injMean3"
#pr_or_dt[[pr_output_name]] <- prepare_rna_protein_dt(pr_output_name, or_res)

pr_output_name <- "zScore_gt5"
pr_or_dt[[pr_output_name]] <- prepare_rna_protein_dt(pr_output_name, or_res, "Z-Scores")


pr_or_dt <- rbindlist(pr_or_dt, fill=TRUE)

# Run stats
results <- pr_or_dt %>%
  group_by(method, protein_direction) %>%
  summarise(get_fisher_stats(cur_data()), .groups = "drop") 


x_axis_breaks <- c(1, 10, 100, 1000)


p <- ggplot(results, aes(x = odds_ratio, y = method, color = method)) +
  # geom_point(size = 3) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high),
                  position = position_dodge(.9)) +
  
  geom_text(aes(y=method, x=ci_high-10, label=p_value),
            hjust=0.1, vjust=-0.1, size=9/.pt)+
  
  scale_x_log10(breaks = x_axis_breaks, ) +
  facet_wrap(~protein_direction, ncol = 1) +
  theme_bw(base_size=12) +
  labs(x = "Odds ratio of shared RNA–protein outliers", y = NULL) + 
  annotation_logticks(sides="b") +
  geom_vline(xintercept = 1, linetype='dashed', color='grey') + 
  coord_cartesian(xlim = c(1, 1000)) +
  theme(
    legend.position="none",
    legend.title = element_blank(),
  )

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/odds_ratio_RNA_pr_outliers.png") ,
    width = 8, height = 6, units = "in", res = 600)
p
dev.off()

