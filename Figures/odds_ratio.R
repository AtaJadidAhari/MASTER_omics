library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

source('~/Scripts/Preprocess/VEP_SORanking.R')
source("~/Scripts/ProteinExpression/load_pr_data.R")




snv_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique.tsv")
indel_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_indel_all_aggregated_unique.tsv")

vep_res_combinded <- rbind(snv_vep_res, indel_vep_res, fill=TRUE)


vep_res_combinded[Consequence_most_severe %in% stop_variants, Consequence_most_severe := "stop (VEP)"]
vep_res_combinded[Consequence_most_severe %in% splice_region_variants, Consequence_most_severe := "splice-region (VEP)"]
vep_res_combinded[Consequence_most_severe %in% splice_site_variants, Consequence_most_severe := "splice-site (VEP)"]
vep_res_combinded[Consequence_most_severe %in% inframe_variants, Consequence_most_severe := "inframe (VEP)"]
vep_res_combinded[Consequence_most_severe %in% up_down_stream_vars, Consequence_most_severe := "up/downstream gene variant (VEP)"]
vep_res_combinded[Consequence_most_severe == "missense_variant", Consequence_most_severe := "VEP_missense_wo_am_pred"]



merge_pr_var_res <- function(pr_output_name){
  pr_res <- load_pr_data(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_", pr_output_name, "/protrider_summary.csv"))
  
  
  
  merged_df_vars <- left_join(pr_res, vep_res_combinded[, c("Gene", "sampleID", "Consequence", "IMPACT", "Consequence_most_severe")], by=c("geneID_short" = "Gene", "sampleID"))
  
  merged_df_vars[,`Rare variant` := FALSE]  
  merged_df_vars[!is.na(Consequence), `Rare variant` := TRUE]  
  
  
  merged_df_vars[, method := pr_output_name]
  
  merged_df_vars[ , pr_expresssion := "Non-outlier"]
  merged_df_vars[zScore < 0 & aberrant == TRUE, pr_expresssion := "Underexpression"]
  merged_df_vars[zScore > 0 & aberrant == TRUE, pr_expresssion := "Overexpression"]
  
  # merged_df_vars <- merged_df_vars[pr_expresssion != "Overexpression"]
  
  merged_df_vars[, is_underexpressed := pr_expresssion == "Underexpression"]
  
  merged_df_vars[is.na(Consequence_most_severe), Consequence_most_severe := ""]
  return (merged_df_vars)
  
}

pr_dt <- list()
pr_output_name <- "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"
pr_dt[[pr_output_name]] <- merge_pr_var_res(pr_output_name)

#pr_output_name <- "t_gs_injMean3"
#pr_dt[[pr_output_name]] <- load_pr_data(pr_output_name)

pr_output_name <- "zScore_gt5"
pr_dt[[pr_output_name]] <- merge_pr_var_res(pr_output_name)


pr_dt <- rbindlist(pr_dt, fill=TRUE)

#pr_dt[method=="gaussian_gs_injMean3", method := "PROTRIDER(gaussian)"]
#pr_dt[method=="t_gs_injMean3", method := "PROTRIDER(t)"]

pr_dt[method=="cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA", method := "PROTRIDER"]

pr_dt[method=="zScore_gt5", method := "zScores\n(|zScore| > 5)"]

# remove na vraiants
pr_dt <- pr_dt[!is.na(Consequence_most_severe)]

no_show_consequence <- c("non_coding_transcript_exon_variant", "", "stop_retained_variant",
                         "non_coding_transcript_exon_variant",
                         "incomplete_terminal_codon_variant", "protein_altering_variant")


# Get unique combinations of consequence and method
combo_list <- unique(pr_dt[!Consequence_most_severe %in% no_show_consequence, .(Consequence_most_severe, method)])

# Initialize results list
results_list <- list()

# Loop over each consequence × method pair
for (i in seq_len(nrow(combo_list))) {
  
  cons <- combo_list$Consequence_most_severe[i]
  meth <- combo_list$method[i]
  
  # Group 1: Variants with this consequence + method
  group_with <- pr_dt[pr_expresssion != "Overexpression" & Consequence_most_severe == cons & method == meth]
  a <- sum(group_with$is_underexpressed)                  # Underexpression
  b <- nrow(group_with) - a                               # Not underexpression
  
  # Group 2: All other variants (everything else)
  group_without <- pr_dt[!(Consequence_most_severe == cons & method == meth)]
  c <- sum(group_without$is_underexpressed)               # Underexpression
  d <- nrow(group_without[pr_expresssion != "Overexpression"]) - c  # Not underexpression
  

  # Only run Fisher’s test if both groups have data
  if ((a + b > 0) && (c + d > 0)) {
    table_matrix <- matrix(c(a, b, c, d), nrow = 2)
    test_result <- fisher.test(table_matrix)
    
    results_list[[i]] <- data.table(
      consequence = cons,
      method = meth,
      odds_ratio = as.numeric(test_result$estimate),
      lower_CI = test_result$conf.int[1],
      upper_CI = test_result$conf.int[2]
    )
  }
}

results <- rbindlist(results_list, fill=TRUE)


results[, consequence := factor(consequence)]
p <- ggplot(results, aes(x = odds_ratio, y = consequence, color = method)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  scale_x_log10() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    x = "Odds ratio (Underexpression vs. Non-outlier)",
    y = NULL,
    color = "Method"
  ) +
  theme_minimal(base_size = 12)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/outlier_odds_ratio_pr_rare_variant_protrider_new.png") ,
    width = 7, height = 6, units = "in", res = 600)
p
dev.off()



# Initialize results list
up_results_list <- list()

# Loop over each consequence × method pair
for (i in seq_len(nrow(combo_list))) {
  
  cons <- combo_list$Consequence_most_severe[i]
  meth <- combo_list$method[i]
  
  # Group 1: Variants with this consequence + method
  group_with <- pr_dt[pr_expresssion != "Underexpression" & Consequence_most_severe == cons & method == meth]
  a <- nrow(group_with[pr_expresssion == "Overexpression"])                  # Underexpression
  b <- nrow(group_with) - a                               # Not underexpression
  
  # Group 2: All other variants (everything else)
  group_without <- pr_dt[!(Consequence_most_severe == cons & method == meth)]
  c <-  nrow(group_without[pr_expresssion == "Overexpression"])                  # Underexpression
  d <- nrow(group_without[pr_expresssion != "Underexpression"]) - c              # Not Overexpression
  
  
  # Only run Fisher’s test if both groups have data
  if ((a + b > 0) && (c + d > 0)) {
    table_matrix <- matrix(c(a, b, c, d), nrow = 2)
    test_result <- fisher.test(table_matrix)
    
    up_results_list[[i]] <- data.table(
      consequence = cons,
      method = meth,
      odds_ratio = as.numeric(test_result$estimate),
      lower_CI = test_result$conf.int[1],
      upper_CI = test_result$conf.int[2]
    )
  }
}

up_results <- rbindlist(up_results_list, fill=TRUE)

print(head(up_results))

up_results[, consequence := factor(consequence)]
p <- ggplot(up_results, aes(x = odds_ratio, y = consequence, color = method)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  scale_x_log10() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    x = "Odds ratio (Overexpression vs. Non-outlier)",
    y = NULL,
    color = "Method"
  ) +
  theme_minimal(base_size = 12)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/outlier_odds_ratio_up_pr_rare_variant_protrider_new.png") ,
    width = 7, height = 6, units = "in", res = 600)
p
dev.off()


