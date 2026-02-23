library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

source('~/Scripts/Preprocess/VEP_SORanking.R')



snv_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_HIGH_aggregated/vep_res_rare_snv_all_aggregated_unique.tsv")
indel_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_HIGH_aggregated/vep_res_rare_indel_all_aggregated_unique.tsv")


vep_res_combinded <- rbind(snv_vep_res, indel_vep_res)

vep_res_combinded[Consequence_most_severe %in% stop_variants, Consequence_most_severe := "stop (VEP)"]
vep_res_combinded[Consequence_most_severe %in% splice_region_variants, Consequence_most_severe := "splice-region (VEP)"]
vep_res_combinded[Consequence_most_severe %in% splice_site_variants, Consequence_most_severe := "splice-site (VEP)"]
vep_res_combinded[Consequence_most_severe %in% inframe_variants, Consequence_most_severe := "inframe (VEP)"]
vep_res_combinded[Consequence_most_severe %in% up_down_stream_vars, Consequence_most_severe := "up/downstream gene variant (VEP)"]
vep_res_combinded[Consequence_most_severe == "missense_variant", Consequence_most_severe := "VEP_missense_wo_am_pred"]



or_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv")
or_res[, geneID:= sub("\\..*", "", geneID)]

merged_df_vars <- left_join(or_res, vep_res_combinded[, c("Gene", "sampleID", "Consequence", "IMPACT", "Consequence_most_severe")], by=c("geneID" = "Gene", "sampleID"))


merged_df_vars[ , expression_direction := "Non-outlier"]
merged_df_vars[zScore >  0 & aberrant == TRUE, expression_direction := "Overexpression"]
merged_df_vars[zScore < 0 & aberrant == TRUE, expression_direction := "Underexpression"]

merged_df_vars[, method := "OUTRIDER"]


merged_df_vars[, total_counts := .N, by = .(method, expression_direction)]



plot_dt <- merged_df_vars[, .N, by = .(method, expression_direction, Consequence_most_severe)][
  , prop := N / sum(N), by = .(method, expression_direction)]

bar_totals <- plot_dt[, .(total_N = sum(N)), by = .(method, expression_direction)]


plot_dt <- plot_dt[!is.na(plot_dt$Consequence_most_severe), ]


no_show_consequence <- c("non_coding_transcript_exon_variant", 
                         "incomplete_terminal_codon_variant", "protein_altering_variant")


variant_colors <- c(
  "stop (VEP)" = "darkred", # "aquamarine3",  
  "frameshift_variant" = "salmon2", # "lightskyblue3",  
  "start_lost" = "darkgreen",
  "splice-site (VEP)" = "darkolivegreen2", 
  "splice_region (VEP)" = "darkolivegreen4", 
  "splice-region (VEP)" = "darkolivegreen4", 
  "intron_variant" = "aquamarine3",
  "stop_retained_variant" = "yellow3",
  # "AbExp" = "yellow1",
  # "VEP_missense" = "salmon1",   
  "pathogenic" = "darkblue", #"darkred",   
  "ambiguous" = "skyblue2", #"salmon3",   
  "benign" = "lightblue1", #"orange2", 
  "VEP_missense_wo_am_pred" = "orange1",
  "5_prime_UTR_variant" = "chocolate4",
  "3_prime_UTR_variant" = "chocolate3",
  # "non_coding (VEP)" = "moccasin",  
  "synonymous_variant" = "darkgrey",
  "no rare variant" = "white",
  "up/downstream gene variant (VEP)"= "purple",
  "inframe (VEP)"= "pink"
)

bar_totals[, method_label := paste0(method, "\nn=", total_N)]

plot_dt <- left_join(plot_dt, bar_totals, by=c("method", "expression_direction"))

plot_dt[, method_label := as.factor(method_label)]

p_n <- ggplot(plot_dt[!Consequence_most_severe %in% no_show_consequence], aes(x=prop, y=method, fill = Consequence_most_severe)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = bar_totals,
            aes(y = method, x = 0.55, label = paste0("n = ", total_N)),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ expression_direction, ncol = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  # scale_y_discrete(labels = c("t_gs_injMean3" = "t_gs_injMean3\nn=8548746", 
  #                             "gaussian_gs_injMean3\nn=8541295",
  #                             "t_gs_injMean3\nn=171",
  #                             "gaussian_gs_injMean3\nn=2802",
  #                             "t_gs_injMean3" ="t_gs_injMean3\nn=690",
  #                             "gaussian_gs_injMean3\nn=5510")) +
  labs(y="", x="Proportion of outliers with\nrare variants") + 
  scale_fill_manual(values = variant_colors) + 
  theme_bw(base_size=10)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/outlier_proportion_rare_variant_outrider.png") ,
    width = 10, height = 6, units = "in", res = 600)
p_n
dev.off()




plot_dt <- merged_df_vars[, .N, by = .(method, expression_direction, Consequence_most_severe)][
  , prop := N / sum(N), by = .(method, expression_direction)]

bar_totals <- plot_dt[, .(total_N = sum(N)), by = .(method, expression_direction)]


plot_dt <- plot_dt[!is.na(plot_dt$Consequence_most_severe), ]




no_show_consequence <- c("non_coding_transcript_exon_variant", 
                         "incomplete_terminal_codon_variant", "protein_altering_variant", 
                         "intron_variant", "up/downstream gene variant (VEP)")

bar_totals[, method_label := paste0(method, "\nn=", total_N)]

plot_dt <- left_join(plot_dt, bar_totals, by=c("method", "expression_direction"))

plot_dt[, method_label := as.factor(method_label)]

p_n <- ggplot(plot_dt[!Consequence_most_severe %in% no_show_consequence], aes(x=prop, y=method, fill = Consequence_most_severe)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = bar_totals,
            aes(y = method, x = 0.095, label = paste0("n = ", total_N)),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ expression_direction, ncol = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.1)) +
  labs(y="", x="Proportion of outliers with\nrare variants") + 
  scale_fill_manual(values = variant_colors) + 
  theme_bw(base_size=10)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/outlier_proportion_rare_variant_outrider_noIntron.png") ,
    width = 10, height = 6, units = "in", res = 600)
p_n
dev.off()



