library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

source('~/Scripts/Preprocess/VEP_SORanking.R')
source("~/Scripts/ProteinExpression/load_pr_data.R")


gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))

snv_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_PROTRIDER_aggregated.tsv")
indel_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_indel_PROTRIDER_aggregated.tsv")

vep_res_combinded <- rbind(snv_vep_res, indel_vep_res)

vep_res_combinded <- vep_res_combinded[, IMPACT := factor(IMPACT, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"), ordered = TRUE)]

vep_res_combinded[, Consequence_most_severe := sapply(strsplit(Consequence, ","), function(x) {
  x <- trimws(x)  # remove extra spaces
  ranked <- x[order(severity_rank[x])]
  ranked[1]  # most severe
})]
setorder(vep_res_combinded, Gene, sampleID, IMPACT)

vep_res_combinded <- unique(vep_res_combinded, by = c("Gene", "sampleID"))

pr_output_name <- "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"

pr_res <-load_pr_data(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_", pr_output_name, "/protrider_summary.csv"))
pr_res[, geneID := gene_id]

pr_res <- pr_res %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))

setnames(pr_res, "aberrant", "PROTEIN_outlier")



protrider_proteins <- unique(pr_res$geneID_short)
protrider_samples <- unique(pr_res$sampleID)


or_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv")
or_res <- or_res[sampleID %in% protrider_samples]
or_res[, geneID:= sub("\\..*", "", geneID)]

or_genes <- unique(or_res[, geneID])

all_genes <- intersect(or_genes, protrider_proteins)

merged_df <- full_join(or_res[geneID %in% all_genes], pr_res[geneID_short %in% all_genes], by=c("sampleID", "geneID"="geneID_short"))

merged_df[, annotation := "Non outlier"]
merged_df[aberrant == TRUE & PROTEIN_outlier == TRUE, annotation := "RNA + protein outlier"]
merged_df[(is.na(aberrant) | aberrant == FALSE) & PROTEIN_outlier == TRUE, annotation := "Protein outlier"]
merged_df[aberrant == TRUE & (is.na(PROTEIN_outlier) | PROTEIN_outlier == FALSE), annotation := "RNA outlier"]




merged_df_vars <- left_join(merged_df, vep_res_combinded[, c("Gene", "sampleID", "Consequence", "IMPACT", "Consequence_most_severe")], by=c("geneID" = "Gene", "sampleID"))

merged_df_vars[,`Rare variant` := FALSE]
merged_df_vars[!is.na(Consequence), `Rare variant` := TRUE]


setnames(merged_df_vars, "zScore.y", "PROTEIN_ZSCORE")
setnames(merged_df_vars, "zScore.x", "zScore")

pr_output_name <- "PROTRIDER"

# Subset the three annotated groups fully
# subset_df <- merged_df_vars %>%
#   filter(annotation != "Non_outlier" | (annotation == "Non_outlier" & row_number() %in% sample(which(annotation == "Non_outlier"), 10000)))

cor_label <- cor(merged_df_vars$zScore.x, merged_df_vars$zScore.y, method="pearson", use = "complete.obs")
range_limit <- range(c(merged_df_vars$zScore.x, merged_df_vars$zScore.y), na.rm = TRUE)

p1 <- ggplot(merged_df_vars, aes(x = zScore.x, y = zScore.y, color = annotation)) +
  geom_point(alpha = 0.7, size = 3, aes(shape=`Rare variant`)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal(base_size = 10) +
  labs(
    x = "RNA z-score",
    y = "Protein z-score",
    color = "Outlier class"
  ) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) +
  scale_color_manual(values=c("grey", "lightblue", "red", "lightgreen")) +
  annotate("text",
           x = min(range_limit)+5,
           y = max(range_limit),
           label = paste0("r=", round(cor_label, 3)),
           size = 5) +
  coord_fixed() +
  xlim(range_limit) +
  ylim(range_limit)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "_rna_pr_zScores.png") ,
    width = 6, height = 6, units = "in", res = 600)
p1
dev.off()


cor_label <- cor(merged_df_vars$zScore, merged_df_vars$PROTEIN_ZSCORE, method="pearson", use = "complete.obs")
range_limit <- range(c(merged_df_vars$zScore.x, merged_df_vars$zScore.y), na.rm = TRUE)

p1.5 <- ggplot(merged_df_vars, aes(x = zScore, y = PROTEIN_ZSCORE)) +
  geom_point(alpha = 0.7, size = 3, aes(shape=`Rare variant`)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal(base_size = 10) +
  labs(
    x = "RNA z-score",
    y = "Protein z-score",
    color = "Outlier class"
  ) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) +
  annotate("text",
           x = min(range_limit)+5,
           y = max(range_limit),
           label = paste0("r=", round(cor_label, 3)),
           size = 5) +
  coord_fixed() +
  xlim(range_limit) +
  ylim(range_limit)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "_rna_pr_zScores_noAnnotation.png") ,
    width = 6, height = 6, units = "in", res = 600)
p1.5
dev.off()


cor_labels <- merged_df_vars %>%
  group_by(annotation) %>%
  summarize(cor = cor(zScore, PROTEIN_ZSCORE, method="pearson", use = "complete.obs")) %>%
  mutate(label = paste("r =", round(cor, 2)))

p2 <- ggplot(merged_df_vars, aes(x = zScore, y = PROTEIN_ZSCORE)) +
  geom_point(alpha = 0.7, size = 3, aes(shape=`Rare variant`)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ annotation) +
  theme_bw(base_size = 10) +
  labs(
    x = "RNA z-score",
    y = "Protein z-score",
    color = "Outlier class"
  ) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16))  +
  geom_text(
    data = cor_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE,
    size = 5
  ) +
  coord_fixed() +
  xlim(range_limit) +
  ylim(range_limit)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "_rna_pr_zScores_facetted.png") ,
    width = 6, height = 6, units = "in", res = 600)
p2
dev.off()

missense_df <- merged_df_vars[grepl("missense", Consequence)]
missense_corr <- cor(missense_df$zScore, missense_df$PROTEIN_ZSCORE, method="pearson", use = "complete.obs")
range_limit <- range(c(missense_df$zScore, missense_df$PROTEIN_ZSCORE), na.rm = TRUE)


p3 <- ggplot(missense_df, aes(x = zScore, y = PROTEIN_ZSCORE, color = annotation)) +
  geom_point(alpha = 0.7, size = 3, aes(shape=`Rare variant`)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal(base_size = 10) +
  labs(
    x = "RNA z-score",
    y = "Protein z-score",
    color = "Outlier class"
  ) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) +
  scale_color_manual(values=c("grey", "lightblue", "red", "lightgreen")) +
  annotate("text", x = min(range_limit) + 5, y = max(range_limit) - 5,
           label = paste0("r = ", round(missense_corr, 3)), size = 5) +
  coord_fixed() +
  xlim(range_limit) +
  ylim(range_limit)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "_rna_pr_zScores_missense.png") ,
    width = 6, height = 6, units = "in", res = 600)
p3
dev.off()



#### proportion plots ####

merged_df_vars[ , pr_expresssion := "Non-outlier"]
merged_df_vars[PROTEIN_ZSCORE > 0 & PROTEIN_outlier == TRUE, pr_expresssion := "Overexpression"]
merged_df_vars[PROTEIN_ZSCORE < 0 & PROTEIN_outlier == TRUE, pr_expresssion := "Underexpression"]

merged_df_vars[, method:= pr_output_name]

plot_dt <- merged_df_vars[, .N, by = .(method, pr_expresssion, Consequence_most_severe)][
  , prop := N / sum(N), by = .(method, pr_expresssion)]

plot_dt <- plot_dt[!is.na(plot_dt$Consequence_most_severe), ]


p_n <- ggplot(plot_dt, aes(y=prop, x = method, fill = Consequence_most_severe)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +                # horizontal bars
  facet_wrap(~ pr_expresssion, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size=10)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "proportion_rare_variant_protrider.png") ,
    width = 6, height = 6, units = "in", res = 600)
p_n
dev.off()


  

#### load size-factor_normalized zScores ####
protein_sf_zScores <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/results/sf_zScores.csv")  
rna_sf_zScores <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_sf_zScores_all.tsv", sep="\t")

protrider_proteins <- unique(protein_sf_zScores$geneID_short)
protrider_samples <- unique(protein_sf_zScores$sampleID)
setnames(protein_sf_zScores, "zScore", "PROTEIN_ZSCORE")

rna_sf_zScores <- rna_sf_zScores[sampleID %in% protrider_samples]
rna_sf_zScores[, geneID:= sub("\\..*", "", geneID)]

or_genes <- unique(rna_sf_zScores[, geneID])

all_genes <- intersect(or_genes, protrider_proteins)

merged_df <- full_join(rna_sf_zScores[geneID %in% all_genes], protein_sf_zScores[geneID_short %in% all_genes], by=c("sampleID", "geneID"="geneID_short"))



merged_df_vars <- left_join(merged_df, vep_res_combinded[, c("Gene", "sampleID", "Consequence", "IMPACT")], by=c("geneID" = "Gene", "sampleID"))

merged_df_vars[,`Rare variant` := FALSE]  
merged_df_vars[!is.na(Consequence), `Rare variant` := TRUE]  


set.seed(123)  # for reproducibility

# Subset the three annotated groups fully
# subset_df <- merged_df_vars %>%
#   filter(annotation != "Non_outlier" | (annotation == "Non_outlier" & row_number() %in% sample(which(annotation == "Non_outlier"), 10000)))

cor_label <- cor(merged_df_vars$zScore, merged_df_vars$PROTEIN_ZSCORE, method="pearson", use = "complete.obs")
range_limit <- range(c(merged_df_vars$zScore, merged_df_vars$PROTEIN_ZSCORE), na.rm = TRUE)

p4 <- ggplot(merged_df_vars, aes(x = zScore, y = PROTEIN_ZSCORE)) +
  geom_point(alpha = 0.7, size = 3, aes(shape=`Rare variant`)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal(base_size = 10) +
  labs(
    x = "RNA z-score",
    y = "Protein z-score",
    color = "Outlier class"
  ) + 
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) + 
  annotate("text", 
           x = min(range_limit)+5, 
           y = max(range_limit), 
           label = paste0("r=", round(cor_label, 3)), 
           size = 5) +
  coord_fixed() +
  xlim(range_limit) +
  ylim(range_limit)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "_rna_pr_sf_zScores.png") ,
    width = 6, height = 6, units = "in", res = 600)
p4
dev.off()







missense_df <- merged_df_vars[grepl("missense", Consequence)]
missense_corr <- cor(missense_df$zScore, missense_df$PROTEIN_ZSCORE, method="pearson", use = "complete.obs")
range_limit <- range(c(missense_df$zScore, missense_df$PROTEIN_ZSCORE), na.rm = TRUE)


p3 <- ggplot(missense_df, aes(x = zScore, y = PROTEIN_ZSCORE)) +
  geom_point(alpha = 0.7, size = 3, aes(shape=`Rare variant`)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal(base_size = 10) +
  labs(
    x = "RNA z-score",
    y = "Protein z-score",
    color = "Outlier class"
  ) + 
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) + 
  annotate("text", x = min(range_limit) + 5, y = max(range_limit) - 5, 
           label = paste0("r = ", round(missense_corr, 3)), size = 5) +
  coord_fixed() +
  xlim(range_limit) +
  ylim(range_limit)

png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "_rna_pr_sf_zScores_missense.png") ,
    width = 6, height = 6, units = "in", res = 600)
p3
dev.off()


