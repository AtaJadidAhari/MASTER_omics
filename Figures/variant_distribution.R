library(ggplot2)
library(data.table)

drop_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")
benchmark_criteria <- "HIGH"

variant_base_path <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/" 
snv_vep_res <- fread(paste0(variant_base_path, "vep_res_rare_snv_", benchmark_criteria, "_aggregated.tsv"))
indel_vep_res <- fread(paste0(variant_base_path, "vep_res_rare_indel_", benchmark_criteria, "_aggregated.tsv"))

vep_res_combinded <- rbind(snv_vep_res, indel_vep_res)

print(nrow(vep_res_combinded))

# Unique_based_on_criteria
if (benchmark_criteria == "CADD_PHRED"){
  vep_res_combinded <- snv_vep_res[order(-CADD_PHRED)][, .SD[1], by = .(sampleID, Gene)]
} else{
  vep_res_combinded <- vep_res_combinded[, IMPACT := factor(IMPACT, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"), ordered = TRUE)]
  setorder(vep_res_combinded, Gene, sampleID, IMPACT)
}

vep_res_combinded <- unique(vep_res_combinded, by = c("Gene", "sampleID"))
print(nrow(vep_res_combinded))



vep_res_combinded <- merge(vep_res_combinded, drop_sa[, c("pid", "Diag", "seq_type")], by.x="sampleID", by.y="pid")
vep_res_combinded <- vep_res_combinded[Diag != "Unstranded_data",]

total_counts <- ggplot(vep_res_combinded) + 
  geom_bar(stat = "count", aes(Diag)) +
  theme_minimal() +
  labs(y="VEP rare high impact variants") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

variant_counts <- as.data.table(vep_res_combinded %>%
  group_by(Diag, sampleID) %>%
  summarise(n_variants = n(), .groups = "drop"))

ggplot(variant_counts, aes(x = Diag, y = n_variants)) +
  geom_boxplot(outlier.colour = "orange", fill = "lightblue") +
  theme_minimal() +
  scale_y_log10() +
  labs(x = "Diag", y = "VEP rare high impact") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 12)
  )

variant_counts <- as.data.table(vep_res_combinded %>%
  group_by(seq_type, sampleID) %>%
  summarise(n_variants = n(), .groups = "drop"))

ggplot(variant_counts, aes(x = seq_type, y = n_variants)) +
  geom_boxplot(outlier.colour = "orange", fill = "lightblue") +
  theme_minimal() +
  scale_y_log10() +
  labs(x = "Sequence type", y = "VEP rare high impact") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 12)
  )

# variant_counts as before
n_per_group <- variant_counts %>%
  group_by(Diag, seq_type) %>%
  summarise(N = n(), .groups = "drop")


variant_counts <- vep_res_combinded %>%
  group_by(Diag, sampleID, seq_type) %>%
  summarise(n_variants = n(), .groups = "drop")

# Plot boxplot: x = Diag, fill = seqType
ggplot(variant_counts, aes(x = Diag, y = n_variants, fill = seq_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.colour = "orange") +
  geom_text(
    data = n_per_group,
    aes(x = Diag, y = max(variant_counts$n_variants) * 1.05, label = N, group = seq_type),
    position = position_dodge(width = 0.8),
    vjust = 0,
    size = 4
  ) +
  theme_minimal(base_size=12) +
  labs(x = "Diagnosis", y = "Number of Variants") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 12)
  ) +
  scale_fill_manual(values = c("WES" = "skyblue", "WGS" = "orange")) + 
  scale_y_log10()





