library(data.table)
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(grid)
library(qqman)

gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))

pr_output_name <- "gaussian_gs_injMean3"
pr_res <- fread(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_", pr_output_name, "/protrider_summary.csv"))

protein_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protrider_sa.tsv")

pr_res <- left_join(pr_res, protein_sa, by=c("sampleID"="full_name"))
pr_res <- left_join(pr_res, gene_annot_dt, by=c("proteinID"="gene_name"))
pr_res[, pos := (start + end) / 2]

pr_res[Diag == "COAD" & proteinID == "GMFB"]
pr_res_sample <- pr_res[pid == "DVXCP3"]
pr_res_sample[, chr_len := max(pos), by = seqnames]
pr_res_sample[, pos_norm := pos / chr_len]

pr_res_sample[, chr := factor(seqnames, levels = chr_levels)]
# Fixed spacing: shift chromosomes by their index
offset <- 1
pr_res_sample[, chr_index := as.integer(chr)]
pr_res_sample[, cum_pos := pos_norm + (chr_index - 1) * offset]

# Midpoints for axis labels
axis_labels <- pr_res_sample[, .(mid = mean(cum_pos)), by = chr]
axis_labels <- axis_labels[order(factor(chr, levels = chr_levels))]  # ensure order
label_shift <- 0.5  # in "chromosome units"

# Manhattan plot
ggplot(pr_res_sample, aes(x = cum_pos, y = -log10(PROTEIN_PVALUE), color = PROTEIN_outlier)) +
  geom_point(size = 1) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = axis_labels$mid,
    labels = axis_labels$chr
  ) +
  scale_y_continuous(name = "-log10(p-value)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values=c("grey", "firebrick")) +
  labs(title = "Manhattan Plot for Protein outliers")



pr_res_aberrant <- pr_res[PROTEIN_outlier == T]


gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))


rna_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")


rna_sample_stranded <- rna_sa[Diag != "Unstranded_data"]




protein_intensities <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protrider_intensities.tsv")

all_proteins <- unique(protein_intensities$protein_ID)
non_na_proteins <- unique(pr_res$proteinID)


expressed_rna <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv")
expressed_rna <- left_join(expressed_rna, gene_annot_dt, by=c("geneID"="gene_id"))
or_res_aberrant <- expressed_rna[aberrant == T]


expressed_rna_genes <- unique(expressed_rna$gene_name)




listInput <- list(`Expressed genes` = expressed_rna, `All proteins` = all_proteins, 
                  `Non-na proteins` = non_na_proteins)

p <- upset(fromList(listInput), order.by = "freq", text.scale = 1.5)

# p + grid.text("Underexpression outliers",x = 0.65, y=0.95, gp=gpar(fontsize=20))

png("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/common_genes.png", 
    width = 8, height = 6, units = "in", res = 600)
p
dev.off()
    


both_pr_rna_samples <- intersect(stranded_rna_samples, pr_samples)



