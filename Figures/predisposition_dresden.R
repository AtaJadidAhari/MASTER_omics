library(data.table)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)

source("~/Scripts/utils/load_gtf_cgc_dresden.R")

or_res_all <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_outliers.tsv")
fr_res_all_genes <- fread( "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_results.tsv")

or_res_all[, predosposition_gene := FALSE]
or_res_all[hgncSymbol %in% dresden_list, predosposition_gene := TRUE]

cohort_props <- as.data.table(or_res_all %>%
  group_by(Diag) %>%
  summarise(proportion = mean(predosposition_gene, na.rm = TRUE)))

bar_plot <- ggplot(cohort_props, aes(y = Diag, x = proportion)) +
  geom_bar(stat = "identity", fill = "orange") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Proportion", title = "Proportion of Predisposed Outliers")



cohort_totals <- or_res_all %>%
  group_by(Diag) %>%
  summarise(total_outlier_genes = n(), .groups = "drop")

# Step 2: For each gene in dresden_list, count predisposed outlier occurrences per cohort (not distinct)
gene_details <- or_res_all %>%
  filter(hgncSymbol %in% dresden_list, predosposition_gene == TRUE) %>%
  group_by(Diag, hgncSymbol) %>%
  summarise(gene_predisposed_outlier_genes = sum(predosposition_gene == TRUE), .groups = "drop") %>%
  pivot_wider(names_from = hgncSymbol,
              values_from = gene_predisposed_outlier_genes,
              values_fill = 0)

# Step 3: Merge with cohort totals
cohort_with_genes <- left_join(cohort_totals, gene_details, by = "Diag")


# intermediate step: calculate enrichments
# dt is your original data.table
dt_long <- as.data.table(melt(cohort_with_genes, id.vars = c("Diag", "total_outlier_genes"), 
                variable.name = "Gene", value.name = "IsOutlier"))

results <- dt_long[, {
  a <- sum(IsOutlier == 1)   # Outliers for gene in this diag
  b <- sum(IsOutlier == 0)   # Non-outliers for gene in this diag
  
  # For all other diags
  other <- dt_long[Diag != .BY$Diag & Gene == .BY$Gene]
  c <- sum(other$IsOutlier == 1)
  d <- sum(other$IsOutlier == 0)
  
  test <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  
  list(Out_in_diag = a, Not_out_in_diag = b, Out_not_diag = c, Not_out_not_diag = d,
       OddsRatio = test$estimate,
       P = test$p.value)
}, by = .(Diag, Gene)]

results[, P_adj := p.adjust(P, method = "BH")]


# Step 4: Calculate proportions for each gene column (no prefixes)
gene_cols <- setdiff(colnames(cohort_with_genes), c("Diag", "total_outlier_genes"))
for (g in gene_cols) {
  cohort_with_genes[[g]] <- cohort_with_genes[[g]] / cohort_with_genes$total_outlier_genes
}

# Step 5: Prepare matrix for heatmap
heatmap_mat <- as.matrix(cohort_with_genes[, gene_cols])
rownames(heatmap_mat) <- cohort_with_genes$Diag

# Step 6: Custom color palette for better visibility
color_palette <- colorRampPalette(c("white", "orange", "red"))(100)

# Step 7: Plot the heatmap
p2 <- pheatmap(heatmap_mat,
         color = color_palette,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         main = "Proportion of Predisposed Outlier Events per Cohort and Gene")

rownames(cohort_props) <- cohort_props$Diag
# cohort_props <- cohort_props[rownames(heatmap_mat), ]

# Define color mapping for the heatmap
col_fun = colorRamp2(c(0, max(heatmap_mat)), c("white", "red"))

# Make barplot annotation using anno_barplot
bar_anno <- rowAnnotation(
  " " = anno_barplot(
    cohort_props$proportion,
    gp = gpar(fill = "orange"),
    width = unit(2, "cm"),
    axis_param = list(direction = "reverse") # bars to the right
  )
)

# Plot heatmap + annotation
Heatmap(
  heatmap_mat,
  name = "Proportion of outlier predisposition\ngenes across outliers ",
  col = col_fun,
  right_annotation = bar_anno,
  cluster_rows = TRUE, # or FALSE for ordered Diag
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 8),        # Reduce row labels font size
  column_names_gp = gpar(fontsize = 8)     # Reduce column labels font size
)

or_res_all_predosposition <- or_res_all[hgncSymbol %in% dresden_list]


or_res_all_predosposition[, expression_direction := "Underexpression"]
or_res_all_predosposition[zScore > 0, expression_direction := "Overexpression"]
table(or_res_all_predosposition$expression_direction)
table(or_res_all_predosposition$Diag)


or_res_all_predosposition_extended <- or_res_all[hgncSymbol %in% extended_dresden_list]
or_res_all_predosposition_extended[, expression_direction := "Underexpression"]
or_res_all_predosposition_extended[zScore > 0, expression_direction := "Overexpression"]
table(or_res_all_predosposition_extended$expression_direction)

or_res_all_predosposition_extended


fr_res_all_predosposition <- fr_res_all_genes[hgncSymbol %in% dresden_list]
table(fr_res_all_predosposition$Diag)


fr_res_all_predosposition_extended <- fr_res_all_genes[hgncSymbol %in% extended_dresden_list]

table(fr_res_all_predosposition_extended$Diag)


cohort_props <- as.data.table(fr_res_all_genes %>%
  group_by(Diag) %>%
  summarise(proportion = mean(predosposition_gene, na.rm = TRUE)))

bar_plot <- ggplot(cohort_props, aes(y = Diag, x = proportion)) +
  geom_bar(stat = "identity", fill = "orange") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Proportion", title = "Proportion of Predisposed Outliers")



cohort_totals <- fr_res_all_genes %>%
  group_by(Diag) %>%
  summarise(total_outlier_genes = n(), .groups = "drop")

# Step 2: For each gene in dresden_list, count predisposed outlier occurrences per cohort (not distinct)
gene_details <- fr_res_all_genes %>%
  filter(hgncSymbol %in% dresden_list, predosposition_gene == TRUE) %>%
  group_by(Diag, hgncSymbol) %>%
  summarise(gene_predisposed_outlier_genes = sum(predosposition_gene == TRUE), .groups = "drop") %>%
  pivot_wider(names_from = hgncSymbol,
              values_from = gene_predisposed_outlier_genes,
              values_fill = 0)

# Step 3: Merge with cohort totals
cohort_with_genes <- left_join(cohort_totals, gene_details, by = "Diag")


# intermediate step: calculate enrichments
# dt is your original data.table
dt_long <- as.data.table(melt(cohort_with_genes, id.vars = c("Diag", "total_outlier_genes"), 
                variable.name = "Gene", value.name = "IsOutlier"))

results <- dt_long[, {
  a <- sum(IsOutlier == 1)   # Outliers for gene in this diag
  b <- sum(IsOutlier == 0)   # Non-outliers for gene in this diag
  
  # For all other diags
  other <- dt_long[Diag != .BY$Diag & Gene == .BY$Gene]
  c <- sum(other$IsOutlier == 1)
  d <- sum(other$IsOutlier == 0)
  
  test <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  
  list(Out_in_diag = a, Not_out_in_diag = b, Out_not_diag = c, Not_out_not_diag = d,
       OddsRatio = test$estimate,
       P = test$p.value)
}, by = .(Diag, Gene)]

results[, P_adj := p.adjust(P, method = "BH")]


# Step 4: Calculate proportions for each gene column (no prefixes)
gene_cols <- setdiff(colnames(cohort_with_genes), c("Diag", "total_outlier_genes"))
for (g in gene_cols) {
  cohort_with_genes[[g]] <- cohort_with_genes[[g]] / cohort_with_genes$total_outlier_genes
}

# Step 5: Prepare matrix for heatmap
heatmap_mat <- as.matrix(cohort_with_genes[, gene_cols])
rownames(heatmap_mat) <- cohort_with_genes$Diag

# Step 6: Custom color palette for better visibility
color_palette <- colorRampPalette(c("white", "orange", "red"))(100)

# Step 7: Plot the heatmap
p2 <- pheatmap(heatmap_mat,
         color = color_palette,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         main = "Proportion of Predisposed Outlier Events per Cohort and Gene")

rownames(cohort_props) <- cohort_props$Diag
# cohort_props <- cohort_props[rownames(heatmap_mat), ]

# Define color mapping for the heatmap
col_fun = colorRamp2(c(0, max(heatmap_mat)), c("white", "red"))

# Make barplot annotation using anno_barplot
bar_anno <- rowAnnotation(
  " " = anno_barplot(
    cohort_props$proportion,
    gp = gpar(fill = "orange"),
    width = unit(2, "cm"),
    axis_param = list(direction = "reverse") # bars to the right
  )
)

# Plot heatmap + annotation
Heatmap(
  heatmap_mat,
  name = "Proportion of outlier predisposition\ngenes across outliers ",
  col = col_fun,
  right_annotation = bar_anno,
  cluster_rows = TRUE, # or FALSE for ordered Diag
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 8),        # Reduce row labels font size
  column_names_gp = gpar(fontsize = 8)     # Reduce column labels font size
)

