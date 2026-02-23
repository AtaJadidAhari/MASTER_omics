library(pheatmap)
library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(ggplot2)




gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))

cgc <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/Cosmic_CancerGeneCensus_v102_GRCh37.tsv")
cgc_tsg <- cgc[grepl("TSG", ROLE_IN_CANCER), ]
cgc_tsg <- inner_join(gene_annot_dt, cgc_tsg, by=c("gene_name"="GENE_SYMBOL"))
# loosing 3 genes

cgc_oncogene <- cgc[grepl("oncogene", ROLE_IN_CANCER), ]
cgc_oncogene <- inner_join(gene_annot_dt, cgc_oncogene, by=c("gene_name"="GENE_SYMBOL"))


pr_output_name <- "gaussian_gs_injMean3"
pr_res <- fread(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_", pr_output_name, "/protrider_summary.csv"))


pr_res_aberrant <- pr_res[PROTEIN_outlier == T]

# merge with sa to get diags
protein_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protrider_sa.tsv")

# protein_sa <- left_join(rna_sa, protein_sa, by="pid")

pr_res <- left_join(pr_res, protein_sa, by=c("sampleID"="full_name"))





create_fisher_dt <- function(aberrant_dt, expression_direction, gene_list, omics_type="PROTRIDER", cluster_col="Diag"){
  
  if (omics_type == "PROTRIDER"){
    omic_id <- "proteinID"
    zscore_col <- "PROTEIN_ZSCORE"
    outlier_col <- "PROTEIN_outlier"
  } else{
    omic_id <- "geneID"
    zscore_col <- "zScore"
    outlier_col <- "aberrant"
  }
  
  
  #aberrant_events_dt_filtered <- aberrant_dt[get(omic_id) %in% gene_list]
  aberrant_events_dt_filtered <- aberrant_dt
  if (expression_direction == "Underexpression"){
    aberrant_events_dt_filtered <- aberrant_events_dt_filtered[get(zscore_col) < 0]
  } else if(expression_direction == "Overexpression"){
    aberrant_events_dt_filtered <- aberrant_events_dt_filtered[get(zscore_col) > 0]
  }
  
  setnames(aberrant_events_dt_filtered, old = omic_id, new = "proteinID")  
  setnames(aberrant_events_dt_filtered, old = cluster_col, new = "cluster_column")  
  
  
  A_counts <- aberrant_events_dt_filtered[
    get(outlier_col) == TRUE,
    uniqueN(sampleID),
    by = c("proteinID", "cluster_column")
  ]
  setnames(A_counts, "V1", "A")
  # B: Samples without (aberrant) Gene G in Disease D
  # This is the total number of samples in a disease minus those with the aberrant gene.
  # First, get total unique samples per disease.
  total_samples_per_disease <- aberrant_events_dt_filtered[, uniqueN(sampleID), by = cluster_column]
  setnames(total_samples_per_disease, "V1", "TotalSamplesInDisease")
  
  # Create a full grid of all gene-disease combinations to ensure all B, C, D are covered
  all_gene_disease_combinations <- CJ(proteinID = unique(aberrant_events_dt_filtered[, proteinID]),
                                      cluster_column = unique(aberrant_events_dt_filtered$cluster_column))
  # Merge A_counts into the full grid
  plot_data <- merge(all_gene_disease_combinations, A_counts, by = c("proteinID", "cluster_column"), all.left = TRUE)
  plot_data[is.na(A), A := 0] # Set A to 0 if no aberrant samples found for that combo
  
  # Merge total samples per disease
  plot_data <- merge(plot_data, total_samples_per_disease, by = "cluster_column", all.left = TRUE)
  
  # Calculate B: Samples in Disease D that do NOT have aberrant Gene G
  plot_data[, B := TotalSamplesInDisease - A]
  # C: Samples with (aberrant) Gene G NOT in Disease D
  # We need total aberrant samples for each gene across ALL diseases first.
  total_aberrant_gene_counts <- aberrant_events_dt_filtered[
    get(outlier_col) == TRUE,
    .(TotalAberrantGeneOverall = uniqueN(sampleID)),
    by = "proteinID"
  ]
  # setnames(total_aberrant_gene_counts, "V1", "TotalAberrantGeneOverall")
  
  # Merge total_aberrant_gene_counts
  plot_data <- merge(plot_data, total_aberrant_gene_counts, by = "proteinID", all.left = TRUE)
  
  plot_data[, C := TotalAberrantGeneOverall - A]
  # D: Samples without (aberrant) Gene G NOT in Disease D
  # Total samples overall for each gene:
  total_samples_overall_per_gene <- aberrant_events_dt_filtered[
    ,
    .(TotalGeneSamplesOverall = uniqueN(sampleID)),
    by = "proteinID"
  ]
  # setnames(total_samples_overall_per_gene, "V1", "TotalGeneSamplesOverall")
  plot_data <- merge(plot_data, total_samples_overall_per_gene, by = "proteinID", all.left = TRUE)
  
  # Samples NOT in Disease D for a specific gene:
  # Total samples for that gene MINUS samples for that gene IN the current Disease_entity
  samples_not_in_disease_for_gene <- aberrant_events_dt_filtered[
    ,
    .(SamplesInDiseaseForGene = uniqueN(sampleID)),
    by = c("proteinID", "cluster_column")
  ]
  # setnames(samples_not_in_disease_for_gene, "V1", "SamplesInDiseaseForGene")
  plot_data <- merge(plot_data, samples_not_in_disease_for_gene, by = c("proteinID", "cluster_column"), all.left = TRUE)
  plot_data[is.na(SamplesInDiseaseForGene), SamplesInDiseaseForGene := 0]
  
  plot_data[, D := (TotalGeneSamplesOverall - SamplesInDiseaseForGene) - C]
  
  
  # --- 3. Calculate Odds Ratio ---
  epsilon <- 0.5 # For stable calculation with zero counts
  
  plot_data[, odds_ratio := {
    or_val = ((A + epsilon) / (B + epsilon)) / ((C + epsilon) / (D + epsilon))
    replace(or_val, is.infinite(or_val) | is.nan(or_val), NA) # Handle Inf/NaN
  }, by = .(proteinID, cluster_column)]
  
  # Cap odds ratio at 1 if less than 1, as per the plot's visual representation
  plot_data[odds_ratio < 1, odds_ratio := 1]
  
}

prepare_plot_data <- function(plot_data){
  
  # --- 4. Prepare data for plotting (categorize Odds Ratio) ---
  my_colors <- c(
    "1" = "#FFFFFF",        # White or very light grey for 1
    "1-10" = "#D9D1E8",     # Light purple for 1-10
    "10-100" = "#A38ECF",   # Medium purple for 10-100
    "100-1000" = "#6D4FB6", # Darker purple for 100-1000
    "ge1000" = "#441E7F"    # Darkest purple for >=1000
  )
  
  # Define breaks and labels for the legend
  or_breaks <- c(1, 10, 100, 1000, Inf)
  or_labels <- c("1-10", "10-100", "100-1000", "ge1000") # Exclude the '1' as a range start point
  
  plot_data[, odds_ratio_category := cut(odds_ratio,
                                         breaks = c(0, or_breaks), # Start from 0 to capture OR=1
                                         labels = c("1", or_labels),
                                         right = FALSE, # [lower, upper)
                                         include.lowest = TRUE)] # Include 1 in the "1" category
  
  plot_data$cluster_column <- factor(plot_data$cluster_column)
  plot_data$proteinID <- factor(plot_data$proteinID)
  
  # Filter out NA levels that might occur if some genes/diseases from the plot are not in your data
  plot_data <- plot_data[!is.na(cluster_column) & !is.na(proteinID)]
  
  
  
  
  odds_ratio_matrix <- as.data.table(dcast(plot_data, cluster_column ~ proteinID, value.var = "odds_ratio"))
  odds_ratio_matrix[is.na(odds_ratio_matrix)] <- 0
  
  # Set row names for the matrix (Disease_entity)
  odds_ratio_matrix_rownames <- odds_ratio_matrix$cluster_column
  
  odds_ratio_matrix[, cluster_column := NULL] # Remove the ID column
  odds_ratio_matrix <- as.matrix(odds_ratio_matrix)
  rownames(odds_ratio_matrix) <- odds_ratio_matrix_rownames
  return (odds_ratio_matrix)
  
}





pr_res <- pr_res[Diag!="Unstranded_data"]
plot_data <- create_fisher_dt(pr_res, "Overexpression", cgc_oncogene$gene_name, omics_type="PROTRIDER")
plot_data <- create_fisher_dt(pr_res, "Overexpression", unique(pr_res$proteinID), omics_type="PROTRIDER")

odds_ratio_matrix <- prepare_plot_data(plot_data)



idx <- which(odds_ratio_matrix > 1, arr.ind = TRUE)

# Build the list of (row_name, col_name, value)
pr_dt <- data.table(
  row_name = rownames(odds_ratio_matrix)[idx[, "row"]],
  col_name = colnames(odds_ratio_matrix)[idx[, "col"]],
  odds_pr_overexpression = odds_ratio_matrix[idx]
)
pr_dt[, gene_cancer := paste0(row_name, "_", col_name)]

pr_dt <- pr_dt[order(-odds_pr_overexpression)]
genes_of_interest <- union(pr_dt[1:100, col_name], cgc_oncogene$gene_name)

# --- 5. Create the Heatmap ---
cols_to_keep <- apply(odds_ratio_matrix, 2, function(x) any(x > 1))

# Subset the matrix
mat_filtered <- odds_ratio_matrix[, intersect(colnames(odds_ratio_matrix), genes_of_interest), drop = FALSE]

breaks_pheatmap <- c(0, 1, 10, 100, 1000, max(odds_ratio_matrix, na.rm = TRUE) + 1) # Define breaks
colors_pheatmap <- colorRampPalette(c("#FFFFFF", "#D9D1E8", "#A38ECF", "#6D4FB6", "#441E7F"))(length(breaks_pheatmap) - 1)


annotation_col <- data.frame(CGC_status = factor(
  ifelse(colnames(mat_filtered) %in% cgc_oncogene$gene_name, "Known oncogene", "Unkown gene"),
  levels = c("Known oncogene", "Unkown gene")
))
rownames(annotation_col) <- colnames(mat_filtered)

annotation_colors <- list(
  CGC_status = c("Known oncogene" = "red", "Unkown gene" = "lightblue")
)

p_heatmap <- pheatmap(mat_filtered,
         color = colors_pheatmap,
         breaks = breaks_pheatmap,
         cluster_rows = TRUE, # Enable clustering for rows (Disease entity)
         cluster_cols = TRUE, # Enable clustering for columns (CGC cancer gene)
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Odds Ratio Heatmap (Clustered)",
         fontsize_row = 9,
         fontsize_col = 8,
         angle_col = 45, # Angle column labels for readability
         # The legend is automatically generated, customize if needed.
         legend_breaks = c(1, 10, 100, 1000), # Explicit legend breaks
         legend_labels = c("1", "10", "100", "\u22651000"), # Explicit legend labels
         width = 10, height = 7,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors

) 

p_heatmap


plot_data <- create_fisher_dt(pr_res, "Underexpression", cgc_tsg$gene_name, omics_type="PROTRIDER")
plot_data <- create_fisher_dt(pr_res, "Underexpression", unique(pr_res$proteinID), omics_type="PROTRIDER")

odds_ratio_matrix <- prepare_plot_data(plot_data)

idx <- which(odds_ratio_matrix > 1, arr.ind = TRUE)

# Build the list of (row_name, col_name, value)
pr_dt <- data.table(
  row_name = rownames(odds_ratio_matrix)[idx[, "row"]],
  col_name = colnames(odds_ratio_matrix)[idx[, "col"]],
  odds_pr_overexpression = odds_ratio_matrix[idx]
)
pr_dt[, gene_cancer := paste0(row_name, "_", col_name)]

pr_dt <- pr_dt[order(-odds_pr_overexpression)]
genes_of_interest <- union(pr_dt[1:30, col_name], cgc_tsg$gene_name)

# --- 5. Create the Heatmap ---
cols_to_keep <- apply(odds_ratio_matrix, 2, function(x) any(x > 1))

# Subset the matrix
mat_filtered <- odds_ratio_matrix[, intersect(colnames(odds_ratio_matrix), genes_of_interest), drop = FALSE]

breaks_pheatmap <- c(0, 1, 10, 100, 1000, max(odds_ratio_matrix, na.rm = TRUE) + 1) # Define breaks
colors_pheatmap <- colorRampPalette(c("#FFFFFF", "#D9D1E8", "#A38ECF", "#6D4FB6", "#441E7F"))(length(breaks_pheatmap) - 1)


annotation_col <- data.frame(CGC_status = factor(
  ifelse(colnames(mat_filtered) %in% cgc_tsg$gene_name, "Known TSG", "Unkown gene"),
  levels = c("Known TSG", "Unkown gene")
))
rownames(annotation_col) <- colnames(mat_filtered)

annotation_colors <- list(
  CGC_status = c("Known TSG" = "red", "Unkown gene" = "lightblue")
)

p_heatmap <- pheatmap(mat_filtered,
                      color = colors_pheatmap,
                      breaks = breaks_pheatmap,
                      cluster_rows = TRUE, # Enable clustering for rows (Disease entity)
                      cluster_cols = TRUE, # Enable clustering for columns (CGC cancer gene)
                      show_rownames = TRUE,
                      show_colnames = TRUE,
                      main = "Odds Ratio Heatmap (Clustered)",
                      fontsize_row = 9,
                      fontsize_col = 8,
                      angle_col = 45, # Angle column labels for readability
                      # The legend is automatically generated, customize if needed.
                      legend_breaks = c(1, 10, 100, 1000), # Explicit legend breaks
                      legend_labels = c("1", "10", "100", "\u22651000"), # Explicit legend labels
                      width = 10, height = 7,
                      annotation_col = annotation_col,
                      annotation_colors = annotation_colors
                      
) 

p_heatmap



#### expression outliers ####


or_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv")
or_res[, geneID:= sub("\\..*", "", geneID)]

rna_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")
or_res <- left_join(or_res,rna_sa[, c("pid", "Oncotree Code", "SEX", "ICD10 Code", "TISSUE", "Tumorzellgehalt (Bioinformatik)")], by=c("sampleID" = "pid"))

or_res <- left_join(or_res, gene_annot_dt, by=c("geneID" = "geneID_short"))
setnames(or_res, "geneID", "geneID_short")
setnames(or_res, "gene_name", "geneID")

plot_data <- create_fisher_dt(or_res, "Overexpression", cgc_oncogene$gene_name, omics_type="OUTRIDER")
plot_data <- plot_data[TotalAberrantGeneOverall > 1]
odds_ratio_matrix <- prepare_plot_data(plot_data)

idx <- which(odds_ratio_matrix > 1, arr.ind = TRUE)

# Build the list of (row_name, col_name, value)
or_dt <- data.table(
  row_name = rownames(odds_ratio_matrix)[idx[, "row"]],
  col_name = colnames(odds_ratio_matrix)[idx[, "col"]],
  odds_pr_overexpression = odds_ratio_matrix[idx]
)
or_dt[, gene_cancer := paste0(row_name, "_", col_name)]
or_dt <- or_dt[order(-odds_pr_overexpression)]
genes_of_interest <- union(or_dt[1:30, col_name], cgc_tsg$gene_name)

breaks <- c(0, 1, 10, 100, max(odds_ratio_matrix, na.rm = TRUE) + 1)
colors <- c("#FFFFFF", "#D9D1E8", "#A38ECF", "#6D4FB6", "#441E7F")
col_fun <- colorRamp2(breaks, colors, space = "RGB")
p_heatmap <-  Heatmap(
  odds_ratio_matrix,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 90,
  name = "Odds Ratio",
  heatmap_legend_param = list(
    at = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "≥1000"),
    color_bar = "discrete",
    legend_direction = "horizontal"
  ),
  column_title = "CGC oncogene gene",
  column_names_gp = gpar(fontsize = 8)
)
png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/or_overexpression_heatmap.png") ,
    width = 12, height = 6, units = "in", res = 600)
p_heatmap
dev.off()




overexpression_intersect <- intersect(or_dt$gene_cancer, pr_dt$gene_cancer)
only_pr <- setdiff(pr_dt$gene_cancer, overexpression_intersect)
only_or <- setdiff(or_dt$gene_cancer, overexpression_intersect)

area1 <- length(unique(pr_dt$gene_cancer))
area2 <- length(unique(or_dt$gene_cancer))
cross <- length(intersect(pr_dt$gene_cancer, or_dt$gene_cancer))

# Open file device (e.g. PNG)
png("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/pr_or_odds_over_venn_diagram.png", width = 800, height = 800, res = 150)

# Draw the plot
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross,
  category = c("Protein\noverexpression", "Gene overexpression"),
  fill = c("skyblue", "tomato"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  ind = FALSE
)

# Render and close the file
grid.draw(venn.plot)
dev.off()



plot_data <- create_fisher_dt(or_res, "Underexpression", cgc_oncogene$gene_name, omics_type="OUTRIDER")
odds_ratio_matrix <- prepare_plot_data(plot_data)

idx <- which(odds_ratio_matrix > 1, arr.ind = TRUE)

# Build the list of (row_name, col_name, value)
or_un_dt <- data.table(
  row_name = rownames(odds_ratio_matrix)[idx[, "row"]],
  col_name = colnames(odds_ratio_matrix)[idx[, "col"]],
  odds_pr_underexpression = odds_ratio_matrix[idx]
)
or_un_dt[, gene_cancer := paste0(row_name, "_", col_name)]


plot_data <- create_fisher_dt(pr_res, "Underexpression", cgc_tsg$gene_name, omics_type="PROTRIDER")
odds_ratio_matrix <- prepare_plot_data(plot_data)

idx <- which(odds_ratio_matrix > 1, arr.ind = TRUE)

# Build the list of (row_name, col_name, value)
pr_un_dt <- data.table(
  row_name = rownames(odds_ratio_matrix)[idx[, "row"]],
  col_name = colnames(odds_ratio_matrix)[idx[, "col"]],
  odds_pr_underexpression = odds_ratio_matrix[idx]
)
pr_un_dt[, gene_cancer := paste0(row_name, "_", col_name)]


underexpression_intersect <- intersect(or_un_dt$gene_cancer, pr_un_dt$gene_cancer)
only_pr <- setdiff(pr_un_dt$gene_cancer, underexpression_intersect)
only_or <- setdiff(or_un_dt$gene_cancer, underexpression_intersect)

area1 <- length(unique(only_pr))
area2 <- length(unique(only_or))
cross <- length(underexpression_intersect)

# Open file device (e.g. PNG)
png("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/pr_or_odds_un_venn_diagram.png", width = 800, height = 800, res = 150)

# Draw the plot
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross,
  category = c("Protein\nunderexpression", "Gene underexpression"),
  fill = c("skyblue", "tomato"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  ind = FALSE
)

# Render and close the file
grid.draw(venn.plot)
dev.off()


pr_res_aberrant <- pr_res[PROTEIN_outlier == T]
or_res_aberrant <- or_res[aberrant == T]

pr_res_aberrant[, gene_sample := paste0(proteinID, "_", pid)]
or_res_aberrant[, gene_sample := paste0(geneID, "_", sampleID)]



underexpression_intersect <- intersect(or_res_aberrant[zScore < 0, gene_sample], pr_res_aberrant[PROTEIN_ZSCORE < 0, gene_sample])
only_pr <- setdiff(pr_res_aberrant[PROTEIN_ZSCORE < 0, gene_sample], underexpression_intersect)
only_or <- setdiff(or_res_aberrant[zScore < 0, gene_sample], underexpression_intersect)

area1 <- length(unique(only_pr))
area2 <- length(unique(only_or))
cross <- length(underexpression_intersect)
# Open file device (e.g. PNG)
png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "/pr_or_outliers_un_venn_diagram.png"), width = 800, height = 800, res = 150)

# Draw the plot
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross,
  category = c("Protein\nunderexpression", "Gene underexpression"),
  fill = c("skyblue", "tomato"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  ind = FALSE
)

# Render and close the file
grid.draw(venn.plot)
dev.off()



underexpression_intersect <- intersect(or_res_aberrant[zScore > 0, gene_sample], pr_res_aberrant[PROTEIN_ZSCORE > 0, gene_sample])
only_pr <- setdiff(pr_res_aberrant[PROTEIN_ZSCORE > 0, gene_sample], underexpression_intersect)
only_or <- setdiff(or_res_aberrant[zScore > 0, gene_sample], underexpression_intersect)

area1 <- length(unique(only_pr))
area2 <- length(unique(only_or))
cross <- length(underexpression_intersect)
# Open file device (e.g. PNG)
png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "/pr_or_outliers_up_venn_diagram.png"), width = 800, height = 800, res = 150)

# Draw the plot
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross,
  category = c("Protein\nunderexpression", "Gene underexpression"),
  fill = c("skyblue", "tomato"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  ind = FALSE
)

# Render and close the file
grid.draw(venn.plot)
dev.off()




underexpression_intersect <- intersect(or_res_aberrant[sampleID %in% pr_res_aberrant$pid & zScore < 0, gene_sample], pr_res_aberrant[PROTEIN_ZSCORE < 0, gene_sample])
only_pr <- setdiff(pr_res_aberrant[PROTEIN_ZSCORE < 0, gene_sample], underexpression_intersect)
only_or <- setdiff(or_res_aberrant[sampleID %in% pr_res_aberrant$pid & zScore < 0, gene_sample], underexpression_intersect)

area1 <- length(unique(only_pr))
area2 <- length(unique(only_or))
cross <- length(underexpression_intersect)
# Open file device (e.g. PNG)
png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "/pr_or_outliers_un_venn_diagram_onlyPR_samples.png"), width = 800, height = 800, res = 150)

# Draw the plot
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross,
  category = c("Protein\nunderexpression", "Gene underexpression"),
  fill = c("skyblue", "tomato"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  ind = FALSE
)

# Render and close the file
grid.draw(venn.plot)
dev.off()





underexpression_intersect <- intersect(or_res_aberrant[sampleID %in% pr_res_aberrant$pid & zScore > 0, gene_sample], pr_res_aberrant[PROTEIN_ZSCORE > 0, gene_sample])
only_pr <- setdiff(pr_res_aberrant[PROTEIN_ZSCORE > 0, gene_sample], underexpression_intersect)
only_or <- setdiff(or_res_aberrant[sampleID %in% pr_res_aberrant$pid & zScore > 0, gene_sample], underexpression_intersect)

area1 <- length(unique(only_pr))
area2 <- length(unique(only_or))
cross <- length(underexpression_intersect)
# Open file device (e.g. PNG)
png(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202507/", pr_output_name, "/pr_or_outliers_up_venn_diagram_onlyPR_samples.png"), width = 800, height = 800, res = 150)

# Draw the plot
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross,
  category = c("Protein\nunderexpression", "Gene underexpression"),
  fill = c("skyblue", "tomato"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  ind = FALSE
)

# Render and close the file
grid.draw(venn.plot)
dev.off()



#### activation outliers ####
mu <- 5
theta <-100
activation_res <- fread(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/nb_act/activation_res_protein_coding_genes_sig_theta", theta, "_mu", mu, "_BH.tsv"))



rna_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")
activation_res <- left_join(activation_res,rna_sa[, c("pid", "Oncotree Code", "SEX", "ICD10 Code", "TISSUE", "Tumorzellgehalt (Bioinformatik)")], by=c("sampleID" = "pid"))

activation_res <- left_join(activation_res, gene_annot_dt, by=c("geneID" = "gene_id"))
setnames(activation_res, "geneID", "geneID_long")
setnames(activation_res, "gene_name", "geneID")
activation_res[, Diag := DROP_GROUP]
activation_res[, aberrant := TRUE]


plot_data <- create_fisher_dt(activation_res, "Activation", cgc_oncogene$gene_name, omics_type="OUTRIDER")
odds_ratio_matrix <- prepare_plot_data(plot_data)

idx <- which(odds_ratio_matrix > 1, arr.ind = TRUE)

# Build the list of (row_name, col_name, value)
or_dt <- data.table(
  row_name = rownames(odds_ratio_matrix)[idx[, "row"]],
  col_name = colnames(odds_ratio_matrix)[idx[, "col"]],
  odds_pr_overexpression = odds_ratio_matrix[idx]
)
or_dt[, gene_cancer := paste0(row_name, "_", col_name)]


