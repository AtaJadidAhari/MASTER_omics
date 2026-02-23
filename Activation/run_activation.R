library(reshape2)
library(dplyr)
library(data.table)
library(OUTRIDER)

source("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/envs/git_packages/gene_activation/scripts/function/nb_act.R")
source("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/envs/git_packages/gene_activation/scripts/function/help_functions.R")


gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))

activation_res_all <- list()
rarely_expressed_genes <- list()

dorp_groups <- c("ACC","Biliary_Tract",  "Bone","Bowel", "Breast", "COAD", "Esophagus_Stomach", "Pancreas", "CNS_Brain", 
                 "Head_and_Neck", "LMS", "Lung", "Other", "Skin", "Soft_Tissue", "SYNS")

protein_coding_genes <- gene_annot_dt[gene_type == "protein_coding"]

mu <- 1
theta <- 100

for (drop_group in dorp_groups){
  ods_unfitted <- readRDS(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/", drop_group, "/ods_unfitted.Rds"))
  pr_coding_ods <- ods_unfitted[rownames(ods_unfitted) %in% protein_coding_genes$gene_id, ]
  pr_coding_ods <- estimateSizeFactors(pr_coding_ods)
  
  rare_exp <- filterRarely(pr_coding_ods)
  rarely_expressed_genes <- union(rarely_expressed_genes, rare_exp$rarely_exp_genes)
  
  # Run NB-act on the subset of rarely expressed genes and adjust pval for multiple testing
  results <- nb_act(pr_coding_ods, rare_exp$rarely_exp_rpkm, adj = "BH", threshold = mu, theta = theta)
  
  
  res_padjust <- as.data.table(results$padj)
  res_padjust[, "geneID" := rownames(results$padj)]
  
  
  res_padjust <- as.data.table(melt(res_padjust, id.vars="geneID"))
  setnames(res_padjust, "variable", "sampleID")
  setnames(res_padjust, "value", "padjust")
  res_padjust[, DROP_GROUP := drop_group]
  res_padjust <- res_padjust[padjust < 0.05]
  activation_res_all[[drop_group]] <- res_padjust
  
}

activation_res_all <- rbindlist(activation_res_all)
rarely_expressed_genes <- unlist(rarely_expressed_genes)

fwrite(activation_res_all, paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/nb_act/activation_res_protein_coding_genes_sig_theta", theta, "_mu", mu, "_BH.tsv"), sep="\t")  
fwrite(as.data.table(rarely_expressed_genes), paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/nb_act/rarely_expressed_protein_coding_genes_", theta, "_mu", mu, "_BH.tsv"), sep="\t")
