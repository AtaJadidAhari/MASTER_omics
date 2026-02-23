library(stringr)
# merge with sa to get diags
protein_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protrider_sa.tsv")

load_pr_data <- function(pr_output){
  pr_res <- fread(pr_output)
  # merge with sa to get diags
  pr_res <- left_join(pr_res, protein_sa, by=c("sampleID"="pid"))
  # if (!"full_name" %in% colnames(pr_res)){
  #   pr_res <- pr_res %>%
  #     mutate(full_name = sampleID) %>% 
  #     mutate(sampleID = sapply(str_split(sampleID, "-"), `[`, 2))
  # }
  if (!"geneID_short" %in% colnames(pr_res)){
    if (!grepl("ENSG",pr_res[, proteinID][1])){
      gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
      gene_annot_dt <- gene_annot_dt %>%
        mutate(geneID_short = sub("\\..*", "", gene_id))
      pr_res <- left_join(pr_res, gene_annot_dt, by=c("proteinID"= "gene_name"))
    } else{
      pr_res[, "geneID_short":= proteinID]
    }
  }

  pr_res[, "geneID":= geneID_short]
  
  # drop na proteins --> TODO: we should check to see why there is NAs
  pr_res <- pr_res[geneID != ""]
  return (pr_res)
  
}


