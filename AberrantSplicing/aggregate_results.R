library(data.table)
library(stringr)
library(dplyr)


gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))


fr_res_dir <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/"
fr_res_files <- list.dirs(fr_res_dir, full.names = TRUE, recursive = FALSE)


# fr_res_agg_all <- c()
# fr_res_agg <- c()

# for (drop_group in fr_res_files){
  
#   drop_group_name <- tail(unlist(strsplit(drop_group, "/")), 1)
  
#   if (drop_group_name == "Unstranded_data"){
#     next
#   }
  
#   fr_res_cohort <- fread(paste0(drop_group, "/results.tsv"))
#   fr_res_cohort[, Diag := drop_group_name]
#   fr_res_agg <- rbind(fr_res_agg, fr_res_cohort, fill=T)
#   fr_res_cohort_all <- fread(paste0(drop_group, "/results_gene_all.tsv"))
#   fr_res_cohort_all[, Diag := drop_group_name]
#   fr_res_agg_all <- rbind(fr_res_agg_all, fr_res_cohort_all, fill=T)
  
  
# }

# fwrite(fr_res_agg, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_results.tsv", sep="\t")
# fwrite(fr_res_agg_all, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_results_gene_all.tsv", sep="\t")



fraser_res_dir <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/"
fr_res_files <- list.dirs(fraser_res_dir, full.names = TRUE, recursive = FALSE)

fr_res_junc <- NULL
for (x in fr_res_files) {
  drop_group_name <- tail(unlist(strsplit(x, "/")), 1)
  
  if (drop_group_name == "Unstranded_data"){
    next
  }
  fr_res_temp <- fread(paste0(x, "/results_per_junction.tsv"))
  
  fr_res_temp <- merge(fr_res_temp, gene_annot_dt[, .(gene_name, gene_type, gene_id, geneID_short)], 
                       by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  #fr_res_temp <- fr_res_temp[gene_type=='protein_coding', ]
  setnames(fr_res_temp, "gene_id", "geneID")
  
  if (exists("fr_res_junc")) {
    fr_res_junc <- rbind(fr_res_junc, fr_res_temp)
  } else {
    fr_res_junc <- fr_res_temp
  }
}
fr_res_junc[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]
fr_res_junc[, junction_id := paste0("Intron ", seqnames, ": ", 
                                    format(start, nsmall=1, big.mark=","), "-", 
                                    format(end, nsmall=1, big.mark=",")), 
            by = rownames(fr_res_junc)]


fwrite(fr_res_junc, paste0(fraser_res_dir, "/aggregated_results_per_junction.tsv"), sep="\t")







