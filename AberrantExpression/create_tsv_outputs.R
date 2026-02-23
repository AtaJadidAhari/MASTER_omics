library(data.table)
library(OUTRIDER)

or_res_dir <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider"
or_res_files <- list.dirs(or_res_dir, full.names = TRUE, recursive = FALSE)

or_res_agg <- c()
for (drop_group in or_res_files){
  
  drop_group_name <- tail(unlist(strsplit(drop_group, "/")), 1)
  
  or_res_cohort <- as.data.table(readRDS(paste0(drop_group, "/OUTRIDER_results_all.Rds")))
  or_res_cohort[, Diag := drop_group_name]
  
  fwrite(or_res_cohort, paste0(drop_group, "/OUTRIDER_results_all.tsv"), sep="\t")
  
}

