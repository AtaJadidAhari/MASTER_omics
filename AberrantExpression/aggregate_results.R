
or_res_dir <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider"
or_res_files <- list.dirs(or_res_dir, full.names = TRUE, recursive = FALSE)

or_res_agg_all <- c()
or_res_agg <- c()
zscores_res_agg_all <- c()

for (drop_group in or_res_files){
  
  drop_group_name <- tail(unlist(strsplit(drop_group, "/")), 1)
  
  if (drop_group_name == "Unstranded_data"){
    next
  }

  # or_res_cohort <- fread(paste0(drop_group, "/OUTRIDER_results.tsv"))
  # or_res_cohort[, Diag := drop_group_name]
  # or_res_agg <- rbind(or_res_agg, or_res_cohort)
  
  # or_res_cohort_all <- fread(paste0(drop_group, "/OUTRIDER_results_all.tsv"))
  # or_res_cohort_all[, Diag := drop_group_name]
  # or_res_agg_all <- rbind(or_res_agg_all, or_res_cohort_all)
  
  zscores_res_agg <-  fread(paste0(drop_group, "/zScores.tsv"))
  zscores_res_agg_all <- rbind(zscores_res_agg, zscores_res_agg_all)
    
}
  
fwrite(or_res_agg, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_outliers.tsv", sep="\t")
fwrite(or_res_agg_all, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_or_res_all.tsv", sep="\t")
fwrite(zscores_res_agg_all, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/aggregated_sf_zScores_all.tsv", sep="\t")





