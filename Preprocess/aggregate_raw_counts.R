or_raw_counts_path <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/exported_counts/"

raw_counts <- list.dirs(or_raw_counts_path, full.names = TRUE, recursive = FALSE)

all_counts <- NULL
for (i in seq(1, length(raw_counts))){
  current_counts <- fread(paste0(raw_counts[i], "/geneCounts.tsv.gz"))
  if (i > 1){
    all_counts <- merge(all_counts, current_counts, by="geneID")
  }
  else{
    all_counts <- current_counts
  }
}



fwrite(all_counts, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aggregated_counts/geneCounts.tsv.gz", sep="\t")
