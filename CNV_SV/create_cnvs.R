
library(MultiAssayExperiment)
library(RaggedExperiment)
library(data.table)


drop_sa <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv")
drop_sa[, master_pid := gsub( "_", ".", nct_pid)]

load("/omics/odcf/analysis/OE0246_projects/datamaster/dataMASTER_release/object/dataMASTER_241104123708.RData")
load("/omics/odcf/analysis/OE0246_projects/datamaster/ext_data/annotations/gencode19_gns_lite.RData")

myobj = dataMASTER[,,c("cnv", "cnv_germline", "sv")] # class: MAE


# 1. Get all unique column names from the MEA
all_colnames <- unlist(colnames(myobj))

# 2. Initialize an empty character vector to store matches
matched_cols <- character()

# 3. Loop through each ID in your sa table
for (pid in drop_sa$master_pid) {
  # We use fixed = TRUE because your IDs have dots (.) 
  # which would otherwise be treated as wildcards in regex
  hits <- grep(pid, all_colnames, value = TRUE, fixed = TRUE)
  
  # Append new hits to our master list
  matched_cols <- c(matched_cols, hits)
}

# 4. Remove duplicates (in case one column matched multiple PIDs)
matched_cols <- unique(matched_cols)

cnv_germline  = getWithColData(myobj, "cnv_germline")
actual_cnv_cols <- colnames(cnv_germline)

# 2. Filter your matched_cols to only include those present in CNV
valid_cnv_matches <- intersect(matched_cols, actual_cnv_cols)

# 3. Now subset safely
cnv_germline <- cnv_germline[, valid_cnv_matches, ]
cnv_germline


cnv_germline@assays[[1]]$Gene = as.character(NA)
sv = as(cnv_germline, "GRangesList")
 
n = length(sv)
sv = lapply(1:n, function(i) {
  svPID = sv[[i]]
  # cat(sprintf("\r%s/%s", i, n))
  idx = findOverlaps(svPID, gencode19_gns_lite)
  gns = tapply(idx@to, idx@from, function(i) {
    paste(gencode19_gns_lite$gene_name[i], collapse=",")
  })
  svPID$Gene[as.numeric(names(gns))] = unname(gns)
  svPID
})
sv_re = RaggedExperiment(sv, colData=colData(cnv_germline))

myobj[["cnv_germline"]] = sv_re
germline_cnv_df <- as.data.frame(myobj[["cnv_germline"]])
setDT(germline_cnv_df)
fwrite(germline_cnv_df, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv.tsv", sep="\t")




cnv  = getWithColData(myobj, "cnv")
actual_cnv_cols <- colnames(cnv)

# 2. Filter your matched_cols to only include those present in CNV
valid_cnv_matches <- intersect(matched_cols, actual_cnv_cols)

# 3. Now subset safely
cnv <- cnv[, valid_cnv_matches, ]
cnv


cnv@assays[[1]]$Gene = as.character(NA)
sv = as(cnv, "GRangesList")
 
n = length(sv)
sv = lapply(1:n, function(i) {
  svPID = sv[[i]]
  # cat(sprintf("\r%s/%s", i, n))
  idx = findOverlaps(svPID, gencode19_gns_lite)
  gns = tapply(idx@to, idx@from, function(i) {
    paste(gencode19_gns_lite$gene_name[i], collapse=",")
  })
  svPID$Gene[as.numeric(names(gns))] = unname(gns)
  svPID
})
sv_re = RaggedExperiment(sv, colData=colData(cnv))

myobj[["cnv"]] = sv_re

cnv_df <- as.data.frame(myobj[["cnv"]])
setDT(cnv_df)

fwrite(cnv_df, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv.tsv", sep="\t")

