library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

source('~/Scripts/Preprocess/VEP_SORanking.R')



snv_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_HIGH_aggregated/vep_res_rare_snv_all_aggregated.tsv")
# indel_vep_res <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_HIGH_aggregated/vep_res_rare_indel_all_aggregated.tsv")


indel_vep_res <- c()
vep_res_combinded <- rbind(snv_vep_res, indel_vep_res)

print(nrow(vep_res_combinded))
vep_res_combinded <- vep_res_combinded[, IMPACT := factor(IMPACT, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"), ordered = TRUE)]

vep_res_combinded[, Consequence_most_severe := sapply(strsplit(Consequence, ","), function(x) {
  x <- trimws(x)  # remove extra spaces
  ranked <- x[order(severity_rank[x])]
  ranked[1]  # most severe
})]
setorder(vep_res_combinded, Gene, sampleID, IMPACT)

vep_res_combinded <- unique(vep_res_combinded, by = c("Gene", "sampleID"))

vep_res_combinded[Consequence_most_severe == "missense_variant" & am_class != "-", Consequence_most_severe := am_class]

print(nrow(vep_res_combinded))
fwrite(vep_res_combinded, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_HIGH_aggregated/vep_res_rare_snv_all_aggregated_unique_2.tsv", sep="\t")
