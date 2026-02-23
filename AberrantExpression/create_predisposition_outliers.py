import pandas as pd
from pathlib import Path
import numpy as np
import os
from statsmodels.stats.proportion import proportion_confint
from scipy.stats import beta


import sys
sys.path.append("/home/a379i/Scripts/")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *


sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")

py_or_res_all = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/protrider_summary.csv")

py_or_res_all = pd.merge(py_or_res_all, sa, left_on="sampleID", right_on="pid")
py_or_res_all["geneID_short"] = py_or_res_all["geneID"].str.split(".").str[0]
py_or_res_all["expressed_in_cohort"] = False
py_or_res_all["gene_sample"] = py_or_res_all["geneID"] + "_" + py_or_res_all["sampleID"]






variant_base_path = Path("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/")
benchmark_criteria = "all_aggregated_unique"

# Read tables
snv_vep_res = pd.read_csv(
    variant_base_path / f"vep_res_rare_snv_{benchmark_criteria}_variant_type.tsv",
    sep="\t"
)
indel_vep_res = pd.read_csv(
    variant_base_path / f"vep_res_rare_indel_{benchmark_criteria}_variant_type.tsv",
    sep="\t"
)

# Combine
vep_res_combined_all = pd.concat([snv_vep_res, indel_vep_res], ignore_index=True)

print(len(vep_res_combined_all))

# Unique_based_on_criteria
if benchmark_criteria == "CADD_PHRED":
    # Sort decreasing by CADD_PHRED
    vep_res_combined = (
        snv_vep_res.sort_values("CADD_PHRED", ascending=False)
        .groupby(["sampleID", "Gene"])
        .head(1)
        .reset_index(drop=True)
    )
# elif benchmark_criteria == "VEP_splice":
#   vep_res_combinded = vep_res_combinded[vep_res_combinded["Consequence"].str.contains("splic", na=False)]
#   vep_res_combinded <- vep_res_combinded[, IMPACT := factor(IMPACT, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"), ordered = TRUE)]
#   setorder(vep_res_combinded, Gene, sampleID, IMPACT)
else:
    # Set IMPACT as ordered category and sort
    impact_order = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
    vep_res_combined_all["IMPACT"] = pd.Categorical(
        vep_res_combined_all["IMPACT"],
        categories=impact_order,
        ordered=True
    )

    vep_res_combined_all = (
        vep_res_combined_all.sort_values(["Gene", "sampleID", "IMPACT"])
    )

# Remove duplicates by Gene + sampleID, keeping first (like data.table)
vep_res_combined_all = (
    vep_res_combined_all.drop_duplicates(subset=["Gene", "sampleID"])
)
print(len(vep_res_combined_all))




py_or_res_aberrant = py_or_res_all[py_or_res_all['padjust'] <= 0.05]
py_or_res_aberrant = pd.merge(py_or_res_aberrant, vep_res_combined_all, left_on=["geneID_short", "sampleID"], right_on= ["Gene", "sampleID"], how="left")
py_or_res_aberrant["IMPACT_bool"] = py_or_res_aberrant["IMPACT"] == "HIGH"

py_or_res_aberrant_predisposed_all = pd.merge(py_or_res_aberrant, dresden_dt_cgc[["gene_name", "gene_type", "geneID", "ROLE_IN_CANCER", "predisposition_gene"]], on="geneID", how="left")
py_or_res_aberrant_predisposed_all["expression_direction"] = "underexpression"
py_or_res_aberrant_predisposed_all.loc[py_or_res_aberrant_predisposed_all["zScore"] >= 0, "expression_direction"] = "Overexpression"
py_or_res_aberrant_predisposed_all[["sampleID", "geneID", "zScore", "pValue", "padjust", "Oncotree Code", "IMPACT", "Consequence", "gene_name", "ROLE_IN_CANCER", "expression_direction", "nct_pid", "predisposition_gene", "ANNOTATION_control", '#Uploaded_variation', "Location", '#CHROM', 'POS', 'ID', "indel_vcf", "snv_vcf"]].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/gene_expression_outliers.tsv", sep="\t", index=None)


py_or_res_outliers_predisposed = py_or_res_aberrant_predisposed_all[py_or_res_aberrant_predisposed_all["predisposition_gene"] == True]
py_or_res_outliers_predisposed[["sampleID", "geneID", "zScore", "pValue", "padjust", "Oncotree Code", "IMPACT", "Consequence", "gene_name", "ROLE_IN_CANCER", "expression_direction", "nct_pid", "predisposition_gene", "ANNOTATION_control", '#Uploaded_variation', "Location", '#CHROM', 'POS', 'ID',  "indel_vcf", "snv_vcf"]].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/gene_expression_predisposition_outliers.tsv", sep="\t", index=None)




#### Protein outliers

pr_output_name = "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"
pr_res_all = load_pr_data("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/protrider_summary.csv")


pr_res_all = pd.merge(pr_res_all, vep_res_combined_all, left_on=["geneID_short", "sampleID"], right_on= ["Gene", "sampleID"], how="left")
pr_res_all["IMPACT_bool"] = pr_res_all["IMPACT"] == "HIGH"

pr_res_aberrant = pr_res_all[pr_res_all["padjust"] <= 0.1]
pr_res_aberrant = pd.merge(pr_res_aberrant, sa, left_on="sampleID", right_on="pid", how="left")


pr_res_aberrant_predisposed_all = pd.merge(pr_res_aberrant, dresden_dt_cgc[["gene_name", "gene_type", "geneID_short", "ROLE_IN_CANCER", "predisposition_gene"]], right_on="geneID_short", left_on="geneID", how="left")
pr_res_aberrant_predisposed_all["expression_direction"] = "underexpression"
pr_res_aberrant_predisposed_all.loc[pr_res_aberrant_predisposed_all["zScore"] >= 0, "expression_direction"] = "Overexpression"
pr_res_aberrant_predisposed_all["geneID"] = pr_res_aberrant_predisposed_all["gene_id"]

pr_res_aberrant_predisposed_all[[ "nct_pid","proteinID", "geneID", "zScore", "pValue", "padjust", "Oncotree Code", "IMPACT", "Consequence", "ROLE_IN_CANCER", "expression_direction", "predisposition_gene", "sampleID", "ANNOTATION_control",'#Uploaded_variation', "Location", '#CHROM', 'POS', 'ID', "indel_vcf", "snv_vcf"]].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/protein_expression_outliers.tsv", sep="\t", index=None)


pr_res_aberrant_predisposed = pr_res_aberrant_predisposed_all[pr_res_aberrant_predisposed_all["predisposition_gene"] == True]
pr_res_aberrant_predisposed[[ "nct_pid","proteinID", "geneID", "zScore", "pValue", "padjust", "Oncotree Code", "IMPACT", "Consequence", "ROLE_IN_CANCER", "expression_direction", "predisposition_gene", "sampleID", "ANNOTATION_control", '#Uploaded_variation', "Location", '#CHROM', 'POS', 'ID', "indel_vcf", "snv_vcf"]].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/protein_expression_predisposition_outliers.tsv", sep="\t", index=None)




#### FRASER results
fr_res_gene = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_outliers.tsv", sep="\t")
fr_res_gene = pd.merge(fr_res_gene, sa, left_on="sampleID", right_on="pid")
fr_res_gene = fr_res_gene.merge(gene_annot_dt[["gene_name", "geneID_short", "gene_id"]], how="left", left_on="hgncSymbol", right_on="gene_name")

fr_res_gene = fr_res_gene.rename(columns={"gene_id":"geneID"})



fr_res_gene = pd.merge(fr_res_gene, vep_res_combined_all, left_on=["geneID_short", "sampleID"], right_on= ["Gene", "sampleID"], how="left")
fr_res_gene_aberrant_predisposed_all = pd.merge(fr_res_gene, dresden_dt_cgc[["gene_type", "geneID_short", "ROLE_IN_CANCER", "predisposition_gene"]], on="geneID_short", how="left")

cols = ["nct_pid", "Oncotree Code", "geneID"] + list(fr_res_gene_aberrant_predisposed_all.columns[:12]) +["annotatedJunction", "potentialImpact", "causesFrameshift", "IMPACT", "Consequence", "ROLE_IN_CANCER", "predisposition_gene", "RNA_BAM_FILE", '#Uploaded_variation', "Location", '#CHROM', 'POS', 'ID', 'ANNOTATION_control', "indel_vcf", "snv_vcf"]


fr_res_gene_aberrant_predisposed_all[cols].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/gene_splicing_outliers.tsv", sep="\t", index=None)
fr_res_gene_aberrant_predisposed = fr_res_gene_aberrant_predisposed_all[fr_res_gene_aberrant_predisposed_all["predisposition_gene"] == True]
fr_res_gene_aberrant_predisposed[cols].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/gene_splicing_predisposition_outliers.tsv", sep="\t", index=None)




### create shared protein-rna outliers

shared_outliers = pr_res_aberrant_predisposed_all.merge(py_or_res_aberrant_predisposed_all, on=["geneID", "sampleID", "nct_pid", "Oncotree Code", "IMPACT", "Consequence", "predisposition_gene", "ROLE_IN_CANCER", "ANNOTATION_control", "#Uploaded_variation"
                                        ,"Location", "#CHROM", "POS", "ID", "indel_vcf", "snv_vcf"]).drop("gene_name", axis=1)
shared_outliers = shared_outliers.rename(columns={"zScore_x": "protein_zScore", "pValue_x": "protein_pValue", "padjust_x": "protein_padjust", "expression_direction_x": "protein_expression_direction", 
                                                  "zScore_y": "RNA_zScore", "pValue_y": "RNA_pValue", "padjust_y": "RNA_padjust", "expression_direction_y": "RNA_expression_direction"})
cols = list(shared_outliers.columns[:7]) + [shared_outliers.columns[10]] + list(shared_outliers.columns[21:]) + list(shared_outliers.columns[7:10]) + list(shared_outliers.columns[11:21])
shared_outliers[cols].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/shared_outliers/shared_rna_protein_outliers.tsv", sep="\t")

shared_outliers[shared_outliers["predisposition_gene"].notna()].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/shared_outliers/shared_rna_protein_predisposition_outliers.tsv", sep="\t")



