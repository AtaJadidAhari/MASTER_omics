import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode
import polars as pl
import os
import polars.selectors as cs

import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *
from utils.util_functions import *

# Concatenate with a separator only if both exist
def join_types(row):
    parts = [p for p in [row["Type_x"], row["Type_y"]] if p != ""]
    return " | ".join(parts) if parts else "No CNV"



cnv_germline = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv_exploded_filtered.tsv", sep="\t")
cnv_germline = cnv_germline[cnv_germline["Confidence"] == "HIGH"]
cnv_germline = cnv_germline[cnv_germline["sampleID"].isin(sa["pid"])]

cnv_germline = cnv_germline.merge(gene_annot_dt[["geneID_short", "gene_name"]], left_on="Gene", right_on="gene_name")
cnv_germline = cnv_germline.drop_duplicates(subset=["sampleID", "geneID_short"])

print("cnv_germline shape", cnv_germline.shape)

cnv = ( pl.scan_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv_exploded_filtered.tsv", separator="\t", infer_schema_length=1000000)
        .select(["sampleID", "Gene", "Type", "width"])
        .with_columns(
        pl.col("sampleID").cast(pl.String)
        )
        .filter(pl.col("Type") != "CNN").filter(pl.col("sampleID").is_in(sa["pid"]))
        .collect(engine="streaming")
      ).to_pandas()

cnv = cnv[((cnv["Type"] == "AMP") | (cnv["Type"] == "HDEL")) & (cnv["width"] >= 10000)]


cnv = cnv.merge(gene_annot_dt[["geneID_short", "gene_name"]], left_on="Gene", right_on="gene_name")

cnv = cnv.drop_duplicates(subset=["sampleID", "geneID_short"])

print("cnv shape", cnv.shape)


drop_groups = list(sa[sa["Diag"] != "Unstranded_data"]["Diag"].unique())
drop_run = "python"

for group in drop_groups:
    print(group, "here")
    # python_outrider
    if drop_run == "python":
        or_res_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/cohort_wise/" + group + "/oht_lr_0_0001_epoc200_cpu/protrider_summary.csv"
        print(or_res_path)
        cohort_wise_res = pd.read_csv(or_res_path)
        cohort_wise_res = cohort_wise_res.rename(columns={"PROTEIN_PVALUE": "pValue", "PROTEIN_ZSCORE": "zScore", "PROTEIN_PADJ": "padjust", "proteinID": "geneID"})
        output_path = f"{or_res_path.replace('protrider_summary.csv', 'or_res_cnv.parquet')}"

    else:
        or_res_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/" + group + "/OUTRIDER_results_all.tsv"
        cohort_wise_res = pd.read_csv(or_res_path, sep="\t")
        output_path = f"{or_res_path.replace('.tsv', '_cnv.parquet')}"

    cohort_wise_res = cohort_wise_res.sort_values("pValue")

    cohort_wise_res = pd.merge(cohort_wise_res, gene_annot_dt, left_on="geneID", right_on="gene_id", how="left")

    
    cohort_wise_res = cohort_wise_res.drop_duplicates(subset=["sampleID", "geneID_short"])

    cohort_wise_res = pd.merge(
        cohort_wise_res, 
        cnv[["geneID_short", "sampleID", "Type"]], 
        on=['geneID_short', 'sampleID'],
        how='left'
    )

    cohort_wise_res = pd.merge(
        cohort_wise_res, 
        cnv_germline[["geneID_short", "sampleID", "Type"]], 
        on=['geneID_short', 'sampleID'], 
        how='left'
    )
    # Replace NaN with empty strings before concatenating

    cohort_wise_res["Type_x"] = [f"Somatic_{x}" if pd.notna(x) and x != "" else "" for x in cohort_wise_res["Type_x"]]

    # Clean Germline labels
    cohort_wise_res["Type_y"] = [f"Germline_{y}" if pd.notna(y) and y != "" else "" for y in cohort_wise_res["Type_y"]]


    cohort_wise_res["CNV"] = cohort_wise_res.apply(join_types, axis=1)

    cohort_wise_res.to_parquet(output_path, index=False)
