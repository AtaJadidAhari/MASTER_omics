import pandas as pd
import polars as pl
import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *
import sys


# Concatenate with a separator only if both exist
def join_types(row):
    parts = [p for p in [row["Type_x"], row["Type_y"]] if p != ""]
    return " | ".join(parts) if parts else "No CNV"



cnv_germline = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv_exploded_filtered.tsv", sep="\t")
cnv_germline = cnv_germline[cnv_germline["Confidence"] == "HIGH"]

cnv = ( pl.scan_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv_exploded_filtered.tsv", separator="\t", infer_schema_length=1000000)
        .filter(pl.col("Type") != "CNN")
        .collect(engine="streaming")
      ).to_pandas()

cnv = cnv[((cnv["Type"] == "AMP") | (cnv["Type"] == "HDEL")) & (cnv["width"] >= 10000)]


cnv = cnv.merge(gene_annot_dt[["geneID_short", "gene_name"]], left_on="Gene", right_on="gene_name")
cnv_germline = cnv_germline.merge(gene_annot_dt[["geneID_short", "gene_name"]], left_on="Gene", right_on="gene_name")

cnv = cnv.drop_duplicates(subset=["sampleID", "geneID_short"])
cnv_germline = cnv_germline.drop_duplicates(subset=["sampleID", "geneID_short"])



# py_or_res_all = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust.parquet")
py_or_res_all = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/zscores/or_variants.csv")

cols = py_or_res_all.columns
final_cols = cols.to_list() + ["CNV"]
print(len(py_or_res_all))
# # Replace NaN with empty strings before concatenating


merged_df = pd.merge(
    py_or_res_all, 
    cnv, 
    on=['geneID_short', 'sampleID'], 
    how='left'
)
merged_df = pd.merge(
    merged_df, 
    cnv_germline, 
    on=['geneID_short', 'sampleID'], 
    how='left'
)
merged_df[ (merged_df["zScore"] > 0)]["Type_x"].value_counts()




merged_df["Type_x"] = [f"Somatic_{x}" if pd.notna(x) and x != "" else "" for x in merged_df["Type_x"]]

# Clean Germline labels
merged_df["Type_y"] = [f"Germline_{y}" if pd.notna(y) and y != "" else "" for y in merged_df["Type_y"]]



merged_df["CNV"] = merged_df.apply(join_types, axis=1)



merged_df[final_cols].to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/zscores/or_variants_cnv.parquet", index=False)

# merged_df[final_cols].to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust_cnv.parquet", index=False)


print(len(merged_df))
sys.exit()


# pr_output_name = "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"
# # pr_output_name = "zScore_gt3"

# pr_res_all = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants_predisppadjust.parquet")
# # pr_res_all = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants.csv")

# pr_cols = pr_res_all.columns.to_list()
# pr_final_cols = pr_cols + ["CNV"]

# print(len(pr_res_all))

# pr_merged_df = pd.merge(
#     pr_res_all, 
#     cnv, 
#     on=['geneID_short', 'sampleID'],
#     how='left'
# )

# pr_merged_df = pd.merge(
#     pr_merged_df, 
#     cnv_germline, 
#     on=['geneID_short', 'sampleID'], 
#     how='left'
# )
# # Replace NaN with empty strings before concatenating


# pr_merged_df["Type_x"] = [f"Somatic_{x}" if pd.notna(x) and x != "" else "" for x in pr_merged_df["Type_x"]]

# # Clean Germline labels
# pr_merged_df["Type_y"] = [f"Germline_{y}" if pd.notna(y) and y != "" else "" for y in pr_merged_df["Type_y"]]


# pr_merged_df["CNV"] = pr_merged_df.apply(join_types, axis=1)
# pr_merged_df[pr_final_cols].to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants_predisppadjust_cnv.parquet", index=False)

# print(len(pr_merged_df))




