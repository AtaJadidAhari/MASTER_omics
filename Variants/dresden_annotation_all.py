"""
Main file to create the variants.
starts reading VEP runs from FRANK: /omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/
Merges them with AbExp and promoterAI, AbSplice, 
Merges them with Annotations from MASTER (/omics/odcf/analysis/hipo/hipo_021/GermlineAnalysis/data_object_master_germline/sequencing/) to get ACMG classifications and CADD --> wtf does the varianst from Franks don't have them ... --> then why am I doing one thing twice ...



"""


import pandas as pd
import polars as pl
import os
import io
import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
import itables
from itables import init_notebook_mode
import polars as pl
init_notebook_mode(all_interactive=True)
import itables.options as opt
opt.lengthMenu = [10, 20, 50, 100]
opt.pageLength = 20   # in newer versions; otherwise use show(..., pageLength=20)

from scipy.stats import mannwhitneyu



import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from utils.util_functions import *
from ProteinExpression.load_pr_data import *


import math


sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")


cpt1 = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/CPT1_hg19_vars_all_protein_coding_genes.parquet")
cpt1 = cpt1[list(cpt1.columns[8:]) + ["REF", "ALT"]]

abexp = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/damaging_max_abexp.parquet").select(["chrom", "start", "end", "ref", "alt", "gene", "abexp_v1.1"]).collect(engine="streaming").to_pandas()
abexp["seqnames"] = abexp["chrom"].str.split("chr").str[1]
abexp = abexp.drop(columns=["chrom"])


promoter_ai = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/promoterAI_tss500_hg19.parquet")
                .filter(abs(pl.col("promoterAI")) >= 0.1)
                .collect(engine="streaming").to_pandas()
              )
df_sorted = promoter_ai.sort_values(
    by="promoterAI",
    key=lambda x: x.abs(),
    ascending=False
)
promoter_ai = df_sorted.drop_duplicates(
    subset=["chrom", "pos", "ref", "alt", "gene_id", "strand"],
    keep="first"
)

promoter_ai = promoter_ai[["chrom", "hg19_start", "ref", "alt", "gene_id", "promoterAI"]]


promoter_ai["chrom"] = promoter_ai["chrom"].str[3:]

promoter_ai['hg19_start'] = promoter_ai['hg19_start'].astype("Int64")
abexp['end'] = abexp['end'].astype("Int64")





absplice_predisp = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/all_scores_gt0_05_max_abSplice_snvs.parquet")
        .select(["chrom", "hg19_end", "ref", "alt", "gene_id", "pangolin_score", "AbSplice2_max"])
        .filter((pl.col("pangolin_score") >= 0.1) | (pl.col("AbSplice2_max") >= 0.1))
        .collect(engine="streaming")
      ).to_pandas()
absplice_predisp.shape


import re
import gzip

def get_csq_headers(vcf_path):
    # Compile a regex to look for the CSQ INFO definition and capture the format string
    csq_pattern = re.compile(r'##INFO=<ID=CSQ,.*Description=".*Format:\s*([^"]+)"')
    open_func = gzip.open if vcf_path.endswith(".gz") else open

    with open_func(vcf_path, "rt", encoding="utf-8") as f:
        for line in f:
            # Stop searching if we accidentally hit the actual variants table
            if not line.startswith("#"):
                break

            match = csq_pattern.search(line)
            if match:
                # Extract the pipe-separated string and split it into a list
                format_string = match.group(1)
                return [field.strip() for field in format_string.split("|")]

    raise ValueError("Could not find the ##INFO=<ID=CSQ...> field in the VCF header.")



# variant_type = "germline"
variant_type = "somatic"

vcf_base_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/0_germline_vcfs"
csq_headers = get_csq_headers("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/0_germline_vcfs/RUQHG9_germline.vcf.gz")


if variant_type == "somatic":
    vcf_base_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/1_somatic_vcfs"
    csq_headers = get_csq_headers("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/1_somatic_vcfs/7QXG3R_somatic.vcf.gz")



chunks = 700
current_chunk = int(sys.argv[1]) # 0-indexed, up to chunks - 1
print("analysiing:", current_chunk)

files = os.listdir(vcf_base_path)
total_files = len(files)

# Calculate chunk size rounding UP so no files are dropped
chunk_length = math.ceil(total_files / chunks)



current_chunk_start = current_chunk * chunk_length
# Cap the end index at the total number of files
current_chunk_end = min((current_chunk + 1) * chunk_length, total_files)

# Slice the list for the current chunk
chunk_files = files[current_chunk_start:current_chunk_end]

print(len(chunk_files))


all_dfs = []

for vep_res in chunk_files:
    sample_id = vep_res.split("_")[0]
    if not sample_id in sa["pid"].values:
        continue
    # 1. Read the VCF file, skipping the metadata lines (starting with ##)
    # We find where the column header ('#CHROM') starts.
    vcf_path = vcf_base_path + vep_res
    
    df = pl.read_csv(
            f"{vcf_base_path}/{vep_res}",
            separator="\t",
            comment_prefix="##",
            has_header=False,
        )
    df.columns = df.row(0)
    df = df.slice(1).to_pandas()
    
    # 2. A clean Python function to turn an INFO string into a flat dictionary
    def parse_info_to_dict(info_str):
        if not info_str or info_str == ".":
            return {}
    
        parsed = {}
        for item in info_str.split(";"):
            if not item:
                continue
            if "=" in item:
                k, v = item.split("=", 1)
                parsed[k] = v
            else:
                # It's a flag (e.g., ';VALIDATED;') -> save it as True
                parsed[item] = True
        return parsed
    
    
    # 3. Map the function and create the new columns dataframe
    info_columns_df = pd.DataFrame(df["INFO"].apply(parse_info_to_dict).tolist())
    
    # 4. Drop the old INFO column and horizontally concatenate the new ones
    pandas_df = pd.concat([df.drop(columns=["INFO"]), info_columns_df], axis=1)
    
    pandas_df = pandas_df.assign(CSQ=pandas_df["CSQ"].str.split(",")).explode("CSQ")
    
    # 6. Split the pipe (|) fields into individual columns using your csq_headers
    # We fill missing or empty CSQ values with empty strings to avoid errors during split
    csq_split_df = pd.DataFrame(
        pandas_df["CSQ"].fillna("").str.split("|").tolist(),
        index=pandas_df.index,  # Keeps alignment with the exploded rows
    )
    
    # If the VCF has trailing pipes or extra fields, truncate or pad to match csq_headers length
    csq_split_df = csq_split_df.iloc[:, : len(csq_headers)]
    csq_split_df.columns = csq_headers[: len(csq_split_df.columns)]
    
    # 7. Drop the temporary CSQ string column and merge the newly created fields
    final_df = pd.concat([pandas_df.drop(columns=["CSQ"]), csq_split_df], axis=1)

    cols_to_fix = [
    "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_AG", 
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"
    ]

    for col in cols_to_fix:
        final_df[col] = pd.to_numeric(final_df[col], errors='coerce')


    final_df["max_spliceai_score"] = final_df[["SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_AG", 
        "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"
    ]].max(axis=1)

    
    final_df = final_df.rename(columns={final_df.columns[8]: "sample_blood", final_df.columns[9]: "sample_plasma"})
    
    res = final_df.copy()
    res['POS'] = res['POS'].astype("Int64")
    

    # merge with promoterAI
    merged_vars = res.merge(promoter_ai, left_on=["#CHROM", "POS", "REF", "ALT", "Gene"], right_on=["chrom", "hg19_start", "ref", "alt", "gene_id"], how="left").drop(columns=[ "chrom", "hg19_start", "ref", "alt", "gene_id"])      
    
    # merge with AbExp
    merged_vars = merged_vars.merge(abexp, left_on=["POS", "#CHROM", "REF", "ALT", "Gene"], right_on=["end", "seqnames", "ref", "alt", "gene"], how="left").drop(columns=["end", "seqnames", "ref", "alt", "gene"])

    merged_vars["#CHROM"] = "chr" + merged_vars["#CHROM"] 

    # merge with spliceAI
    merged_vars = merged_vars.merge(absplice_predisp, left_on=["#CHROM", "POS", "REF", "ALT", "Gene"], right_on=["chrom", "hg19_end", "ref", "alt", "gene_id"], how="left").drop(columns=["chrom", "hg19_end", "ref", "alt", "gene_id"])
    
    # merge with cpt1
    merged_vars = merged_vars.merge(cpt1, left_on=["#CHROM", "POS", "REF", "ALT", "Gene"], right_on=["hg19_chrom", "hg19_pos", "REF", "ALT", "geneID_short"], how="left").drop(columns=["hg19_chrom", "hg19_pos", "geneID_short"])

    
    all_dfs.append(merged_vars[(merged_vars["IMPACT"].str.contains("HIGH")) |  (merged_vars["Consequence"].str.contains("missense")) | (merged_vars["max_spliceai_score"] >= 0.1) | (merged_vars["AbSplice2_max"] >= 0.1) | (merged_vars["pangolin_score"] >= 0.1) | 
         (abs(merged_vars["abexp_v1.1"]) >= 0.1) |  (abs(merged_vars["promoterAI"]) >= 0.1)])

merged_vars = pd.concat(all_dfs)



merged_vars = merged_vars.rename(columns={"#Uploaded_variation": f"{variant_type}_#Uploaded_variation", "Consequence": f"{variant_type}_Consequence", "IMPACT": f"{variant_type}_IMPACT", "am_pathogenicity": f"{variant_type}_am_pathogenicity","LoF": f"{variant_type}_am_LoF", "max_spliceai_score": f"{variant_type}_max_spliceai_score", "pangolin_score": f"{variant_type}_pangolin_score", "AbSplice2_max": f"{variant_type}_AbSplice2_max", "promoterAI": f"{variant_type}_promoterAI", "abexp_v1.1": f"{variant_type}_abexp_v1.1", "CPT1_score": f"{variant_type}_CPT1_score" })



merged_vars.to_parquet(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/{variant_type}/aggregated_chunks/chunk{current_chunk}_{variant_type}_damaging_aggregated_all.parquet", index=None)


dresden_germline_all = merged_vars.merge(sa[["pid", "nct_pid"]], how="left", left_on="PID", right_on="pid")

annotated_vars = (pl.scan_parquet(
                  "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA.parquet")
                  .filter((pl.col("POS").is_in(dresden_germline_all["POS"])))
                  .collect(engine="streaming")
                 )



annotated_vars = annotated_vars.to_pandas().drop_duplicates(subset = ["#CHROM", "POS", "REF", "ALT", "HUGO_Symbol"])

annotated_vars["vcf_path"] = "."

print(dresden_germline_all.columns)
pids = dresden_germline_all["nct_pid"].unique()



all_dfs = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/l_lp_vars.parquet")


all_dfs = pd.concat((annotated_vars[["#CHROM", "POS", "REF", "ALT", "ACMG_Classification","CharGer_Classification", "ClinVar_Pathogenicity", "HUGO_Symbol", "CADD_PHRED", "vcf_path"]], all_dfs[["#CHROM", "POS", "REF", "ALT", "ACMG_Classification","CharGer_Classification", "ClinVar_Pathogenicity"]], annotated_vars[["#CHROM", "POS", "REF", "ALT", "ACMG_Classification","CharGer_Classification", "ClinVar_Pathogenicity", "HUGO_Symbol", "CADD_PHRED", "vcf_path"]]))

all_dfs = all_dfs.drop_duplicates(subset = ["#CHROM", "POS", "REF", "ALT", "HUGO_Symbol"])




dresden_germline_all["#CHROM"] = dresden_germline_all["#CHROM"].str.split("chr").str[1]

dresden_germline_all = dresden_germline_all.merge(all_dfs[["#CHROM", "POS", "REF", "ALT", "ACMG_Classification","CharGer_Classification", "ClinVar_Pathogenicity", "HUGO_Symbol", "CADD_PHRED"]], how="left", left_on=["#CHROM", "POS", "REF", "ALT", "SYMBOL"], right_on = ["#CHROM", "POS", "REF", "ALT", "HUGO_Symbol"]).rename(columns={"ACMG_Classification": f"{variant_type}_ACMG_Classification", "CharGer_Classification": f"{variant_type}_CharGer_Classification", "ClinVar_Pathogenicity": f"{variant_type}_ClinVar_Pathogenicity", "CADD_PHRED": f"{variant_type}_CADD_PHRED"})





dresden_germline_all.to_parquet(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/{variant_type}/aggregated_chunks/chunk{current_chunk}_{variant_type}_aggregated_damaging_acmg.parquet", index=None)
