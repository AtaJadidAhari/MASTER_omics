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


# vcf_base_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/0_germline_vcfs"
vcf_base_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/1_somatic_vcfs"



sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")


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


# Use it on your file

csq_headers = get_csq_headers("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/1_somatic_vcfs/7QXG3R_somatic.vcf.gz")
# csq_headers = get_csq_headers("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/0_germline_vcfs/RUQHG9_germline.vcf.gz")

# csq_headers


import polars as pl
all_dfs = []

for vep_res in os.listdir(vcf_base_path):
    sample_id = vep_res.split("_")[0]
    # if not sample_id in sa["pid"].values:
    #     continue
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
    final_df = final_df.rename(columns={final_df.columns[8]: "sample_blood", final_df.columns[9]: "sample_plasma"})
    all_dfs.append(final_df[(final_df["Gene"].isin(genes_of_interest["geneID_short"])) | (final_df["SYMBOL"].isin(genes_of_interest["gene_name"]))])

res = pd.concat(all_dfs)


valid_positions = (
    res["POS"]
    .dropna()
    .unique()
    .astype(int)
    .tolist()
)


abexp = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/predisp_cgc_max_abexp.parquet").filter((pl.col("start").is_in(valid_positions)) | (pl.col("end").is_in(valid_positions))).collect(engine="streaming").to_pandas()

abexp["seqnames"] = abexp["chrom"].str.split("chr").str[1]


promoter_ai = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/promoterAI_tss500_hg19.parquet")
                .filter((pl.col("gene_id").is_in(genes_of_interest["geneID_short"])))
                .filter(pl.col("hg19_start").is_in(valid_positions))
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
promoter_ai["chrom"] = promoter_ai["chrom"].str[3:]


res['POS'] = res['POS'].astype(float)
promoter_ai['hg19_start'] = promoter_ai['hg19_start'].astype(float)
abexp['end'] = abexp['end'].astype(float)


merged_vars = res.merge(promoter_ai, left_on=["#CHROM", "POS", "REF", "ALT", "Gene"], right_on=["chrom", "hg19_start", "ref", "alt", "gene_id"], how="left").drop(columns=[ "chrom"])      

merged_vars = merged_vars.merge(abexp, left_on=["POS", "#CHROM", "REF", "ALT", "Gene"], right_on=["end", "seqnames", "ref", "alt", "gene"], how="left")


variant_type = "somatic"


merged_vars = merged_vars.rename(columns={"#Uploaded_variation": f"{variant_type}_#Uploaded_variation", "Consequence": f"{variant_type}_Consequence", "IMPACT": f"{variant_type}_IMPACT", "am_pathogenicity": f"{variant_type}_am_pathogenicity","LoF": f"{variant_type}_am_LoF", "max_spliceai_score": f"{variant_type}_max_spliceai_score", "pangolin_score": f"{variant_type}_pangolin_score", "AbSplice2_max": f"{variant_type}_AbSplice2_max", "promoterAI": f"{variant_type}_promoterAI", "abexp_v1.1": f"{variant_type}_abexp_v1.1" })



merged_vars["Variant"] = (
    merged_vars["seqnames"] + ":" + 
    merged_vars["POS"].astype(str) + ":" + 
    merged_vars["REF"] + ":" + 
    merged_vars["ALT"]
)


# merged_vars.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/germline_aggregated_predisp_cgc.parquet", index=None)

merged_vars.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/somatic_aggregated_predisp_cgc.parquet", index=None)
