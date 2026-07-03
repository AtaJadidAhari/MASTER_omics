import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode
import polars as pl
import os
import glob
from mizani.formatters import custom_format

init_notebook_mode(all_interactive=True)
from scipy.stats import mannwhitneyu


import sys
sys.path.append("/home/a379i/Scripts/")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *
from utils.util_functions import *



germline_snvs  = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA.parquet").collect(engine="streaming").to_pandas()

germline_vep_res = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA_vep.parquet",).collect(engine="streaming").to_pandas()

germline_snvs[["REF_trimmed", "ALT_trimmed"]] = germline_snvs.apply(left_align_trim_vcf, axis=1)

germline_snvs.loc[(germline_snvs["REF"].str.len() > 1) | (germline_snvs["ALT"].str.len() > 1), "ALT"] = germline_snvs.loc[(germline_snvs["REF"].str.len() > 1) | (germline_snvs["ALT"].str.len() > 1), "ALT_trimmed"]
germline_snvs.loc[(germline_snvs["REF"].str.len() > 1) | (germline_snvs["ALT"].str.len() > 1), "REF"] = germline_snvs.loc[(germline_snvs["REF"].str.len() > 1) | (germline_snvs["ALT"].str.len() > 1), "REF_trimmed"]


germline_vep_res["chrom"] = germline_vep_res["Location"].str.split(":").str[0]
germline_vep_res["POS"] = 0
germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] == "SNV", "POS"] = germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] == "SNV", "Location"].str.split(":").str[1].astype(int)
germline_vep_res["POS"] = germline_vep_res["POS"] + 1

germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] != "SNV", "POS"] = germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] != "SNV", "#Uploaded_variation"].str.split(":").str[1].astype(int)



germline_snvs = germline_snvs.merge(germline_vep_res, left_on=["#CHROM",  "POS", "REF", "ALT", "HUGO_Symbol"], right_on=["chrom", "POS", "REF_ALLELE", "Allele", "SYMBOL"], how="left") 

germline_snvs.drop(columns=["INFO", "__index_level_0___x", "REF_trimmed", "ALT_trimmed", "__index_level_0___y", "chrom"]).to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA_merged_vep.parquet", index=None)