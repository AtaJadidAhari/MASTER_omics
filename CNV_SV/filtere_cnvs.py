import pandas as pd

import pandas as pd
from itables import init_notebook_mode
import polars as pl
import os
import numpy as np
init_notebook_mode(all_interactive=True)
import plotnine as pn



import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *


sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")
sample_id = "2DWXNH"
sa[sa["pid"] == sample_id]

sa["key"] = sa["pid"] +  "." + sa["nct_pid"].str.split("_").str[1]



cnv_germline = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv_exploded.tsv", sep="\t")
cnv = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv_exploded.tsv", sep="\t")
cnv = cnv[(cnv["Type"].notna()) & (cnv["Type"] != "CNN")]




cnv["key"] = cnv["group_name"].str[4:]
cnv_germline["key"] = cnv_germline["group_name"].str[4:]

cnv = cnv.merge(sa[["key"]], on="key", how="inner")
cnv_germline = cnv_germline.merge(sa[["key"]], on="key", how="inner")


cnv.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv_exploded_filtered.tsv", sep="\t", index=None)
cnv_germline.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv_exploded_filtered.tsv", sep="\t", index=None)
