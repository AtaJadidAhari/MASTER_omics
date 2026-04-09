import pandas as pd

import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

cnv_germline = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv.tsv", sep="\t")
cnv = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv.tsv", sep="\t")
cnv_germline.loc[:, 'Gene'] = cnv_germline['Gene'].str.split(',')
cnv.loc[:, 'Gene'] = cnv['Gene'].str.split(',')

cnv_germline = cnv_germline.explode('Gene').reset_index(drop=True)
cnv = cnv.explode('Gene').reset_index(drop=True).drop_duplicates(subset=["group_name", "Gene"])
cnv_germline = cnv_germline.explode('Gene').reset_index(drop=True).drop_duplicates(subset=["group_name", "Gene"])

cnv_germline["sampleID"] = cnv_germline["group_name"].str.split(".").str[1]
cnv["sampleID"] = cnv["group_name"].str.split(".").str[1]



cnv.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv_exploded.tsv", sep="\t", index=False)
cnv_germline.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv_exploded.tsv", sep="\t", index=False)


cnv_predisp = cnv[cnv["Gene"].isin(extended_dresden_dt["gene_name"])]
cnv_germline_predisp = cnv_germline[cnv_germline["Gene"].isin(extended_dresden_dt["gene_name"])]


cnv_predisp.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/cnv_predisp_exploded.tsv", sep="\t", index=False)
cnv_germline_predisp.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_cnv_predisp_exploded.tsv", sep="\t", index=False)
