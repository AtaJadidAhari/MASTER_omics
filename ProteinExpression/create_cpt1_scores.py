import pandas as pd
import polars as pl
from liftover import get_lifter


import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *
from Variants.add_missing_acmg import *
 

cpt1 = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/CPT1_all_protein_coding_genes.parquet").collect().to_pandas()
eve_db = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/grch38_popEVE_ukbb_20250715.vcf.gz", sep="\t")

df_parsed = (
    eve_db["INFO"]
    .str.split(";")
    .explode()
    .str.extract(r"(?P<key>[^=]+)=(?P<value>.*)")
    .reset_index()
    .pivot(index="index", columns="key", values="value")
)

df_result = eve_db.join(df_parsed)

converter = get_lifter('hg38', 'hg19', one_based=True)
hg38_start = []
hg38_chrom = []

for chrom, pos in zip(df_result["#CHROM"], df_result["POS"]):
    mappings = converter[str(chrom)][pos - 1]

    if mappings:
        c, p, strand = mappings[0]
        hg38_chrom.append(c)
        # REMOVED str(): Keep as integers
        hg38_start.append(p + 1)  
    else:
        hg38_chrom.append(None)
        hg38_start.append(None)

# REMOVED the [str(i) for i in ...] list comprehensions

df_result["hg19_chrom"] = hg38_chrom

# Convert to nullable Int64 so Polars sees them as integers, not objects/strings
df_result["hg19_pos"] = hg38_start



df_result["#CHROM"] = (
    df_result["#CHROM"]
    .astype(str)
    .str.replace(r"\.0$", "", regex=True)
)

# 2. Perform the merge
df_result = df_result.merge(
    cpt1, 
    left_on=["gene", "mutant"], 
    right_on=["uniprot_id", "mutant"], 
    how="right"
)
df_result["POS"] = df_result["POS"].astype("Int64")
df_result["hg19_pos"] = df_result["hg19_pos"].astype("Int64")


df_result = df_result.merge(gene_annot_dt[["gene_name", "geneID_short"]], left_on="uniprot_id", right_on="gene_name").drop(columns="gene_name")


df_result.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/CPT1_hg19_vars_all_protein_coding_genes.parquet", index=None)