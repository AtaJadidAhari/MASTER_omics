import pandas as pd
import polars as pl
from liftover import get_lifter
import plotnine as pn

from itables import init_notebook_mode
init_notebook_mode(all_interactive=True)
import itables.options as opt
opt.lengthMenu = [10, 20, 50, 100]
opt.pageLength = 20   # in newer versions; otherwise use show(..., pageLength=20)

import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *
from Variants.add_missing_acmg import *



def make_mutant(row):
    aa = row["Amino_acids"]
    pos = row["Protein_position"]

    if pd.isna(aa) or pd.isna(pos) or "/" not in str(aa):
        return pd.NA

    ref, alt = str(aa).split("/")
    
    # Protein_position can sometimes be "15-16"; keep only simple positions
    pos = str(pos).split("-")[0]

    return f"{ref}{pos}{alt}"





cpt1_all = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/CPT1_all_protein_coding_genes.parquet")
# cpt1_all

# somatic variants

# somatic = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/somatic_aggregated_predisp_cgc.parquet")

# somatic = somatic[somatic.columns[:-11]]

# somatic["protein_id"] = somatic["SWISSPROT"].str.split("_").str[0]


# somatic["mutant"] = somatic.apply(make_mutant, axis=1)
# somatic_cpt = somatic.merge(cpt1_all, how="left", right_on=["mutant", "uniprot_id"], left_on=["mutant", "protein_id"])

# somatic_cpt = somatic_cpt.rename(columns= {"CPT1_score": "somatic_CPT1_score"})
# somatic_cpt.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/somatic_aggregated_predisp_cgc.parquet", index=None)


# germline
germline = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/germline_aggregated_predisp_cgc.parquet")
germline = germline[germline.columns[:-12]]


germline["protein_id"] = germline["SWISSPROT"].str.split("_").str[0]
germline["mutant"] = germline.apply(make_mutant, axis=1)

germline_cpt = germline.merge(cpt1_all, how="left", right_on=["mutant", "uniprot_id"], left_on=["mutant", "protein_id"])

germline_cpt = germline_cpt.rename(columns= {"CPT1_score": "germline_CPT1_score"})
germline_cpt.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/germline_aggregated_predisp_cgc.parquet", index=None)
