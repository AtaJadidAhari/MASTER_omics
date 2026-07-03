import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode
import polars as pl
import os
from mizani.formatters import custom_format

init_notebook_mode(all_interactive=True)
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt



import sys
sys.path.append("/home/a379i/Scripts/")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *
from utils.util_functions import *


somatic_snvs = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/somatic_snvs_vep_annotated.parquet")



sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")





germline_snvs  = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA.parquet").collect(engine="streaming").to_pandas()

germline_snvs = germline_snvs.rename(columns={"#CHROM": "seqnames", "REF": "Ref", "ALT": "Alt"})
print(germline_snvs.shape)

germline_snvs["chrom"] = "chr" + germline_snvs["seqnames"]
germline_snvs.shape

germline_snvs = germline_snvs.merge(gene_annot_dt[["gene_name", "geneID_short"]], left_on="HUGO_Symbol", right_on="gene_name")

valid_positions = (
    germline_snvs["POS"]
    .dropna()
    .unique()
    .astype(int)
    .tolist()
)

germline_snvs[["VEP_Ref", "VEP_Alt"]] = germline_snvs.apply(
    to_vep_alleles, axis=1
)

germline_snvs["Location"] = germline_snvs.apply(
    lambda row: compute_vep_location(
        {
            "#CHROM": row["seqnames"],
            "POS": row["POS"],
            "REF": row["VEP_Ref"],
            "ALT": row["VEP_Alt"],
        }
    ),
    axis=1,
)






germline_vep_res = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA_vep.parquet",).collect(engine="streaming").to_pandas()


germline_vep_res["chrom"] = germline_vep_res["Location"].str.split(":").str[0]
germline_vep_res["POS"] = 0
germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] == "SNV", "POS"] = germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] == "SNV", "Location"].str.split(":").str[1].astype(int)
# germline_vep_res["POS"] = germline_vep_res["POS"] + 1
germline_snvs = germline_snvs.merge(germline_vep_res[["Location", "chrom", "VARIANT_CLASS", "POS", "Allele", "REF_ALLELE", "IMPACT", "Consequence", "Gene", "#Uploaded_variation", "am_pathogenicity", "am_class", "LoF", "CADD_PHRED", "CADD_RAW","existing_InFrame_oORFs",  "existing_OutOfFrame_oORFs","existing_uORFs", "five_prime_UTR_variant_annotation", "five_prime_UTR_variant_consequence", "max_spliceai_score"]], left_on=["Location", "VEP_Ref", "VEP_Alt", "geneID_short"], right_on=["Location", "REF_ALLELE", "Allele", "Gene"], how="left") 

print(germline_snvs.head())

germline_snvs["seqnames"] = "chr" + germline_snvs["seqnames"] 
print(germline_snvs["VARIANT_CLASS"].value_counts())


absplice_predisp = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/predisp_cgc_max_abSplice_snvs_hg19.parquet")
        .filter((pl.col("hg19_end").is_in(valid_positions)))
        .collect(engine="streaming")
      ).to_pandas()
absplice_predisp.shape




print(len(germline_snvs))
merged_vars = germline_snvs.merge(absplice_predisp, left_on=["seqnames", "POS_x", "Ref", "Alt"], right_on=["chrom", "hg19_end", "ref", "alt"], how="left")
print(len(merged_vars))
absplice_samples = merged_vars[(merged_vars["AbSplice2_max"] >= 0.05) | (merged_vars["pangolin_score"] >= 0.2) | (merged_vars["IMPACT"].str.contains("HIGH")) | (merged_vars["spliceAI_DS_gt_04"].notna()) | (merged_vars["max_spliceai_score"] >= 0.2) ]["sampleID"]
# merged_vars.loc[merged_vars["AbSplice2_max"].isna(), "AbSplice2_max"] = 0

print(merged_vars["VARIANT_CLASS"].value_counts())


all_res = []

for cohort in sa[sa["Diag"] != "Unstranded_data"]["Diag"].unique():
    
    fr_res_new =( pl.scan_csv(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/{cohort}/results_per_junction.tsv", 
                                schema_overrides={
                                    "deltaPsi": pl.Float64,
                                    "seqnames": pl.String,
                                }, 
                                null_values=["NA"],
                                 separator="\t")
                    .select("sampleID", "hgncSymbol", "pValue", "padjust", 'seqnames', 'start', 'end', 'strand', 'deltaPsi', 'counts', 'totalCounts', 'nonsplitCounts')
                    .filter(pl.col("sampleID").is_in(absplice_samples)).collect(engine="streaming")
                    ).rename({"start": "junction_start", "end": "junction_end"}).to_pandas()
    # fr_res_new[(fr_res_new["padjust"] <= 0.05) & (fr_res_new["padjust_Genes_to_test_on_all_samples"].notna())]
    fr_res_new = fr_res_new.merge(sa, left_on="sampleID", right_on="pid")
    

    # fr_res_new["key"] = fr_res_new["sampleID"] +  "." + fr_res_new["nct_pid"].str.split("_").str[1]
    predisp_absplice = merged_vars.merge(fr_res_new, right_on=["sampleID", "hgncSymbol"], left_on=["sampleID", "HUGO_Symbol"], how="inner")
    all_res.append(predisp_absplice[(predisp_absplice["padjust"] <= 0.1) & (predisp_absplice["deltaPsi"].abs() >= 0.1)])
all_res = pd.concat(all_res)


all_res["Variant"] = (
    all_res["seqnames_x"] + ":" + 
    all_res["POS_x"].astype(str) + ":" + 
    all_res["Ref"] + ":" + 
    all_res["Alt"]
)

columns = ["sampleID", "Variant", "junction_start" , "IMPACT", "five_prime_UTR_variant_consequence", "Consequence", "junction_end", "deltaPsi", "max_spliceai_score", "pangolin_score", "AbSplice2_max", "Oncotree Code", "hgncSymbol", "Oncotree Text", "VEP_Most_Severe_Consequence", "VEP_IMPACT", "Diag", "Tumorzellgehalt (Bioinformatik)",  "padjust", "pValue", "ClinVar_Pathogenicity", "ACMG_Classification", "CharGer_Classification",
                 "score_category", "strand", "Tumor_VAF", "Control_VAF", "FILTER" , "counts", "totalCounts", "nonsplitCounts", "nct_pid", "spliceAI_DS_gt_04"]
#   
print(all_res.shape, "herere")
all_res.drop_duplicates(subset=["sampleID", "hgncSymbol"]) # total number of outliers with germline variant

print(all_res.shape, "after duplication")
all_res_splice = all_res[(all_res["AbSplice2_max"] >= 0.05) | (all_res["pangolin_score"] >= 0.2) |(all_res["Consequence"].isin(["splice_acceptor_variant", "splice_donor_variant"])) | (all_res["max_spliceai_score"] >= 0.2)]
# & (all_res_splice["POS"] >= all_res_splice["junction_start"] - 1000) & (all_res_splice["POS"] <= all_res_splice["junction_end"] + 1000)]




all_res_splice = all_res_splice.reset_index(drop=True)
all_res_splice.to_parquet("~/splicing_vus.parquet", index=None)

all_res.to_parquet("~/all_splicing.parquet", index=None)

(merged_vars[["sampleID"]]
 .drop_duplicates()
 .to_parquet("~/germline_variant_samples.parquet", index=False))