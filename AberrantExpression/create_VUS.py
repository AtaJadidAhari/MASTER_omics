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

from utils.util_functions import *
from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

somatic_snvs = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/somatic_snvs_vep_annotated.parquet")


sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")



or_res_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_interesting_genes_padjust_cnv.parquet"
needed_cols = ["sampleID", "zScore", "pValue", "padjust", "IMPACT", "geneID",  "geneID_short",
               "#Uploaded_variation_snv", "IMPACT_snv", "ANNOTATION_control_snv", "Consequence_snv",
               "Location_indel", "IMPACT_indel", "ANNOTATION_control_indel", "Consequence_indel", "promoterAI_snv",
               "padjust_predisp", "padjust_predisp_extended", "CNV", "padjust_genes_of_interest" ]

py_or_res_all = (
    pl.read_parquet(or_res_path, use_pyarrow=True)
    .filter(
        (pl.col("padjust_genes_of_interest") <= 0.1)
    )
    .select(needed_cols)
    .to_pandas()
)

py_or_res_all["Method"] = "OUTRIDER"

# py_or_res_all = py_or_res_all[(py_or_res_all["padjust_predisp_extended"] <= 0.05) | (py_or_res_all["padjust"] <= 0.05)]


py_or_res_all = pd.merge(py_or_res_all, sa[["pid", "Diag", "seq_type", "Oncotree Code", "Oncotree Text", "nct_pid", "Tumorzellgehalt (Bioinformatik)"]], left_on="sampleID", right_on="pid")
py_or_res_all = pd.merge(py_or_res_all, cgc[[ "ROLE_IN_CANCER", "geneID"]], on="geneID", how="left")

py_or_res_all = py_or_res_all.sort_values("pValue")
py_or_res_aberrant = py_or_res_all.drop_duplicates(subset=["sampleID", "geneID_short"])

# py_or_res_aberrant["key"] = py_or_res_aberrant["sampleID"] +  "." + py_or_res_aberrant["nct_pid"].str.split("_").str[1]
py_or_res_aberrant = py_or_res_aberrant[(py_or_res_aberrant["padjust_predisp_extended"].notna()) | (py_or_res_aberrant["geneID_short"].isin(cgc["geneID_short"]))]
sample_ids = py_or_res_aberrant["sampleID"].unique()
py_or_res_aberrant = py_or_res_aberrant.merge(somatic_snvs[["somatic_snv_#Uploaded_variation", "sampleID", "MASTER_annotated_gene", "vep_Gene", "somatic_snv_IMPACT", "somatic_snv_Consequence", "somatic_snv_max_spliceai_score", "somatic_snv_am_pathogenicity", "somatic_snv_AbSplice2_max", "somatic_snv_promoterAI", "somatic_snv_abexp_v1.1"]] , left_on=["geneID_short", "sampleID"], right_on=["vep_Gene", "sampleID"], how="left")



pr_output_name = "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"

pr_res_all = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants_interesting_genes_padjust_cnv.parquet").filter((pl.col("padjust_genes_of_interest") <= 0.1)).collect(engine="streaming").to_pandas()
pr_res_all = pd.merge(pr_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID_short", "ROLE_IN_CANCER", "predisposition_gene"]], right_on="geneID_short", left_on="geneID", how="left")
pr_res_all = pd.merge(pr_res_all, sa, left_on="sampleID", right_on="pid")
pr_res_aberrant = pr_res_all[(pr_res_all["padjust_genes_of_interest"] <= 0.1)]
pr_res_aberrant = pr_res_aberrant.merge(somatic_snvs[["somatic_snv_#Uploaded_variation", "sampleID", "MASTER_annotated_gene", "vep_Gene", "somatic_snv_IMPACT", "somatic_snv_Consequence", "somatic_snv_max_spliceai_score", "somatic_snv_am_pathogenicity", "somatic_snv_AbSplice2_max", "somatic_snv_promoterAI", "somatic_snv_abexp_v1.1"]] , left_on=["geneID", "sampleID"], right_on=["vep_Gene", "sampleID"], how="left")

pr_samples = pr_res_aberrant["sampleID"].unique()




germline_snvs  = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA.parquet").filter((pl.col("sampleID").is_in(sample_ids)) | (pl.col("sampleID").is_in(pr_samples))).collect(engine="streaming").to_pandas()

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
germline_snvs = germline_snvs.merge(germline_vep_res[["Location", "VARIANT_CLASS", "chrom", "POS", "Allele", "REF_ALLELE", "IMPACT", "Consequence", "Gene", "#Uploaded_variation", "am_pathogenicity", "am_class", "LoF", "CADD_PHRED", "CADD_RAW","existing_InFrame_oORFs",  "existing_OutOfFrame_oORFs","existing_uORFs", "five_prime_UTR_variant_annotation", "five_prime_UTR_variant_consequence", "max_spliceai_score"]], left_on=["Location", "VEP_Ref", "VEP_Alt", "geneID_short"], right_on=["Location", "REF_ALLELE", "Allele", "Gene"], how="left") 


print(germline_snvs["VARIANT_CLASS"].value_counts())



abexp = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/predisp_cgc_max_abexp.parquet").filter((pl.col("start").is_in(valid_positions)) | (pl.col("end").is_in(valid_positions))).collect(engine="streaming").to_pandas()

# abexp = abexp.sort_values("abs_val", ascending=False).drop_duplicates(subset=["chrom", "start", "end", "ref", "alt"])


# abexp = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/predisp_max_abexp.parquet").filter(pl.col("chrom") == "chr4").collect(engine="streaming").to_pandas()


abexp["seqnames"] = abexp["chrom"].str.split("chr").str[1]

promoter_ai = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/promoterAI_tss500_hg19.parquet")
                .filter((pl.col("gene_id").is_in(extended_dresden_dt["geneID_short"])) | (pl.col("gene_id").is_in(cgc[cgc["geneID_short"].notna()]["geneID_short"])))
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


merged_vars = germline_snvs.merge(promoter_ai, left_on=["seqnames", "POS_x", "Ref", "Alt", "geneID_short"], right_on=["chrom", "hg19_start", "ref", "alt", "gene_id"], how="left").drop(columns=["chrom_x", "chrom"])      

merged_vars = merged_vars.merge(abexp, left_on=["POS_x", "seqnames", "Ref", "Alt", "geneID_short"], right_on=["end", "seqnames", "ref", "alt", "gene"], how="left")

merged_vars = merged_vars.rename(columns={"#Uploaded_variation": "germline_#Uploaded_variation", "Consequence": "germline_Consequence", "IMPACT": "germline_IMPACT", "am_pathogenicity": "germline_am_pathogenicity","LoF": "germline_am_LoF", "max_spliceai_score": "germilne_max_spliceai_score", "pangolin_score": "germline_pangolin_score", "AbSplice2_max": "germline_AbSplice2_max", "promoterAI": "germline_promoterAI", "abexp_v1.1": "germline_abexp_v1.1" })


merged_vars["Variant"] = (
    merged_vars["seqnames"] + ":" + 
    merged_vars["POS_x"].astype(str) + ":" + 
    merged_vars["Ref"] + ":" + 
    merged_vars["Alt"]
)

print(len(py_or_res_aberrant[py_or_res_aberrant['zScore']<=0]), "gene underexpression outleiers")


predisp_absplice = merged_vars.merge(py_or_res_aberrant, right_on=["sampleID", "geneID_short"], left_on=["sampleID", "gene_id"], how="left")
# predisp_absplice = predisp_absplice.drop_duplicates(subset=["pid", "start", "end", "Ref", "Alt"])

predisp_absplice["Outlier status"] = False
predisp_absplice.loc[predisp_absplice["geneID_short_y"].notna(), "Outlier status"] = True

print(predisp_absplice["VARIANT_CLASS"].value_counts())


predisp_absplice[predisp_absplice["geneID_short_y"].notna()]["VARIANT_CLASS"].value_counts()



smallest_columns = ["sampleID", "Variant", "germline_promoterAI", "germline_abexp_v1.1", "germline_IMPACT", "five_prime_UTR_variant_consequence", "germline_Consequence", "somatic_snv_promoterAI", "somatic_snv_abexp_v1.1", "somatic_snv_IMPACT", "somatic_snv_Consequence", "VEP_IMPACT", "VEP_Most_Severe_Consequence", "HGVSg", "HGVSc", "CharGer_Classification" , "ClinVar_Pathogenicity" , "zScore", "HUGO_Symbol", "Oncotree Code", "Oncotree Text", "IMPACT_snv", "ROLE_IN_CANCER", "CNV", "Tumorzellgehalt (Bioinformatik)"]



predisp_snvs_underexpression = predisp_absplice[(predisp_absplice["zScore"] <= 0) & (predisp_absplice["Outlier status"] == True)]

final_gene_expression = predisp_snvs_underexpression[
    (~predisp_snvs_underexpression["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])) &
    ( (predisp_snvs_underexpression["germline_abexp_v1.1"] <= -0.1) |
    (predisp_snvs_underexpression["germline_promoterAI"] <= -0.1) |
      ((predisp_snvs_underexpression["five_prime_UTR_variant_consequence"] != "-") & (predisp_snvs_underexpression["five_prime_UTR_variant_consequence"].notna())) |
    (predisp_snvs_underexpression["germline_IMPACT"].str.contains("HIGH")))]

# (predisp_snvs_underexpression["CharGer_Classification"] != "Pathogenic") &


final_gene_expression["VUS"] = "RNA_outlier"

final_gene_expression
# predisp_snvs_underexpression[smallest_columns]


final_gene_expression = final_gene_expression.reset_index(drop=True)

predisp_snvs_underexpression.to_parquet("~/gene_exp_vus.parquet", index=None)


print(len(pr_res_aberrant[pr_res_aberrant['zScore']<=0]), "protein underexpression outleiers")


pr_vars = merged_vars.merge(pr_res_aberrant, left_on=["sampleID", "HUGO_Symbol"],  right_on=["sampleID", "proteinID"],how="left")
pr_vars["Outlier status"] = False
pr_vars.loc[pr_vars["proteinID"].notna(), "Outlier status"] = True


# smallest_columns = ["sampleID", "Variant", "promoterAI", "abexp_v1.1", "vep_impact", "vep_consequence", "zScore", "proteinID", "Oncotree Code",  "IMPACT_snv", "#Uploaded_variation_snv", "Consequence_snv", "CharGer_Classification", "control_zygosity", "CNV", "Tumorzellgehalt (Bioinformatik)"]


final_pr_vus = pr_vars[(pr_vars["CharGer_Classification"] != "Pathogenic") & (pr_vars["proteinID"].isin(AD_inheritence["Approved symbol (HGNC)"])) & (pr_vars["Outlier status"] == True) & (pr_vars["zScore"] <= 0)  & 
                        ((pr_vars["germline_abexp_v1.1"] <= -0.1) | (pr_vars["germline_IMPACT"].str.contains("HIGH")) | (pr_vars["germline_Consequence"].str.contains("missense")) |  
                        (pr_vars["germline_promoterAI"] <= -0.1) | ((pr_vars["five_prime_UTR_variant_consequence"] != "-") & (pr_vars["five_prime_UTR_variant_consequence"].notna())))]
final_pr_vus["VUS"] = "Protein_outlier"


final_pr_vus = final_pr_vus.reset_index(drop=True)

pr_vars.to_parquet("~/protein_exp_vus.parquet", index=None)
