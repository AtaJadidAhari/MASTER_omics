import pandas as pd

import sys
sys.path.append("/home/a379i/Scripts/")

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

#pr_output_name = "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"
pr_output_name = "sf_zScores"

pr_res_all = load_pr_data("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/protrider_summary.csv")

print(pr_res_all.shape)


impact_rank = {
    "HIGH": 3,
    "MODERATE": 2,
    "LOW": 1,
    "MODIFIER": 0
}

def collapse_variants(df, impact_col="IMPACT"):
    return (
        df
        .assign(impact_rank=df[impact_col].map(impact_rank))
        .sort_values("impact_rank", ascending=False)
        .drop_duplicates(subset=["Gene", "sampleID"])
        .drop(columns="impact_rank")
    )
    
snv_variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38_promoterAI.tsv", sep="\t")
snv_variants["variant_type"] = "snv"

indel_variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_indel_all_aggregated_unique_variant_type.tsv", sep="\t")
indel_variants["variant_type"] = "indel"

snv_collapsed = collapse_variants(snv_variants)[
    ["Gene", "sampleID", "IMPACT", "Consequence", "ANNOTATION_control", "#CHROM", "POS", "am_class", "promoterAI", "Allele", "ref", "#Uploaded_variation"]
].rename(columns={
    "IMPACT": "IMPACT_snv",
    "Consequence": "Consequence_snv",
    "ANNOTATION_control": "ANNOTATION_control_snv",
    "promoterAI": "promoterAI_snv",
    "#CHROM": "#CHROM_snv",
    "POS": "POS_snv",
    "Allele": "Allele_snv",
    "ref": "ref_snv",
    "am_class": "am_class_snv",
    "#Uploaded_variation": "#Uploaded_variation_snv"
})

indel_collapsed = collapse_variants(indel_variants)[
    ["Gene", "sampleID", "IMPACT", "Consequence", "am_class", "ANNOTATION_control", "#CHROM", "POS", "Allele", "#Uploaded_variation"]
].rename(columns={
    "IMPACT": "IMPACT_indel",
    "Consequence": "Consequence_indel",
    "ANNOTATION_control": "ANNOTATION_control_indel",
    "#CHROM": "#CHROM_indel",
    "POS": "POS_indel",
    "Allele": "Allele_indel",
    "am_class": "am_class_indel",
    "#Uploaded_variation": "#Uploaded_variation_indel"
})


pr_res_all = pr_res_all.merge(
    snv_collapsed,
    left_on=["geneID_short", "sampleID"],
    right_on=["Gene", "sampleID"],
    how="left"
)

print(pr_res_all.shape)

pr_res_all = pr_res_all.merge(
    indel_collapsed,
    left_on=["geneID_short", "sampleID"],
    right_on=["Gene", "sampleID"],
    how="left",
    suffixes=("", "_indel")
)

print(pr_res_all.shape)

pr_res_all["variant_type"] = (
    pr_res_all["IMPACT_snv"].notna().map({True: "snv", False: ""}) +
    pr_res_all["IMPACT_indel"].notna().map({True: ",indel", False: ""})
).str.strip(",")


pr_res_all["IMPACT"] = (
    pr_res_all[["IMPACT_snv", "IMPACT_indel"]]
    .apply(
        lambda x: max(
            x.dropna(),
            key=lambda v: impact_rank.get(v, -1),
            default=pd.NA
        ),
        axis=1
    )
)

if "padjust" in pr_res_all.columns:
    pr_res_all['Outlier status'] = "Non-outlier"
    pr_res_all.loc[(pr_res_all['zScore'] < 0) & (pr_res_all['padjust'] <= 0.1), 'Outlier status'] = "Underexpression"
    pr_res_all.loc[(pr_res_all['zScore'] > 0) & (pr_res_all['padjust'] <= 0.1), 'Outlier status'] = "Overexpression"

print(pr_res_all.shape)



#pr_res_aberrant = pd.merge(pr_res_aberrant, sa, left_on="sampleID", right_on="pid", how="left")
pr_res_all.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants.csv", index=None)