import pandas as pd
import plotnine as pn
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


outlier_res = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/pyoutrider_zscores/aggregated_sf_zScores_all.parquet"
# outlier_res = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/protrider_summary.parquet"
if outlier_res.split(".")[-1] == "parquet":
    py_or_res_all = pd.read_parquet(outlier_res)
else:
    py_or_res_all = pd.read_csv(outlier_res, sep="\t")

# py_or_res_all = py_or_res_all.rename(columns={"sf_zScore": "zScore"})

print(py_or_res_all.shape)
py_or_res_all["geneID_short"] = py_or_res_all["geneID"].str.split(".").str[0]

py_or_res_all = py_or_res_all.merge(
    snv_collapsed,
    left_on=["geneID_short", "sampleID"],
    right_on=["Gene", "sampleID"],
    how="left"
)

print(py_or_res_all.shape)

py_or_res_all = py_or_res_all.merge(
    indel_collapsed,
    left_on=["geneID_short", "sampleID"],
    right_on=["Gene", "sampleID"],
    how="left",
    suffixes=("", "_indel")
).drop(columns="Gene_indel")

print(py_or_res_all.shape)

py_or_res_all["variant_type"] = (
    py_or_res_all["IMPACT_snv"].notna().map({True: "snv", False: ""}) +
    py_or_res_all["IMPACT_indel"].notna().map({True: ",indel", False: ""})
).str.strip(",")


py_or_res_all["IMPACT"] = (
    py_or_res_all[["IMPACT_snv", "IMPACT_indel"]]
    .apply(
        lambda x: max(
            x.dropna(),
            key=lambda v: impact_rank.get(v, -1),
            default=pd.NA
        ),
        axis=1
    )
)

if "padjust" in py_or_res_all.columns: 
    py_or_res_all['Outlier status'] = "Non-outlier"
    py_or_res_all.loc[(py_or_res_all['zScore'] < 0) & (py_or_res_all['padjust'] <= 0.05), 'Outlier status'] = "Underexpression"
    py_or_res_all.loc[(py_or_res_all['zScore'] > 0) & (py_or_res_all['padjust'] <= 0.05), 'Outlier status'] = "Overexpression"

print(py_or_res_all.shape)


py_or_res_all.to_csv(f"{"/".join(outlier_res.split("/")[:-1])}/or_variants.parquet", index=None)

# joined = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/rna_protein_outliers.parquet")

# joined["Outlier status"] = "Non outlier"
# joined.loc[(joined["RNA_aberrant"] == False) & (joined["Protein_aberrant"] == True), "Outlier status"] = "Protein"
# joined.loc[(joined["RNA_aberrant"] == True) & (joined["Protein_aberrant"] == True), "Outlier status"] = "RNA & Protein"
# joined.loc[(joined["RNA_aberrant"] == True) & (joined["Protein_aberrant"] == False), "Outlier status"] = "RNA"


# joined["VUS"] = "No rare variant"
# joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == True), "VUS"] = "RNA & protein variant"
# joined.loc[(joined["VUS_x"] == False) & (joined["VUS_y"] == True), "VUS"] = "Protein variant"
# joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == False), "VUS"] = "RNA variant"


# p = (pn.ggplot(joined) + 
#     pn.geom_point(pn.aes(x="RNA_zScore", y="Protein_zScore", color="Outlier status", shape="VUS"), alpha=0.5) +
#     pn.scale_size_manual(
#         values={
#             "No rare variant": 1,
#             "RNA variant": 3,
#             "Protein variant": 3,
#             "RNA & protein variant": 3
#         }
#     ) + 
#      pn.scale_color_manual(
#         values={
#             "Non outlier": "grey",
#             "Protein": "lightblue",
#             "RNA": "firebrick",
#             "RNA & Protein": "lightgreen"
#         }
#     ) + 
#      pn.theme_bw()+
#     pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
#      pn.coord_cartesian(ylim=[-25, 22], xlim=[-25, 22])
    
# )

# p.save("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/promoter_ai/or_vs_pr_zScores.png", dpi=600, width=8, height=6, units="in")


