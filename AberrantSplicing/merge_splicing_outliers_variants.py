import pandas as pd

import sys
sys.path.append("/home/a379i/Scripts/")

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

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
    ["Gene", "sampleID", "IMPACT", "Consequence", "ANNOTATION_control", "promoterAI", "Location", "Allele", "ref", "#Uploaded_variation"]
].rename(columns={
    "IMPACT": "IMPACT_snv",
    "Consequence": "Consequence_snv",
    "ANNOTATION_control": "ANNOTATION_control_snv",
    "promoterAI": "promoterAI_snv",
    "#Uploaded_variation": "#Uploaded_variation_snv"
})

indel_collapsed = collapse_variants(indel_variants)[
    ["Gene", "sampleID", "IMPACT", "Consequence", "ANNOTATION_control",  "Location", "Allele"]
].rename(columns={
    "IMPACT": "IMPACT_indel",
    "Consequence": "Consequence_indel",
    "ANNOTATION_control": "ANNOTATION_control_indel",
})

fr_res_gene = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_results_gene_all.tsv", sep="\t")

fr_res_gene = fr_res_gene.merge(gene_annot_dt[["gene_name", "geneID_short", "gene_id"]], how="left", left_on="hgncSymbol", right_on="gene_name")

fr_res_gene = fr_res_gene.rename(columns={"gene_id":"geneID"})


print(fr_res_gene.shape)

fr_res_gene = fr_res_gene.merge(
    snv_collapsed,
    left_on=["geneID_short", "sampleID"],
    right_on=["Gene", "sampleID"],
    how="left"
)





fr_res_gene = fr_res_gene.merge(
    indel_collapsed,
    left_on=["geneID_short", "sampleID"],
    right_on=["Gene", "sampleID"],
    how="left",
    suffixes=("", "_indel")
)

print(fr_res_gene.shape)

fr_res_gene["variant_type"] = (
    fr_res_gene["IMPACT_snv"].notna().map({True: "snv", False: ""}) +
    fr_res_gene["IMPACT_indel"].notna().map({True: ",indel", False: ""})
).str.strip(",")


fr_res_gene["IMPACT"] = (
    fr_res_gene[["IMPACT_snv", "IMPACT_indel"]]
    .apply(
        lambda x: max(
            x.dropna(),
            key=lambda v: impact_rank.get(v, -1),
            default=pd.NA
        ),
        axis=1
    )
)

print(fr_res_gene.shape)



#pr_res_aberrant = pd.merge(pr_res_aberrant, sa, left_on="sampleID", right_on="pid", how="left")
fr_res_gene.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_results_gene_all_variants.tsv", index=None, sep="\t")