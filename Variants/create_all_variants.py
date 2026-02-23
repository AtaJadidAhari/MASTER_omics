# import pandas as pd
# import sys
# sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

# from Preprocess.VEP_SORanking import *

# # Read input
# snv_vep_res = pd.read_csv(
#     "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated.tsv",
#     sep="\t"
# )

# # drop all variants that do not match to any gene
# snv_vep_res = snv_vep_res[snv_vep_res["Gene"] != '-']

# # IMPACT as ordered categorical
# impact_levels = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
# snv_vep_res["IMPACT"] = pd.Categorical(
#     snv_vep_res["IMPACT"],
#     categories=impact_levels,
#     ordered=True
# )

# # Compute most severe consequence
# def most_severe_consequence(consequence, severity_rank):
#     parts = [x.strip() for x in consequence.split(",")]
#     parts_sorted = sorted(parts, key=lambda x: severity_rank.get(x, float("inf")))
#     return parts_sorted[0]

# snv_vep_res["Consequence_most_severe"] = snv_vep_res["Consequence"].apply(
#     lambda x: most_severe_consequence(x, severity_rank)
# )

# # Sort like setorder(Gene, sampleID, IMPACT)
# snv_vep_res = snv_vep_res.sort_values(
#     by=["Gene", "sampleID", "IMPACT"]
# )

# # Replace missense_variant with am_class when applicable
# mask = (
#     (snv_vep_res["Consequence_most_severe"] == "missense_variant") &
#     (snv_vep_res["am_class"] != "-")
# )

# snv_vep_res.loc[mask, "Consequence_most_severe"] = snv_vep_res.loc[mask, "am_class"]




# # Keep first row per Gene + sampleID (after sorting)
# vep_res_unique = snv_vep_res.drop_duplicates(
#     subset=["Gene", "sampleID"],
#     keep="first"
# )

# # Print number of rows
# print(len(vep_res_unique))

# snv_variant_type = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/variant_types/snv_variant_types.tsv", sep="\t")

# vep_res_unique = pd.merge(vep_res_unique, snv_variant_type, on=["sampleID", "Location"], how="left")
# vep_res_unique.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique.tsv", sep="\t", index=False)



import polars as pl
import sys
sys.path.append("/home/a379i/Scripts")

from Preprocess.VEP_SORanking import *

# -------------------------
# Read input
# -------------------------
snv_vep_res = pl.read_csv(
    "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated.tsv",
    separator="\t"
)

# -------------------------
# Drop variants without gene
# -------------------------
snv_vep_res = snv_vep_res.filter(pl.col("Gene") != "-")

# -------------------------
# IMPACT as ordered categorical
# -------------------------
impact_rank = {
    "HIGH": 0,
    "MODERATE": 1,
    "LOW": 2,
    "MODIFIER": 3,
}

snv_vep_res = snv_vep_res.with_columns(
    pl.col("IMPACT")
    .map_elements(lambda x: impact_rank.get(x, 99))
    .alias("IMPACT_rank")
)

snv_vep_res = snv_vep_res.sort(
    ["Gene", "sampleID", "IMPACT_rank"]
)

snv_vep_res = snv_vep_res.drop("IMPACT_rank")

# -------------------------
# Compute most severe consequence
# -------------------------
def most_severe_consequence(consequence: str) -> str:
    parts = [x.strip() for x in consequence.split(",")]
    parts_sorted = sorted(parts, key=lambda x: severity_rank.get(x, float("inf")))
    return parts_sorted[0]

snv_vep_res = snv_vep_res.with_columns(
    pl.col("Consequence")
    .map_elements(most_severe_consequence)
    .alias("Consequence_most_severe")
)

# -------------------------
# Sort like setorder(Gene, sampleID, IMPACT)
# -------------------------




# -------------------------
# Replace missense_variant with am_class when applicable
# -------------------------
snv_vep_res = snv_vep_res.with_columns(
    pl.when(
        (pl.col("Consequence_most_severe") == "missense_variant") &
        (pl.col("am_class") != "-")
    )
    .then(pl.col("am_class"))
    .otherwise(pl.col("Consequence_most_severe"))
    .alias("Consequence_most_severe")
)

# -------------------------
# Keep first row per Gene + sampleID
# (after sorting)
# -------------------------
vep_res_unique = (
    snv_vep_res
    .group_by(["Gene", "sampleID"], maintain_order=True)
    .agg(pl.all().first())
)

# -------------------------
# Print number of rows
# -------------------------
print(vep_res_unique.height)

# -------------------------
# Merge variant types
# -------------------------
# snv_variant_type = pl.read_csv(
#     "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/variant_types/snv_variant_types.tsv",
#     separator="\t"
# )

# vep_res_unique = vep_res_unique.join(
#     snv_variant_type,
#     on=["sampleID", "Location"],
#     how="left"
# )

# -------------------------
# Write output
# -------------------------
vep_res_unique.write_csv(
    "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique.tsv",
    separator="\t"
)

