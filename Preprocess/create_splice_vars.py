import pandas as pd
from pathlib import Path

variant_base_path = Path("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/")

indel_vep_res = pd.read_csv(
    variant_base_path / f"vep_res_rare_snv_all_aggregated_unique.tsv",
    sep="\t"
)

indel_vep_res = indel_vep_res[indel_vep_res["Consequence"].str.contains("splic", na=False)]
# Set ordered categorical levels for IMPACT
impact_order = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
indel_vep_res["IMPACT"] = pd.Categorical(
        indel_vep_res["IMPACT"],
        categories=impact_order,
        ordered=True
)

# Sort by Gene, sampleID, IMPACT
indel_vep_res = indel_vep_res.sort_values(["Gene", "sampleID", "IMPACT"])
indel_vep_res.shape

indel_vep_res.to_csv(variant_base_path / f"vep_res_rare_snv_VEP_splice_aggregated.tsv",
    sep="\t")
