variant_base_path = Path("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/")
benchmark_criteria = "HIGH"

# Read tables
snv_vep_res = pd.read_csv(
    variant_base_path / f"vep_res_rare_snv_{benchmark_criteria}_aggregated.tsv",
    sep="\t"
)
indel_vep_res = pd.read_csv(
    variant_base_path / f"vep_res_rare_indel_{benchmark_criteria}_aggregated.tsv",
    sep="\t"
)

# Combine
vep_res_combined = pd.concat([snv_vep_res, indel_vep_res], ignore_index=True)

print(len(vep_res_combined))

# Unique_based_on_criteria
if benchmark_criteria == "CADD_PHRED":
    # Sort decreasing by CADD_PHRED
    vep_res_combined = (
        snv_vep_res.sort_values("CADD_PHRED", ascending=False)
        .groupby(["sampleID", "Gene"])
        .head(1)
        .reset_index(drop=True)
    )
else:
    # Set IMPACT as ordered category and sort
    impact_order = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
    vep_res_combined["IMPACT"] = pd.Categorical(
        vep_res_combined["IMPACT"],
        categories=impact_order,
        ordered=True
    )

    vep_res_combined = (
        vep_res_combined.sort_values(["Gene", "sampleID", "IMPACT"])
    )

# Remove duplicates by Gene + sampleID, keeping first (like data.table)
vep_res_combined = (
    vep_res_combined.drop_duplicates(subset=["Gene", "sampleID"])
)
print(len(vep_res_combined))