import pandas as pd

snv_variant_type = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/variant_types/snv_variant_types.tsv", sep="\t")

vep_res_high = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique.tsv", sep="\t")
vep_res_high = pd.merge(vep_res_high, snv_variant_type, on=["sampleID", "Location"], how="left")
vep_res_high.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique_variant_type.tsv", sep="\t", index=False)



# indel_variant_type = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/variant_types/indel_variant_types.tsv", sep="\t")

# vep_res_high = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_indel_all_aggregated_unique.tsv", sep="\t")
# vep_res_high = pd.merge(vep_res_high, indel_variant_type, on=["sampleID", "Location"], how="left")
# vep_res_high.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_indel_all_aggregated_unique_variant_type.tsv", sep="\t", index=False)


# vep_res_high = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_VEP_splice_aggregated.tsv", sep="\t")
# vep_res_high = pd.merge(vep_res_high, snv_variant_type, on=["sampleID", "Location"], how="left")
# vep_res_high.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_VEP_splice_aggregated_variant_type.tsv", sep="\t")
