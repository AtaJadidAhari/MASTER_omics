from liftover import get_lifter
import pandas as pd



hg19_variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique_variant_type.tsv", sep="\t")
print(hg19_variants.shape, " rows in hg19 variants")


converter = get_lifter('hg19', 'hg38', one_based=True)
hg38_chrom = []
hg38_pos = []

for chrom, pos in zip(hg19_variants["#CHROM"], hg19_variants["POS"]):
    #chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    mappings = converter[str(chrom)][pos - 1]  # convert to 0-based

    if mappings:
        c, p, strand = mappings[0]
        hg38_chrom.append(c)
        hg38_pos.append(str(p + 1))  # back to 1-based
    else:
        hg38_chrom.append(None)
        hg38_pos.append(None)
hg38_pos = [str(i) for i in hg38_pos]
hg19_variants["hg38_chrom"] = hg38_chrom
hg19_variants["hg38_pos"] = hg38_pos
print(hg19_variants.shape, "hg38 variants")

hg19_variants.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38.tsv", sep="\t", index=None)