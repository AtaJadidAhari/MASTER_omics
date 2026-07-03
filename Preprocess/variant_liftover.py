from liftover import get_lifter
import pandas as pd
import numpy as np
import sys


# var_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type"
# # hg19_variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique_variant_type.tsv", sep="\t")
# # print(hg19_variants.shape, " rows in hg19 variants")

# # germline SNVS

# var_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_snv"
# hg19_variants = pd.read_csv(f"{var_path}.tsv", sep="\t")

# converter = get_lifter('hg19', 'hg38', one_based=True)
# hg38_chrom = []
# hg38_pos = []

# for chrom, pos in zip(hg19_variants["seqnames"], hg19_variants["POS"]):
#     #chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
#     mappings = converter[str(chrom)][pos - 1]  # convert to 0-based

#     if mappings:
#         c, p, strand = mappings[0]
#         hg38_chrom.append(c)
#         hg38_pos.append(str(p + 1))  # back to 1-based
#     else:
#         hg38_chrom.append(None)
#         hg38_pos.append(None)
# hg38_pos = [str(i) for i in hg38_pos]
# hg19_variants["hg38_chrom"] = hg38_chrom
# hg19_variants["hg38_pos"] = hg38_pos
# print(hg19_variants.shape, "hg38 variants")

# hg19_variants.to_csv(f"{var_path}_hg38.tsv", sep="\t", index=None)
# hg19_variants.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38.tsv", sep="\t", index=None)


# hg19_variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique_variant_type.tsv", sep="\t")
# print(hg19_variants.shape, " rows in hg19 variants")



# var_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/promoterAI_tss500.parquet"
# hg19_variants = pd.read_parquet(var_path)

# converter = get_lifter('hg38', 'hg19', one_based=True)
# hg38_chrom = []
# hg38_start = []
# hg38_end = []

# for chrom, pos in zip(hg19_variants["chrom"], hg19_variants["pos"]):
#     mappings = converter[str(chrom)][pos - 1]

#     if mappings:
#         c, p, strand = mappings[0]
#         hg38_chrom.append(c)
#         # REMOVED str(): Keep as integers
#         hg38_start.append(p + 1)  
#         hg38_end.append(p + 2)    
#     else:
#         hg38_chrom.append(None)
#         hg38_start.append(None)
#         hg38_end.append(None)

# hg19_variants["hg19_chrom"] = hg38_chrom

# # Convert to nullable Int64 so Polars sees them as integers, not objects/strings
# hg19_variants["hg19_start"] = hg38_start
# hg19_variants["hg19_end"] = hg38_end

# print(hg19_variants.shape, "hg19 variants")

# hg19_variants.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/promoterAI_tss500_hg19.parquet")

# sys.exit()



# germline absplice

var_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/predisp_cgc_max_abSplice_snvs.parquet"
hg19_variants = pd.read_parquet(var_path)

converter = get_lifter('hg38', 'hg19', one_based=True)
hg38_chrom = []
hg38_start = []
hg38_end = []

for chrom, pos in zip(hg19_variants["chrom"], hg19_variants["start"]):
    mappings = converter[str(chrom)][pos - 1]

    if mappings:
        c, p, strand = mappings[0]
        hg38_chrom.append(c)
        # REMOVED str(): Keep as integers
        hg38_start.append(p + 1)  
        hg38_end.append(p + 2)    
    else:
        hg38_chrom.append(None)
        hg38_start.append(None)
        hg38_end.append(None)

# REMOVED the [str(i) for i in ...] list comprehensions

hg19_variants["hg19_chrom"] = hg38_chrom

# Convert to nullable Int64 so Polars sees them as integers, not objects/strings
hg19_variants["hg19_start"] = hg38_start
hg19_variants["hg19_end"] = hg38_end

print(hg19_variants.shape, "hg19 variants")

hg19_variants.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/predisp_cgc_max_abSplice_snvs_hg19.parquet")



# # germline indel
# var_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_indel"
# # germline SNVS
# hg19_variants = pd.read_csv(f"{var_path}.tsv.gz", sep="\t")

# converter = get_lifter('hg19', 'hg38', one_based=True)
# hg38_chrom = []
# hg38_start = []
# hg38_end = []


# for i, row in hg19_variants.iterrows():
#     #chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
#     mappings_start = converter[str(row["seqnames"])][row["start"] - 1]  # convert to 0-based
#     mappings_end = converter[str(row["seqnames"])][row["end"] - 1]  # convert to 0-based
#     if mappings_start:
#         c, p, strand = mappings_start[0]
#         hg38_start.append(str(p + 1))  # back to 1-based
#     else:
#         hg38_start.append(None)
#         hg38_chrom.append(None)
#         hg38_end.append(None)
#         continue

#     if mappings_start:
#         c, p, strand = mappings_start[0]
#         hg38_end.append(str(p + 1))  # back to 1-based
#     else:
#         hg38_end.append(None)
#         hg38_chrom.append(None)
#         hg38_start.append(None)
#         continue

#     hg38_chrom.append(c)

        
# hg38_start = [str(i) for i in hg38_start]
# hg38_end = [str(i) for i in hg38_end]
# hg19_variants["hg38_chrom"] = hg38_chrom
# hg19_variants["hg38_start"] = hg38_start
# hg19_variants["hg38_end"] = hg38_start

# print(hg19_variants.shape, "hg38 variants")

# hg19_variants.to_csv(f"{var_path}_hg38.tsv.gz", sep="\t", compression="gzip", index=None)