

# %%
import pandas as pd
import os
import gzip


# %%
import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

# %%
from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

# %%
sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")
samples = sa.pid


# %%
genes = dresden_dt_cgc["gene_name"].values
output = []
all_genes = []
# for sample in sa.pid:
#     vcf_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/standard_vcf_snv/snvs_H021-{sample}.standard_snv.vcf.gz"
#     if not os.path.exists(vcf_path):
#         vcf_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/standard_vcf_snv/snvs_P021-{sample}.standard_snv.vcf.gz"
#         if not os.path.exists(vcf_path):
#             print('no vcf for', sample)
#         continue
#     with gzip.open(vcf_path, "rt") as fin, open(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/predisposition_variants/snvs_P021-{sample}.standard_snv_predisposition.vcf.gz", "w") as fout:
#         for line in fin:
#             if line.startswith("#"):
#                 fout.write(line)
#                 continue
#             if line.split("\t")[6] == "PASS":
#                 info = line.split("\t")[7]
#                 if "GENE=" not in info:
#                             continue

# %% [markdown]
#                 gene_field = info.split("GENE=")[1].split(";")[0]
#                 gene_names = gene_field.split(",")

# %% [markdown]
#                 if any(gene in genes for gene in gene_names):
#                     fout.write(line)

# %%

# %%
for idx, row in sa.iterrows():
    sample = row.pid
    # Try both file paths
    vcf_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/standard_vcf_indel/indel_H021-{sample}.standard_indel.vcf.gz"
    if not os.path.exists(vcf_path):
        vcf_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/standard_vcf_indel/indel_P021-{sample}.standard_indel.vcf.gz"
        if not os.path.exists(vcf_path):
            print('No VCF for', sample)
            continue
    seq_type = row.seq_type

    temp_vcf = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/predisposition_variants/{seq_type}/indel/indel_{sample}.temp.vcf"
    final_vcf = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/predisposition_variants/{seq_type}/indel/{sample}.standard_indel_predisposition.vcf.gz"

    # Filter and write to temporary VCF
    with gzip.open(vcf_path, "rt") as fin, open(temp_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.strip().split("\t")
            info = fields[7]
            if "GENE=" not in info:
                continue
            gene_field = info.split("GENE=")[1].split(";")[0]
            gene_names = gene_field.split(",")
            if any(gene in genes for gene in gene_names):
                fout.write(line)

    # Compress with bgzip
    os.system(f"bgzip -f -c {temp_vcf} > {final_vcf}")
    os.system(f"tabix -p vcf {final_vcf}")
    # Remove temporary file
    os.remove(temp_vcf)

