import pandas as pd
from itables import init_notebook_mode
import gzip
import os
import numpy as np
init_notebook_mode(all_interactive=True)



import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *


# py_or_res_all = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants.csv")
# #py_or_res_all = pd.merge(py_or_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID", "ROLE_IN_CANCER", "predisposition_gene"]], on="geneID", how="left")

# py_or_res_all["chrom_snv"] =  py_or_res_all["Location_snv"].str.split(":").str[0]
# py_or_res_all["pos_snv"] = py_or_res_all["Location_snv"].str.split(":").str[1]

#samples = np.unique(py_or_res_all.sampleID)

sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")

samples = np.unique(sa.pid)
len(samples)


## load header
input_vcf_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_indel_rare/"
sampleID = "9PJRHA"

complete_path = f"{input_vcf_path}indel_H021-{sampleID}.vep_res_indel_rare_filtered.vcf.gz"


with gzip.open(complete_path, 'rt', encoding='utf-8') as f:
    for line in f:
        if line.startswith('#Uploaded_variation'):
            header = line.strip().split('\t')
            break


dfs = []
for sampleID in samples:
    #print(sampleID)
    complete_path = f"{input_vcf_path}indel_H021-{sampleID}.vep_res_indel_rare_filtered.vcf.gz"
    
    # Read the file into a DataFrame
    try:
        sample_vcf = pd.read_csv(complete_path, sep='\t', comment='#', names=header)
    except:
        try:
            complete_path = f"{input_vcf_path}indel_P021-{sampleID}.vep_res_indel_rare_filtered.vcf.gz"
            sample_vcf = pd.read_csv(complete_path, sep='\t', comment='#', names=header)
        except:
            print(sampleID , " not done yet")
            continue

    sample_vcf["sampleID"] = sampleID
    sample_vcf["chrom_indel"] =  sample_vcf["Location"].str.split(":").str[0]
    sample_vcf["pos_indel"] = sample_vcf["Location"].str.split(":").str[1]
    
    sample_vcf = sample_vcf.rename(columns={"#Uploaded_variation": "#Uploaded_variation_indel"})
    print(sample_vcf.shape, sampleID)
    
    #merged = py_or_res_all.merge(sample_vcf[["#Uploaded_variation_indel", "Location", "chrom_indel", "pos_indel", "REF_ALLELE", "Gene", "Allele", "Feature", "sampleID", "Existing_variation"]], left_on=["Location_indel", "Gene", "Allele_indel", "sampleID"], right_on=["Location", "Gene", "Allele", "sampleID"], )
    
    
    impact_levels = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
    sample_vcf["IMPACT"] = pd.Categorical(
        sample_vcf["IMPACT"],
        categories=impact_levels,
        ordered=True
    )
    
    # Sort like setorder(Gene, sampleID, IMPACT)
    sample_vcf = sample_vcf.sort_values(
        by=["Gene", "sampleID", "IMPACT"]
    )
    
    
    # Keep first row per Gene + sampleID (after sorting)
    sample_vcf = sample_vcf.drop_duplicates(
        subset=["Gene", "chrom_indel", "pos_indel", "REF_ALLELE", "Allele"],
        keep="first"
    ).drop(columns="Feature")
    
    dfs.append(sample_vcf[["#Uploaded_variation_indel", "chrom_indel", "pos_indel", "REF_ALLELE", "Allele", "Location", "Gene", "Existing_variation"]].rename(columns={"Allele": "alt_indel", "REF_ALLELE": "ref_indel", "Location": "Location_indel", "Gene": "Gene_indel"}))


final_df = pd.concat(dfs)
final_df = final_df.drop_duplicates(subset=["chrom_indel", "pos_indel", "alt_indel", "ref_indel", "Gene_indel"])

final_df.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/unique_indes.csv")


