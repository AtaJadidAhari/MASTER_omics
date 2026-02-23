import pandas as pd
import sys
import os
import pickle

import sys

sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from Preprocess.VEP_SORanking import *



var_type = sys.argv[1]   # "snv" or "indel"
print(var_type)

def get_most_severe(consequence_str):
    """Replicates your R sapply(strsplit(...)) logic"""
    if pd.isna(consequence_str) or consequence_str == "":
        return None
    # Split by comma, trim whitespace
    terms = [t.strip() for t in str(consequence_str).split(',')]
    # Return the term with the lowest rank (most severe)
    return min(terms, key=lambda x: severity_rank.get(x, 999))

def resolve_ref(row):
    var = row["#Uploaded_variation"]
    if not isinstance(var, str) or not var.startswith("rs"):
        return row["ref"]
    return get_ref_from_rsid_and_alt_hg19(var, row["Allele"])


CACHE_FILE = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/alt_vs_rs.pkl"
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, "rb") as f:
        _rs_cache = pickle.load(f)

snv_vep_res = pd.read_csv(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_{var_type}_all_aggregated.tsv", sep="\t")

snv_vep_res['Consequence_most_severe'] = snv_vep_res['Consequence'].apply(get_most_severe)

print(snv_vep_res.shape, "initial vars")


try:
    # add variant type
    snv_variant_type = pd.read_csv(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/variant_types/{var_type}_variant_types.tsv", sep="\t")
    snv_vep_res = pd.merge(snv_vep_res, snv_variant_type, on=["sampleID", "Location"], how="left")
except:
    pass
print(snv_vep_res.shape, "after adding varinat type")

# 1. Define Impact Priority (HIGH > MODERATE > LOW > MODIFIER)
# This ensures that when we pick 'one' transcript, we pick the most important one.
impact_map = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
snv_vep_res['impact_rank'] = snv_vep_res['IMPACT'].map(impact_map)

snv_vep_res['germline_rank'] = (snv_vep_res['ANNOTATION_control'] == 'germline').astype(int)

# 2. Update consequence if needed (replicating your R logic)
mask = (snv_vep_res['Consequence_most_severe'] == 'missense_variant') & (snv_vep_res['am_class'] != '-')
snv_vep_res.loc[mask, 'Consequence_most_severe'] = snv_vep_res['am_class']

# 3. Sort the data
# We sort by Gene/Sample/Location, then by Impact Rank (descending)
# This puts the "Best" transcript for each specific variant at the top.
snv_vep_res = snv_vep_res[snv_vep_res["Gene"] != "-"]
snv_vep_res = snv_vep_res.sort_values(
    by=['sampleID', 'Gene', '#Uploaded_variation', 'impact_rank', "germline_rank"], 
    ascending=[True, True, True, False, False]
)

# 4. Deduplicate
# IMPORTANT: Include 'Location' (or 'Start') so we keep multiple variants per gene!
df_final = snv_vep_res.drop_duplicates(
    subset=['sampleID', 'Gene', '#Uploaded_variation']
)

print(df_final.shape, "after uniquing for #Uploaded_variation")
# Clean up temporary column
df_final = df_final.drop(columns=['impact_rank', 'germline_rank'])


try:
    df_final.loc[df_final['#Uploaded_variation'].str.startswith("rs"), "ref"] = variants.loc[variants['ID'].str.startswith("rs")].apply(resolve_ref, axis=1)
    
    df_final[df_final["ref"].notna()].shape
except:
    pass
print(df_final.shape, "after adding ref")



df_final.to_csv(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_{var_type}_all_aggregated_unique_variants.tsv", sep="\t", index=False)