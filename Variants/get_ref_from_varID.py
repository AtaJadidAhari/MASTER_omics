import requests
import pandas as pd
import pickle
import os


variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38.tsv", sep="\t")

variants["ref"] = variants["#Uploaded_variation"].str.split("_").str[2].str.split("/").str[0]

variants[variants["ref"].notna()].shape

CACHE_FILE = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/alt_vs_rs.pkl"
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, "rb") as f:
        _rs_cache = pickle.load(f)
else:
    _rs_cache = {}

_new_since_checkpoint = 0
CHECKPOINT_EVERY = 1000

def save_cache():
    with open(CACHE_FILE, "wb") as f:
        pickle.dump(_rs_cache, f)
    

def get_ref_from_rsid_and_alt_hg19(rsid: str, alt: str):
    global _new_since_checkpoint
    global CHECKPOINT_EVERY
    # return cached result if available
    if rsid in _rs_cache:
        return _rs_cache[rsid].get(alt)

    url = f"https://grch37.rest.ensembl.org/variation/human/{rsid}"

    try:
        r = requests.get(
            url,
            headers={"Content-Type": "application/json"},
            timeout=10
        )
        r.raise_for_status()
        data = r.json()
    except requests.RequestException:
        # cache failure to avoid retrying forever
        _rs_cache[rsid] = {}
        return None

    alt_to_ref = {}

    for m in data.get("mappings", []):
        alleles = m["allele_string"].split("/")
        ref = alleles[0]
        for a in alleles[1:]:
            alt_to_ref[a] = ref

    _rs_cache[rsid] = alt_to_ref

    _new_since_checkpoint += 1

    if _new_since_checkpoint >= CHECKPOINT_EVERY:
        save_cache()
        _new_since_checkpoint = 0
    
    return alt_to_ref.get(alt)



def resolve_ref(row):
    var = row["#Uploaded_variation"]
    if not isinstance(var, str) or not var.startswith("rs"):
        return row["ref"]
    return get_ref_from_rsid_and_alt_hg19(var, row["Allele"])


# apply to whole dataframe
variants.loc[variants['ID'].str.startswith("rs"), "ref"] = variants.loc[variants['ID'].str.startswith("rs")].apply(resolve_ref, axis=1)

variants[variants["ref"].notna()].shape


variants.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38_ref.tsv")

with open(CACHE_FILE, "wb") as f:
    pickle.dump(_rs_cache, f)


