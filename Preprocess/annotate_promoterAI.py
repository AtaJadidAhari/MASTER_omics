import pandas as pd
import pickle as pkl
from collections import defaultdict


variants= pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38.tsv", sep="\t")
#variants = variants[variants["sampleID"] == "11FX5L"]

print(variants.shape)

variants["hg38_pos"] = variants["hg38_pos"].astype("Int64")

with open('/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/alt_vs_rs.pkl', 'rb') as handle:
    b = pkl.load(handle)
rows = []

for var_id, mapping in b.items():
    ref_to_alts = defaultdict(list)
    for alt, ref in mapping.items():
        ref_to_alts[ref].append(alt)

    for ref, alts in ref_to_alts.items():
        rows.append({
            "var_id": var_id,
            "ref": ref,
            "alts": alts
        })

df = pd.DataFrame(rows)
df_exploded = df.explode("alts").rename(columns={"alts": "alt"})

variants = variants.merge(df_exploded, left_on=["ID", "Allele"], right_on=["var_id", "alt"], how="left").drop(columns="alt")
variants.loc[~variants["ref"].notna(),"ref"] = variants.loc[~variants["ref"].notna()]["#Uploaded_variation"].str.split("_").str[2].str.split("/").str[0]

#variants["ref"] = variants["#Uploaded_variation"].str.split("_").str[2].str.split("/").str[0]

promoter_ai = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/promoterAI_tss500.parquet")
promoter_ai["Feature"] = promoter_ai["transcript_id"].str.split(".").str[0]

# some genes have multiple transcript. take onky the most sever transcript
promoter_ai["promoterAI"] = pd.to_numeric(promoter_ai["promoterAI"])
df_sorted = promoter_ai.sort_values(
    by="promoterAI",
    key=lambda x: x.abs(),
    ascending=False
)
promoter_ai = df_sorted.drop_duplicates(
    subset=["chrom", "pos", "ref", "alt", "gene_id", "strand"],
    keep="first"
)

print(promoter_ai.shape)



merged = pd.merge(variants, promoter_ai[["pos", "chrom", "alt", "ref", "promoterAI", "strand", "gene_id"]], how="left", left_on = ["hg38_chrom", "hg38_pos", "Allele", "ref", "Gene", "STRAND"], right_on=["chrom", "pos", "alt", "ref", "gene_id", "strand"])

merged = merged.drop(["pos", "chrom", "alt"], axis=1)

print(merged.shape)


merged.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38_promoterAI.tsv", sep="\t", index=None)