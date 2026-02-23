import pandas as pd
import os
import polars as pl
import numpy as np
import sys


def compute_location(row):
    chrom = row["#CHROM"]
    pos = int(row["POS"])
    ref = row["REF"]
    alt = row["ALT"]

    len_ref = len(ref)
    len_alt = len(alt)

    # --- SNV ---
    if len_ref == 1 and len_alt == 1:
        return f"{chrom}:{pos}"

    # --- Deletion (REF > ALT, ALT is anchor only) ---
    if len_ref > 1 and len_alt == 1:
        start = pos + 1
        end = pos + len_ref - 1
        if start == end:
            return f"{chrom}:{start}"
        return f"{chrom}:{start}-{end}"

    # --- Insertion (ALT > REF, REF is anchor only) ---
    if len_ref == 1 and len_alt > 1:
        # VEP shows insertions as POS-POS+1
        return f"{chrom}:{pos}-{pos+1}"

    # --- Complex Indel ---
    span = max(len_ref, len_alt)
    start = pos
    end = pos + span - 1
    if start == end:
        return f"{chrom}:{start}"
    return f"{chrom}:{start}-{end}"

def extract_vcf_data_universal(row):
    # 1. Mandatory columns: index 7 is INFO, index 8 is FORMAT
    info_str = row.iloc[7]
    format_str = row.iloc[8]
    
    # 2. Sample data: index 9 is the FIRST sample column
    # (The one with the long path: /icgc/dkfzlsdf/...)
    sample_str = str(row.iloc[9])
    
    # Parse INFO into a dictionary
    info_parts = dict(item.split("=") for item in info_str.split(";") if "=" in item)
    
    # Parse FORMAT and SAMPLE into a dictionary
    fmt_keys = format_str.split(':')
    sample_values = sample_str.split(':')
    sample_data = dict(zip(fmt_keys, sample_values))
    
    nr, nv = 0.0, 0.0
    
    try:
        # LOGIC A: Check for BCFtools/Samtools style (DP4 in INFO)
        if 'DP4' in info_parts:
            dp4 = [float(x) for x in info_parts['DP4'].split(',')]
            nr = sum(dp4)
            nv = dp4[2] + dp4[3]
        
        # LOGIC B: Check for Platypus/Standard style (NV/NR in FORMAT)
        elif 'NV' in sample_data and 'NR' in sample_data:
            nr = float(sample_data['NR'])
            nv = float(sample_data['NV'])
            
        # LOGIC C: Check for basic DP/DV style
        elif 'DP' in sample_data and 'DV' in sample_data:
            nr = float(sample_data['DP'])
            nv = float(sample_data['DV'])

        vaf = nv / nr if nr > 0 else 0.0
        return nr, nv, vaf
        
    except (ValueError, TypeError, IndexError):
        return 0.0, 0.0, 0.0


def extract_single_vaf(row, sample_id):
    # Split the FORMAT definition (e.g., GT:GL:GOF:GQ:NR:NV)
    fmt_keys = row['FORMAT'].split(':')
    # Split the actual sample values (e.g., 0/1:-33.6:15:99:18:7)
    sample_values = str(row[sample_id]).split(':')
    
    # Map them together
    data = dict(zip(fmt_keys, sample_values))
    
    try:
        # Platypus uses NR for total depth and NV for variant depth
        nr = float(data.get('NR', 0))
        nv = float(data.get('NV', 0))
        return nr, nv, nv / nr if nr > 0 else 0.0
    except (ValueError, TypeError):
        return 0.0, 0, 0
        
chunk_idx = int(sys.argv[1])
total_chunks = int(sys.argv[2])
print(chunk_idx, total_chunks)

var_type = "snv"
#variants_unique = pl.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_snv_all_aggregated_unique_variant_type.tsv", separator="\t")
variants_unique = pl.read_csv(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_{var_type}_all_aggregated_unique_variants.tsv", separator="\t", infer_schema_length=1000000)
print(variants_unique.shape)

sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")


vep_res_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_{var_type}_rare/"

keys = ["#Uploaded_variation", "sampleID", "Location", "Allele", "Feature"]
vep_dfs = []
i = 0

all_files = sorted([f for f in os.listdir(vep_res_path) if not f.endswith("warnings.txt")])
chunk_size = int(np.ceil(len(all_files) / total_chunks))
my_files = all_files[chunk_idx * chunk_size : (chunk_idx + 1) * chunk_size]

sample_ids = []

for vep_res in my_files:
    i += 1
    # read in vep_res
    if not vep_res.endswith("warnings.txt"):
        sample_id = vep_res.split("-")[1].split(".")[0]
        sample_ids.append(sample_id)
        try:
            orig_vcf_path = sa.loc[sa["pid"] == sample_id, var_type + "_vcf"].values[0]
            if var_type == "indel":
                orig_vcf_path = sa.loc[sa["pid"] == sample_id, var_type + "_vcf"].values[0].split("/")[-1].split('.')[0]
                orig_vcf_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/normalized_vcf_indel_rare/{orig_vcf_path}.normalized_standard_indel_rare_filtered.vcf.gz"
        except:
            continue

        original_vcf = pl.read_csv(
            orig_vcf_path,
            separator="\t",
            comment_prefix="##",
            has_header=False,
        )
        original_vcf.columns = original_vcf.row(0)
        original_vcf = original_vcf.slice(1)

        if var_type == "snv":
            original_vcf = original_vcf.with_columns(
                pl.nth(11)
                  .str.extract(r"VAF=([^;]+)", 1)
                  .cast(pl.Float64, strict=False)
                  .alias("VAF")
            )
            pdf = original_vcf.to_pandas()
            pdf[['NR', 'NV', 'VAF_samtools']] = pdf.apply(
                extract_vcf_data_universal, 
                axis=1, 
                result_type='expand'
            )
            original_vcf = pl.from_pandas(pdf)   
            
        else:
            pdf = original_vcf.to_pandas()
            pdf["VAF_samtools"] = 0
            pdf[['NR', 'NV', 'VAF']] = pdf.apply(
                extract_single_vaf, 
                axis=1, 
                args=(sample_id,), 
                result_type='expand'
            )
            original_vcf = pl.from_pandas(pdf)   
            
        original_vcf = original_vcf.with_columns(
            pl.col("ALT")  
            .str.split(",")
        ).explode("ALT")

        
        df = pl.read_csv(
            vep_res_path + vep_res,
            separator="\t",
            comment_prefix="##",
            has_header=False,
        )
        


        df.columns = df.row(0)
        df = df.slice(1)


    
        
        df = df.unique(subset=["#Uploaded_variation", "Allele", "Feature"])
        
        df_small = df.select([
            "#Uploaded_variation",
            "Location",
            "Allele",
            "Feature",
            "gnomADe_AF",
            "gnomADg_AF",
        ]).with_columns([
            pl.col("gnomADe_AF").cast(pl.Float64, strict=False).fill_null(0.0),
            pl.col("gnomADg_AF").cast(pl.Float64, strict=False).fill_null(0.0),
            pl.lit(sample_id).alias("sampleID")
        ])
        df_small = df_small.with_columns([
        pl.col("Location").str.split(":").list.get(0).alias("#CHROM"),
        pl.col("Location").str.split(":").list.get(1).alias("POS"),
        ])
        
        original_vcf = original_vcf.with_columns(
            Location = pl.struct(["#CHROM", "POS", "REF", "ALT"])
                .map_elements(lambda x: compute_location(x), return_dtype=pl.String)
            )
        
        df_small = df_small.join(
            original_vcf.select(pl.col("Location", "VAF", "NR", "NV", "VAF_samtools")),
            how="left",
            on="Location",
        )
        vep_dfs.append(df_small)


vep_all = pl.concat(vep_dfs, rechunk=True)

variants_unique =  variants_unique.filter(pl.col("sampleID").is_in(sample_ids))
# merge once with variants
variants_unique = variants_unique.join(
    vep_all,
    on=keys,
    how="left"
)

variants_unique = variants_unique.with_columns(
    pl.col("gnomADg_AF")
      .cast(pl.Float64, strict=False)
      .fill_null(0.0)
)

variants_unique = variants_unique.with_columns(
    pl.col("gnomADe_AF")
      .cast(pl.Float64, strict=False)
      .fill_null(0.0)
)

print(variants_unique.shape)
output_file = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/snv_chunk/vep_res_rare_{var_type}_chunk_{chunk_idx}_all_aggregated_unique_variants_VAF.tsv"
variants_unique.write_csv(output_file, separator="\t")
