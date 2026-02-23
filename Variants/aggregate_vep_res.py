import pandas as pd
import os
import glob
import polars as pl

def create_all_vars_unique(vep_res_files, base_path):
    vep_list = []
    
    # 1. Define Severity Mapping
    # Higher number = More severe
    impact_map = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
    
    cols_to_keep = [
                    "#Uploaded_variation", "Location", "Allele", "Gene",
                    "sampleID", "IMPACT", "Consequence", "SYMBOL", "PHENO", "STRAND", "Feature_type", "Feature",
                    "am_class", "am_pathogenicity", "LoF", "LoF_filter", "AF", "gnomADe_AF", "gnomADg_AF",
                    "SpliceAI_pred", "CADD_PHRED","CADD_RAW"
    ]

    for i, vep_file in enumerate(vep_res_files):
        if i % 500 == 2:
            print(f"Processing file {i}...")

        file_path = os.path.join(base_path, vep_file)
        
        try:
            # Read VEP TSV (skipping ## header lines)
            df = pl.read_csv(file_path, separator="\t", comment_prefix="##", infer_schema_length=100000)
            df = df.to_pandas()
            
            # Extract Sample ID
            sample_id = vep_file.split(".")[0].split("-")[1]
            df["sampleID"] = sample_id

            # 2. Assign Rank for Sorting
            # We map the IMPACT string to a numeric value
            df['impact_rank'] = df['IMPACT'].map(impact_map).fillna(0)

            # 3. Sort by Variant and Impact
            # We sort descending so the '4' (HIGH) comes before '1' (MODIFIER)
            df = df.sort_values(
                by=["Gene", "#Uploaded_variation", "impact_rank"], 
                ascending=[True, True, False]
            )

            # 4. Deduplicate: Keep only the first (most severe) row per variant
            df_unique = df.drop_duplicates(subset=["#Uploaded_variation", "Gene"])

            # 5. Filter columns and store
            existing_cols = [c for c in cols_to_keep if c in df_unique.columns]
            vep_list.append(df_unique[existing_cols])

        except Exception as e:
            print(f"Error in {vep_file}: {e}")

    # Combine all samples into one massive table
    if not vep_list:
        return pd.DataFrame()
        
    return pd.concat(vep_list, ignore_index=True)

# Usage
vep_indel_rare_dir = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_indel_rare/"
vep_files = [f for f in os.listdir(vep_indel_rare_dir) if f.endswith(".gz")]
final_df = create_all_vars_unique(vep_files, vep_indel_rare_dir)

var_type = "indel"
final_df.to_csv(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated/vep_res_rare_{var_type}_all_aggregated.tsv", sep="\t")

