
# import pandas as pd
# import polars as pl
# import os
# import io
# import pandas as pd
# import plotnine as pn
# import numpy as np
# from scipy.stats import beta
# import itables
# from itables import init_notebook_mode
# import polars as pl
# init_notebook_mode(all_interactive=True)
# import itables.options as opt

# from scipy.stats import mannwhitneyu



# import sys
# sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

# from utils.load_gtf_cgc_dresden import *
# from utils.util_functions import *
# from ProteinExpression.load_pr_data import *



# sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")




# pids = sa["nct_pid"].unique()
# all_dfs = []

# for pid in pids:
#     found = False
#     for seq_type in ["whole_genome_sequencing", "exon_sequencing"]:
        
#         base_path = f"/omics/odcf/analysis/hipo/hipo_021/GermlineAnalysis/data_object_master_germline/sequencing/{seq_type}/results_per_pid/"
#         pid_parts = pid.split("_")
#         base_sample = f"H021-{pid_parts[0]}"
#         variant_suffix = pid_parts[1]
        
        
#         try:
#             sample_path = (
#                     f"{base_sample}/germline_smallVariants.{variant_suffix}.blood/v5/"
#                     f"clinvar_2023-09/smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                 )
                
#             if not os.path.exists(f"{base_path}/{sample_path}"):
#                 # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                 # (Skipping '02' since we already checked it)
#                 suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                
#                 found = False
#                 for suffix in suffixes:
#                     # Swap out the .buffy_coat02/ with the current suffix choice
#                     potential_path = sample_path.replace(".blood/", f".blood{suffix}/")
#                     if os.path.exists(f"{base_path}/{potential_path}"):
#                         print(potential_path)
#                         sample_path = potential_path
#                         found = True
#                         break
                        
#             df = pl.scan_csv(f"{base_path}/{sample_path}",
#                              separator="\t", comment_prefix="##",  schema_overrides={
#                                                 "QUAL": pl.String,
#                                                 "#CHROM": pl.String,
#                                                 "RNA_VARIANT_AF": pl.Float64,
#                                             }, 
#                                 null_values=["NA", "."],
#             ).filter((pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"])) | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]) ) | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))).collect(engine="streaming").to_pandas()
#             if df.empty:
#                 # Create a single row filled with NaN for all existing columns
#                 df = pd.DataFrame(index=[0], columns=df.columns)
#             df.columns.values[9:11] = ["sample", "tumor"]
#             df["sampleID"] = pid.split("_")[0]
#             df["vcf_path"] = f"{base_path}/{sample_path}"
#             all_dfs.append(df)
#             print("first try", pid.split("_")[0])
#             break
#         except:
#             try:
#                 sample_path = (
#                     f"{base_sample}/germline_smallVariants.{variant_suffix}.buffy-coat/v5/"
#                     f"clinvar_2023-09/smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                 )
                
#                 if not os.path.exists(f"{base_path}/{sample_path}"):
#                     # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                     # (Skipping '02' since we already checked it)
#                     suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                    
#                     found = False
#                     for suffix in suffixes:
#                         # Swap out the .buffy_coat02/ with the current suffix choice
#                         potential_path = sample_path.replace(".buffy-coat/", f".buffy-coat{suffix}/")
#                         if os.path.exists(potential_path):
#                             sample_path = potential_path
#                             found = True
#                             break
#                     if found == False:
#                         sample_path = sample_path.replace(".buffy-coat", f".buffy_coat")
#                         for suffix in suffixes:
#                             # Swap out the .buffy_coat02/ with the current suffix choice
#                             potential_path = sample_path.replace(".buffy_coat/", f".buffy_coat{suffix}/")
#                             if os.path.exists(f"{base_path}/{potential_path}"):
#                                 sample_path = potential_path
#                                 found = True
#                                 break
#                 df = pl.scan_csv(f"{base_path}/{sample_path}",
#                                  separator="\t", comment_prefix="##",  schema_overrides={
#                                                     "QUAL": pl.String,
#                                                     "#CHROM": pl.String,
#                                                     "RNA_VARIANT_AF": pl.Float64,
#                                                 }, 
#                                     null_values=["NA", "."],
#                 ).filter((pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"])) | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]) ) | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))).collect(engine="streaming").to_pandas()
#                 if df.empty:
#                     # Create a single row filled with NaN for all existing columns
#                     df = pd.DataFrame(index=[0], columns=df.columns)
#                 df.columns.values[9:11] = ["sample", "tumor"]
#                 df["sampleID"] = pid.split("_")[0]
#                 df["vcf_path"] = f"{base_path}/{sample_path}"
#                 all_dfs.append(df)
#                 print("second try", pid.split("_")[0])
#                 break
            
#             except:
#                 try:
#                     sample_path = (
#                             f"{base_sample}/germline_smallVariants.{variant_suffix}.blood/v4.0.2/"
#                             f"smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                         )
                        
#                     if not os.path.exists(f"{base_path}/{sample_path}"):
#                         # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                         # (Skipping '02' since we already checked it)
#                         suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                        
#                         found = False
#                         for suffix in suffixes:
#                             # Swap out the .buffy_coat02/ with the current suffix choice
#                             potential_path = sample_path.replace(".blood/", f".blood{suffix}/")
#                             if os.path.exists(f"{base_path}/{potential_path}"):
#                                 print(potential_path)
#                                 sample_path = potential_path
#                                 found = True
#                                 break
                                
#                     df = pl.scan_csv(f"{base_path}/{sample_path}",
#                                      separator="\t", comment_prefix="##",  schema_overrides={
#                                                         "QUAL": pl.String,
#                                                         "#CHROM": pl.String,
#                                                         "RNA_VARIANT_AF": pl.Float64,
#                                                     }, 
#                                         null_values=["NA", "."],
#                     ).filter((pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"])) | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]) ) | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))).collect(engine="streaming").to_pandas()
#                     if df.empty:
#                         # Create a single row filled with NaN for all existing columns
#                         df = pd.DataFrame(index=[0], columns=df.columns)
#                     df.columns.values[9:11] = ["sample", "tumor"]
#                     df["sampleID"] = pid.split("_")[0]
#                     df["vcf_path"] = f"{base_path}/{sample_path}"
#                     all_dfs.append(df)
#                     break
                   
#                 except:
#                         try:
#                             sample_path = (
#                                     f"{base_sample}/germline_smallVariants.{variant_suffix}.buffy-coat/v4.0.2/"
#                                     f"smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                                 )
                                
#                             if not os.path.exists(f"{base_path}/{sample_path}"):
#                                 # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                                 # (Skipping '02' since we already checked it)
#                                 suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                                
#                                 found = False
#                                 for suffix in suffixes:
#                                     # Swap out the .buffy_coat02/ with the current suffix choice
#                                     potential_path = sample_path.replace(".buffy-coat/", f".buffy-coat{suffix}/")
#                                     if os.path.exists(f"{base_path}/{potential_path}"):
#                                         print(potential_path)
#                                         sample_path = potential_path
#                                         found = True
#                                         break
#                                 if found == False:
#                                     sample_path = sample_path.replace(".buffy-coat", f".buffy_coat")
#                                     for suffix in suffixes:
#                                         # Swap out the .buffy_coat02/ with the current suffix choice
#                                         potential_path = sample_path.replace(".buffy_coat/", f".buffy_coat{suffix}/")
#                                         if os.path.exists(f"{base_path}/{potential_path}"):
#                                             sample_path = potential_path
#                                             found = True
#                                             break
                                        
#                             df = pl.scan_csv(f"{base_path}/{sample_path}",
#                                              separator="\t", comment_prefix="##",  schema_overrides={
#                                                                 "QUAL": pl.String,
#                                                                 "#CHROM": pl.String,
#                                                                 "RNA_VARIANT_AF": pl.Float64,
#                                                             }, 
#                                                 null_values=["NA", "."],
#                             ).filter((pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"])) | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]) ) | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))).collect(engine="streaming").to_pandas()
#                             if df.empty:
#                                 # Create a single row filled with NaN for all existing columns
#                                 df = pd.DataFrame(index=[0], columns=df.columns)
#                             df.columns.values[9:11] = ["sample", "tumor"]
#                             df["sampleID"] = pid.split("_")[0]
#                             df["vcf_path"] = f"{base_path}/{sample_path}"
#                             all_dfs.append(df)
#                             break
                       
#                         except:
#                             try:
#                                 sample_path = (
#                                     f"{base_sample}/germline_smallVariants.{variant_suffix}.blood/v3.3.1/"
#                                     f"smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                                 )
                                
#                                 if not os.path.exists(f"{base_path}/{sample_path}"):
#                                     # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                                     # (Skipping '02' since we already checked it)
#                                     suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                                    
#                                     found = False
#                                     for suffix in suffixes:
#                                         # Swap out the .buffy_coat02/ with the current suffix choice
#                                         potential_path = sample_path.replace(".blood/", f".blood{suffix}/")
#                                         if os.path.exists(potential_path):
#                                             sample_path = potential_path
#                                             found = True
#                                             break
#                                 df = pl.scan_csv(f"{base_path}/{sample_path}",
#                                                  separator="\t", comment_prefix="##",  schema_overrides={
#                                                                     "QUAL": pl.String,
#                                                                     "#CHROM": pl.String,
#                                                                     "RNA_VARIANT_AF": pl.Float64,
#                                                                 }, 
#                                                     null_values=["NA", "."],
#                                 ).filter((pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"])) | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]) ) | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))).collect(engine="streaming").to_pandas()
#                                 if df.empty:
#                                         # Create a single row filled with NaN for all existing columns
#                                         df = pd.DataFrame(index=[0], columns=df.columns)
#                                 df.columns.values[9:11] = ["sample", "tumor"]
#                                 df["sampleID"] = pid.split("_")[0]
#                                 df["vcf_path"] = f"{base_path}/{sample_path}"
#                                 all_dfs.append(df)
#                                 print("third try", pid.split("_")[0])
#                                 break
#                             except:
#                                 try:
#                                     sample_path = (
#                                         f"{base_sample}/germline_smallVariants.{variant_suffix}.buffy-coat/v3.3.1/"
#                                         f"smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                                     )
                                    
#                                     if not os.path.exists(f"{base_path}/{sample_path}"):
#                                         # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                                         # (Skipping '02' since we already checked it)
#                                         suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                                        
#                                         found = False
#                                         for suffix in suffixes:
#                                             # Swap out the .buffy_coat02/ with the current suffix choice
#                                             potential_path = sample_path.replace(".buffy-coat/", f".buffy-coat{suffix}/")
#                                             if os.path.exists(potential_path):
#                                                 sample_path = potential_path
#                                                 found = True
#                                                 break
#                                         if found == False:
#                                             sample_path = sample_path.replace(".buffy-coat", f".buffy_coat")
#                                             for suffix in suffixes:
#                                                 # Swap out the .buffy_coat02/ with the current suffix choice
#                                                 potential_path = sample_path.replace(".buffy_coat/", f".buffy_coat{suffix}/")
#                                                 if os.path.exists(f"{base_path}/{potential_path}"):
#                                                     sample_path = potential_path
#                                                     found = True
#                                                     break
#                                     df = pl.scan_csv(f"{base_path}/{sample_path}",
#                                                      separator="\t", comment_prefix="##",  schema_overrides={
#                                                                         "QUAL": pl.String,
#                                                                         "#CHROM": pl.String,
#                                                                         "RNA_VARIANT_AF": pl.Float64,
#                                                                     }, 
#                                                         null_values=["NA", "."],
#                                     ).filter((pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"])) | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]) ) | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))).collect(engine="streaming").to_pandas()
#                                     if df.empty:
#                                         # Create a single row filled with NaN for all existing columns
#                                         df = pd.DataFrame(index=[0], columns=df.columns)
#                                     df.columns.values[9:11] = ["sample", "tumor"]
#                                     df["sampleID"] = pid.split("_")[0]
#                                     df["vcf_path"] = f"{base_path}/{sample_path}"
#                                     all_dfs.append(df)
#                                     print("third try", pid.split("_")[0])
#                                     break
#                                 except:
#                                     try:
#                                         sample_path = (
#                                                 f"{base_sample}/germline_smallVariants.{variant_suffix}.blood/v3.0.0/"
#                                                 f"smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                                             )
                                            
#                                         if not os.path.exists(f"{base_path}/{sample_path}"):
#                                             # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                                             # (Skipping '02' since we already checked it)
#                                             suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                                            
#                                             found = False
#                                             for suffix in suffixes:
#                                                 # Swap out the .buffy_coat02/ with the current suffix choice
#                                                 potential_path = sample_path.replace(".blood/", f".blood{suffix}/")
#                                                 if os.path.exists(f"{base_path}/{potential_path}"):
#                                                     sample_path = potential_path
#                                                     found = True
#                                                     break
#                                         df = (
#                                             pl.scan_csv(
#                                                 f"{base_path}/{sample_path}",
#                                                 separator="\t",
#                                                 comment_prefix="##",
#                                                 schema_overrides={
#                                                     "QUAL": pl.String,
#                                                     "#CHROM": pl.String,
#                                                     "RNA_VARIANT_AF": pl.Float64,
#                                                 },
#                                                 null_values=["NA", "."],
#                                             )
#                                             .filter(
#                                                 (pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"]))
#                                                 | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]))
#                                                 | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))
#                                             )
#                                             .collect(engine="streaming")
#                                             .to_pandas()
#                                         )
#                                         if df.empty:
#                                             # Create a single row filled with NaN for all existing columns
#                                             df = pd.DataFrame(index=[0], columns=df.columns)
#                                         df.columns.values[9:11] = ["sample", "tumor"]
#                                         df["sampleID"] = pid.split("_")[0]
#                                         df["vcf_path"] = f"{base_path}/{sample_path}"
#                                         all_dfs.append(df)
#                                         break
#                                     except:
#                                         try:
                                            
#                                             sample_path = (
#                                                     f"{base_sample}/germline_smallVariants.{variant_suffix}.buffy-coat/v3.0.0/"
#                                                     f"smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
#                                                 )
                                                
#                                             if not os.path.exists(f"{base_path}/{sample_path}"):
#                                                 # This list generates: ['', '01', '03', '04', '05', '06', '07', '08'] 
#                                                 # (Skipping '02' since we already checked it)
#                                                 suffixes = [""] + [f"{i:02d}" for i in range(1, 9)]
                                                
#                                                 found = False
#                                                 for suffix in suffixes:
#                                                     # Swap out the .buffy_coat02/ with the current suffix choice
#                                                     potential_path = sample_path.replace(".buffy-coat/", f".buffy-coat{suffix}/")
#                                                     if os.path.exists(f"{base_path}/{potential_path}"):
#                                                         sample_path = potential_path
#                                                         found = True
#                                                         break
#                                                 if found == False:
#                                                     sample_path = sample_path.replace(".buffy-coat", f".buffy_coat")
#                                                     for suffix in suffixes:
#                                                         # Swap out the .buffy_coat02/ with the current suffix choice
#                                                         potential_path = sample_path.replace(".buffy_coat/", f".buffy_coat{suffix}/")
#                                                         if os.path.exists(f"{base_path}/{potential_path}"):
#                                                             sample_path = potential_path
#                                                             found = True
#                                                             break
                                                        
#                                             df = pl.scan_csv(f"{base_path}/{sample_path}",
#                                                              separator="\t", comment_prefix="##",  schema_overrides={
#                                                                                 "QUAL": pl.String,
#                                                                                 "#CHROM": pl.String,
#                                                                                 "RNA_VARIANT_AF": pl.Float64,
#                                                                             }, 
#                                                                 null_values=["NA", "."],
#                                             ).filter((pl.col("ACMG_Classification").is_in(["Pathogenic", "Likely Pathogenic"])) | (pl.col("CharGer_Classification").is_in(["Pathogenic", "Likely Pathogenic"]) ) | (pl.col("ClinVar_Pathogenicity").is_in(["Pathogenic", "Likely Pathogenic"]))).collect(engine="streaming").to_pandas()
#                                             if df.empty:
#                                                 # Create a single row filled with NaN for all existing columns
#                                                 df = pd.DataFrame(index=[0], columns=df.columns)
#                                             df.columns.values[9:11] = ["sample", "tumor"]
#                                             df["sampleID"] = pid.split("_")[0]
#                                             df["vcf_path"] = f"{base_path}/{sample_path}"
#                                             all_dfs.append(df)
#                                             break
#                                         except:
#                                             if seq_type == "exon_sequencing":
#                                                 print("no vcf found for this sample", pid)
#                                             else:
#                                                 print("checking exon for", pid.split("_")[0])
#                                             continue

# all_dfs = pd.concat(all_dfs)
# # 
# all_dfs['CLNILinMASTER'] = all_dfs['CLNILinMASTER'].astype(str)

# all_dfs.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/l_lp_vars.parquet", index=None)
import glob
import os
import re
import pandas as pd
import polars as pl

sa = pd.read_csv(
    "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv",
    sep="\t",
)
pids = sa["nct_pid"].unique()

all_dfs = []


def natural_sort_key(path):
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r"(\d+)", path)
    ]


for pid in pids:
    pid_parts = pid.split("_")
    base_sample = f"H021-{pid_parts[0]}"
    variant_suffix = pid_parts[1]  # e.g., 'metastasis'
    found_vcf = False

    for seq_type in ["whole_genome_sequencing", "exon_sequencing"]:
        base_path = f"/omics/odcf/analysis/hipo/hipo_021/GermlineAnalysis/data_object_master_germline/sequencing/{seq_type}/results_per_pid"

        # Match any folder starting with 'germline_smallVariants.<variant_suffix>'
        # e.g., 'germline_smallVariants.metastasis02.blood'
        pattern = os.path.join(
            base_path,
            base_sample,
            f"germline_smallVariants.{variant_suffix}*",
            "v*",
            "**",
            f"smallVariants_{base_sample}.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf*",
        )

        matching_files = glob.glob(pattern, recursive=True)

        if matching_files:
            # Sort to prefer the highest version number first
            matching_files.sort(key=natural_sort_key, reverse=True)
            sample_path = matching_files[0]
            found_vcf = True

            try:
                df = (
                    pl.scan_csv(
                        sample_path,
                        separator="\t",
                        comment_prefix="##",
                        schema_overrides={
                            "QUAL": pl.String,
                            "#CHROM": pl.String,
                            "RNA_VARIANT_AF": pl.Float64,
                        },
                        null_values=["NA", "."],
                    )
                    .filter(
                        (
                            pl.col("ACMG_Classification").is_in(
                                ["Pathogenic", "Likely Pathogenic"]
                            )
                        )
                        | (
                            pl.col("CharGer_Classification").is_in(
                                ["Pathogenic", "Likely Pathogenic", "Likely Benign", "Benign"]
                            )
                        )
                        | (
                            pl.col("ClinVar_Pathogenicity").is_in(
                                ["Pathogenic", "Likely Pathogenic"]
                            )
                        )
                    )
                    .collect(engine="streaming")
                    .to_pandas()
                )

                if df.empty:
                    df = pd.DataFrame(index=[0], columns=df.columns)

                df.columns.values[9:11] = ["sample", "tumor"]
                df["sampleID"] = pid_parts[0]
                df["vcf_path"] = sample_path
                all_dfs.append(df)
                print(
                    f"Found VCF for {pid_parts[0]} in {seq_type} -> {sample_path}"
                )
                break  # Exit seq_type loop once file is retrieved

            except Exception as e:
                print(f"Error processing VCF {sample_path}: {e}")
                found_vcf = False

    if not found_vcf:
        print(f"No VCF found for sample: {pid}")

if all_dfs:
    final_df = pd.concat(all_dfs, ignore_index=True)
    if "CLNILinMASTER" in final_df.columns:
        final_df["CLNILinMASTER"] = final_df["CLNILinMASTER"].astype(str)

    final_df.to_parquet(
        "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/l_lp_vars.parquet",
        index=False,
    )