import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode
import polars as pl
import os
import polars.selectors as cs


init_notebook_mode(all_interactive=True)

import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *


sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")

## load variants

dresden_somatic = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/somatic_damaging_aggregated_all.parquet")
                    .select(["Gene", "PID", "somatic_IMPACT", "somatic_Consequence", "somatic_abexp_v1.1", "somatic_promoterAI", "somatic_max_spliceai_score"])
                    .collect(engine="streaming")
                    
                  ).to_pandas()
# dresden_somatic

print("after reading somatic")


dresden_germline = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/germline_damaging_aggregated_all.parquet")
                    .select(["Gene", "PID", "germline_IMPACT", "germline_Consequence", "germline_abexp_v1.1", "germline_promoterAI", "germline_max_spliceai_score"])
                    .collect(engine="streaming")        
                  ).to_pandas()
# dresden_germline


# py_or_res_all = pl.scan_parquet(or_res_path).select(needed_cols).filter((pl.col("geneID_short").is_in(genes_of_interest[(genes_of_interest["cgc"].isna()) |(genes_of_interest["ROLE_IN_CANCER"].str.contains("TSG"))]["geneID_short"]))).collect().to_pandas()

or_res_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_interesting_genes_padjust_cnv.parquet"

print("before reading py or res")

py_or_res_all = pl.scan_parquet(or_res_path).select(["sampleID", "pValue", "geneID", "geneID_short", "zScore", "CNV"]).collect().to_pandas()



py_or_res_all["Method"] = "OUTRIDER"



py_or_res_all = pd.merge(py_or_res_all, gene_annot_dt[["gene_id", "gene_type"]], left_on="geneID", right_on="gene_id", how="left")

py_or_res_all = py_or_res_all[py_or_res_all["gene_type"] == "protein_coding"]

print("before soring pValue")
py_or_res_all = py_or_res_all.sort_values("pValue")
py_or_res_all = py_or_res_all.drop_duplicates(subset=["sampleID", "geneID_short"])

print("before merging variants")

# py_or_res_aberrant["key"] = py_or_res_aberrant["sampleID"] +  "." + py_or_res_aberrant["nct_pid"].str.split("_").str[1]
py_or_res_all = py_or_res_all.merge(dresden_germline[["Gene", "PID", "germline_IMPACT", "germline_Consequence", "germline_abexp_v1.1", "germline_promoterAI", "germline_max_spliceai_score"]]  , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")
py_or_res_all = py_or_res_all.merge(dresden_somatic[["Gene", "PID", "somatic_IMPACT", "somatic_Consequence", "somatic_abexp_v1.1", "somatic_promoterAI", "somatic_max_spliceai_score"]] , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")

py_or_res_all = py_or_res_all.drop_duplicates(subset=["sampleID", "geneID_short"])



py_or_res_all["germline_somatic_IMPACT"] = py_or_res_all["germline_IMPACT"]
py_or_res_all.loc[py_or_res_all["somatic_IMPACT"].str.contains("HIGH"), "germline_somatic_IMPACT"] = "HIGH"

def get_auc_metrics(df, expression_direction="underexpression", log_scale=True):
    """
    Calculates Area Under the Curve (AUC) and Mean Proportion for variant enrichment.
    If log_scale=True, the AUC represents the area as seen on a log-x plot.
    The 'mean_proportion' represents the average height of the curve over that range.
    """
    # Filter for valid ranks and remove NaNs to prevent calculation errors
    mask = (df[f"{expression_direction}_rank"] > 100) & \
           (df[f"{expression_direction}_rank"].notna())
    
    subset = df[mask].sort_values(f"{expression_direction}_rank")
    
    if subset.empty:
        return {}

    # Define x-coordinates
    x_linear = subset[f"{expression_direction}_rank"].values
    x = np.log10(x_linear) if log_scale else x_linear
    
    # Calculate the span of the x-axis (width)
    x_range = x.max() - x.min()
    
    results = {}
    for variant_type in ["CNV", "VEP"]:
        y_col = f"{expression_direction}_proportions_{variant_type}"
        
        if y_col in subset.columns:
            y = subset[y_col].values
            
            # 1. Calculate Raw AUC (Trapezoidal Rule)
            area = np.trapz(y, x)
            
            # 2. Calculate Mean Proportion (Average Height)
            # This normalizes the area by the width of the plot
            mean_prop = area / x_range if x_range > 0 else 0
            
            results[variant_type] = {
                "raw_auc": area,
                "mean_proportion": mean_prop
            }
            
    return results

def calculate_proportions(df, expression_col="CNV", variant_type="germline", omic_type = "gene"):
    # Work on a copy to prevent SettingWithCopyWarning and index fragmentation
    df = df.copy()
    
    for expression_direction in ["overexpression", "underexpression"]:
        vus_col = f"VUS_{expression_direction}_{expression_col}"
        df[vus_col] = False

        if expression_direction == "underexpression":
            if (df["Method"] == "sf-normalized zScores").any():
                df = df.sort_values("zScore", ascending=True) 
        elif expression_direction == "overexpression":
            if (df["Method"] == "sf-normalized zScores").any():
                df = df.sort_values("zScore", ascending=False)

        # 1. Define VUS based on CNV/VEP type
        if expression_col == "CNV":
            if expression_direction == "underexpression":
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=True) 
                if variant_type == "germline":
                    df.loc[df[expression_col].str.contains("Germline_DEL", na=False), vus_col] = True
                else:
                    df.loc[df[expression_col].str.contains("Somatic_HDEL", na=False), vus_col] = True
            else:
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=False)
                df.loc[df[expression_col].str.contains("AMP|DUP", na=False), vus_col] = True
                
        elif expression_col == "VEP":
            df.loc[(df[f"{variant_type}_IMPACT"].str.contains("HIGH")), vus_col] = True
            if omic_type == "protein":
                df.loc[(df[f"{variant_type}_Consequence"].str.contains("missense")), vus_col] = True
            
        elif expression_col == "CNV_VEP":
            # Combined logic: either high impact VEP or relevant CNV
            if expression_direction == "underexpression":
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=True)
                df.loc[df["CNV"].str.contains("DEL", na=False) | (df["IMPACT"] == "HIGH"), vus_col] = True
            else:
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=False)
                df.loc[df["CNV"].str.contains("AMP|DUP", na=False) | (df["IMPACT"] == "HIGH"), vus_col] = True
        
        else: # Fallback for other custom expression_cols
            if expression_direction == "underexpression":
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=True)   
                df.loc[df[expression_col].str.contains("DEL", na=False), vus_col] = True
            else:
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=False)
                df.loc[df[expression_col].str.contains("AMP|DUP", na=False), vus_col] = True
                df.loc[(df["IMPACT"] == "HIGH"), vus_col] = True

        # 2. Subset based on zScore direction for proportion calculation
        if expression_direction == "underexpression":
            mask = df['zScore'] < 0
        else:
            mask = df['zScore'] > 0
            
        # Create a temporary subset to calculate running statistics
        # This ensures we don't have indexing conflicts with the main DF
        subset = df.loc[mask].copy()
        
        if len(subset) > 0:
            # Calculate ranks and cumulative counts
            u_rank = np.arange(1, len(subset) + 1)
            u_vus = (subset[vus_col] == True).to_numpy()
            u_x = np.cumsum(u_vus)
            
            # Vectorized Clopper–Pearson interval calculation
            alpha = 0.05
            # beta.ppf handles the statistics for the CI
            ci_min = beta.ppf(alpha/2, u_x, u_rank - u_x + 1)
            ci_max = beta.ppf(1 - alpha/2, u_x + 1, u_rank - u_x)

            # Map the calculated values back to the original 'df' using the mask
            df.loc[mask, f"{expression_direction}_rank"] = u_rank
            df.loc[mask, f"{expression_direction}_proportions_{expression_col}_{variant_type}"] = u_x / u_rank
            df.loc[mask, f"{expression_direction}_ci_min_{expression_col}_{variant_type}"] = ci_min
            df.loc[mask, f"{expression_direction}_ci_max_{expression_col}_{variant_type}"] = ci_max

    return df


def plot_variant_proportions(plot_data, plot_title = "OUTRIDER", expression_direction = "underexpression", variant_type="germline"):
    if expression_direction == "underexpression":
        plot_data = plot_data[plot_data["zScore"] < 0]
        CNV_label = "Deletion CNV"
        
    else:
        plot_data = plot_data[plot_data["zScore"] > 0]
        CNV_label = "AMP/DUP CNV"

    return (
        pn.ggplot(plot_data[
                                (plot_data[f"{expression_direction}_rank"] > 100) & 
                                (plot_data[f"{expression_direction}_rank"] < 1e10)])
            
            # --- Line 1 (Deletion) ---
            # We map 'color' to a string label. Plotnine sees this as a category.
            + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                y=f"{expression_direction}_proportions_CNV", 
                                color=f"'{CNV_label}'")) 
            + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                    ymin=f"{expression_direction}_ci_min_CNV", ymax=f"{expression_direction}_ci_max_CNV", 
                                    fill=f"'{CNV_label}'"), 
                            alpha=0.2, outline_type='none')
            
            # --- Line 2 (VEP) ---
            + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                y=f"{expression_direction}_proportions_VEP_{variant_type}", 
                                color="'VEP HIGH impact'"))
            + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                    ymin=f"{expression_direction}_ci_min_VEP_{variant_type}", ymax=f"{expression_direction}_ci_max_VEP_{variant_type}", 
                                    fill="'VEP HIGH impact'"), 
                            alpha=0.2, outline_type='none')
            
            # --- Assigning Colors Manually ---
            # We use a dictionary to ensure 'Deletion CNV' is blue and 'VEP' is red
            + pn.scale_color_manual(name="Variant Type", 
                                    values={CNV_label: "blue", "VEP HIGH impact": "red"})
            + pn.scale_fill_manual(name="Variant Type", 
                                values={CNV_label: "blue", "VEP HIGH impact": "red"})
            
            + pn.scale_x_log10()
            + pn.annotation_logticks(sides="b")
            + pn.labs(
                y=f"Proportion of {expression_direction} outliers\nwith rare variants",
                x="Outlier rank cutoff",
                title = plot_title
            )
            + pn.theme_bw(base_size=12)
    )


def plot_variant_proportions(plot_data, plot_title = "OUTRIDER", expression_direction = "underexpression"):
    if expression_direction == "underexpression":
        plot_data = plot_data[plot_data["zScore"] < 0]
        CNV_label = "Deletion CNV"
        
    else:
        plot_data = plot_data[plot_data["zScore"] > 0]
        CNV_label = "AMP/DUP CNV"

    return (
        pn.ggplot(plot_data[
                                (plot_data[f"{expression_direction}_rank"] > 100) & 
                                (plot_data[f"{expression_direction}_rank"] < 1e10)])
            
            # --- Line 1 (Deletion) ---
            # We map 'color' to a string label. Plotnine sees this as a category.
            + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                y=f"{expression_direction}_proportions_CNV", 
                                color=f"'{CNV_label}'")) 
            + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                    ymin=f"{expression_direction}_ci_min_CNV", ymax=f"{expression_direction}_ci_max_CNV", 
                                    fill=f"'{CNV_label}'"), 
                            alpha=0.2, outline_type='none')
            
            # --- Line 2 (VEP) ---
            + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                y=f"{expression_direction}_proportions_VEP", 
                                color="'VEP HIGH impact'"))
            + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                    ymin=f"{expression_direction}_ci_min_VEP", ymax=f"{expression_direction}_ci_max_VEP", 
                                    fill="'VEP HIGH impact'"), 
                            alpha=0.2, outline_type='none')
            
            # --- Assigning Colors Manually ---
            # We use a dictionary to ensure 'Deletion CNV' is blue and 'VEP' is red
            + pn.scale_color_manual(name="Variant Type", 
                                    values={CNV_label: "blue", "VEP HIGH impact": "red"})
            + pn.scale_fill_manual(name="Variant Type", 
                                values={CNV_label: "blue", "VEP HIGH impact": "red"})
            
            + pn.scale_x_log10()
            + pn.annotation_logticks(sides="b")
            + pn.labs(
                y=f"Proportion of {expression_direction} outliers\nwith rare variants",
                x="Outlier rank cutoff",
                title = plot_title
            )
            + pn.theme_bw(base_size=12)
    )



# py_or_res_all = calculate_proportions(py_or_res_all, expression_col="CNV_VEP")
# py_or_res_all = calculate_proportions(py_or_res_all, expression_col="VEP")

# print(py_or_res_all.columns)

# # underexpression outrider
# p = plot_variant_proportions(
#     py_or_res_all, plot_title="OUTRIDER (across all genes)", expression_direction="underexpression")
# p.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_underexpression_outliers.png", width=6, height=4, dpi=300)

print("before loading gene_zscores")

gene_zscore = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/pyoutrider_zscores/or_variants_cnv.parquet").select(["sampleID", "pvalue", "geneID_short", "zScore", "CNV"]).collect().to_pandas()

gene_zscore = pd.merge(gene_zscore, gene_annot_dt[["geneID_short", "gene_type"]], on="geneID_short", how="left")



# gene_zscore = pd.merge(gene_zscore, sa[["pid", "Diag", "seq_type", "Oncotree Code"]], left_on="sampleID", right_on="pid")
gene_zscore["Method"] = "sf-normalized zScores"
gene_zscore = gene_zscore.drop_duplicates(subset=["geneID_short", "sampleID"])

gene_zscore = gene_zscore[gene_zscore["gene_type"] == "protein_coding"]

# gene_zscore = gene_zscore.drop(columns=["IMPACT"])



# py_or_res_aberrant["key"] = py_or_res_aberrant["sampleID"] +  "." + py_or_res_aberrant["nct_pid"].str.split("_").str[1]
gene_zscore = gene_zscore.merge(dresden_germline[["Gene", "PID", "germline_IMPACT", "germline_Consequence", "germline_abexp_v1.1", "germline_promoterAI", "germline_max_spliceai_score"]] , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")
gene_zscore = gene_zscore.merge(dresden_somatic[["Gene", "PID", "somatic_IMPACT", "somatic_Consequence", "somatic_abexp_v1.1", "somatic_promoterAI", "somatic_max_spliceai_score"]] , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")


gene_zscore = gene_zscore.drop_duplicates(subset=["geneID_short", "sampleID"])


gene_zscore["germline_somatic_IMPACT"] = gene_zscore["germline_IMPACT"]
gene_zscore.loc[gene_zscore["somatic_IMPACT"].str.contains("HIGH"), "germline_somatic_IMPACT"] = "HIGH"




# gene_zscore = calculate_proportions(gene_zscore, "CNV")
# gene_zscore = calculate_proportions(gene_zscore, "VEP")
# gene_zscore = calculate_proportions(gene_zscore, "CNV_VEP")


# # underexpression zscore
# p2 = plot_variant_proportions(
#     gene_zscore, plot_title="sf-normalized zScores (across all genes)", expression_direction="underexpression")
# p2.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_underexpression_outliers_sf_normalized_zscores.png", width=6, height=4, dpi=300)



# # overexpression outrider
# p3 = plot_variant_proportions(
#     py_or_res_all, plot_title="OUTRIDER (across all genes)", expression_direction="overexpression")
# p3.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_overexpression_outliers.png", width=6, height=4, dpi=300)


# # overexpression zscore
# p3 = plot_variant_proportions(
#     gene_zscore, plot_title="sf-normalized zScores (across all genes)", expression_direction="overexpression")
# p3.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_overexpression_outliers_sf_normalized_zscores.png", width=6, height=4, dpi=300)


# # predispostions
# py_or_res_predisp = py_or_res_all[py_or_res_all["padjust_predisp_extended"].notna()]
# py_or_res_predisp = py_or_res_predisp.sort_values("padjust_predisp_extended")

# gene_zscore_predisp = gene_zscore[gene_zscore["geneID_short"].isin(extended_dresden_dt["geneID_short"])]


# gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "CNV_VEP")
# gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "CNV")
# gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "VEP")



# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV")
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "VEP")
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV_VEP")


# p4 = plot_variant_proportions(
#     py_or_res_predisp, plot_title="OUTRIDER (only predispostion genes)", expression_direction="underexpression")
# p4.save("/home/a379i/Scripts/AberrantExpression/proportions_predisposition_genes_underexpression_outliers.png", width=6, height=4, dpi=300)


# p5 = plot_variant_proportions(
#     py_or_res_predisp, plot_title="OUTRIDER (only predispostion genes)", expression_direction="overexpression")
# p5.save("/home/a379i/Scripts/AberrantExpression/proportions_predisposition_genes_overexpression_outliers.png", width=6, height=4, dpi=300)



# zScores vs OUTRIDER: both CNV and VEP_HIGH

# py_or_res_all = calculate_proportions(py_or_res_all, "CNV_VEP")
# gene_zscore = calculate_proportions(gene_zscore, "CNV_VEP")



# plot_dt = pd.concat([py_or_res_all, gene_zscore])



# for expression_direction in ["underexpression", "overexpression"]:

#     if expression_direction == "underexpression":
#             plot_data = plot_dt[plot_dt["zScore"] < 0]
#             CNV_label = "Deletion CNV"
#             CNV_label = "AMP/DUP CNV"
#     else:
#         plot_data = plot_dt[plot_dt["zScore"] > 0]


#     p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
#                                     (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
#                 # --- Line 1 (Deletion) ---
#                 # We map 'color' to a string label. Plotnine sees this as a category.
#                 + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
#                                     y=f"{expression_direction}_proportions_CNV_VEP", 
#                                     color="Method")) 
#                 + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
#                                         ymin=f"{expression_direction}_ci_min_CNV_VEP", ymax=f"{expression_direction}_ci_max_CNV_VEP", 
#                                         fill="Method"), 
#                                 alpha=0.2, outline_type='none')
                
            
#                 + pn.scale_x_log10()
#                 + pn.annotation_logticks(sides="b")
#                 + pn.labs(
#                     y=f"Proportion of {expression_direction} outliers\nwith variants (Rare SNV/indel and CNV)",
#                     x="Outlier rank cutoff",
#                     title="Across all genes"
#                 )
#                 + pn.theme_bw(base_size=12)
#     )

#     p.save(f"/home/a379i/Scripts/AberrantExpression/proportions_{expression_direction}_or_vs_zScore.png", width=6, height=4, dpi=300)

#     p = ( pn.ggplot(plot_data[(plot_data["Method"] == "OUTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
#                                     (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
#                 # --- Line 1 (Deletion) ---
#                 # We map 'color' to a string label. Plotnine sees this as a category.
#                 + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
#                                     y=f"{expression_direction}_proportions_CNV_VEP", 
#                                     )) 
#                 + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
#                                         ymin=f"{expression_direction}_ci_min_CNV_VEP", ymax=f"{expression_direction}_ci_max_CNV_VEP"), 
#                                 alpha=0.2, outline_type='none')
                
            
#                 + pn.scale_x_log10()
#                 + pn.annotation_logticks(sides="b")
#                 + pn.labs(
#                     y=f"Proportion of {expression_direction} outliers\nwith variants (Rare SNV/indel and CNV)",
#                     x="Outlier rank cutoff",
#                     title="Across all genes"
#                 )
#                 + pn.theme_bw(base_size=12)
#     )

#     p.save(f"/home/a379i/Scripts/AberrantExpression/proportions_{expression_direction}_vep_cnv.png", width=6, height=4, dpi=300)

# zScores vs OUTRIDER: both CNV and VEP_HIGH


# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV_VEP")
# gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "CNV_VEP")


plot_df = calculate_proportions(py_or_res_all, "VEP", "germline_somatic")
plot_df["Variant_type"] = "Small germline or somatic"

plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_germline_somatic": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_germline_somatic": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_germline_somatic": "underexpression_ci_min"})


plot_df2 = calculate_proportions(gene_zscore, "VEP", "germline_somatic")
plot_df2["Variant_type"] = "Small germline or somatic"

plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_germline_somatic": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_germline_somatic": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_germline_somatic": "underexpression_ci_min"})


plot_df2["Method"] = "sf-normalized\nzScores"



plot_dt = pd.concat([plot_df, plot_df2])

plot_dt["facet_label"] = f"VEP high impact variants (n={py_or_res_all[py_or_res_all["germline_somatic_IMPACT"] == "HIGH"].shape[0]})"


#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["underexpression"]:
    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                + pn.scale_fill_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
               + pn.scale_color_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of gene {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across expressed protein coding genes"
                )
                + pn.theme_bw(base_size=8)
                + pn.theme(
                     strip_text = pn.element_text(size =7),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/gene_expression_proportions/all_genes/proportions_all_genes_{expression_direction}_vep_small_somatic_germline_or_vs_zScore.png", width=4, height=3.6, dpi=300)



py_or_res_all["expression_damaging_IMPACT"] = "False"

py_or_res_all.loc[(py_or_res_all["germline_IMPACT"] == "HIGH") | (py_or_res_all["somatic_IMPACT"] == "HIGH") | (py_or_res_all["somatic_promoterAI"] <= -0.1) | (py_or_res_all["germline_promoterAI"] <= -0.1)
    | (py_or_res_all["germline_max_spliceai_score"] >= 0.1) | (py_or_res_all["somatic_max_spliceai_score"] >= 0.1)  | (py_or_res_all["germline_abexp_v1.1"] <= -0.1) | (py_or_res_all["somatic_abexp_v1.1"] <= -0.1), "expression_damaging_IMPACT"] = "HIGH"


plot_df = calculate_proportions(py_or_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline or somatic"


gene_zscore["expression_damaging_IMPACT"] = "False"

gene_zscore.loc[(gene_zscore["germline_IMPACT"] == "HIGH") | (gene_zscore["somatic_IMPACT"] == "HIGH") | (gene_zscore["somatic_promoterAI"] <= -0.1) | (gene_zscore["germline_promoterAI"] <= -0.1)
    | (gene_zscore["germline_max_spliceai_score"] >= 0.1) | (gene_zscore["somatic_max_spliceai_score"] >= 0.1)  | (gene_zscore["germline_abexp_v1.1"] <= -0.1) | (gene_zscore["somatic_abexp_v1.1"] <= -0.1), "expression_damaging_IMPACT"] = "HIGH"


plot_df2 = calculate_proportions(gene_zscore, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline or somatic"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})


plot_df2["Method"] = "sf-normalized\nzScores"

plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Expression disrupting variants\nVEP high impact, AbExp ≤ -0.1,\npromoterAI ≤ -0.1, SpliceAI ≥ 0.1 (n={py_or_res_all[py_or_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["underexpression"]:
    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of gene {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across expressed protein coding genes"
                )
                + pn.scale_fill_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.theme_bw(base_size=8)
                + pn.theme(
                     strip_text = pn.element_text(size =7),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/gene_expression_proportions/all_genes/proportions_germline_somatic_all_genes_{expression_direction}_expression_disrupting_or_vs_zScore.png", width=4, height=3.6, dpi=300)



py_or_res_all["expression_damaging_IMPACT"] = "False"

py_or_res_all.loc[(py_or_res_all["germline_IMPACT"] == "HIGH") | (py_or_res_all["somatic_promoterAI"] <= -0.1) | (py_or_res_all["germline_promoterAI"] <= -0.1)
    | (py_or_res_all["germline_max_spliceai_score"] >= 0.1) | (py_or_res_all["somatic_max_spliceai_score"] >= 0.1)  | (py_or_res_all["germline_abexp_v1.1"] <= -0.1) | (py_or_res_all["somatic_abexp_v1.1"] <= -0.1), "expression_damaging_IMPACT"] = "HIGH"


plot_df = calculate_proportions(py_or_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"


gene_zscore["expression_damaging_IMPACT"] = "False"

gene_zscore.loc[(gene_zscore["germline_IMPACT"] == "HIGH") | (gene_zscore["somatic_promoterAI"] <= -0.1) | (gene_zscore["germline_promoterAI"] <= -0.1)
    | (gene_zscore["germline_max_spliceai_score"] >= 0.1) | (gene_zscore["somatic_max_spliceai_score"] >= 0.1)  | (gene_zscore["germline_abexp_v1.1"] <= -0.1) | (gene_zscore["somatic_abexp_v1.1"] <= -0.1), "expression_damaging_IMPACT"] = "HIGH"


plot_df2 = calculate_proportions(gene_zscore, "VEP", "expression_damaging")

plot_df2["Variant_type"] = "Small germline"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})


plot_df2["Method"] = "sf-normalized\nzScores"

plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Expression disrupting variants\nVEP high impact, AbExp ≤ -0.1,\npromoterAI ≤ -0.1, SpliceAI ≥ 0.1 (n={py_or_res_all[py_or_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["underexpression"]:
    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of gene {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across expressed protein coding genes"
                )
                + pn.scale_fill_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.theme_bw(base_size=8)
                + pn.theme(
                     strip_text = pn.element_text(size =6),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/gene_expression_proportions/all_genes/proportions_germline_all_genes_{expression_direction}_expression_disrupting_or_vs_zScore.png", width=4, height=3.6, dpi=300)




py_or_res_all["expression_damaging_IMPACT"] = "False"

py_or_res_all.loc[(py_or_res_all["germline_promoterAI"] >= 0.1) | (py_or_res_all["somatic_promoterAI"] >= 0.1) | (py_or_res_all["CNV"].isin(["Somatic_AMP", "Germline_DUP"])) , "expression_damaging_IMPACT"] = "HIGH"




plot_df = calculate_proportions(py_or_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"

gene_zscore["expression_damaging_IMPACT"] = "False"

gene_zscore.loc[(gene_zscore["germline_promoterAI"] >= 0.1) | (gene_zscore["somatic_promoterAI"] >= 0.1) | (gene_zscore["CNV"].isin(["Somatic_AMP", "Germline_DUP"])) , "expression_damaging_IMPACT"] = "HIGH"

plot_df2 = calculate_proportions(gene_zscore, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline"


plot_df2 = plot_df2.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Gene expression disrupting variants\npromoterAI ≥ 0.1, Somatic AMP\nGermline DUP (n={py_or_res_all[py_or_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"




#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["overexpression"]:
    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels or CNV",
                    x="Outlier rank cutoff",
                    title="Across expressed protein coding genes"
                )
                + pn.scale_fill_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.theme_bw(base_size=8)
                + pn.theme(
                     strip_text = pn.element_text(size =7),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/gene_expression_proportions/all_genes/proportions_all_genes_{expression_direction}_or_vs_zScore_cnv.png",  width=4, height=3.6, dpi=300)




py_or_res_all["expression_damaging_IMPACT"] = "False"

py_or_res_all.loc[(py_or_res_all["germline_IMPACT"] == "HIGH") | (py_or_res_all["somatic_IMPACT"] == "HIGH") | (py_or_res_all["somatic_promoterAI"] <= -0.1) | (py_or_res_all["germline_promoterAI"] <= -0.1)
    | (py_or_res_all["germline_max_spliceai_score"] >= 0.1) | (py_or_res_all["somatic_max_spliceai_score"] >= 0.1)  | (py_or_res_all["germline_abexp_v1.1"] <= -0.1) | (py_or_res_all["somatic_abexp_v1.1"] <= -0.1) | (py_or_res_all["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])), "expression_damaging_IMPACT"] = "HIGH"


plot_df = calculate_proportions(py_or_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline or somatic"


gene_zscore["expression_damaging_IMPACT"] = "False"

gene_zscore.loc[(gene_zscore["germline_IMPACT"] == "HIGH") | (gene_zscore["somatic_IMPACT"] == "HIGH") | (gene_zscore["somatic_promoterAI"] <= -0.1) | (gene_zscore["germline_promoterAI"] <= -0.1)
    | (gene_zscore["germline_max_spliceai_score"] >= 0.1) | (gene_zscore["somatic_max_spliceai_score"] >= 0.1)  | (gene_zscore["germline_abexp_v1.1"] <= -0.1) | (gene_zscore["somatic_abexp_v1.1"] <= -0.1) | (gene_zscore["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])), "expression_damaging_IMPACT"] = "HIGH"


plot_df2 = calculate_proportions(gene_zscore, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline or somatic"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})


plot_df2["Method"] = "sf-normalized\nzScores"

plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Expression disrupting variants\nVEP high impact, AbExp ≤ -0.1,\npromoterAI ≤ -0.1, SpliceAI ≥ 0.1 \n Deletion CNV (n={py_or_res_all[py_or_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["underexpression"]:
    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of gene {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across expressed protein coding genes"
                )
                + pn.scale_fill_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.theme_bw(base_size=8)
                + pn.theme(
                     strip_text = pn.element_text(size =7),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/gene_expression_proportions/all_genes/proportions_germline_somatic_cnv_all_genes_{expression_direction}_expression_disrupting_or_vs_zScore.png", width=4, height=3.6, dpi=300)







py_or_res_all["expression_damaging_IMPACT"] = "False"

py_or_res_all.loc[(py_or_res_all["germline_promoterAI"] >= 0.1) | (py_or_res_all["somatic_promoterAI"] >= 0.1)  , "expression_damaging_IMPACT"] = "HIGH"




plot_df = calculate_proportions(py_or_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"

gene_zscore["expression_damaging_IMPACT"] = "False"

gene_zscore.loc[(gene_zscore["germline_promoterAI"] >= 0.1) | (gene_zscore["somatic_promoterAI"] >= 0.1) , "expression_damaging_IMPACT"] = "HIGH"

plot_df2 = calculate_proportions(gene_zscore, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline"


plot_df2 = plot_df2.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Gene expression disrupting variants\npromoterAI ≥ 0.1 (n={py_or_res_all[py_or_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"




#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["overexpression"]:
    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels or CNV",
                    x="Outlier rank cutoff",
                    title="Across expressed protein coding genes"
                )
                + pn.scale_fill_manual(values={
                        "OUTRIDER": "#f8766d", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                       "OUTRIDER": "#f8766d",  
                        "sf-normalized\nzScores": "grey"
                })
                + pn.theme_bw(base_size=8)
                + pn.theme(
                     strip_text = pn.element_text(size =7),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/gene_expression_proportions/all_genes/proportions_all_genes_{expression_direction}_or_vs_zScore_promoterAI_gt_01.png",  width=4, height=3.6, dpi=300)

