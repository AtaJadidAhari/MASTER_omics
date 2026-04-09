import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode
import polars as pl
import os

init_notebook_mode(all_interactive=True)

import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")


or_res_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust_cnv.parquet"
needed_cols = ["sampleID", "zScore", "pValue", "padjust", "IMPACT", "geneID",  "geneID_short",
               "#Uploaded_variation_snv", "IMPACT_snv", "ANNOTATION_control_snv", "Consequence_snv", "promoterAI_snv",
               "Location_indel", "IMPACT_indel", "ANNOTATION_control_indel", "Consequence_indel",
               "padjust_predisp", "padjust_predisp_extended", "CNV" ]
py_or_res_all = pd.read_parquet(or_res_path, columns=needed_cols)
py_or_res_all["Method"] = "OUTRIDER"

# py_or_res_all = py_or_res_all[py_or_res_all["padjust_predisp_extended"].notna()]
py_or_res_all = pd.merge(py_or_res_all, sa[["pid", "Diag", "seq_type", "Oncotree Code"]], left_on="sampleID", right_on="pid")
py_or_res_all = pd.merge(py_or_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID", "ROLE_IN_CANCER"]], on="geneID", how="left")

py_or_res_all = py_or_res_all.sort_values("pValue")
py_or_res_all = py_or_res_all.drop_duplicates(subset=["sampleID", "geneID_short"])




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

def calculate_proportions(df, expression_col="CNV"):
    # Work on a copy to prevent SettingWithCopyWarning and index fragmentation
    df = df.copy()
    
    for expression_direction in ["overexpression", "underexpression"]:
        vus_col = f"VUS_{expression_direction}_{expression_col}"
        df[vus_col] = False

        # 1. Define VUS based on CNV/VEP type
        if expression_col == "CNV":
            if expression_direction == "underexpression":
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=True) 
                df.loc[df[expression_col].str.contains("DEL", na=False), vus_col] = True
            else:
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=False)
                df.loc[df[expression_col].str.contains("AMP|DUP", na=False), vus_col] = True
                
        elif expression_col == "VEP":
            df.loc[(df["IMPACT"] == "HIGH"), vus_col] = True
            
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
            df.loc[mask, f"{expression_direction}_proportions_{expression_col}"] = u_x / u_rank
            df.loc[mask, f"{expression_direction}_ci_min_{expression_col}"] = ci_min
            df.loc[mask, f"{expression_direction}_ci_max_{expression_col}"] = ci_max

    return df


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


py_or_res_all = calculate_proportions(py_or_res_all, expression_col="CNV")
py_or_res_all = calculate_proportions(py_or_res_all, expression_col="VEP")

print(py_or_res_all.columns)

# underexpression outrider
p = plot_variant_proportions(
    py_or_res_all, plot_title="OUTRIDER (across all genes)", expression_direction="underexpression")
p.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_underexpression_outliers.png", width=6, height=4, dpi=300)



gene_zscore = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/zscores/or_variants_cnv.parquet",)
gene_zscore = pd.merge(gene_zscore, sa[["pid", "Diag", "seq_type", "Oncotree Code"]], left_on="sampleID", right_on="pid")
gene_zscore["Method"] = "sf-normalized zScores"
gene_zscore = gene_zscore.drop_duplicates(subset=["geneID_short", "sampleID"])


gene_zscore = calculate_proportions(gene_zscore, "CNV")
gene_zscore = calculate_proportions(gene_zscore, "VEP")


# underexpression zscore
p2 = plot_variant_proportions(
    gene_zscore, plot_title="sf-normalized zScores (across all genes)", expression_direction="underexpression")
p2.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_underexpression_outliers_sf_normalized_zscores.png", width=6, height=4, dpi=300)



# overexpression outrider
p3 = plot_variant_proportions(
    py_or_res_all, plot_title="OUTRIDER (across all genes)", expression_direction="overexpression")
p3.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_overexpression_outliers.png", width=6, height=4, dpi=300)


# overexpression zscore
p3 = plot_variant_proportions(
    gene_zscore, plot_title="sf-normalized zScores (across all genes)", expression_direction="overexpression")
p3.save("/home/a379i/Scripts/AberrantExpression/proportions_all_genes_overexpression_outliers_sf_normalized_zscores.png", width=6, height=4, dpi=300)


# predispostions
py_or_res_predisp = py_or_res_all[py_or_res_all["padjust_predisp_extended"].notna()]
py_or_res_predisp = py_or_res_predisp.sort_values("padjust_predisp_extended")

gene_zscore_predisp = gene_zscore[gene_zscore["geneID_short"].isin(extended_dresden_dt["geneID_short"])]


gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "CNV_VEP")
gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "CNV")
gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "VEP")



py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV")
py_or_res_predisp = calculate_proportions(py_or_res_predisp, "VEP")
py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV_VEP")


p4 = plot_variant_proportions(
    py_or_res_predisp, plot_title="OUTRIDER (only predispostion genes)", expression_direction="underexpression")
p4.save("/home/a379i/Scripts/AberrantExpression/proportions_predisposition_genes_underexpression_outliers.png", width=6, height=4, dpi=300)


p5 = plot_variant_proportions(
    py_or_res_predisp, plot_title="OUTRIDER (only predispostion genes)", expression_direction="overexpression")
p5.save("/home/a379i/Scripts/AberrantExpression/proportions_predisposition_genes_overexpression_outliers.png", width=6, height=4, dpi=300)



# zScores vs OUTRIDER: both CNV and VEP_HIGH

py_or_res_all = calculate_proportions(py_or_res_all, "CNV_VEP")
gene_zscore = calculate_proportions(gene_zscore, "CNV_VEP")



plot_dt = pd.concat([py_or_res_all, gene_zscore])



for expression_direction in ["underexpression", "overexpression"]:

    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP/DUP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions_CNV_VEP", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min_CNV_VEP", ymax=f"{expression_direction}_ci_max_CNV_VEP", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
            
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of {expression_direction} outliers\nwith variants (Rare SNV/indel and CNV)",
                    x="Outlier rank cutoff",
                    title="Across all genes"
                )
                + pn.theme_bw(base_size=12)
    )

    p.save(f"/home/a379i/Scripts/AberrantExpression/proportions_{expression_direction}_or_vs_zScore.png", width=6, height=4, dpi=300)

    p = ( pn.ggplot(plot_data[(plot_data["Method"] == "OUTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions_CNV_VEP", 
                                    )) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min_CNV_VEP", ymax=f"{expression_direction}_ci_max_CNV_VEP"), 
                                alpha=0.2, outline_type='none')
                
            
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of {expression_direction} outliers\nwith variants (Rare SNV/indel and CNV)",
                    x="Outlier rank cutoff",
                    title="Across all genes"
                )
                + pn.theme_bw(base_size=12)
    )

    p.save(f"/home/a379i/Scripts/AberrnatExpression/proportions_{expression_direction}_vep_cnv.png", width=6, height=4, dpi=300)

# zScores vs OUTRIDER: both CNV and VEP_HIGH

plot_dt = pd.concat([gene_zscore_predisp, py_or_res_predisp])



for expression_direction in ["underexpression", "overexpression"]:

    if expression_direction == "underexpression":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP/DUP CNV"
    else:
        plot_data = plot_dt[plot_dt["zScore"] > 0]


    p = ( pn.ggplot(plot_data[(plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions_CNV_VEP", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min_CNV_VEP", ymax=f"{expression_direction}_ci_max_CNV_VEP", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
            
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of {expression_direction} outliers\nwith variants (Rare SNV/indel and CNV)",
                    x="Outlier rank cutoff",
                    title="Across predisposistion genes"
                )
                + pn.theme_bw(base_size=12)
    )

    p.save(f"/home/a379i/Scripts/AberrantExpression/proportions_predisp_{expression_direction}_or_vs_zScore.png", width=6, height=4, dpi=300)

    p = ( pn.ggplot(plot_data[(plot_data["Method"] == "OUTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions_CNV_VEP", 
                                    )) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min_CNV_VEP", ymax=f"{expression_direction}_ci_max_CNV_VEP"), 
                                alpha=0.2, outline_type='none')
                
            
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of {expression_direction} outliers\nwith variants (Rare SNV/indel and CNV)",
                    x="Outlier rank cutoff",
                    title="Across all genes"
                )
                + pn.theme_bw(base_size=12)
    )

    p.save(f"/home/a379i/Scripts/AberrnatExpression/proportions_predisp{expression_direction}_vep_cnv.png", width=6, height=4, dpi=300)