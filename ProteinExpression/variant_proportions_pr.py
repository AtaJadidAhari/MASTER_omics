# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     sync_on_save: true
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
# ---

# %%
import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode
import polars as pl
import os
import polars.selectors as cs


# %%
init_notebook_mode(all_interactive=True)

# %%
import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

# %%
from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

# %%
sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")



## load variants

dresden_somatic = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/somatic_aggregated_predisp_cgc.parquet")
                    .with_columns(
                            cs.starts_with("SpliceAI_pred_DS_").cast(pl.Float64, strict=False)
                    )
                   .filter(
                           (pl.col("somatic_abexp_v1.1").abs() >= 0.1) | (pl.col("somatic_promoterAI").abs() >= 0.1) | (pl.col("somatic_IMPACT").str.contains("HIGH")) |
                            (pl.col("SpliceAI_pred_DS_AL") >= 0.1) |  
                               (pl.col("SpliceAI_pred_DS_AG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DL") >= 0.1)
                          )
                    .collect(engine="streaming")
                    
                  ).to_pandas().drop(columns=["ref", "ref"])
# dresden_somatic



dresden_germline = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/germline_aggregated_predisp_cgc_tinda.parquet")
                    .with_columns(
                            cs.starts_with("SpliceAI_pred_DS_").cast(pl.Float64, strict=False)
                    ).filter((pl.col("PID").is_in(sa["pid"])) & (
                               (pl.col("germline_abexp_v1.1").abs() >= 0.1) | (pl.col("germline_promoterAI").abs() >= 0.1) | (pl.col("germline_IMPACT").str.contains("HIGH")) | 
                              (pl.col("SpliceAI_pred_DS_AL") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_AG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DL") >= 0.1)
                            )
                          )
                    .collect(engine="streaming")        
                  ).to_pandas()


dresden_germline["germline_max_spliceai_score"] = dresden_germline[["SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_AG", 
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"
]].max(axis=1)


dresden_somatic["somatic_max_spliceai_score"] = dresden_somatic[["SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_AG", 
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"
]].max(axis=1)



# %%
pr_output_name = "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"

# %%
# pr_res_all = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants.csv")
pr_res_all = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants_predisppadjust_cnv.parquet").filter(pl.col("geneID_short").is_in(genes_of_interest["geneID_short"])).collect().to_pandas()

# %%
pr_res_all = pd.merge(pr_res_all, sa, left_on="sampleID", right_on="pid")

# %%
pr_res_all = pr_res_all.sort_values("pValue")
pr_res_all = pr_res_all.drop_duplicates(subset=["proteinID", "sampleID"])
print(len(pr_res_all), " protrider length")

# %%
pr_res_all["Method"] = "PROTRIDER"
print(pr_res_all.shape)



# py_or_res_aberrant["key"] = py_or_res_aberrant["sampleID"] +  "." + py_or_res_aberrant["nct_pid"].str.split("_").str[1]
pr_res_all = pr_res_all.merge(dresden_germline , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")
pr_res_all = pr_res_all.merge(dresden_somatic , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")

pr_res_all["germline_somatic_IMPACT"] = pr_res_all["germline_IMPACT"]
pr_res_all.loc[pr_res_all["somatic_IMPACT"].str.contains("HIGH"), "germline_somatic_IMPACT"] = "HIGH"
pr_res_all.loc[pr_res_all["somatic_Consequence"].str.contains("missense"), "germline_somatic_IMPACT"] = "HIGH"
pr_res_all.loc[pr_res_all["germline_Consequence"].str.contains("missense"), "germline_somatic_IMPACT"] = "HIGH"




# %%
pr_output_name_2 = "sf_zScores"

# %%
pr_res_all_zScores = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name_2 + "/pr_variants_predisppadjust_cnv.parquet").filter(pl.col("geneID_short").is_in(genes_of_interest["geneID_short"])).collect().to_pandas()

pr_res_all_zScores = pd.merge(pr_res_all_zScores, sa, left_on="sampleID", right_on="pid")

# %%
pr_res_all_zScores = pr_res_all_zScores.drop_duplicates(subset=["proteinID", "sampleID"])
pr_res_all_zScores["Method"] = "sf-normalized zScores"
print(len(pr_res_all_zScores), " zscores length")


pr_res_all_zScores = pr_res_all_zScores.merge(dresden_germline , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")
pr_res_all_zScores = pr_res_all_zScores.merge(dresden_somatic , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")




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
    

# %%
def plot_variant_proportions(plot_data, auc_dict, plot_title = "OUTRIDER", expression_direction = "underexpression"):

    
    if expression_direction == "underexpression":
        direction_mask = plot_data["zScore"] < 0
        plot_data = plot_data[plot_data["zScore"] < 0]
        CNV_label = "Deletion CNV"
        
    else:
        direction_mask = plot_data["zScore"] > 0
        plot_data = plot_data[plot_data["zScore"] > 0]
        CNV_label = "AMP/DUP CNV"

    # Apply filtering for the plot range (same as in ggplot call)
    plot_subset = plot_data[
        direction_mask & 
        (plot_data[f"{expression_direction}_rank"] > 100) & 
        (plot_data[f"{expression_direction}_rank"] < 1e10)
    ].copy()

    # 2. Calculate Dynamic X and Y
    # Get the furthest point on the right
    max_x = plot_subset[f"{expression_direction}_rank"].max()
    
    # Find the absolute highest point in the plot (top of the CI ribbons)
    ci_cols = [f"{expression_direction}_ci_max_CNV", f"{expression_direction}_ci_max_VEP"]
    # We take the max of those columns, then add a 5% buffer of the total y-span
    actual_y_max = plot_subset[ci_cols].max().max()
    dynamic_y = actual_y_max + (actual_y_max * 0.05)
    
    cnv_text = f"{CNV_label} \n (AUC: {auc_dict['CNV']['mean_proportion']:.1%})"
    vep_text = f"VEP \n (AUC:{auc_dict['VEP']['mean_proportion']:.1%})"
    annotation_text = f"{cnv_text}\n{vep_text}"

    max_rank = plot_data[f"{expression_direction}_rank"].max()

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
                                    values={CNV_label: "blue", "VEP HIGH impact": "red"},
                                   labels={CNV_label: cnv_text, "VEP HIGH impact": vep_text})
            + pn.scale_fill_manual(name="Variant Type", 
                                values={CNV_label: "blue", "VEP HIGH impact": "red"},
                                  labels={CNV_label: cnv_text, "VEP HIGH impact": vep_text})
            
            + pn.scale_x_log10()
            + pn.annotation_logticks(sides="b")
            + pn.labs(
                y=f"Proportion of {expression_direction} outliers\nwith rare variants",
                x="Outlier rank cutoff",
                title = plot_title
            )
            + pn.theme_bw(base_size=12)
    )


# # %%
# pr_res_all = calculate_proportions(pr_res_all, expression_col="CNV")
# pr_res_all = calculate_proportions(pr_res_all, expression_col="VEP")


# # %%
# pr_res_all_zScores = calculate_proportions(pr_res_all_zScores, "CNV")
# pr_res_all_zScores = calculate_proportions(pr_res_all_zScores, "VEP")


# protrider_auc_underexpression = get_auc_metrics(pr_res_all)

# zscore_auc_underexpression = get_auc_metrics(pr_res_all_zScores)

# protrider_auc_overexpression = get_auc_metrics(pr_res_all, expression_direction="overexpression")

# zscore_auc_overexpression = get_auc_metrics(pr_res_all_zScores, expression_direction="overexpression")

# # %%
# print(pr_res_all.shape, "len protrider after adding proportions")

# # %%
# p = plot_variant_proportions(
#     pr_res_all, protrider_auc_underexpression, plot_title="PTORDIER (across all genes)", expression_direction="underexpression")
# p.save("/home/a379i/Scripts/ProteinExpression/proportions_all_genes_underexpression_outliers.png", width=6, height=4, dpi=300)

# # %%






# # %%
# p2 = plot_variant_proportions(
#     pr_res_all_zScores, zscore_auc_underexpression, plot_title="sf-normalized zScores (across all genes)", expression_direction="underexpression")
# p2.save("/home/a379i/Scripts/ProteinExpression/proportions_all_genes_underexpression_outliers_sf_normalized_zscores.png", width=6, height=4, dpi=300)

# # %%

# # %%
# # overexpression
# p3 = plot_variant_proportions(
#     pr_res_all, zscore_auc_overexpression, plot_title="PROTRIDER (across all genes)", expression_direction="overexpression")
# p3.save("/home/a379i/Scripts/ProteinExpression/proportions_all_genes_overexpression_outliers.png", width=6, height=4, dpi=300)


# # %%
# # predispostions
# py_or_res_predisp = pr_res_all[pr_res_all["padjust_predisp_extended"].notna()]
# py_or_res_predisp = py_or_res_predisp.sort_values("padjust_predisp_extended")

# # %%
# gene_zscore_predisp = pr_res_all_zScores[pr_res_all_zScores["proteinID"].isin(extended_dresden_list)]


# # %%
# gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "CNV_VEP")
# gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "CNV")
# gene_zscore_predisp = calculate_proportions(gene_zscore_predisp, "VEP")

# # %%

# # %%
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV")
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "VEP")
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV_VEP")



# protrider_auc_underexpression_predisp = get_auc_metrics(py_or_res_predisp)

# zscore_auc_underexpression_predisp = get_auc_metrics(gene_zscore_predisp)

# protrider_auc_overexpression_predisp = get_auc_metrics(py_or_res_predisp, expression_direction="overexpression")

# zscore_auc_overexpression_predisp = get_auc_metrics(gene_zscore_predisp, expression_direction="overexpression")

# # %%
# p4 = plot_variant_proportions(
#     py_or_res_predisp, protrider_auc_underexpression_predisp, plot_title="OUTRIDER (only predispostion genes)", expression_direction="underexpression")
# p4.save("/home/a379i/Scripts/ProteinExpression/proportions_predisposition_genes_underexpression_outliers.png", width=6, height=4, dpi=300)


# # %%
# p5 = plot_variant_proportions(
#     py_or_res_predisp, protrider_auc_overexpression_predisp, plot_title="OUTRIDER (only predisposition genes)", expression_direction="overexpression")
# p5.save("/home/a379i/Scripts/ProteinExpression/proportions_predisposition_genes_overexpression_outliers.png", width=6, height=4, dpi=300)

# # %%


# %% [markdown]
# zScores vs OUTRIDER: both CNV and VEP_HIGH


pr_res_all_zScores["germline_somatic_IMPACT"] = pr_res_all_zScores["germline_IMPACT"]
pr_res_all_zScores.loc[pr_res_all_zScores["somatic_IMPACT"].str.contains("HIGH"), "germline_somatic_IMPACT"] = "HIGH"

pr_res_all_zScores.loc[pr_res_all_zScores["somatic_Consequence"].str.contains("missense"), "germline_somatic_IMPACT"] = "HIGH"
pr_res_all_zScores.loc[pr_res_all_zScores["germline_Consequence"].str.contains("missense"), "germline_somatic_IMPACT"] = "HIGH"



plot_df = calculate_proportions(pr_res_all, "VEP", "germline_somatic")
plot_df["Variant_type"] = "Small germline or somatic"

plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_germline_somatic": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_germline_somatic": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_germline_somatic": "underexpression_ci_min"})


plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "germline_somatic")
plot_df2["Variant_type"] = "Small germline or somatic"

plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_germline_somatic": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_germline_somatic": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_germline_somatic": "underexpression_ci_min"})


plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])

plot_dt["facet_label"] = f"VEP high impact/missense variants \n (n={pr_res_all[pr_res_all["germline_somatic_IMPACT"] == "HIGH"].shape[0]})"

print(plot_dt["facet_label"])

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
                    y=f"Proportion of protein {expression_direction} outliers\nwith rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_cgc_predisp_{expression_direction}_vep_small_somatic_germline_or_vs_zScore.png", width=4, height=3.6, dpi=300)

    p = ( pn.ggplot(plot_data[(plot_data["Method"] == "PROTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    )) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of protein {expression_direction} outliers\nwith rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
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

    p.save(f"/home/a379i/Scripts/ProteinExpression/proportions_cgc_predisp_{expression_direction}_vep_small_somatic_germline.png", width=4, height=3.6, dpi=300)


pr_res_all["expression_damaging_IMPACT"] = "False"

pr_res_all.loc[(pr_res_all["germline_IMPACT"] == "HIGH") | (pr_res_all["somatic_promoterAI"] <= -0.1) | (pr_res_all["germline_promoterAI"] <= -0.1)
    | (pr_res_all["germline_max_spliceai_score"] >= 0.1) | (pr_res_all["somatic_max_spliceai_score"] >= 0.1)  | (pr_res_all["germline_abexp_v1.1"] <= -0.1) | (pr_res_all["somatic_abexp_v1.1"] <= -0.1), "expression_damaging_IMPACT"] = "HIGH"


pr_res_all.loc[pr_res_all["somatic_Consequence"].str.contains("missense"), "expression_damaging_IMPACT"] = "HIGH"
pr_res_all.loc[pr_res_all["germline_Consequence"].str.contains("missense"), "expression_damaging_IMPACT"] = "HIGH"



plot_df = calculate_proportions(pr_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline or somatic"


pr_res_all_zScores["expression_damaging_IMPACT"] = "False"


pr_res_all_zScores.loc[(pr_res_all_zScores["germline_IMPACT"] == "HIGH") | (pr_res_all_zScores["somatic_promoterAI"] <= -0.1) | (pr_res_all_zScores["germline_promoterAI"] <= -0.1)
    | (pr_res_all_zScores["germline_max_spliceai_score"] >= 0.1) | (pr_res_all_zScores["somatic_max_spliceai_score"] >= 0.1)  | (pr_res_all_zScores["germline_abexp_v1.1"] <= -0.1) | (pr_res_all_zScores["somatic_abexp_v1.1"] <= -0.1), "expression_damaging_IMPACT"] = "HIGH"


pr_res_all_zScores.loc[pr_res_all_zScores["somatic_Consequence"].str.contains("missense"), "expression_damaging_IMPACT"] = "HIGH"
pr_res_all_zScores.loc[pr_res_all_zScores["germline_Consequence"].str.contains("missense"), "expression_damaging_IMPACT"] = "HIGH"


plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline or somatic"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Protein disrupting variants \n missense, VEP high impact, AbExp ≤ -0.1,\npromoterAI ≤ -0.1, SpliceAI ≥ 0.1 (n={pr_res_all[pr_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


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
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_tsg_predisp_{expression_direction}_expression_disrupting_or_vs_zScore.png",  width=4, height=3.6, dpi=300)

    p = ( pn.ggplot(plot_data[(plot_data["Method"] == "PROTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    )) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_tsg_predisp_{expression_direction}_expression_disrupting.png",  width=4, height=3.6, dpi=300)



## 


pr_res_all["expression_damaging_IMPACT"] = "False"

pr_res_all.loc[pr_res_all["germline_IMPACT"].notna(), "expression_damaging_IMPACT"] = "HIGH"

pr_res_all.loc[pr_res_all["germline_Consequence"].str.contains("missense"), "expression_damaging_IMPACT"] = "HIGH"



plot_df = calculate_proportions(pr_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"


pr_res_all_zScores["expression_damaging_IMPACT"] = "False"

pr_res_all_zScores.loc[((pr_res_all_zScores["germline_IMPACT"].notna()) & (pr_res_all_zScores["germline_IMPACT"] == "HIGH")), "expression_damaging_IMPACT"] = "HIGH"
pr_res_all_zScores.loc[pr_res_all_zScores["germline_Consequence"].str.contains("missense"), "expression_damaging_IMPACT"] = "HIGH"


plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Protein disrupting variants \n missense, VEP high impact, AbExp ≤ -0.1,\npromoterAI ≤ -0.1, SpliceAI ≥ 0.1 (n={pr_res_all[pr_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


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
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_tsg_predisp_{expression_direction}_expression_disrupting_smallGermline_or_vs_zScore.png",  width=4, height=3.6, dpi=300)

    p = ( pn.ggplot(plot_data[(plot_data["Method"] == "PROTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    )) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_tsg_predisp_{expression_direction}_expression_disrupting_smallGermline.png",  width=4, height=3.6, dpi=300)



### alphamissense


pr_res_all["expression_damaging_IMPACT"] = "False"

pr_res_all.loc[(pr_res_all["germline_IMPACT"].notna() & (pr_res_all["am_class_x"] == "pathogenic")), "expression_damaging_IMPACT"] = "HIGH"
pr_res_all.loc[(pr_res_all["somatic_IMPACT"].notna() & (pr_res_all["am_class_y"] == "pathogenic")), "expression_damaging_IMPACT"] = "HIGH"




plot_df = calculate_proportions(pr_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"


pr_res_all_zScores["expression_damaging_IMPACT"] = "False"

pr_res_all_zScores.loc[(pr_res_all_zScores["germline_IMPACT"].notna()) & (pr_res_all_zScores["am_class_x"] == "pathogenic"), "expression_damaging_IMPACT"] = "HIGH"
pr_res_all_zScores.loc[(pr_res_all_zScores["somatic_IMPACT"].notna()) & (pr_res_all_zScores["am_class_y"] == "pathogenic"), "expression_damaging_IMPACT"] = "HIGH"

plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline or somatic"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"alphamissense pathogenic variants \n (n={pr_res_all[pr_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


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
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_germline_tsg_predisp_{expression_direction}_alphamissense_or_vs_zScore.png",  width=4, height=3.6, dpi=300)

    p = ( pn.ggplot(plot_data[(plot_data["Method"] == "PROTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    )) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/roportions_germline_tsg_predisp_{expression_direction}_alphamissense.png",  width=4, height=3.6, dpi=300)



##### CPT1 


pr_res_all["expression_damaging_IMPACT"] = "False"

pr_res_all["germline_CPT1_score"] = pd.to_numeric(
    pr_res_all["germline_CPT1_score"],
    errors="coerce"
)
pr_res_all["somatic_CPT1_score"] = pd.to_numeric(
    pr_res_all["somatic_CPT1_score"],
    errors="coerce"
)
pr_res_all.loc[pr_res_all["germline_IMPACT"].notna() & (pr_res_all["germline_CPT1_score"] >= 0.2), "expression_damaging_IMPACT"] = "HIGH"
pr_res_all.loc[pr_res_all["somatic_IMPACT"].notna() & (pr_res_all["somatic_CPT1_score"] >= 0.2), "expression_damaging_IMPACT"] = "HIGH"




plot_df = calculate_proportions(pr_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"



pr_res_all_zScores["germline_CPT1_score"] = pd.to_numeric(
    pr_res_all_zScores["germline_CPT1_score"],
    errors="coerce"
)
pr_res_all_zScores["somatic_CPT1_score"] = pd.to_numeric(
    pr_res_all_zScores["somatic_CPT1_score"],
    errors="coerce"
)
pr_res_all_zScores["expression_damaging_IMPACT"] = "False"

pr_res_all_zScores.loc[(pr_res_all_zScores["germline_IMPACT"].notna()) & (pr_res_all_zScores["germline_CPT1_score"] >= 0.2), "expression_damaging_IMPACT"] = "HIGH"
pr_res_all_zScores.loc[(pr_res_all_zScores["somatic_IMPACT"].notna()) & (pr_res_all_zScores["somatic_CPT1_score"] >= 0.2), "expression_damaging_IMPACT"] = "HIGH"

plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline or somatic"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"CPT1 ≥ 0.5 variants \n (n={pr_res_all[pr_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


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
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_germline_tsg_predisp_{expression_direction}_cpt1_gt0_2_or_vs_zScore.png",  width=4, height=3.6, dpi=300)

    p = ( pn.ggplot(plot_data[(plot_data["Method"] == "PROTRIDER") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions", 
                                    )) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min", ymax=f"{expression_direction}_ci_max"), 
                                alpha=0.2, outline_type='none')
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_germline_tsg_predisp_{expression_direction}_cpt1_gt0_2.png",  width=4, height=3.6, dpi=300)





### overexpression
## 


dresden_somatic = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/somatic_aggregated_predisp_cgc.parquet")
                    .with_columns(
                            cs.starts_with("SpliceAI_pred_DS_").cast(pl.Float64, strict=False)
                    )
                   .filter((pl.col("somatic_abexp_v1.1") >= 0.1) | (pl.col("somatic_promoterAI") >= 0.1)
                          )
                    .collect(engine="streaming")
                  ).to_pandas().drop(columns=["ref", "ref"])
# dresden_somatic

dresden_germline = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/germline_aggregated_predisp_cgc_tinda.parquet")
                    .with_columns(
                            cs.starts_with("SpliceAI_pred_DS_").cast(pl.Float64, strict=False)
                    ).filter((pl.col("PID").is_in(sa["pid"])) & (
                               (pl.col("germline_abexp_v1.1") >= 0.1) | (pl.col("germline_promoterAI") >= 0.1))
                          )
                    .collect(engine="streaming")        
                  ).to_pandas()
# dresden_germline




pr_res_all["expression_damaging_IMPACT"] = "False"

pr_res_all.loc[(pr_res_all["germline_promoterAI"] >= 0.1) | (pr_res_all["germline_abexp_v1.1"] >= 0.1) | (pr_res_all["somatic_abexp_v1.1"] >= 0.1) | (pr_res_all["somatic_promoterAI"] >= 0.1) , "expression_damaging_IMPACT"] = "HIGH"




plot_df = calculate_proportions(pr_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"


pr_res_all_zScores["expression_damaging_IMPACT"] = "False"

pr_res_all_zScores.loc[(pr_res_all_zScores["germline_promoterAI"] >= 0.1) | (pr_res_all_zScores["germline_abexp_v1.1"] >= 0.1) | (pr_res_all_zScores["somatic_abexp_v1.1"] >= 0.1) | (pr_res_all_zScores["somatic_promoterAI"] >= 0.1) , "expression_damaging_IMPACT"] = "HIGH"

plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline"


plot_df2 = plot_df2.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Protein disrupting variants \n AbExp ≥ 0.1,\npromoterAI ≥ 0.1(n={pr_res_all[pr_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"




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
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
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

    p.save(f"/home/a379i/plots/protein_expression_proportion_plots/predisp_tsg/proportions_cgc_predisp_{expression_direction}_or_vs_zScore_somatic_germline.png",  width=4, height=3.6, dpi=300)








pr_res_all["expression_damaging_IMPACT"] = "False"

pr_res_all.loc[(pr_res_all["germline_promoterAI"] >= 0.1) | (pr_res_all["germline_abexp_v1.1"] >= 0.1) | (pr_res_all["somatic_promoterAI"] >= 0.1) | (pr_res_all["CNV"].isin(["Somatic_AMP", "Germline_DUP"])) , "expression_damaging_IMPACT"] = "HIGH"




plot_df = calculate_proportions(pr_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline"


pr_res_all_zScores["expression_damaging_IMPACT"] = "False"

pr_res_all_zScores.loc[(pr_res_all_zScores["germline_promoterAI"] >= 0.1) | (pr_res_all_zScores["germline_abexp_v1.1"] >= 0.1) | (pr_res_all_zScores["somatic_promoterAI"] >= 0.1) | (pr_res_all["CNV"].isin(["Somatic_AMP", "Germline_DUP"])) , "expression_damaging_IMPACT"] = "HIGH"

plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline"


plot_df2 = plot_df2.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"overexpression_proportions_VEP_expression_damaging": "overexpression_proportions"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_max_VEP_expression_damaging": "overexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"overexpression_ci_min_VEP_expression_damaging": "overexpression_ci_min"})

plot_df2["Method"] = "sf-normalized\nzScores"


plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Protein expression disrupting variants\npromoterAI ≥ 0.1, Somatic AMP\nGermline DUP (n={pr_res_all[pr_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"




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
                    #title="Across tier 3 genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.theme_bw(base_size=17)
                + pn.theme(
                     strip_text = pn.element_text(size=16),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/ASHG_poster/protein_proportions_cgc_predisp_{expression_direction}_or_vs_zScore_somatic_germline_cnv.png",  width=8, height=5.4, dpi=300)






pr_res_all["expression_damaging_IMPACT"] = "False"

pr_res_all.loc[(pr_res_all["germline_IMPACT"] == "HIGH") | (pr_res_all["somatic_promoterAI"] <= -0.1) | (pr_res_all["germline_promoterAI"] <= -0.1)
    | (pr_res_all["germline_max_spliceai_score"] >= 0.1) | (pr_res_all["somatic_max_spliceai_score"] >= 0.1)  | (pr_res_all["germline_abexp_v1.1"] <= -0.1) | (pr_res_all["somatic_abexp_v1.1"] <= -0.1) | (pr_res_all["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])), "expression_damaging_IMPACT"] = "HIGH"


plot_df = calculate_proportions(pr_res_all, "VEP", "expression_damaging")
plot_df["Variant_type"] = "Small germline or somatic"


pr_res_all_zScores["expression_damaging_IMPACT"] = "False"

pr_res_all_zScores.loc[(pr_res_all_zScores["germline_IMPACT"] == "HIGH") | (pr_res_all_zScores["somatic_promoterAI"] <= -0.1) | (pr_res_all_zScores["germline_promoterAI"] <= -0.1)
    | (pr_res_all_zScores["germline_max_spliceai_score"] >= 0.1) | (pr_res_all_zScores["somatic_max_spliceai_score"] >= 0.1)  | (pr_res_all_zScores["germline_abexp_v1.1"] <= -0.1) | (pr_res_all_zScores["somatic_abexp_v1.1"] <= -0.1) | (pr_res_all_zScores["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])), "expression_damaging_IMPACT"] = "HIGH"


plot_df2 = calculate_proportions(pr_res_all_zScores, "VEP", "expression_damaging")
plot_df2["Variant_type"] = "Small germline or somatic"


plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})



plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})


plot_df2["Method"] = "sf-normalized\nzScores"

plot_dt = pd.concat([plot_df, plot_df2])


plot_dt["facet_label"] = f"Protein expression disrupting variants\nVEP high impact, AbExp ≤ -0.1,\npromoterAI ≤ -0.1, SpliceAI ≥ 0.1 \n Deletion CNV (n={pr_res_all[pr_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"


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
                    y=f"Proportion of protein {expression_direction}\noutliers w/ rare SNV/indels or CNV",
                    x="Outlier rank cutoff",
                    # title="Across expressed protein coding genes"
                )
                + pn.scale_fill_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.scale_color_manual(values={
                        "PROTRIDER": "#4A7BB0", 
                        "sf-normalized\nzScores": "grey"
                })
                + pn.theme_bw(base_size=17)
                + pn.theme(
                     strip_text = pn.element_text(size=16),

                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    #figure_size=(8, 6),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/ASHG_poster/protein_proportions_germline_somatic_cnv_predisp_cgc_{expression_direction}_expression_disrupting_or_vs_zScore.png", width=8, height=5.4, dpi=300)