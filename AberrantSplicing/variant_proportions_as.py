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

# %%
init_notebook_mode(all_interactive=True)

# %%
import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

# %%
from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *

# %%
from scipy.integrate import trapezoid



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
        if expression_direction == "underexpression" and omic_type != "splicing":
            mask = df['zScore'] < 0
        elif omic_type != "splicing":
            mask = df['zScore'] > 0
        else:
            mask = df.index.notna()
            
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
        
    elif expression_direction == "overexpression":
        direction_mask = plot_data["zScore"] > 0
        plot_data = plot_data[plot_data["zScore"] > 0]
        CNV_label = "AMP/DUP CNV"
    else:
        direction_mask = plot_data["hgncSymbol"] == plot_data["hgncSymbol"]
        plot_data = plot_data.copy()
        CNV_label = ""

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



all_res = []
# padjust_Genes_with_rare_splice_variants
for cohort in sa[sa["Diag"] != "Unstranded_data"]["Diag"].unique():
    # if cohort == "Other":
    #     splicing_res_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/{cohort}/results_per_junction.tsv"
    # else:
    #     splicing_res_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/{cohort}/results_per_junction_dna_first.tsv"
    #| (pl.col("padjust_Genes_with_rare_splice_variants") <= 0.1)
    splicing_res_path = f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/{cohort}/results_all.tsv.gz"
    fr_res_new =( pl.scan_csv(splicing_res_path, 
                                schema_overrides={
                                    "deltaPsi": pl.Float64,
                                    "seqnames": pl.String,
                                    "padjust_Genes_with_rare_splice_variants": pl.Float64,
                                }, 
                                null_values=["NA"],
                                 separator="\t")
                    #.select("sampleID", "hgncSymbol", "pValue", "padjust", 'seqnames', 'start', 'end', 'strand', 'deltaPsi', 'counts', 'totalCounts', 'nonsplitCounts')
                    .filter(pl.col("pValue") < 0.5)
                    .filter(pl.col("hgncSymbol").is_in(genes_of_interest["gene_name"]))
                    .collect(engine="streaming")
                    ).rename({"start": "junction_start", "end": "junction_end"}).to_pandas()
    # fr_res_new[(fr_res_new["padjust"] <= 0.05) & (fr_res_new["padjust_Genes_to_test_on_all_samples"].notna())]
    fr_res_new = fr_res_new.merge(sa, left_on="sampleID", right_on="pid")
    all_res.append(fr_res_new)
all_res = pd.concat(all_res)

print(all_res.shape, "all splicing outliers")
print(all_res["hgncSymbol"].unique())

splicing_outliers_merged = all_res.merge(dresden_germline,  right_on=["PID", "SYMBOL"], left_on=["sampleID", "hgncSymbol"], how="left")
splicing_outliers_merged = splicing_outliers_merged.merge(dresden_somatic,  right_on=["PID", "SYMBOL"], left_on=["sampleID", "hgncSymbol"], how="left")



splicing_outliers_merged["germline_somatic_IMPACT"] = "False"
splicing_outliers_merged.loc[(splicing_outliers_merged["germline_IMPACT"].notna()) | (splicing_outliers_merged["somatic_IMPACT"].notna()), "germline_somatic_IMPACT"] = "HIGH"


splicing_outliers_merged["CNV"] = "No CNV"

splicing_outliers_merged = splicing_outliers_merged.drop_duplicates(subset=["sampleID", "junction_start", "junction_end", "hgncSymbol"])
splicing_outliers_merged = splicing_outliers_merged.sort_values("pValue")
splicing_outliers_merged["Method"] = "FRASER"

print(splicing_outliers_merged.shape, "splicng outliers after merging with variants")


# %%

# %%

# %%
# pr_res_all = calculate_proportions(pr_res_all, expression_col="CNV")

omic_type = "splicing"
plot_df = calculate_proportions(splicing_outliers_merged, "VEP", "germline_somatic",  omic_type=omic_type)
plot_df["Variant_type"] = "Small germline or somatic"

plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_germline_somatic": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_germline_somatic": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_germline_somatic": "underexpression_ci_min"})

# %%
plot_dt = plot_df.copy()

var_num_df = plot_dt.copy()
var_num_df = var_num_df.drop_duplicates(subset=["sampleID", "hgncSymbol"])
print(var_num_df.shape)

plot_dt["facet_label"] = f"Splice disrupting variants \n VEP splice-related, SpliceAI ≥ 0.1,\nPangolin ≥ 0.1, AbSplice2 ≥ 0.1 (n={var_num_df[var_num_df["germline_somatic_IMPACT"] == "HIGH"].shape[0]})"

#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["underexpression"]:
    if expression_direction == "underexpression" and omic_type != "splicing":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    elif omic_type != "splicing":
        plot_data = plot_dt[plot_dt["zScore"] > 0]
    else:
        plot_data = plot_dt.copy()


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
                        "FRASER": "#fc8d59", 
                    })
                  + pn.scale_color_manual(values={
                        "FRASER": "#fc8d59", 
                    })
                
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of {omic_type} outliers\nwith rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.theme_bw(base_size=17) 
                + pn.theme(
                    strip_text = pn.element_text(size=16),
                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    # figure_size=(5, 5),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/ASHG_poster/splicing_proportions_cgc_predisp_vep_small_somatic_germline_fraser.png", width=8, height=5.4, dpi=300)





splicing_outliers_merged["germline_somatic_IMPACT"] = "False"
splicing_outliers_merged.loc[(splicing_outliers_merged["germline_IMPACT"].notna()), "germline_somatic_IMPACT"] = "HIGH"


splicing_outliers_merged["CNV"] = "No CNV"

splicing_outliers_merged = splicing_outliers_merged.drop_duplicates(subset=["sampleID", "junction_start", "junction_end", "hgncSymbol"])
splicing_outliers_merged = splicing_outliers_merged.sort_values("pValue")
splicing_outliers_merged["Method"] = "FRASER"

print(splicing_outliers_merged.shape, "splicng outliers after merging with variants")




omic_type = "splicing"
plot_df = calculate_proportions(splicing_outliers_merged, "VEP", "germline_somatic",  omic_type=omic_type)
plot_df["Variant_type"] = "Small germline or somatic"

plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_germline_somatic": "underexpression_proportions"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_germline_somatic": "underexpression_ci_max"})
plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_germline_somatic": "underexpression_ci_min"})

# %%
plot_dt = plot_df.copy()

var_num_df = plot_dt.copy()
var_num_df = var_num_df.drop_duplicates(subset=["sampleID", "hgncSymbol"])
print(var_num_df.shape)

plot_dt["facet_label"] = f"Splice disrupting variants \n VEP splice-related, SpliceAI ≥ 0.1,\nPangolin ≥ 0.1, AbSplice2 ≥ 0.1 (n={var_num_df[var_num_df["germline_somatic_IMPACT"] == "HIGH"].shape[0]})"

#for expression_direction in ["underexpression", "overexpression"]:
for expression_direction in ["underexpression"]:
    if expression_direction == "underexpression" and omic_type != "splicing":
            plot_data = plot_dt[plot_dt["zScore"] < 0]
            CNV_label = "Deletion CNV"
            CNV_label = "AMP CNV"
    elif omic_type != "splicing":
        plot_data = plot_dt[plot_dt["zScore"] > 0]
    else:
        plot_data = plot_dt.copy()


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
                        "FRASER": "#fc8d59", 
                    })
                  + pn.scale_color_manual(values={
                        "FRASER": "#fc8d59", 
                    })
                
                
                + pn.facet_grid(". ~ facet_label")
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of {omic_type} outliers\nwith rare SNV/indels",
                    x="Outlier rank cutoff",
                    title="Across tier 3 genes"
                )
                + pn.theme_bw(base_size=17) 
                + pn.theme(
                    strip_text = pn.element_text(size=16),
                    # Center multi-line y-axis title vertically
                    axis_title_y=pn.element_text(ha='center'),
                    # figure_size=(5, 5),
                    plot_margin=0.05,
                )
    )

    p.save(f"/home/a379i/plots/ASHG_poster/splicing_proportions_cgc_predisp_vep_small_germline_fraser.png", width=8, height=5.4, dpi=300)



