# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: protrider_env
#     language: python
#     name: python3
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
import gc
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


# %%
import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *
from utils.util_functions import *


cohort = sys.argv[1]


# %%
def calculate_proportions(df, expression_col="CNV", variant_type="germline", omic_type = "gene"):
    # Work on a copy to prevent SettingWithCopyWarning and index fragmentation

    for expression_direction in ["overexpression", "underexpression"]:
        if expression_direction == "underexpression":
            if (df["Method"] == "sf-normalized zScores").any():
                df = df.sort_values("zScore", ascending=True) 
        elif expression_direction == "overexpression":
            if (df["Method"] == "sf-normalized zScores").any():
                df = df.sort_values("zScore", ascending=False)
        
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
dresden_somatic = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/somatic/somatic_aggregated_predisp_cgc.parquet")
                    .with_columns(
                            cs.starts_with("SpliceAI_pred_DS_").cast(pl.Float64, strict=False)
                    )
                   .filter((pl.col("somatic_abexp_v1.1").abs() >= 0.1) | (pl.col("somatic_promoterAI").abs() >= 0.1) | (pl.col("somatic_IMPACT").str.contains("HIGH")) |
                              (pl.col("SpliceAI_pred_DS_AL") >= 0.1) |  
                               (pl.col("SpliceAI_pred_DS_AG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DL") >= 0.1)
                          )
                    .collect(engine="streaming")
                    
                  ).to_pandas()
# dresden_somatic

dresden_somatic["somatic_max_spliceai_score"] = dresden_somatic[["SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_AG", 
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"
]].max(axis=1)

dresden_somatic = (
    dresden_somatic
    .drop_duplicates(subset=["Gene", "PID"])
)
                  

dresden_germline = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/dresden_annotations/germline/germline_aggregated_predisp_cgc_tinda.parquet")
                    .with_columns(
                            cs.starts_with("SpliceAI_pred_DS_").cast(pl.Float64, strict=False)
                    ).filter((pl.col("PID").is_in(sa["pid"])) & (
                               (pl.col("germline_abexp_v1.1") <= -0.1) | (pl.col("germline_promoterAI") <= -0.1) | (pl.col("germline_IMPACT").str.contains("HIGH")) | 
                              (pl.col("SpliceAI_pred_DS_AL") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_AG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DG") >= 0.1) | 
                                (pl.col("SpliceAI_pred_DS_DL") >= 0.1)
                            )
                          )
                    .collect(engine="streaming")        
                  ).to_pandas()
# dresden_germline


dresden_germline["germline_max_spliceai_score"] = dresden_germline[["SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_AG", 
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"
]].max(axis=1)

dresden_germline = (
    dresden_germline
    .drop_duplicates(subset=["Gene", "PID"])
)




# %%
cohorts = list(sa[sa["Diag"] != "Unstranded_data"]["Diag"].unique())

# %%
# cohort = "ACC"

# %%
for i in range(1):
    cohort = cohort
    cohort_samples = sa[sa["Diag"] == cohort]
    print(cohort)
    or_res_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust_cnv.parquet"
    needed_cols = ["sampleID", "zScore", "pValue", "padjust", "IMPACT", "geneID",  "geneID_short",
                "padjust_predisp", "padjust_predisp_extended", "CNV",] #"padjust_genes_of_interest"]
    # py_or_res_all = pl.scan_parquet(or_res_path).select(needed_cols).filter((pl.col("geneID_short").is_in(genes_of_interest[(genes_of_interest["cgc"].isna()) |(genes_of_interest["ROLE_IN_CANCER"].str.contains("TSG"))]["geneID_short"]))).collect().to_pandas()

    or_res_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust_cnv.parquet"

    py_or_res_all = pl.scan_parquet(or_res_path).select(needed_cols).filter(pl.col("sampleID").is_in(cohort_samples["pid"])).collect().to_pandas()


    py_or_res_all["Method"] = "OUTRIDER_all_cohorts"

    # py_or_res_all = py_or_res_all[(py_or_res_all["padjust_predisp_extended"] <= 0.05) | (py_or_res_all["padjust"] <= 0.05)]


    py_or_res_all = pd.merge(py_or_res_all, sa[["pid", "Diag", "seq_type", "Oncotree Code", "Oncotree Text", "nct_pid", "Tumorzellgehalt (Bioinformatik)"]], left_on="sampleID", right_on="pid")
    py_or_res_all = pd.merge(py_or_res_all, gene_annot_dt.drop(columns="geneID_short"), left_on="geneID", right_on="gene_id", how="left")

    py_or_res_all = py_or_res_all.sort_values("pValue")
    py_or_res_all = py_or_res_all.drop_duplicates(subset=["sampleID", "geneID_short"])

    # py_or_res_aberrant["key"] = py_or_res_aberrant["sampleID"] +  "." + py_or_res_aberrant["nct_pid"].str.split("_").str[1]
    py_or_res_all = py_or_res_all.merge(dresden_germline , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")
    py_or_res_all = py_or_res_all.merge(dresden_somatic , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")

    cohort_wise_res = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/cohort_wise/" + cohort + "/oht_lr_0_0001_epoc200_cpu/or_res_cnv.parquet")


    cohort_wise_res["Method"] = "OUTRIDER_cohort_wise_py"

    # py_or_res_all = py_or_res_all[(py_or_res_all["padjust_predisp_extended"] <= 0.05) | (py_or_res_all["padjust"] <= 0.05)]


    cohort_wise_res = pd.merge(cohort_wise_res, sa[["pid", "Diag", "seq_type", "Oncotree Code", "Oncotree Text", "nct_pid", "Tumorzellgehalt (Bioinformatik)"]], left_on="sampleID", right_on="pid")
    # cohort_wise_res = pd.merge(cohort_wise_res, gene_annot_dt, left_on="proteinID", right_on="gene_id", how="left")

    
    
    if "pValue" in cohort_wise_res.columns:
        cohort_wise_res = cohort_wise_res.sort_values("pValue")
    else:
        cohort_wise_res = cohort_wise_res.sort_values("PROTEIN_PVALUE")
    
    if "PROTEIN_ZSCORE" in cohort_wise_res.columns:
        cohort_wise_res = cohort_wise_res.rename(columns={"PROTEIN_ZSCORE": "zScore"})

    cohort_wise_res = cohort_wise_res.drop_duplicates(subset=["sampleID", "geneID_short"])

    # py_or_res_aberrant["key"] = py_or_res_aberrant["sampleID"] +  "." + py_or_res_aberrant["nct_pid"].str.split("_").str[1]
    cohort_wise_res = cohort_wise_res.merge(dresden_germline , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")
    cohort_wise_res = cohort_wise_res.merge(dresden_somatic , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")

    # cohort_wise_res = cohort_wise_res.rename(columns={"PROTEIN_ZSCORE": "zScore"})


    r_or_res = pd.read_parquet(f"/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/{cohort}/OUTRIDER_results_all_cnv.parquet")
    r_or_res["Method"] = "OUTRIDER_cohort_wise_R"

    # py_or_res_all = py_or_res_all[(py_or_res_all["padjust_predisp_extended"] <= 0.05) | (py_or_res_all["padjust"] <= 0.05)]


    r_or_res = pd.merge(r_or_res, sa[["pid", "Diag", "seq_type", "Oncotree Code", "Oncotree Text", "nct_pid", "Tumorzellgehalt (Bioinformatik)"]], left_on="sampleID", right_on="pid")
    # r_or_res = pd.merge(r_or_res, gene_annot_dt, left_on="proteinID", right_on="gene_id", how="left")

    r_or_res = r_or_res.sort_values("pValue")
    r_or_res = r_or_res.drop_duplicates(subset=["sampleID", "geneID_short"])

    # py_or_res_aberrant["key"] = py_or_res_aberrant["sampleID"] +  "." + py_or_res_aberrant["nct_pid"].str.split("_").str[1]
    r_or_res = r_or_res.merge(dresden_germline , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")
    r_or_res = r_or_res.merge(dresden_somatic , left_on=["geneID_short", "sampleID"], right_on=["Gene", "PID"], how="left")


    # venn digaram of outliers detected in cohort wise vs all cohorts
    cohort_wise_sig = cohort_wise_res[(cohort_wise_res["PROTEIN_PADJ"] <= 0.1)]
    all_cohorts_sig = py_or_res_all[(py_or_res_all["padjust"] <= 0.1)]

    if "geneID" in cohort_wise_sig.columns:
        cohort_wise_sig["gene_sample"] = cohort_wise_sig["geneID"] + "_" + cohort_wise_sig["sampleID"]
    else:
        cohort_wise_sig["gene_sample"] = cohort_wise_sig["proteinID"] + "_" + cohort_wise_sig["sampleID"]

    all_cohorts_sig["gene_sample"] = all_cohorts_sig["geneID"] + "_" + all_cohorts_sig["sampleID"]

    set_cohort = set(cohort_wise_sig["gene_sample"])  # adjust column name
    set_all = set(all_cohorts_sig["gene_sample"])

    plt.figure(figsize=(6, 4))
    venn2([set_cohort, set_all], set_labels=("cohort_wise_sig", "all_cohorts_sig"))
    plt.title(f"{cohort} gene expression outliers", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"/home/a379i/plots/{cohort}_or_outliers_venn.png", dpi=150)


    ## all cohorts
    py_or_res_all["expression_damaging_IMPACT"] = "False"

    py_or_res_all.loc[(py_or_res_all["germline_IMPACT"] == "HIGH") | (py_or_res_all["somatic_IMPACT"] == "HIGH") | (py_or_res_all["somatic_promoterAI"] <= -0.1) | (py_or_res_all["germline_promoterAI"] <= -0.1)
        | (py_or_res_all["germline_max_spliceai_score"] >= 0.1) | (py_or_res_all["somatic_max_spliceai_score"] >= 0.1)  | (py_or_res_all["germline_abexp_v1.1"] <= -0.1) | (py_or_res_all["somatic_abexp_v1.1"] <= -0.1) | (py_or_res_all["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])), "expression_damaging_IMPACT"] = "HIGH"


    plot_df = calculate_proportions(py_or_res_all, "VEP", "expression_damaging")
    plot_df["Variant_type"] = "Small germline or somatic"

    # cohort_wise python
    cohort_wise_res["expression_damaging_IMPACT"] = "False"

    cohort_wise_res.loc[(cohort_wise_res["germline_IMPACT"] == "HIGH") | (cohort_wise_res["somatic_IMPACT"] == "HIGH") | (cohort_wise_res["somatic_promoterAI"] <= -0.1) | (cohort_wise_res["germline_promoterAI"] <= -0.1)
        | (cohort_wise_res["germline_max_spliceai_score"] >= 0.1) | (cohort_wise_res["somatic_max_spliceai_score"] >= 0.1)  | (cohort_wise_res["germline_abexp_v1.1"] <= -0.1) | (cohort_wise_res["somatic_abexp_v1.1"] <= -0.1) | (cohort_wise_res["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])), "expression_damaging_IMPACT"] = "HIGH"




    plot_df2 = calculate_proportions(cohort_wise_res, "VEP", "expression_damaging")
    plot_df2["Variant_type"] = "Small germline or somatic"


    # cohort wise R
    r_or_res["expression_damaging_IMPACT"] = "False"

    r_or_res.loc[(r_or_res["germline_IMPACT"] == "HIGH") | (r_or_res["somatic_IMPACT"] == "HIGH") | (r_or_res["somatic_promoterAI"] <= -0.1) | (r_or_res["germline_promoterAI"] <= -0.1)
        | (r_or_res["germline_max_spliceai_score"] >= 0.1) | (r_or_res["somatic_max_spliceai_score"] >= 0.1)  | (r_or_res["germline_abexp_v1.1"] <= -0.1) | (r_or_res["somatic_abexp_v1.1"] <= -0.1) | (r_or_res["CNV"].isin(["Somatic_HDEL", "Germline_DEL"])), "expression_damaging_IMPACT"] = "HIGH"


    plot_df3 = calculate_proportions(r_or_res, "VEP", "expression_damaging")
    plot_df3["Variant_type"] = "Small germline or somatic"


        
    plot_df = plot_df.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
    plot_df = plot_df.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
    plot_df = plot_df.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})


    plot_df2 = plot_df2.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
    plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
    plot_df2 = plot_df2.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})


    plot_df3 = plot_df3.rename(columns = {f"underexpression_proportions_VEP_expression_damaging": "underexpression_proportions"})
    plot_df3 = plot_df3.rename(columns = {f"underexpression_ci_max_VEP_expression_damaging": "underexpression_ci_max"})
    plot_df3 = plot_df3.rename(columns = {f"underexpression_ci_min_VEP_expression_damaging": "underexpression_ci_min"})





    plot_dt = pd.concat([plot_df, plot_df2, plot_df3])


    plot_dt["facet_label"] = f"Expression disrupting variants\nVEP high impact, AbExp ≤ -0.1,\npromoterAI ≤ -0.1, SpliceAI ≥ 0.1 \n Deletion CNV (n={py_or_res_all[py_or_res_all["expression_damaging_IMPACT"] == "HIGH"].shape[0]})"



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
                    # + pn.scale_fill_manual(values={
                    #         "OUTRIDER": "#f8766d", 
                    #         "sf-normalized\nzScores": "grey"
                    # })
                #    + pn.scale_color_manual(values={
                #             "OUTRIDER": "#f8766d", 
                #             "sf-normalized\nzScores": "grey"
                #     })
                    
                    + pn.facet_grid(". ~ facet_label")
                    + pn.scale_x_log10()
                    + pn.annotation_logticks(sides="b")
                    + pn.labs(
                        y=f"Proportion of gene {expression_direction}\noutliers w/ rare SNV/indels",
                        x="Outlier rank cutoff",
                        title="Across all expressed genes " + cohort
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

        p.save(f"/home/a379i/plots/gene_expression_proportions/cohort_wise_vs_all_cohorts/{cohort}_proportions_all_genes_{expression_direction}_vep_small_somatic_germline_cnv.png", width=4, height=3.6, dpi=300)
        

        
        del (
            cohort_samples,
            py_or_res_all,
            cohort_wise_res,
            plot_df,
            plot_df2,
            plot_dt,
            plot_data,
            p,
        )

        gc.collect()
# %%

