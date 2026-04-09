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
#     name: protrider_env
# ---

# %%
import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode
import polars as pl
import sys
init_notebook_mode(all_interactive=True)


# %%
# !which pip

# %%
import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *


# %%
sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")


# %%
needed_cols = ["sampleID", "zScore", "pValue", "padjust", "IMPACT", "geneID",  "geneID_short",
               "#Uploaded_variation_snv", "IMPACT_snv", "ANNOTATION_control_snv", "Consequence_snv", "promoterAI_snv",
               "Location_indel", "IMPACT_indel", "ANNOTATION_control_indel", "Consequence_indel",
               "padjust_predisp", "padjust_predisp_extended", "CNV"]
py_or_res_all = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust.parquet",
                                columns=needed_cols)
# py_or_res_all = py_or_res_all[py_or_res_all["padjust_predisp_extended"].notna()]
py_or_res_all = pd.merge(py_or_res_all, sa[["pid", "Diag", "seq_type", "Oncotree Code"]], left_on="sampleID", right_on="pid")
py_or_res_all = pd.merge(py_or_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID", "ROLE_IN_CANCER"]], on="geneID", how="left")


# %%
# py_or_res_all = py_or_res_all[py_or_res_all['zScore'] < 0]
# py_or_res_all = py_or_res_all.sort_values("pValue")[:10000000]

# %%
# needed_cols = ["sampleID", "zScore", "geneID",  "geneID_short", "IMPACT"]
# gene_zscore = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/zscores/or_variants.csv",
#                          usecols=needed_cols)
# gene_zscore = gene_zscore[gene_zscore["geneID_short"].isin(extended_dresden_dt["geneID_short"])]
# gene_zscore = pd.merge(gene_zscore, sa[["pid", "Diag", "seq_type", "Oncotree Code"]], left_on="sampleID", right_on="pid")

# gene_zscore = gene_zscore[gene_zscore["zScore"] < 0][:10000000]

# %%
# gene_zscore.shape

# %%
py_or_res_all = py_or_res_all.sort_values("promoterAI_snv")
py_or_res_all["aberrant"] = False
py_or_res_all.loc[py_or_res_all["padjust"] <= 0.05, "aberrant"] = True

# %% [markdown]
# # PrmoterAI

# %%
under_rank = np.arange(1, len(py_or_res_all) + 1)
py_or_res_all['underexpression_rank'] = under_rank

under_vus = (py_or_res_all['aberrant'].to_numpy() == True)
x = np.cumsum(under_vus)

py_or_res_all['underexpression_proportions'] = x / under_rank


# %%
py_or_res_all = py_or_res_all.sort_values("promoterAI_snv", ascending=False)


# %%
py_or_res_all[(py_or_res_all["predisposition_gene"] == True) & (py_or_res_all["zScore"] > 0) & (py_or_res_all["padjust"] < 0.05) & (py_or_res_all["promoterAI_snv"] > 0)]

# %%
over_rank = np.arange(1, len(py_or_res_all) + 1)
py_or_res_all['overexpression_rank'] = over_rank

over_vus = (py_or_res_all['aberrant'].to_numpy() == True)
x = np.cumsum(over_vus)

py_or_res_all['overexpression_proportions'] = x / over_rank


# %%
py_or_res_all["facet_label"] = f"Underexpression outliers (n={py_or_res_all[(py_or_res_all["zScore"] < 0) & (py_or_res_all["aberrant"] == True)].shape[0]})"
p = (
    pn.ggplot(py_or_res_all[(py_or_res_all["promoterAI_snv"].notna()) & (py_or_res_all["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of variants \nwith underexpression outliers",
        x="Promoterai variant rank (sorted from -1 to 1)"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
py_or_res_all["facet_label"] = f"Overexpression outliers (n={py_or_res_all[(py_or_res_all["zScore"] > 0) & (py_or_res_all["aberrant"] == True)].shape[0]})"
p = (
    pn.ggplot(py_or_res_all[(py_or_res_all["promoterAI_snv"].notna()) & (py_or_res_all["overexpression_rank"] < 1e8)], pn.aes(x="overexpression_rank", 
                   y="overexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of variants \nwith overexpression outliers",
        x="Promoterai varaint rank (sorted from 1 to -1)"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# # quick check on wgs vs wes

# %%
p = (pn.ggplot(sa[sa["Diag"] != "Unstranded_data"]) + 
        pn.geom_bar(pn.aes(x="Diag", fill="seq_type")) +
        pn.theme_bw()+
        pn.theme(axis_text_x=pn.element_text(rotation=90))
    )
p

# %%
sa[sa["proteomics"].notna()]["seq_type"].value_counts()

# %%
p = (pn.ggplot(sa[(sa["Diag"] != "Unstranded_data") & (sa["proteomics"].notna())]) + 
        pn.geom_bar(pn.aes(x="Diag", fill="seq_type")) +
        pn.theme_bw()+
        pn.theme(axis_text_x=pn.element_text(rotation=90))
    )
p

# %% [markdown]
# # Transcriptomics

# %% [markdown]
# ### All types of vcfs

# %%

py_or_res_all.loc[:, "VUS_snv"] = False
py_or_res_all.loc[py_or_res_all["promoterAI_snv"] <= -0.1, "VUS_snv"] = True

py_or_res_all.loc[:, "over_VUS_snv"] = False
py_or_res_all.loc[py_or_res_all["promoterAI_snv"] >= 0.1, "over_VUS_snv"] = True

print(py_or_res_all[py_or_res_all['promoterAI_snv'].notna()]["VUS_snv"].value_counts())
py_or_res_all[py_or_res_all['promoterAI_snv'].notna()]["over_VUS_snv"].value_counts()


# %%
## all underexpression outliers
print("all underexpression outliers")
a = py_or_res_all[(py_or_res_all["padjust"] <= 0.05) & (py_or_res_all["zScore"] <= 0)]
print(a.shape) 

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["VUS"] == True]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["IMPACT_snv"] == "HIGH")) |
      ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["IMPACT_indel"] == "HIGH"))]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %%
py_or_res_all[(py_or_res_all["promoterAI_snv"] <= -0.1) & (py_or_res_all["predisposition_gene"] == True)]

# %%
py_or_res_all.loc[:, "VUS"] = False
py_or_res_all.loc[(py_or_res_all["IMPACT"] == "HIGH"), "VUS"] = True


print(py_or_res_all["VUS"].value_counts())

# %%
py_or_res_all['proportions'] = (
    py_or_res_all['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(py_or_res_all) + 1)
)

py_or_res_all['Rank'] = np.arange(1, len(py_or_res_all) + 1)

# Underexpression subset
under_idx = py_or_res_all['zScore'] < 0
under = py_or_res_all.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
py_or_res_all.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()

# %%
py_or_res_all["facet_label"] = f"VEP HIGH impact (n={py_or_res_all[py_or_res_all["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(py_or_res_all[(py_or_res_all["zScore"] < 0) & (py_or_res_all["underexpression_rank"] > 100) & (py_or_res_all["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of underexpression outliers\nwith rare (WGS, WES)  variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %% [markdown]
# ### size factor zScores method

# %%
gene_zscore.loc[:, "VUS"] = False
gene_zscore.loc[gene_zscore["IMPACT"] == "HIGH", "VUS"] = True


print(gene_zscore["VUS"].value_counts())


# %%
gene_zscore['proportions'] = (
    gene_zscore['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(gene_zscore) + 1)
)

gene_zscore['Rank'] = np.arange(1, len(gene_zscore) + 1)

# Underexpression subset
under_idx = gene_zscore['zScore'] < 0
under = gene_zscore.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
gene_zscore.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()

# %%
gene_zscore["Method"] = "Z-scores"
py_or_res_all["Method"] = "pyOUTRIDER"
plot_data = pd.concat((gene_zscore, py_or_res_all))

# %%
plot_data["facet_label"] = f"VEP HIGH impact (n={plot_data[plot_data["VUS"] == True].shape[0]//2})"
p = (
    pn.ggplot(plot_data[(plot_data["zScore"] < 0) & (plot_data["underexpression_rank"] > 100) & (plot_data["underexpression_rank"] < 2e6)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions", color="Method"))
    +   pn.geom_ribbon(
        pn.aes(
            fill = "Method",
            ymin="ci_min",
            ymax="ci_max"
        ),
        outline_type='none',
        alpha=0.2
    )
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        # title="Predisposition genes",
        y="Proportion of underexpression outliers\nwith rare variants",
        x="Outlier rank cutoff"
    )
     + pn.scale_color_manual(values={
        "Z-scores": "#00BFC4",
        "pyOUTRIDER": "#F8766D"
    })
    + pn.scale_fill_manual(values={
        "Z-scores": "#00BFC4",
        "pyOUTRIDER": "#F8766D"
    })

    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %% [markdown]
# ### WES analysis

# %%
py_or_res_all_wes = py_or_res_all[py_or_res_all["seq_type"] == "WES"]
print(py_or_res_all_wes.shape)


# %% [markdown]
# #### promoterAI

# %%
py_or_res_all_wes.loc[:, "VUS"] = False
py_or_res_all_wes.loc[py_or_res_all_wes["promoterAI_snv"] <= -0.1, "VUS"] = True

py_or_res_all_wes.loc[:, "over_VUS"] = False
py_or_res_all_wes.loc[py_or_res_all_wes["promoterAI_snv"] >= 0.1, "over_VUS"] = True

print(py_or_res_all_wes[py_or_res_all_wes['promoterAI_snv'].notna()]["VUS"].value_counts())
py_or_res_all_wes[py_or_res_all_wes['promoterAI_snv'].notna()]["over_VUS"].value_counts()


# %%
py_or_res_all_wes['proportions'] = (
    py_or_res_all_wes['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(py_or_res_all_wes) + 1)
)

py_or_res_all_wes['Rank'] = np.arange(1, len(py_or_res_all_wes) + 1)

# Underexpression subset
under_idx = py_or_res_all_wes['zScore'] < 0
under = py_or_res_all_wes.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
py_or_res_all_wes.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# Overexpression subset
over_idx = py_or_res_all_wes['zScore'] > 0
over = py_or_res_all_wes.loc[over_idx].copy()

over_rank = np.arange(1, len(over) + 1)
over['overexpression_rank'] = over_rank

over_vus = (over['over_VUS'].to_numpy() == True)
x = np.cumsum(over_vus)

over['overexpression_proportions'] = x / over_rank

# Vectorized Clopper–Pearson
alpha = 0.05
over['over_ci_min'] = beta.ppf(alpha/2, x, over_rank - x + 1)
over['over_ci_max'] = beta.ppf(1 - alpha/2, x + 1, over_rank - x)

# Efficient vectorized assignment (NO MERGE)
py_or_res_all_wes.loc[over_idx, [
    "overexpression_proportions",
    "overexpression_rank",
    "over_ci_min",
    "over_ci_max"
]] = over[
    ["overexpression_proportions", "overexpression_rank", "over_ci_min", "over_ci_max"]
].to_numpy()

# %%
## all underexpression outliers
print("all underexpression outliers")
a = py_or_res_all_wes[(py_or_res_all_wes["padjust"] <= 0.05) & (py_or_res_all_wes["Outlier status"] == "Underexpression")]
print(a.shape)

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["VUS"] == True]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["IMPACT_snv"] == "HIGH")) |
      ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["IMPACT_indel"] == "HIGH"))]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %%
c

# %%
py_or_res_all_wes["facet_label"] = f"promoterAI < -0.1 (n={py_or_res_all_wes[py_or_res_all_wes["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(py_or_res_all_wes[(py_or_res_all_wes["zScore"] < 0) & (py_or_res_all_wes["underexpression_rank"] > 100) & (py_or_res_all_wes["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of underexpression outliers\nwith rare WES variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
py_or_res_all_wes["facet_label"] = f"promoterAI > 0.1 (n={py_or_res_all_wes[py_or_res_all_wes["over_VUS"] == True].shape[0]})"

p2 = (
    pn.ggplot(py_or_res_all_wes[(py_or_res_all_wes["zScore"] > 0) & (py_or_res_all_wes["overexpression_rank"] > 100) & (py_or_res_all_wes["overexpression_rank"] < 1e8)], pn.aes(x="overexpression_rank", 
                   y="overexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of oevrexpression outliers\nwith rare WES variant",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p2

# %% [markdown]
# #### VEP high impact

# %%
py_or_res_all_wes.loc[:, "VUS"] = False
py_or_res_all_wes.loc[py_or_res_all_wes["IMPACT"] == "HIGH", "VUS"] = True


print(py_or_res_all_wes["VUS"].value_counts())

# %%
## all underexpression outliers
print("all underexpression outliers")
a = py_or_res_all_wes[(py_or_res_all_wes["padjust"] <= 0.05) & (py_or_res_all_wes["Outlier status"] == "Underexpression")]
print(a.shape)

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["VUS"] == True]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["IMPACT_snv"] == "HIGH")) |
      ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["IMPACT_indel"] == "HIGH"))]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %%
d 

# %%
py_or_res_all_wes['proportions'] = (
    py_or_res_all_wes['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(py_or_res_all_wes) + 1)
)

py_or_res_all_wes['Rank'] = np.arange(1, len(py_or_res_all_wes) + 1)

# Underexpression subset
under_idx = py_or_res_all_wes['zScore'] < 0
under = py_or_res_all_wes.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
py_or_res_all_wes.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# %%
py_or_res_all_wes["facet_label"] = f"VEP HIGH impact (n={py_or_res_all_wes[py_or_res_all_wes["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(py_or_res_all_wes[(py_or_res_all_wes["zScore"] < 0) & (py_or_res_all_wes["underexpression_rank"] > 100) & (py_or_res_all_wes["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of underexpression outliers\nwith rare WES variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %% [markdown]
# ### WGS

# %%
py_or_res_all_wgs = py_or_res_all[py_or_res_all["seq_type"] == "WGS"]


# %% [markdown]
# #### promoterAI
#

# %%
py_or_res_all_wgs.loc[:, "VUS"] = False
py_or_res_all_wgs.loc[py_or_res_all_wgs["promoterAI_snv"] <= -0.1, "VUS"] = True

py_or_res_all_wgs.loc[:, "over_VUS"] = False
py_or_res_all_wgs.loc[py_or_res_all_wgs["promoterAI_snv"] >= 0.1, "over_VUS"] = True

print(py_or_res_all_wgs[py_or_res_all_wgs['promoterAI_snv'].notna()]["VUS"].value_counts())
py_or_res_all_wgs[py_or_res_all_wgs['promoterAI_snv'].notna()]["over_VUS"].value_counts()


# %%
py_or_res_all_wgs['proportions'] = (
    py_or_res_all_wgs['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(py_or_res_all_wgs) + 1)
)

py_or_res_all_wgs['Rank'] = np.arange(1, len(py_or_res_all_wgs) + 1)

# Underexpression subset
under_idx = py_or_res_all_wgs['zScore'] < 0
under = py_or_res_all_wgs.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
py_or_res_all_wgs.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# Overexpression subset
over_idx = py_or_res_all_wgs['zScore'] > 0
over = py_or_res_all_wgs.loc[over_idx].copy()

over_rank = np.arange(1, len(over) + 1)
over['overexpression_rank'] = over_rank

over_vus = (over['over_VUS'].to_numpy() == True)
x = np.cumsum(over_vus)

over['overexpression_proportions'] = x / over_rank

# Vectorized Clopper–Pearson
alpha = 0.05
over['over_ci_min'] = beta.ppf(alpha/2, x, over_rank - x + 1)
over['over_ci_max'] = beta.ppf(1 - alpha/2, x + 1, over_rank - x)

# Efficient vectorized assignment (NO MERGE)
py_or_res_all_wgs.loc[over_idx, [
    "overexpression_proportions",
    "overexpression_rank",
    "over_ci_min",
    "over_ci_max"
]] = over[
    ["overexpression_proportions", "overexpression_rank", "over_ci_min", "over_ci_max"]
].to_numpy()

# %%
py_or_res_all_wgs["facet_label"] = f"promoterAI < -0.1 (n={py_or_res_all_wgs[py_or_res_all_wgs["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(py_or_res_all_wgs[(py_or_res_all_wgs["zScore"] < 0) & (py_or_res_all_wgs["underexpression_rank"] > 100) & (py_or_res_all_wgs["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of underexpression outliers\nwith rare WGS variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
py_or_res_all_wgs["facet_label"] = f"promoterAI > 0.1 (n={py_or_res_all_wgs[py_or_res_all_wgs["over_VUS"] == True].shape[0]})"

p2 = (
    pn.ggplot(py_or_res_all_wgs[(py_or_res_all_wgs["zScore"] > 0) & (py_or_res_all_wgs["overexpression_rank"] > 100) & (py_or_res_all_wgs["overexpression_rank"] < 1e8)], pn.aes(x="overexpression_rank", 
                   y="overexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of oevrexpression outliers\nwith rare WGS variant",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p2

# %% [markdown]
# #### VEP HIHG IMPACT

# %%
py_or_res_all_wgs.loc[:, "VUS"] = False
py_or_res_all_wgs.loc[py_or_res_all_wgs["IMPACT"] == "HIGH", "VUS"] = True


print(py_or_res_all_wgs["VUS"].value_counts())

# %%
## all underexpression outliers
print("all underexpression outliers")
a = py_or_res_all_wgs[(py_or_res_all_wgs["padjust"] <= 0.1) & (py_or_res_all_wgs["Outlier status"] == "Underexpression")]
print(a.shape)

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["VUS"] == True]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["IMPACT_snv"] == "HIGH")) |
      ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["IMPACT_indel"] == "HIGH"))]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %%
py_or_res_all_wgs['proportions'] = (
    py_or_res_all_wgs['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(py_or_res_all_wgs) + 1)
)

py_or_res_all_wgs['Rank'] = np.arange(1, len(py_or_res_all_wgs) + 1)

# Underexpression subset
under_idx = py_or_res_all_wgs['zScore'] < 0
under = py_or_res_all_wgs.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
py_or_res_all_wgs.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# %%
py_or_res_all_wgs["facet_label"] = f"VEP HIGH impact (n={py_or_res_all_wgs[py_or_res_all_wgs["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(py_or_res_all_wgs[(py_or_res_all_wgs["zScore"] < 0) & (py_or_res_all_wgs["underexpression_rank"] > 100) & (py_or_res_all_wgs["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of underexpression outliers\nwith rare WGS variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %% [markdown]
# ### Combined WGS and WES

# %%
### check differences of WGS and WES
combined_res = pd.concat([py_or_res_all_wes, py_or_res_all_wgs])


# %%
combined_res["facet_label"] = f"VEP HIGH impact (n={combined_res[combined_res["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(combined_res[(combined_res["zScore"] < 0) & (combined_res["underexpression_rank"] > 100) & (combined_res["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions", color="seq_type"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of underexpression outliers\nwith rare variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
# Filter
df_filt = py_or_res_all[
    (py_or_res_all["promoterAI"] <= -0.1) &
    (py_or_res_all["padjust"] <= 0.05) &
    (py_or_res_all["zScore"] <= 0)
]

# Count consequences
counts = (
    df_filt["Consequence"]
    .value_counts()
    .reset_index()
)
counts.columns = ["Consequence", "count"]

# Drop singletons
counts = counts[counts["count"] > 1]

# Bar plot
p = (
    pn.ggplot(counts, pn.aes(x="reorder(Consequence, -count)", y="count"))
    + pn.geom_col()
    + pn.theme_bw()
    + pn.theme(
        axis_text_x=pn.element_text(rotation=45, ha="right")
    )
    + pn.labs(
        x="Consequence",
        y="Underexpression Outliers w/ promoterAI <= -0.1"
    )
)

p

# %%
df_filt['IMPACT'].value_counts()

# %%
# Filter
df_filt = py_or_res_all[
    (py_or_res_all["promoterAI"] <= -0.1) &
    (py_or_res_all["padjust"] <= 0.05) & 
    (py_or_res_all["zScore"] <= 0)
]

# Count consequences
counts = (
    df_filt["IMPACT"]
    .value_counts()
    .reset_index()
)
counts.columns = ["IMPACT", "count"]

# Drop singletons
counts = counts[counts["count"] > 1]

# Bar plot
p = (
    pn.ggplot(counts, pn.aes(x="reorder(IMPACT, -count)", y="count"))
    + pn.geom_col(width=0.5)
    + pn.theme_bw()
    + pn.theme(
        axis_text_x=pn.element_text(rotation=45, ha="right")
    )
    + pn.labs(
        x="IMPACT",
        y="Underexpression Outliers w/ promoterAI <= -0.1"
    )
)

p

# %%
py_or_res_aberrant_predisposed_all[
    (py_or_res_aberrant_predisposed_all["promoterAI"] <= -0.1) &
    (py_or_res_aberrant_predisposed_all["padjust"] <= 0.05) & (py_or_res_aberrant_predisposed_all["predisposition_gene"] == True)]["IMPACT"]

# %%
df_filt['ANNOTATION_control'].value_counts()

# %%
# Filter
df_filt = py_or_res_all[
    (py_or_res_all["promoterAI"] <= -0.1) &
    (py_or_res_all["padjust"] <= 0.05) &
    (py_or_res_all["zScore"] <= 0)
]

# Count consequences
counts = (
    df_filt["ANNOTATION_control"]
    .value_counts()
    .reset_index()
)
counts.columns = ["ANNOTATION_control", "count"]

# Drop singletons
counts = counts[counts["count"] > 1]

# Bar plot
p = (
    pn.ggplot(counts, pn.aes(x="reorder(ANNOTATION_control, -count)", y="count"))
    + pn.geom_col(width=0.5)
    + pn.theme_bw()
    + pn.theme(
        axis_text_x=pn.element_text(rotation=45, ha="right")
    )
    + pn.labs(
        x="IMPACT",
        y="Underexpression Outliers w/ promoterAI <= -0.1"
    )
)

p

# %%
py_or_res_aberrant_predisposed_all = pd.merge(py_or_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID", "ROLE_IN_CANCER", "predisposition_gene"]], on="geneID", how="left")


# %%
py_or_res_aberrant_predisposed_all[py_or_res_aberrant_predisposed_all["predisposition_gene"].notna()]

# %%
import gzip
import pandas as pd

path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_snv_rare/snvs_H021-019E.vep_res_snv_rare_filtered.vcf.gz"

with gzip.open(path, "rt") as f:
    for line in f:
        if line.startswith("#Uploaded_variation"):
            header = line.lstrip("#").strip().split("\t")
            break

vep_res_sample = pd.read_csv(
    path,
    sep="\t",
    comment="#",
    names=header,
    compression="gzip"
)


# %%
vep_res_sample["AF"].value_counts()

# %% [markdown]
# # Proteomics

# %%
pr_output_name = "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"

pr_res_all = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants.csv")
pr_res_all = pd.merge(pr_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID_short", "ROLE_IN_CANCER", "predisposition_gene"]], right_on="geneID_short", left_on="geneID", how="left")
pr_res_all = pd.merge(pr_res_all, sa, left_on="sampleID", right_on="pid")

# %%
pr_output_name_2 = "zScore_gt3"

pr_res_all_t = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name_2 + "/pr_variants.csv")
pr_res_all_t = pd.merge(pr_res_all_t, dresden_dt_cgc[["gene_name", "gene_type", "geneID_short", "ROLE_IN_CANCER", "predisposition_gene"]], right_on="geneID_short", left_on="geneID", how="left")
pr_res_all_t = pd.merge(pr_res_all_t, sa, left_on="sampleID", right_on="pid")

# %%

# %% [markdown]
# ## all variants

# %%
protein_disrupting_variants = ["stop_gained", "missense_variant", "stop_lost", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant"]

pr_res_all["VUS"] = False

mask_snv = (
    pr_res_all["Consequence_snv"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

mask_indel = (
    pr_res_all["Consequence_indel"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

pr_res_all.loc[mask_snv, "VUS"] = True
pr_res_all.loc[mask_indel, "VUS"] = True

print(pr_res_all["VUS"].value_counts())

pr_res_all.loc[pr_res_all["IMPACT"] == "HIGH", "VUS"] = True

print(pr_res_all["VUS"].value_counts())



# %%
pr_res_all['proportions'] = (
    pr_res_all['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(pr_res_all) + 1)
)

pr_res_all['Rank'] = np.arange(1, len(pr_res_all) + 1)

# Underexpression subset
under_idx = pr_res_all['zScore'] < 0
under = pr_res_all.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# %%
protein_disrupting_variants = ["stop_gained", "missense_variant", "stop_lost", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant"]

pr_res_all_t["VUS"] = False

mask_snv = (
    pr_res_all_t["Consequence_snv"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

mask_indel = (
    pr_res_all_t["Consequence_indel"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

pr_res_all_t.loc[mask_snv, "VUS"] = True
pr_res_all_t.loc[mask_indel, "VUS"] = True

print(pr_res_all_t["VUS"].value_counts())

pr_res_all_t.loc[pr_res_all_t["IMPACT"] == "HIGH", "VUS"] = True

print(pr_res_all_t["VUS"].value_counts())



# %%
pr_res_all_t['proportions'] = (
    pr_res_all_t['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(pr_res_all_t) + 1)
)

pr_res_all_t['Rank'] = np.arange(1, len(pr_res_all_t) + 1)

# Underexpression subset
under_idx = pr_res_all_t['zScore'] < 0
under = pr_res_all_t.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all_t.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# %%
pr_res_all["Method"] = "PROTRIDER"
pr_res_all_t["Method"] = "Z-scores"
to_plot = pd.concat((pr_res_all[pr_res_all["underexpression_rank"] <= 1e6], pr_res_all_t[pr_res_all_t["underexpression_rank"] <= 1e6]))


# %%
to_plot["facet_label"] = f"Protein disrupting variants (n={to_plot[to_plot["VUS"] == True].shape[0]//2})"
p = (
    pn.ggplot(to_plot[(to_plot["zScore"] < 0) & (to_plot["underexpression_rank"] > 100) & (to_plot["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions", color="Method"))
     +   pn.geom_ribbon(
        pn.aes(
            fill = "Method",
            ymin="ci_min",
            ymax="ci_max"
        ),
        outline_type='none',
        alpha=0.2
    )
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein underexpression outliers\nwith rare variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
pr_res_all["facet_label"] = f"Protein disrupting variants (n={pr_res_all[pr_res_all["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(pr_res_all[(pr_res_all["zScore"] < 0) & (pr_res_all["underexpression_rank"] > 100) & (pr_res_all["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein underexpression outliers\nwith rare variants",
        x="Outlier rank cutoff"
    ) +
    pn.geom_ribbon(
        pn.aes(
            ymin="ci_min",
            ymax="ci_max"
        ),
        alpha=0.2
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%

# %%
print(pr_res_all[(pr_res_all["VUS"] == True) & (pr_res_all["seq_type"] == "WGS")].shape)
pr_res_all[(pr_res_all["VUS"] == True) & (pr_res_all["seq_type"] == "WES")].shape

# %%
## all underexpression outliers
print("all underexpression outliers")
a = pr_res_all[(pr_res_all["padjust"] <= 0.1) & (pr_res_all["Outlier status"] == "Underexpression")]
print(a.shape)

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["IMPACT"] == "HIGH"]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
d = c[c["ANNOTATION_control"].str.contains("germline", na=False)]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %% [markdown]
# ## WES

# %%
pr_res_all_wes = pr_res_all[pr_res_all["seq_type"] == "WES"]
print(pr_res_all_wes.shape)


# %% [markdown]
# ### All types of variants

# %% [markdown]
# ### VEP protein disrupting variants

# %%
protein_disrupting_variants = ["stop_gained", "missense_variant", "stop_lost", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant"]

pr_res_all_wes["VUS"] = False

mask_snv = (
    pr_res_all_wes["Consequence_snv"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

mask_indel = (
    pr_res_all_wes["Consequence_indel"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

pr_res_all_wes.loc[mask_snv, "VUS"] = True
pr_res_all_wes.loc[mask_indel, "VUS"] = True

print(pr_res_all_wes["VUS"].value_counts())

pr_res_all_wes.loc[pr_res_all_wes["IMPACT"] == "HIGH", "VUS"] = True

print(pr_res_all_wes["VUS"].value_counts())

# %%
## all underexpression outliers
print("all underexpression outliers")
a = pr_res_all_wes[(pr_res_all_wes["padjust"] <= 0.1) & (pr_res_all_wes["Outlier status"] == "Underexpression")]
print(a.shape)

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["VUS"] == True]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
# d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["IMPACT_snv"] == "HIGH")) |
#       ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["IMPACT_indel"] == "HIGH"))]
d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["VUS"] == True)) |
       ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["VUS"] == True))]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %%
pr_res_all_wes['proportions'] = (
    pr_res_all_wes['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(pr_res_all_wes) + 1)
)

pr_res_all_wes['Rank'] = np.arange(1, len(pr_res_all_wes) + 1)

# Underexpression subset
under_idx = pr_res_all_wes['zScore'] < 0
under = pr_res_all_wes.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all_wes.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# %%
pr_res_all_wes["facet_label"] = f"Protein disrupting variants (n={pr_res_all_wes[pr_res_all_wes["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(pr_res_all_wes[(pr_res_all_wes["zScore"] < 0) & (pr_res_all_wes["underexpression_rank"] > 100) & (pr_res_all_wes["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein underexpression outliers\nwith rare WES variants",
        x="Outlier rank cutoff"
    ) +
    + pn.geom_ribbon(
        pn.aes(
            ymin="ci_min",
            ymax="ci_max"
        ),
        alpha=0.2
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# ### PromoterAI

# %%
pr_res_all_wes.loc[:, "VUS"] = False
pr_res_all_wes.loc[pr_res_all_wes["promoterAI_snv"] <= -0.1, "VUS"] = True

pr_res_all_wes.loc[:, "over_VUS"] = False
pr_res_all_wes.loc[pr_res_all_wes["promoterAI_snv"] >= 0.1, "over_VUS"] = True

print(pr_res_all_wes[pr_res_all_wes['promoterAI_snv'].notna()]["VUS"].value_counts())
pr_res_all_wes[pr_res_all_wes['promoterAI_snv'].notna()]["over_VUS"].value_counts()


# %%
pr_res_all_wes['proportions'] = (
    pr_res_all_wes['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(pr_res_all_wes) + 1)
)

pr_res_all_wes['Rank'] = np.arange(1, len(pr_res_all_wes) + 1)

# Underexpression subset
under_idx = pr_res_all_wes['zScore'] < 0
under = pr_res_all_wes.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all_wes.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# Overexpression subset
over_idx = pr_res_all_wes['zScore'] > 0
over = pr_res_all_wes.loc[over_idx].copy()

over_rank = np.arange(1, len(over) + 1)
over['overexpression_rank'] = over_rank

over_vus = (over['over_VUS'].to_numpy() == True)
x = np.cumsum(over_vus)

over['overexpression_proportions'] = x / over_rank

# Vectorized Clopper–Pearson
alpha = 0.05
over['over_ci_min'] = beta.ppf(alpha/2, x, over_rank - x + 1)
over['over_ci_max'] = beta.ppf(1 - alpha/2, x + 1, over_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all_wes.loc[over_idx, [
    "overexpression_proportions",
    "overexpression_rank",
    "over_ci_min",
    "over_ci_max"
]] = over[
    ["overexpression_proportions", "overexpression_rank", "over_ci_min", "over_ci_max"]
].to_numpy()

# %%
pr_res_all_wes["facet_label"] = f"promoterAI < -0.1 (n={pr_res_all_wes[pr_res_all_wes["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(pr_res_all_wes[(pr_res_all_wes["zScore"] < 0) & (pr_res_all_wes["underexpression_rank"] > 100) & (pr_res_all_wes["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein underexpression outliers\nwith rare WES variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
pr_res_all_wes["facet_label"] = f"promoterAI > 0.1 (n={pr_res_all_wes[pr_res_all_wes["over_VUS"] == True].shape[0]})"

p2 = (
    pn.ggplot(pr_res_all_wes[(pr_res_all_wes["zScore"] > 0) & (pr_res_all_wes["overexpression_rank"] > 100) & (pr_res_all_wes["overexpression_rank"] < 1e8)], pn.aes(x="overexpression_rank", 
                   y="overexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein oevrexpression outliers\nwith rare WGS variant",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p2

# %% [markdown]
# ## WGS

# %%
pr_res_all_wgs = pr_res_all[pr_res_all["seq_type"] == "WGS"]
print(pr_res_all_wgs.shape)


# %% [markdown]
# ### VEP protein disrupting genes

# %%
protein_disrupting_variants = ["stop_gained", "missense_variant", "stop_lost", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant"]

pr_res_all_wgs["VUS"] = False

mask_snv = (
    pr_res_all_wgs["Consequence_snv"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

mask_indel = (
    pr_res_all_wgs["Consequence_indel"]
    .fillna("")
    .str.split(",")
    .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
)

pr_res_all_wgs.loc[mask_snv, "VUS"] = True
pr_res_all_wgs.loc[mask_indel, "VUS"] = True
pr_res_all_wgs.loc[pr_res_all_wgs["IMPACT"] == "HIGH", "VUS"] = True


print(pr_res_all_wgs["VUS"].value_counts())

# %%
## all underexpression outliers
print("all underexpression outliers")
a = pr_res_all_wgs[(pr_res_all_wgs["padjust"] <= 0.1) & (pr_res_all_wgs["Outlier status"] == "Underexpression")]
print(a.shape)

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["VUS"] == True]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
# d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["IMPACT_snv"] == "HIGH")) |
#       ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["IMPACT_indel"] == "HIGH"))]
d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["VUS"] == True)) |
       ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["VUS"] == True))]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %%
pr_res_all_wgs['proportions'] = (
    pr_res_all_wgs['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(pr_res_all_wgs) + 1)
)

pr_res_all_wgs['Rank'] = np.arange(1, len(pr_res_all_wgs) + 1)

# Underexpression subset
under_idx = pr_res_all_wgs['zScore'] < 0
under = pr_res_all_wgs.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all_wgs.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# %%
pr_res_all_wgs["facet_label"] = f"Protein disrupting variants (n={pr_res_all_wgs[pr_res_all_wgs["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(pr_res_all_wgs[(pr_res_all_wgs["zScore"] < 0) & (pr_res_all_wgs["underexpression_rank"] > 100) & (pr_res_all_wgs["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein underexpression outliers\nwith rare WGS variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %% [markdown]
# ### promoter AI

# %%
pr_res_all_wgs.loc[:, "VUS"] = False
pr_res_all_wgs.loc[pr_res_all_wgs["promoterAI_snv"] <= -0.1, "VUS"] = True

pr_res_all_wgs.loc[:, "over_VUS"] = False
pr_res_all_wgs.loc[pr_res_all_wgs["promoterAI_snv"] >= 0.1, "over_VUS"] = True

print(pr_res_all_wgs[pr_res_all_wgs['promoterAI_snv'].notna()]["VUS"].value_counts())
pr_res_all_wgs[pr_res_all_wgs['promoterAI_snv'].notna()]["over_VUS"].value_counts()


# %%
## all underexpression outliers
print("all underexpression outliers")
a = pr_res_all_wgs[(pr_res_all_wgs["padjust"] <= 0.1) & (pr_res_all_wgs["Outlier status"] == "Underexpression")]
print(a.shape)

## underexpression outliers in predisposition genes
print("underexpression outliers in predisposition genes")
b = a[(a["predisposition_gene"] == True)]
print(b.shape, b.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant
print("underexpression outliers in predisposition genes with supporting variant")
c = b[b["VUS"] == True]
print(c.shape, c.shape[0]/a.shape[0])

## underexpression outliers in predisposition genes with supporting variant in germline
print("underexpression outliers in predisposition genes with supporting variant in germline")
# d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["IMPACT_snv"] == "HIGH")) |
#       ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["IMPACT_indel"] == "HIGH"))]
d = c[((c["ANNOTATION_control_snv"].str.contains("germline", na=False)) & (c["VUS"] == True)) |
       ((c["ANNOTATION_control_indel"].str.contains("germline", na=False)) & (c["VUS"] == True))]
print(d.shape, d.shape[0]/a.shape[0])

# two genes not in CGC --> could be interesting: TRIM37 in ACC and DDX41 in Lung (NSCLC), both stop_gained



# %%
b

# %%
pr_res_all_wgs['proportions'] = (
    pr_res_all_wgs['VUS'].eq(True).fillna(False).cumsum()
    / np.arange(1, len(pr_res_all_wgs) + 1)
)

pr_res_all_wgs['Rank'] = np.arange(1, len(pr_res_all_wgs) + 1)

# Underexpression subset
under_idx = pr_res_all_wgs['zScore'] < 0
under = pr_res_all_wgs.loc[under_idx].copy()

under_rank = np.arange(1, len(under) + 1)
under['underexpression_rank'] = under_rank

under_vus = (under['VUS'].to_numpy() == True)
x = np.cumsum(under_vus)

under['underexpression_proportions'] = x / under_rank

# Vectorized Clopper–Pearson
alpha = 0.05
under['ci_min'] = beta.ppf(alpha/2, x, under_rank - x + 1)
under['ci_max'] = beta.ppf(1 - alpha/2, x + 1, under_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all_wgs.loc[under_idx, [
    "underexpression_proportions",
    "underexpression_rank",
    "ci_min",
    "ci_max"
]] = under[
    ["underexpression_proportions", "underexpression_rank", "ci_min", "ci_max"]
].to_numpy()


# Overexpression subset
over_idx = pr_res_all_wgs['zScore'] > 0
over = pr_res_all_wgs.loc[over_idx].copy()

over_rank = np.arange(1, len(over) + 1)
over['overexpression_rank'] = over_rank

over_vus = (over['over_VUS'].to_numpy() == True)
x = np.cumsum(over_vus)

over['overexpression_proportions'] = x / over_rank

# Vectorized Clopper–Pearson
alpha = 0.05
over['over_ci_min'] = beta.ppf(alpha/2, x, over_rank - x + 1)
over['over_ci_max'] = beta.ppf(1 - alpha/2, x + 1, over_rank - x)

# Efficient vectorized assignment (NO MERGE)
pr_res_all_wgs.loc[over_idx, [
    "overexpression_proportions",
    "overexpression_rank",
    "over_ci_min",
    "over_ci_max"
]] = over[
    ["overexpression_proportions", "overexpression_rank", "over_ci_min", "over_ci_max"]
].to_numpy()

# %%
pr_res_all_wgs["facet_label"] = f"promoterAI < -0.1 (n={pr_res_all_wgs[pr_res_all_wgs["VUS"] == True].shape[0]})"
p = (
    pn.ggplot(pr_res_all_wgs[(pr_res_all_wgs["zScore"] < 0) & (pr_res_all_wgs["underexpression_rank"] > 100) & (pr_res_all_wgs["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein underexpression outliers\nwith rare WGS variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
pr_res_all_wgs["facet_label"] = f"promoterAI > 0.1 (n={pr_res_all_wgs[pr_res_all_wgs["over_VUS"] == True].shape[0]})"

p2 = (
    pn.ggplot(pr_res_all_wgs[(pr_res_all_wgs["zScore"] > 0) & (pr_res_all_wgs["overexpression_rank"] > 100) & (pr_res_all_wgs["overexpression_rank"] < 1e8)], pn.aes(x="overexpression_rank", 
                   y="overexpression_proportions"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of protein oevrexpression outliers\nwith rare WGS variant",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p2

# %% [markdown]
# ## combined WGS and WES

# %%
### check differences of WGS and WES
combined_res = pd.concat([pr_res_all_wes, pr_res_all_wgs])


# %%
combined_res["facet_label"] = f"Protein disrupting (n={combined_res[combined_res["VUS"] == True].shape[0]})\n (VEP HIGH, missense_variant)"
p = (
    pn.ggplot(combined_res[(combined_res["zScore"] < 0) & (combined_res["underexpression_rank"] > 100) & (combined_res["underexpression_rank"] < 1e8)], pn.aes(x="underexpression_rank", 
                   y="underexpression_proportions", color="seq_type"))
    + pn.geom_line()
    + pn.scale_x_log10()
    + pn.annotation_logticks(sides="b") +
    pn.facet_grid(". ~ facet_label") +
    pn.labs(
        y="Proportion of underexpression outliers\nwith rare variants",
        x="Outlier rank cutoff"
    )
    + pn.theme_bw(base_size=12)
    #+ pn.facet_grid(". ~ plot_title")
)
p

# %%
pr_res_all_wes[(pr_res_all_wes["IMPACT"].notna()) & (pr_res_all_wes["VUS"] == False)]

# %% [markdown]
# # Splicing

# %%
fr_res =  pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_outliers_variants.tsv", sep="\t")


# %%
fr_res = pd.merge(fr_res, dresden_dt_cgc[["gene_name", "gene_type", "geneID_short", "ROLE_IN_CANCER", "predisposition_gene"]], on="geneID_short", how="left")
fr_res = pd.merge(fr_res, sa, left_on="sampleID", right_on="pid")


# %% [markdown]
# # shsred RNA-protein

# %%
# joined.to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/joined.csv")
joined = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/rna_protein_outliers.parquet")

# %%
joined.columns

# %%
joined

# %%
joined.loc[joined["padjust_predisp_extended_x"] <= 0.05, "RNA_aberrant"] = True
joined.loc[joined["padjust_predisp_extended_y"] <= 0.1, "Protein_aberrant"] = True


# %%
joined["Outlier status"] = "Non outlier"
joined.loc[(joined["RNA_aberrant"] == False) & (joined["Protein_aberrant"] == True), "Outlier status"] = "Protein"
joined.loc[(joined["RNA_aberrant"] == True) & (joined["Protein_aberrant"] == True), "Outlier status"] = "RNA & Protein"
joined.loc[(joined["RNA_aberrant"] == True) & (joined["Protein_aberrant"] == False), "Outlier status"] = "RNA"



# %%
joined[(joined["Outlier status"] == "RNA & Protein") & (joined["RNA_zScore"] < 0) & (joined["Protein_zScore"] < 0) & 
    ((joined["ANNOTATION_control_snv_x"] == "germline") | (joined["ANNOTATION_control_indel_x"] == "germline"))]

# %%
joined["VUS"] = "No rare variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == True), "VUS"] = "RNA & protein variant"
joined.loc[(joined["VUS_x"] == False) & (joined["VUS_y"] == True), "VUS"] = "Protein variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == False), "VUS"] = "RNA variant"



# %%

# %%
p = (pn.ggplot(joined) + 
    pn.geom_point(pn.aes(x="RNA_zScore", y="Protein_zScore", color="Outlier status", shape="VUS"), alpha=0.5) +
    pn.scale_size_manual(
        values={
            "No rare variant": 1,
            "RNA variant": 3,
            "Protein variant": 3,
            "RNA & protein variant": 3
        }
    ) + 
     pn.scale_color_manual(
        values={
            "Non outlier": "grey",
            "Protein": "lightblue",
            "RNA": "firebrick",
            "RNA & Protein": "lightgreen"
        }
    ) + 
     pn.theme_bw()+
    pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
     pn.coord_cartesian(ylim=[-25, 22], xlim=[-25, 22])
    
)
p

# %%
joined["VUS"] = "No rare variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == True), "VUS"] = "RNA & protein variant"

predisp_var = joined[joined["padjust_predisp_extended_x"].notna()]
mask_missense = (
    predisp_var["Consequence_snv_x"].str.contains("missense", na=False) |
    predisp_var["Consequence_indel_x"].str.contains("missense", na=False)
)

mask_stop = (
    predisp_var["Consequence_snv_x"].str.contains("stop", na=False) |
    predisp_var["Consequence_indel_x"].str.contains("stop", na=False)
)

mask_frameshift = (
    predisp_var["Consequence_snv_x"].str.contains("frameshift_variant", na=False) |
    predisp_var["Consequence_indel_x"].str.contains("frameshift_variant", na=False)
)

mask_splice = (
    predisp_var["Consequence_snv_x"].str.contains("splice_donor_variant|splice_acceptor_variant", na=False) |
    predisp_var["Consequence_indel_x"].str.contains("splice_donor_variant|splice_acceptor_variant", na=False)
)

predisp_var.loc[mask_missense, "VUS"] = "missense"
predisp_var.loc[mask_stop, "VUS"] = "VEP stop"
predisp_var.loc[mask_frameshift, "VUS"] = "VEP frameshift_variant"
predisp_var.loc[mask_splice, "VUS"] = "VEP splice"

predisp_var["VUS"].value_counts()


# %%
predisp_var[(predisp_var["Outlier status"] == "RNA & Protein")]

# %%

p = (pn.ggplot(predisp_var[(predisp_var["VUS"] != "No rare variant") & 
     ((predisp_var["ANNOTATION_control_snv_x"].str.contains("germline")) | 
      (predisp_var["ANNOTATION_control_indel_x"].str.contains("germline"))
     )]
    ) + 
    pn.geom_point(pn.aes(x="RNA_zScore", y="Protein_zScore", color="Outlier status", shape="VUS", alpha="VUS")) +
    pn.scale_size_manual(
        values={
            "No rare variant": 0.5,
            "RNA variant": 3,
            "Protein variant": 3,
            "RNA & protein variant": 3, 
            "missense": 4,
            "VEP splice": 4,
            "VEP frameshift_variant": 4,
            "VEP frameshift_variant": 4,
        }
    ) + 
     pn.scale_alpha_manual(
        values={
            "No rare variant": 0.4,
            "RNA variant": 0.8,
            "Protein variant": 0.8,
            "RNA & protein variant": 0.8, 
            "missense": 0.7,
            "VEP stop": 0.7,
            "VEP frameshift_variant": 0.7,
            "VEP splice": 0.7,

        }
    ) + 
     pn.scale_color_manual(
        values={
            "Non outlier": "grey",
            "Protein": "blue",
            "RNA": "firebrick",
            "RNA & Protein": "lightgreen"
        }
    ) + 
     pn.labs(title="Predisposition genes with rare germline variant") +
     pn.theme_bw()+
    pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
    pn.coord_cartesian(ylim=[-20, 20], xlim=[-20, 20]) 

    
)
p

# %%
interesting =predisp[(predisp["Protein_aberrant"] == True) & (predisp["RNA_aberrant"] == False)]

# %%
interesting["chrom_snv"] =  interesting["Location_snv"].str.split(":").str[0]
interesting["pos_snv"] = interesting["Location_snv"].str.split(":").str[1]

interesting["chrom_indel"] =  interesting["Location_indel_x"].str.split(":").str[0]
interesting["pos_indel"] = interesting["Location_indel_x"].str.split(":").str[1]

interesting = interesting.rename(columns={"pValue_x": "RNA_pValue", "pValue_y": "protein_pValue", 
                            "padjust_x": "RNA_padjust", "padjust_y": "protein_padjust", 
                            "geneID_x": "geneID", "IMPACT_snv_x": "IMPACT_snv", "Consequence_snv_x" : "Consequence_snv", 
                           "ANNOTATION_control_snv_x": "ANNOTATION_control_snv", "Allele_snv": "alt_snv", "IMPACT_indel_x": "IMPACT_indel",
                           "Consequence_indel_x": "Consequence_indel", "ANNOTATION_control_indel_x": "ANNOTATION_control_indel",
                            "Outlier status_x": "RNA_outlier_status", "Outlier status_y": "protein_outlier_status", 
                            "ROLE_IN_CANCER_x": "ROLE_IN_CANCER"
                            },)

# %%
interesting = interesting.merge(sa, left_on="sampleID", right_on="pid")


# %%
cols = ["sampleID", "geneID", "Oncotree Code", "Oncotree Text", "gene_name", "RNA_pValue", "RNA_padjust", "protein_pValue", "protein_padjust", "chrom_snv", "pos_snv", "ref_snv", "alt_snv", 
        "IMPACT_snv", "Consequence_snv", "ANNOTATION_control_snv", "chrom_indel", "pos_indel", "ref_indel", "alt_indel",
         "IMPACT_indel", "Consequence_indel", "ANNOTATION_control_indel", "RNA_outlier_status", "protein_outlier_status", "ROLE_IN_CANCER", "nct_pid"]

interesting[cols]

# %%

# %%
interesting[cols].to_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/protein_nonRNA_predisposition_outliers.tsv", sep="\t", index=None)

# %%
interesting
