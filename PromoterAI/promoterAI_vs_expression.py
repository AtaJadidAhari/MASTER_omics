import pandas as pd
from plotnine import *
from itertools import combinations
from scipy.stats import mannwhitneyu


variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38_promoterAI.tsv", sep="\t")
variants = variants[variants['promoterAI'].notna()]


py_or_res_all = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/protrider_summary.csv")
py_or_res_all["geneID_short"] = py_or_res_all["geneID"].str.split(".").str[0]


py_or_res_all = pd.merge(py_or_res_all, variants[["Gene", "sampleID", "IMPACT", "Consequence", "ANNOTATION_control", "promoterAI"]], left_on=["geneID_short", "sampleID"], right_on= ["Gene", "sampleID"], how="inner")


py_or_res_all['Outlier status'] = "Non-outlier"
py_or_res_all.loc[(py_or_res_all['zScore'] < 0) & (py_or_res_all['padjust'] <= 0.05), 'Outlier status'] = "Underexpression"
py_or_res_all.loc[(py_or_res_all['zScore'] > 0) & (py_or_res_all['padjust'] <= 0.05), 'Outlier status'] = "Overexpression"



order = ["Non-outlier", "Overexpression", "Underexpression"]

groups = {
    k: py_or_res_all.loc[py_or_res_all["Outlier status"] == k, "promoterAI"].values
    for k in order
}

results = []

for g1, g2 in combinations(order, 2):
    stat, p = mannwhitneyu(groups[g1], groups[g2], alternative="two-sided")
    results.append({"g1": g1, "g2": g2, "p": p})

df = pd.DataFrame(results)

from statsmodels.stats.multitest import multipletests

df = pd.DataFrame(results)
df["p_adj"] = multipletests(df["p"], method="fdr_bh")[1]


def p_to_star(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

df["label"] = df["p_adj"].apply(p_to_star)
df


x_pos = {k: i + 1 for i, k in enumerate(order)}

y_min = py_or_res_all["promoterAI"].min()
y_max = py_or_res_all["promoterAI"].max()
step = 0.08 * (y_max - y_min)

df["y"] = y_max + (df.index + 1) * step


p = (
    ggplot(py_or_res_all, aes("Outlier status", "promoterAI"))
    + geom_boxplot(outlier_alpha=0)
    + theme_bw()
)

# Add brackets + labels
for _, r in df.iterrows():
    x1 = x_pos[r["g1"]]
    x2 = x_pos[r["g2"]]
    y  = r["y"]

    p += geom_segment(
        aes(x=x1, xend=x2, y=y, yend=y),
        inherit_aes=False
    )
    p += geom_segment(
        aes(x=x1, xend=x1, y=y, yend=y - step/4),
        inherit_aes=False
    )
    p += geom_segment(
        aes(x=x2, xend=x2, y=y, yend=y - step/4),
        inherit_aes=False
    )
    p += annotate(
        "text",
        x=(x1 + x2) / 2,
        y=y + step/6,
        label=r["label"],
        size=10
    )

p.save("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/promoter_ai/pAI_scores_per_outlier_status.png", dpi=600, width=8, height=5, units="in")
