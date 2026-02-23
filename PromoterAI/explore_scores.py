import pandas as pd
import plotnine as pn


variants = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/vep_res_aggregated_hg38/vep_res_rare_snv_all_aggregated_unique_variant_type_hg38_promoterAI.tsv", sep="\t")

p1 = (pn.ggplot(variants) +
      pn.geom_boxplot(pn.aes(x="IMPACT", y="promoterAI")) +
      pn.theme_bw()
     )
p1.save("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/promoter_ai/pAI_scores_per_impact.png", dpi=600, width=8, height=5, units="in")


p2 = (pn.ggplot(variants) +
    pn.geom_boxplot(pn.aes(x="am_class", y="promoterAI")) +
    pn.theme_bw()
    )
p2.save("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/promoter_ai/pAI_scores_per_am_class.png", dpi=600, width=8, height=5, units="in")


p3 = (pn.ggplot(variants) +
     pn.geom_boxplot(pn.aes(x="Consequence_most_severe", y="promoterAI")) +
     pn.theme_bw() +
     pn.theme(axis_text_x=pn.element_text(rotation=90))
    )
p3.save("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/results/res_202511/promoter_ai/pAI_scores_per_consequence.png", dpi=600, width=8, height=5, units="in")