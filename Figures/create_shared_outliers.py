
import pandas as pd
import plotnine as pn
import numpy as np
from scipy.stats import beta
from itables import init_notebook_mode

init_notebook_mode(all_interactive=True)



import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *
from ProteinExpression.load_pr_data import *


sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")

# py_or_res_all = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust_cnv.parquet")

# #py_or_res_all = pd.merge(py_or_res_all, sa, left_on="sampleID", right_on="pid")
# py_or_res_all = pd.merge(py_or_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID", "ROLE_IN_CANCER", "predisposition_gene"]], on="geneID", how="left")




# py_or_res_all.loc[:, "VUS_snv"] = False
# py_or_res_all.loc[py_or_res_all["promoterAI_snv"] <= -0.1, "VUS_snv"] = True

# py_or_res_all.loc[:, "over_VUS_snv"] = False
# py_or_res_all.loc[py_or_res_all["promoterAI_snv"] >= 0.1, "over_VUS_snv"] = True

# print(py_or_res_all[py_or_res_all['promoterAI_snv'].notna()]["VUS_snv"].value_counts())
# py_or_res_all[py_or_res_all['promoterAI_snv'].notna()]["over_VUS_snv"].value_counts()


# py_or_res_all.loc[:, "VUS"] = False
# py_or_res_all.loc[py_or_res_all["IMPACT"] == "HIGH", "VUS"] = True


# print(py_or_res_all["VUS"].value_counts())


# ### load proteomics

# pr_output_name = "cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA"

# pr_res_all = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_" + pr_output_name + "/pr_variants_predisppadjust_cnv.parquet")
# pr_res_all = pd.merge(pr_res_all, dresden_dt_cgc[["gene_name", "gene_type", "geneID_short", "ROLE_IN_CANCER", "predisposition_gene"]], right_on="geneID_short", left_on="geneID", how="left")
# #pr_res_all = pd.merge(pr_res_all, sa, left_on="sampleID", right_on="pid")




# protein_disrupting_variants = ["stop_gained", "missense_variant", "stop_lost", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant"]

# pr_res_all["VUS"] = False

# mask_snv = (
#     pr_res_all["Consequence_snv"]
#     .fillna("")
#     .str.split(",")
#     .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
# )

# mask_indel = (
#     pr_res_all["Consequence_indel"]
#     .fillna("")
#     .str.split(",")
#     .apply(lambda x: bool(set(protein_disrupting_variants).intersection(x)))
# )

# pr_res_all.loc[mask_snv, "VUS"] = True
# pr_res_all.loc[mask_indel, "VUS"] = True

# print(pr_res_all["VUS"].value_counts())

# pr_res_all.loc[pr_res_all["IMPACT"] == "HIGH", "VUS"] = True

# print(pr_res_all["VUS"].value_counts())




# py_or_res_all["gene_sample"] = py_or_res_all["geneID_short"] + "_" + py_or_res_all["sampleID"]
# py_or_res_outleirs = set(py_or_res_all[py_or_res_all["padjust"] <= 0.05]["gene_sample"])



# pr_res_all["gene_sample"] = pr_res_all["geneID"] + "_" + pr_res_all["sampleID"]
# pr_res_all_outliers = set(pr_res_all[pr_res_all["padjust"] <= 0.1]["gene_sample"])

# rna_genes = set(py_or_res_all["geneID_short"])
# proteins = set(pr_res_all["geneID"])
# proteins_rna = proteins.intersection(rna_genes)
# print(len(rna_genes), len(proteins), len(proteins_rna))

# all_outliers = pr_res_all_outliers.union(py_or_res_outleirs)


# # sub_or = py_or_res_all[py_or_res_all["gene_sample"].isin(all_outliers)]
# sub_or = py_or_res_all[py_or_res_all["geneID_short"].isin(proteins_rna)]
# sub_or.loc[sub_or["padjust_predisp"] <= 0.05, "aberrant"] = True
# sub_or = sub_or.rename(columns={"zScore": "RNA_zScore", "aberrant": "RNA_aberrant"})


# # sub_pr = pr_res_all[pr_res_all["gene_sample"].isin(all_outliers)]
# sub_pr = pr_res_all[pr_res_all["geneID"].isin(proteins_rna)]
# sub_pr.loc[sub_pr["padjust_predisp"] <= 0.1, "aberrant"] = True
# sub_pr = sub_pr.rename(columns={"zScore": "Protein_zScore", "aberrant": "Protein_aberrant"})

# joined = sub_or.merge(sub_pr, left_on=["sampleID", "geneID_short"], right_on=["sampleID", "geneID"], how="outer")


# joined = joined.drop(columns=['PROTEIN_EXPECTED_LOG2INT_x', 'PROTEIN_INT_x','PROTEIN_EXPECTED_LOG2INT_y', 'PROTEIN_INT_y', 'geneID_x'])


# joined.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/rna_protein_outliers.parquet")

joined = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/rna_protein_outliers.parquet")


joined.loc[joined["padjust_predisp_extended_x"] <= 0.05, "RNA_aberrant"] = True
joined.loc[joined["padjust_predisp_extended_y"] <= 0.1, "Protein_aberrant"] = True

joined["Outlier status"] = "Non outlier"
joined.loc[(joined["RNA_aberrant"] == False) & (joined["Protein_aberrant"] == True), "Outlier status"] = "Protein"
joined.loc[(joined["RNA_aberrant"] == True) & (joined["Protein_aberrant"] == True), "Outlier status"] = "RNA & Protein"
joined.loc[(joined["RNA_aberrant"] == True) & (joined["Protein_aberrant"] == False), "Outlier status"] = "RNA"

joined[(joined["Outlier status"] == "RNA & Protein") & (joined["RNA_zScore"] < 0) & (joined["Protein_zScore"] < 0) & 
    ((joined["ANNOTATION_control_snv_x"] == "germline") | (joined["ANNOTATION_control_indel_x"] == "germline"))]

joined["VUS"] = "No rare variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == True), "VUS"] = "RNA & protein variant"
joined.loc[(joined["VUS_x"] == False) & (joined["VUS_y"] == True), "VUS"] = "Protein variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == False), "VUS"] = "RNA variant"


joined["VUS"] = "No rare variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == True), "VUS"] = "RNA & protein variant"

mask_missense = (
    joined["Consequence_snv_x"].str.contains("missense", na=False) |
    joined["Consequence_indel_x"].str.contains("missense", na=False)
)

mask_stop = (
    joined["Consequence_snv_x"].str.contains("stop", na=False) |
    joined["Consequence_indel_x"].str.contains("stop", na=False)
)

mask_frameshift = (
    joined["Consequence_snv_x"].str.contains("frameshift_variant", na=False) |
    joined["Consequence_indel_x"].str.contains("frameshift_variant", na=False)
)

mask_splice = (
    joined["Consequence_snv_x"].str.contains("splice_donor_variant|splice_acceptor_variant", na=False) |
    joined["Consequence_indel_x"].str.contains("splice_donor_variant|splice_acceptor_variant", na=False)
)

joined.loc[mask_missense, "VUS"] = "missense"
joined.loc[mask_stop, "VUS"] = "VEP stop"
joined.loc[mask_frameshift, "VUS"] = "VEP frameshift_variant"
joined.loc[mask_splice, "VUS"] = "VEP splice"

joined.loc[(joined["CNV_x"].str.contains("DEL")) | (joined["CNV_y"].str.contains("DEL")), "VUS"] = "CNV Del"
joined.loc[(joined["CNV_x"].str.contains("AMP")) | (joined["CNV_y"].str.contains("AMP")), "VUS"] = "CNV AMP"
joined.loc[(joined["CNV_x"].str.contains("DUP")) | (joined["CNV_y"].str.contains("DUP")), "VUS"] = "CNV DUP"



joined["VUS"].value_counts()

predisp = joined[joined["padjust_predisp_extended_x"].notna()]



p = (pn.ggplot(predisp[(predisp["VUS"] != "No rare variant") & 
     ((predisp["ANNOTATION_control_snv_x"].str.contains("germline")) | 
      (predisp["ANNOTATION_control_indel_x"].str.contains("germline"))
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
            "CNV Del": 4,
            "CNV AMP": 4,
            "CNV DUP": 4

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
            "CNV Del": 0.7,
            "CNV AMP": 0.7,
            "CNV DUP": 0.7

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
     pn.labs(title="Predisposition genes with rare small germline variant") +
     pn.theme_bw()+
    pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
    pn.coord_cartesian(ylim=[-20, 20], xlim=[-20, 20]) 

    
)
p.save("/home/a379i/rna_pr_zsocres_predisp.png", dpi=600, width=8, height=5, units="in")


p = (pn.ggplot(predisp[(predisp["VUS"] != "No rare variant")]
    ) + 
    pn.geom_point(pn.aes(x="RNA_zScore", y="Protein_zScore", color="Outlier status", shape="VUS", alpha="VUS")) +
    pn.scale_size_manual(
        values={
            "RNA variant": 3,
            "Protein variant": 3,
            "RNA & protein variant": 3, 
            "missense": 4,
            "VEP splice": 4,
            "VEP frameshift_variant": 4,
            "CNV Del": 4,
            "CNV AMP": 4,
            "CNV DUP": 4

        }
    ) + 
     pn.scale_alpha_manual(
        values={
            "RNA variant": 0.8,
            "Protein variant": 0.8,
            "RNA & protein variant": 0.8, 
            "missense": 0.7,
            "VEP stop": 0.7,
            "VEP frameshift_variant": 0.7,
            "VEP splice": 0.7,
            "CNV Del": 0.7,
            "CNV AMP": 0.7,
            "CNV DUP": 0.7

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
     pn.labs(title="Predisposition genes with rare small variant") +
     pn.theme_bw()+
    pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
    pn.coord_cartesian(ylim=[-20, 20], xlim=[-20, 20]) 

    
)
p.save("/home/a379i/rna_pr_zsocres_predisp_nongermline.png", dpi=600, width=8, height=5, units="in")


p = (pn.ggplot(predisp[(predisp["VUS"] != "No rare variant") & ((predisp["RNA_aberrant"] == True ) | ((predisp["Protein_aberrant"] == True )))]
    ) + 
    pn.geom_point(pn.aes(x="RNA_zScore", y="Protein_zScore", color="Outlier status", shape="VUS", alpha="VUS")) +
    pn.scale_size_manual(
        values={
            "RNA variant": 3,
            "Protein variant": 3,
            "RNA & protein variant": 3, 
            "missense": 4,
            "VEP splice": 4,
            "VEP frameshift_variant": 4,
            "CNV Del": 4,
            "CNV AMP": 4,
            "CNV DUP": 4

        }
    ) + 
     pn.scale_alpha_manual(
        values={
            "RNA variant": 0.8,
            "Protein variant": 0.8,
            "RNA & protein variant": 0.8, 
            "missense": 0.7,
            "VEP stop": 0.7,
            "VEP frameshift_variant": 0.7,
            "VEP splice": 0.7,
            "CNV Del": 0.7,
            "CNV AMP": 0.7,
            "CNV DUP": 0.7

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
     pn.labs(title="Predisposition genes with rare small variant") +
     pn.theme_bw()+
    pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
    pn.coord_cartesian(ylim=[-20, 20], xlim=[-20, 20]) 

    
)
p.save("/home/a379i/rna_pr_zsocres_predisp_nongermline_nonoutlier.png", dpi=600, width=8, height=5, units="in")



joined["VUS"] = "No rare variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == True), "VUS"] = "RNA & protein variant"
joined.loc[(joined["VUS_x"] == False) & (joined["VUS_y"] == True), "VUS"] = "Protein variant"
joined.loc[(joined["VUS_x"] == True) & (joined["VUS_y"] == False), "VUS"] = "RNA variant"



p1 = (pn.ggplot(joined) + 
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
p1
p1.save("/home/a379i/rna_pr_zsocres_all_genes.png", dpi=600, width=9, height=6, units="in")


p2 = (pn.ggplot(predisp) + 
    pn.geom_point(pn.aes(x="RNA_zScore", y="Protein_zScore", color="Outlier status", shape="VUS"), alpha=0.6) +
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
     pn.labs("Lmiting to only predisposition genes") +
     pn.theme_bw()+
    pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
     pn.coord_cartesian(ylim=[-25, 22], xlim=[-25, 22])
    
)
p2.save("/home/a379i/rna_pr_predisp_zsocres.png", dpi=600, width=9, height=6, units="in")






extended_predisp = joined[(joined["padjust_predisp_extended_x"].notna()) & (joined["padjust_predisp_extended_y"].notna())]
p3 = (pn.ggplot(extended_predisp[(extended_predisp["Consequence_snv_x"].notna()) & (extended_predisp["Consequence_indel_x"].notna())])  + 
    pn.geom_point(pn.aes(x="RNA_zScore", y="Protein_zScore", color="Outlier status", shape="VUS"), alpha=0.6) +
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
     pn.labs(title="Lmiting to only extended predisposition genes with variants") +
     pn.theme_bw()+
    pn.geom_abline(intercept=0, slope=1, linetype="dashed") +
     pn.coord_cartesian(ylim=[-25, 22], xlim=[-25, 22])
    
)
p3.save("/home/a379i/rna_pr_extended_predisp_variants_zsocres.png", dpi=600, width=8, height=5, units="in")






