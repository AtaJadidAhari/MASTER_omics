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
            area = np.trapezoid(y, x)
            
            # 2. Calculate Mean Proportion (Average Height)
            # This normalizes the area by the width of the plot
            mean_prop = area / x_range if x_range > 0 else 0
            
            results[variant_type] = {
                "raw_auc": area,
                "mean_proportion": mean_prop
            }
            
    return results


# %%
def calculate_proportions(df, expression_col="CNV", expression_directions = ["overexpression", "underexpression"]):
    # Work on a copy to prevent SettingWithCopyWarning and index fragmentation
    df = df.copy()
    
    for expression_direction in expression_directions:
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
        elif expression_direction == "overexpression":
            mask = df['zScore'] > 0
        else: # spcling outlier
            mask = df["hgncSymbol"] == df["hgncSymbol"]
            
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

somatic_snvs = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/somatic_snvs_vep_annotated.parquet")

sa = pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/master_drop_sample_annotation_sizeFactorFiltered_0.1.tsv", sep="\t")



fr_res =  pl.scan_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/aggregated_results_gene_all_variants.tsv", separator="\t").filter(pl.col("geneID_short").is_in(genes_of_interest["geneID_short"])).collect(engine="streaming").to_pandas()

fr_res = pd.merge(fr_res, sa, left_on="sampleID", right_on="pid")
fr_res["Method"] = "FRASER2"

print(len(fr_res), "total_number of outliers")

# fr_res = fr_res[fr_res["geneID_short"].isin(genes_of_interest["geneID_short"])]

print(len(fr_res), "total_number of outliers in genesa of interest")

fr_res = fr_res.merge(somatic_snvs[["somatic_snv_#Uploaded_variation", "sampleID", "MASTER_annotated_gene", "vep_Gene", "somatic_snv_IMPACT", "somatic_snv_Consequence", "somatic_snv_max_spliceai_score", "somatic_snv_am_pathogenicity", "somatic_snv_AbSplice2_max", "somatic_snv_promoterAI", "somatic_snv_abexp_v1.1"]] , left_on=["geneID_short", "sampleID"], right_on=["vep_Gene", "sampleID"], how="left")



gene_annot_pl = pl.from_pandas(gene_annot_dt[["gene_name", "geneID_short"]])

germline_snvs = (
    pl.scan_parquet(
        "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA.parquet"
    )
    # Apply all filters first to reduce data size
    .filter(
        (pl.col("sampleID").is_in(fr_res["sampleID"]))
    )
    # Perform the join within Polars using join_as_coerced
    .join(
        gene_annot_pl.lazy(), left_on="HUGO_Symbol", right_on="gene_name", how="inner"
    )
    .filter(pl.col("geneID_short").is_in(genes_of_interest["geneID_short"]))
    # Collect using the streaming engine and convert to Pandas
    .collect(engine="streaming")
    .to_pandas()
)


germline_snvs = germline_snvs.rename(columns={"#CHROM": "seqnames", "REF": "Ref", "ALT": "Alt"})
print(germline_snvs.shape)

germline_snvs = germline_snvs[(germline_snvs["Ref"].str.len() == 1) & (germline_snvs["Alt"].str.len() == 1)]
germline_snvs["chrom"] = "chr" + germline_snvs["seqnames"]
germline_snvs.shape

valid_positions = (
    germline_snvs["POS"]
    .dropna()
    .unique()
    .astype(int)
    .tolist()
)


germline_vep_res = pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/vcf/cnv/germline_predisp_cgc_all_TINDA_vep.parquet",).collect(engine="streaming").to_pandas()


germline_vep_res["chrom"] = germline_vep_res["Location"].str.split(":").str[0]
germline_vep_res["POS"] = 0
germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] == "SNV", "POS"] = germline_vep_res.loc[germline_vep_res["VARIANT_CLASS"] == "SNV", "Location"].str.split(":").str[1].astype(int)
# germline_vep_res["POS"] = germline_vep_res["POS"] + 1
germline_snvs = germline_snvs.merge(germline_vep_res[["chrom", "POS", "Allele", "REF_ALLELE", "IMPACT", "Consequence", "Gene", "#Uploaded_variation", "am_pathogenicity", "am_class", "LoF", "CADD_PHRED", "CADD_RAW","existing_InFrame_oORFs",  "existing_OutOfFrame_oORFs","existing_uORFs", "five_prime_UTR_variant_annotation", "five_prime_UTR_variant_consequence", "max_spliceai_score"]], left_on=["seqnames",  "POS", "Ref", "Alt", "geneID_short"], right_on=["chrom", "POS", "REF_ALLELE", "Allele", "Gene"], how="left") 


absplice_predisp = (pl.scan_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/predisp_max_abSplice_snvs_hg19.parquet")
        .filter((pl.col("hg19_end").is_in(valid_positions)))
        .collect(engine="streaming")
      ).to_pandas()
absplice_predisp.shape



merged_vars = germline_snvs.merge(absplice_predisp, left_on=["chrom_x", "POS", "Ref", "Alt"], right_on=["chrom", "hg19_end", "ref", "alt"], how="left")
# merged_vars.loc[merged_vars["AbSplice2_max"].isna(), "AbSplice2_max"] = 0






merged_vars["germline_splicing"] = False
merged_vars.loc[(merged_vars["AbSplice2_max"] >= 0.1) | (merged_vars["IMPACT"].str.contains("spli")) | (merged_vars["Consequence"].isin(["splice_acceptor_variant", "splice_donor_variant"])) |     (merged_vars["max_spliceai_score"] >= 0.1) | (merged_vars["pangolin_score"] >= 0.1)  , "germline_splicing"] = True
print(sum(merged_vars["germline_splicing"]), "germline_splicing merged_vars")


fr_res_annotated = fr_res.merge(merged_vars, left_on=["sampleID", "geneID_short"], right_on=["sampleID", "geneID_short"], how="left")


fr_res_annotated["germline_splicing"] = False
fr_res_annotated.loc[(fr_res_annotated["AbSplice2_max"] >= 0.1) | (fr_res_annotated["Consequence"].str.contains("spli")) |(fr_res_annotated["Consequence"].isin(["splice_acceptor_variant", "splice_donor_variant"])) | (fr_res_annotated["max_spliceai_score"] >= 0.1) | (fr_res_annotated["pangolin_score"] >= 0.1)  , "germline_splicing"] = True


fr_res_annotated["somatic_splicing"] = False
fr_res_annotated.loc[(fr_res_annotated["somatic_snv_AbSplice2_max"] >= 0.05) | (fr_res_annotated["somatic_snv_Consequence"].str.contains("spli")) |(fr_res_annotated["somatic_snv_Consequence"].isin(["splice_acceptor_variant", "splice_donor_variant"])) | (fr_res_annotated["somatic_snv_max_spliceai_score"] >= 0.1) , "somatic_splicing"] = True



fr_res_annotated["IMPACT"] = "False"
fr_res_annotated.loc[(fr_res_annotated["germline_splicing"] == True) | (fr_res_annotated["somatic_splicing"] == True), "IMPACT"] = "HIGH"


print(sum(fr_res_annotated["germline_splicing"]))
print(sum(fr_res_annotated["somatic_splicing"]))

print(len(fr_res_annotated[fr_res_annotated["IMPACT"] == "HIGH"]), "fr_res impact")

fr_res_annotated["CNV"] = "No CNV"

fr_res_annotated = fr_res_annotated.sort_values("pValue")

# %%

# %%

# %%
# pr_res_all = calculate_proportions(pr_res_all, expression_col="CNV")
fr_res_annotated = calculate_proportions(fr_res_annotated, expression_col="VEP", expression_directions = ["splicing"])


# %%


# protrider_auc_underexpression = get_auc_metrics(pr_res_all)

protrider_auc_overexpression = get_auc_metrics(fr_res_annotated, expression_direction="splicing")

# %%
print(fr_res_annotated.shape, "len protrider after adding proportions")

# %%
# p = plot_variant_proportions(
#     fr_res_annotated, protrider_auc_overexpression, plot_title="FRASER", expression_direction="splicing")
# p.save("/home/a379i/Scripts/AberrantSplicing/proportions_all_genes_underexpression_outliers.png", width=6, height=4, dpi=300)

# %%


plot_data = fr_res_annotated.copy()
expression_direction = "splicing"
p = ( pn.ggplot(plot_data[(plot_data["Method"] == "FRASER2") & (plot_data[f"{expression_direction}_rank"] > 100) & 
                                    (plot_data[f"{expression_direction}_rank"] < 1e10)])
                
                # --- Line 1 (Deletion) ---
                # We map 'color' to a string label. Plotnine sees this as a category.
                + pn.geom_line(pn.aes(x=f"{expression_direction}_rank", 
                                    y=f"{expression_direction}_proportions_VEP", 
                                    color="Method")) 
                + pn.geom_ribbon(pn.aes(x=f"{expression_direction}_rank", 
                                        ymin=f"{expression_direction}_ci_min_VEP", 
                                        ymax=f"{expression_direction}_ci_max_VEP", 
                                        fill="Method"), 
                                alpha=0.2, outline_type='none')
                
            
                + pn.scale_x_log10()
                + pn.annotation_logticks(sides="b")
                + pn.labs(
                    y=f"Proportion of {expression_direction} outliers with \nrare splice-disrupting variants",
                    x="Outlier rank cutoff",
                    #title="Across all genes"
                )
                + pn.theme_bw(base_size=12)
    )





p.save("/home/a379i/Scripts/AberrantSplicing/proportions_all_genes_underexpression_outliers.png", width=6, height=4, dpi=300)

# %%








# %%
# predispostions
# py_or_res_predisp = pr_res_all[pr_res_all["padjust_predisp_extended"].notna()]
# py_or_res_predisp = py_or_res_predisp.sort_values("padjust_predisp_extended")



# # %%
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV")
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "VEP")
# py_or_res_predisp = calculate_proportions(py_or_res_predisp, "CNV_VEP")



# protrider_auc_underexpression_predisp = get_auc_metrics(py_or_res_predisp)




# # %%
# p4 = plot_variant_proportions(
#     py_or_res_predisp, protrider_auc_underexpression_predisp, plot_title="OUTRIDER (only predispostion genes)", expression_direction="underexpression")
# p4.save("/home/a379i/Scripts/ProteinExpression/proportions_predisposition_genes_underexpression_outliers.png", width=6, height=4, dpi=300)


# # %%
# p5 = plot_variant_proportions(
#     py_or_res_predisp, protrider_auc_overexpression_predisp, plot_title="OUTRIDER (only predisposition genes)", expression_direction="overexpression")
# p5.save("/home/a379i/Scripts/ProteinExpression/proportions_predisposition_genes_overexpression_outliers.png", width=6, height=4, dpi=300)

# # %%


# # %% [markdown]
# # zScores vs OUTRIDER: both CNV and VEP_HIGH

# # %%
# pr_res_all = calculate_proportions(pr_res_all, "CNV_VEP")

# # %%


