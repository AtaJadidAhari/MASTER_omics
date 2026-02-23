import pandas as pd
import re

def load_pr_data(pr_output):
    # Load PR result file
    pr_res = pd.read_csv(pr_output, sep=",")   # equivalent to fread()
    # If geneID_short is missing
    if "geneID_short" not in pr_res.columns:
        
        # Check whether proteinID is ENSG format
        if "proteinID" not in pr_res.columns:
            pr_res = pr_res.rename(columns={"geneID": "proteinID"})
        first_protein = str(pr_res["proteinID"].iloc[0])
        
        if not first_protein.startswith("ENSG"):
            # Load gene annotation file
            gene_annot_file = (
                "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/"
                "drop_master_202502_allGenes/processed_data/preprocess/v19/"
                "gene_name_mapping_v19.tsv"
            )
            gene_annot_dt = pd.read_csv(gene_annot_file, sep=",")

            # Create short gene ID
            gene_annot_dt["geneID_short"] = gene_annot_dt["gene_id"].str.replace(
                r"\..*", "", regex=True
            )

            # Left join: proteinID → gene_name
            pr_res = pr_res.merge(
                gene_annot_dt,
                how="left",
                left_on="proteinID",
                right_on="gene_name"
            )
        else:
            # If proteinID is ENSG already, just copy
            pr_res["geneID_short"] = pr_res["proteinID"]

    # Rename geneID := geneID_short
    pr_res["geneID"] = pr_res["geneID_short"]

    # Drop empty geneID values
    pr_res = pr_res[pr_res["geneID"].astype(str) != ""]

    if "pValue" not in pr_res.columns:
        pr_res = pr_res.rename(columns={"PROTEIN_PVALUE": "pValue", "PROTEIN_PADJ": "padjust", "PROTEIN_ZSCORE": "zScore", "PROTEIN_outlier": "aberrant"})
    pr_res = pr_res.sort_values("pValue")
    
    return pr_res