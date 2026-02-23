import pandas as pd
import plotnine as pn
from pathlib import Path
import numpy as np
import os

base_path = "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider/"
cohorts = os.listdir(base_path)
for cohort in cohorts:
    if not cohort.endswith(".tsv"):
        cohort_res =  pd.read_csv(base_path + cohort + "OUTRIDER_results_all.tsv", sep="\t")
        cohort_expressed_genes = np.unique(coad_res["geneID"])
        cohort_expressed_genes = pd.DataFrame({"geneID": coad_expressed_genes})
        cohort_expressed_genes["geneID_short"] = cohort_expressed_genes["geneID"].str.split(".").str[0]
        cohort_expressed_genes.to_csv(base_path + cohort + "expressed_genes.tsv", sep="\t")
        
        
        