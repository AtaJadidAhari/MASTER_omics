import pandas as pd

pr_or_res_aberrant =  pd.read_csv("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA/pr_variants.csv",)
pr_or_res_aberrant = pr_or_res_aberrant[(pr_or_res_aberrant["padjust"] <= 0.05) | (pr_or_res_aberrant["padjust_predisp"] <= 0.05) | (pr_or_res_aberrant["padjust_predisp_extended"] <= 0.05)]

print(pr_or_res_aberrant.info(verbose=True))

pr_or_res_aberrant.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/protrider_runs/output_cov_gaussian_gs_lr_0_001_epoc2000_noInitPCA/pr_varinats_outliers.parquet", index=None)