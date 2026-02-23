import pandas as pd

py_or_res_aberrant =  pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust.parquet",)
print(py_or_res_aberrant.columns)
py_or_res_aberrant = py_or_res_aberrant[(py_or_res_aberrant["padjust"] <= 0.05) | (py_or_res_aberrant["padjust_predisp"] <= 0.05) | (py_or_res_aberrant["padjust_predisp_extended"] <= 0.05)]


py_or_res_aberrant.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_outliers.parquet", index=None)