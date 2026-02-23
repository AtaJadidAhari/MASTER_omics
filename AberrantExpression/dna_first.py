import pandas as pd
from itables import init_notebook_mode
import plotnine as pn
init_notebook_mode(all_interactive=True)
import numpy as np


import sys
sys.path.append("/home/a379i/Scripts")   # path to folder containing the python file

from utils.load_gtf_cgc_dresden import *


def _false_discovery_control(ps, *, axis=0, method='bh'):
    # Input Validation and Special Cases
    ps = np.asarray(ps)

    ps_in_range = (np.issubdtype(ps.dtype, np.number)
                   and np.all(ps == np.clip(ps, 0, 1)))
    if not ps_in_range:
        raise ValueError("`ps` must include only numbers between 0 and 1.")

    methods = {'bh', 'by'}
    if method.lower() not in methods:
        raise ValueError(f"Unrecognized `method` '{method}'."
                         f"Method must be one of {methods}.")
    method = method.lower()

    if axis is None:
        axis = 0
        ps = ps.ravel()

    axis = np.asarray(axis)[()]
    if not np.issubdtype(axis.dtype, np.integer) or axis.size != 1:
        raise ValueError("`axis` must be an integer or `None`")

    if ps.size <= 1 or ps.shape[axis] <= 1:
        return ps[()]

    ps = np.moveaxis(ps, axis, -1)
    m = ps.shape[-1]

    # Main Algorithm
    # Equivalent to the ideas of [1] and [2], except that this adjusts the
    # p-values as described in [3]. The results are similar to those produced
    # by R's p.adjust.

    # "Let [ps] be the ordered observed p-values..."
    order = np.argsort(ps, axis=-1)
    ps = np.take_along_axis(ps, order, axis=-1)  # this copies ps

    # Equation 1 of [1] rearranged to reject when p is less than specified q
    i = np.arange(1, m + 1)
    ps *= m / i

    # Theorem 1.3 of [2]
    if method == 'by':
        ps *= np.sum(1 / i)

    # accounts for rejecting all null hypotheses i for i < k, where k is
    # defined in Eq. 1 of either [1] or [2]. See [3]. Starting with the index j
    # of the second to last element, we replace element j with element j+1 if
    # the latter is smaller.
    np.minimum.accumulate(ps[..., ::-1], out=ps[..., ::-1], axis=-1)

    # Restore original order of axes and data
    np.put_along_axis(ps, order, values=ps.copy(), axis=-1)
    ps = np.moveaxis(ps, -1, axis)

    return np.clip(ps, 0, 1)


or_res = pd.read_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants.parquet")

or_res_predispostion = or_res[or_res["geneID"].isin(dresden_dt["geneID"])]

or_res_predispostion['padjust_predisp'] = _false_discovery_control(or_res_predispostion['pValue'].values, method='bh')

or_res = or_res.merge(or_res_predispostion[["sampleID", "geneID", "padjust_predisp"]], how="left", on=["sampleID", "geneID"])



or_res_predispostion_extended = or_res[or_res["geneID_short"].isin(extended_dresden_dt["geneID_short"])]
or_res_predispostion_extended['padjust_predisp_extended'] = _false_discovery_control(or_res_predispostion_extended['pValue'].values, method='bh')
or_res = or_res.merge(or_res_predispostion_extended[["sampleID", "geneID", "padjust_predisp_extended"]], how="left", on=["sampleID", "geneID"])



or_res.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_predisppadjust.parquet", index=None)



py_or_res_aberrant = or_res[(or_res["padjust"] <= 0.05) | (or_res["padjust_predisp"] <= 0.05) | (or_res["padjust_predisp_extended"] <= 0.05)]


py_or_res_aberrant.to_parquet("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/py_outrider_runs/all_cohorts/oht_cov_diag_lr_0_0001_epoc200_gpu/or_variants_outliers.parquet", index=None)

