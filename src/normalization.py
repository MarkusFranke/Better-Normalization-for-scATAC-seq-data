import episcanpy as epi
import scipy as sc
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import logging
from scipy import sparse
import numpy as np

logger = logging.getLogger()
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()


def log1pPF(adata):
    # proportional fitting to mean of cell depth
    proportional_fitting = epi.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1pPF_normalization"] = epi.pp.log1p(adata, copy=True).X

    counts = np.ravel(adata.X.sum(1))
    counts_greater_than_zero = counts[counts > 0]
    median_counts = np.median(counts_greater_than_zero, axis=0)
    counts += counts == 0
    counts = counts / median_counts
    adata.obs['pf_size_factors'] = counts

    print("log1pPF_normalization layer added")
    print("pf_size_factors layer added")


def scran(adata):
    # Preliminary clustering for differentiated normalisation
    adata_pp = adata.copy()
    epi.pp.normalize_total(adata_pp)
    epi.pp.log1p(adata_pp)
    epi.pp.pca(adata_pp, n_comps=15)
    epi.pp.neighbors(adata_pp)
    epi.tl.leiden(adata_pp, key_added="groups")

    data_mat = adata.X.T
    # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
    if sc.sparse.issparse(data_mat):
        if data_mat.nnz > 2 ** 31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
    ro.globalenv["data_mat"] = data_mat
    ro.globalenv["input_groups"] = adata_pp.obs["groups"]
    del adata_pp

    ro.r("source('r/normalization.R')")
    size_factors = ro.r("normScran(%s, %s)" % ('data_mat', 'input_groups'))
    del ro.globalenv["data_mat"]
    del ro.globalenv["input_groups"]

    adata.obs["scran_size_factors"] = size_factors
    scran = adata.X / adata.obs["scran_size_factors"].values[:, None]
    adata.layers["scran_normalization"] = sparse.csr_matrix(epi.pp.log1p(scran))
    print("scran_normalization layer added")
    print("scran_size_factors obs added")
