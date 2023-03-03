import episcanpy as epi
import scipy as scipy
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import logging
from scipy import sparse
import numpy as np
from importlib.resources import files
from scipy.sparse import csr_matrix
from multiprocessing.pool import ThreadPool
from sklearn.linear_model import LogisticRegression
from tqdm import tqdm
import scanpy as sc
import anndata as ad

RNORM_SCRIPT = files('epinorm.r').joinpath('normalization.R')

logger = logging.getLogger()
rcb.logger.setLevel(logging.ERROR)


def log1pPF(adata, layer='raw'):
    # proportional fitting to mean of cell depth
    proportional_fitting = sc.pp.normalize_total(adata, layer=layer, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["%s_log1pPF_normalization" % layer] = epi.pp.log1p(proportional_fitting['X'])
    adata.layers["%s_PFlog1pPF_normalization" % layer] = sc.pp.normalize_total(
        adata, target_sum=None, layer="%s_log1pPF_normalization" % layer, inplace=False
    )["X"]

    counts = np.ravel(adata.X.sum(1))
    counts_greater_than_zero = counts[counts > 0]
    median_counts = np.median(counts_greater_than_zero, axis=0)
    counts += counts == 0
    counts = counts / median_counts
    # adata.obs['pf_size_factors'] = counts

    print("%s_PFlog1pPF_normalization layer added" % layer)
    print("%s_log1pPF_normalization layer added" % layer)
    # print("pf_size_factors layer added")


def scran(adata, layer='raw'):
    try:
        ro.pandas2ri.activate()
        anndata2ri.activate()
        # Preliminary clustering for differentiated normalisation
        adata_pp = adata.copy()
        adata_pp.X = adata_pp.layers[layer]
        epi.pp.normalize_total(adata_pp)
        epi.pp.log1p(adata_pp)
        epi.pp.pca(adata_pp, n_comps=15)
        epi.pp.neighbors(adata_pp)
        epi.tl.leiden(adata_pp, key_added="groups")

        data_mat = adata.layers[layer].T
        # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
        if scipy.sparse.issparse(data_mat):
            if data_mat.nnz > 2 ** 31 - 1:
                data_mat = data_mat.tocoo()
            else:
                data_mat = data_mat.tocsc()
        ro.globalenv["data_mat"] = data_mat
        ro.globalenv["input_groups"] = adata_pp.obs["groups"]
        del adata_pp

        ro.r("source('%s')" % RNORM_SCRIPT)
        size_factors = ro.r("normScran(%s, %s)" % ('data_mat', 'input_groups'))
        del ro.globalenv["data_mat"]
        del ro.globalenv["input_groups"]

        # adata.obs["scran_size_factors"] = size_factors
        scran = adata.layers[layer] / size_factors[:, None]
        adata.layers["%s_scran_normalization" % layer] = sparse.csr_matrix(epi.pp.log1p(scran))
        print("%s_scran_normalization layer added" % layer)
        # print("scran_size_factors obs added")
    finally:
        ro.pandas2ri.deactivate()
        anndata2ri.deactivate()


def sctransform(adata):
    try:
        ro.pandas2ri.activate()
        anndata2ri.activate()

        # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
        if scipy.sparse.issparse(adata.X):
            if not adata.X.has_sorted_indices:
                adata.X.sort_indices()

        adata_cp = adata.copy()
        del adata_cp.uns
        ro.globalenv["adata"] = adata_cp

        ro.r("source('%s')" % RNORM_SCRIPT)
        ro.globalenv["res"] = ro.r("normSctransform(adata)")
        norm_x = ro.r("res@assays$SCT@scale.data").T
        bins = ro.r("res@assays$SCT@meta.features").index.values
        bins = np.char.replace(bins.astype('str'), '-', '_')
        adata.var["idx"] = list(range(len(adata.var)))
        idxs = adata.var['idx'][bins]
        del adata.var['idx']
        tmp = csr_matrix(adata.shape)
        tmp[:, idxs] = norm_x
        adata.layers["scTransform_normalization"] = tmp
        del ro.globalenv["adata"]
        del ro.globalenv["res"]

        print("scTransform_normalization layer added")
    finally:
        ro.pandas2ri.deactivate()
        anndata2ri.deactivate()


def binary_residuals(adata):
    if 'binarized' not in adata.layers.keys():
        raise ValueError('the adata object should contain a binarized layer')
    current_layer = adata.X
    adata.X = adata.layers['binarized']

    if sparse.issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()

    def process_column(i):
        return np.asarray(np.squeeze(adata.X[:, i].toarray()))

    with ThreadPool() as p:
        y = list(tqdm(p.imap(process_column, range(adata.X.shape[1])), total=adata.X.shape[1]))
    y = np.array(y)

    X = np.asarray(np.sum(adata.X, axis=1)).astype(float)

    def process_column_2(i):
        if not np.any(y[i]):
            return np.zeros(len(y[i]))
        model = LogisticRegression(solver='liblinear', warm_start=True).fit(X=X, y=y[i])
        return y[i] - model.predict_proba(X)[:, 1]

    with ThreadPool() as p:
        X_res = np.array(list(tqdm(p.imap(process_column_2, range(adata.X.shape[1])), total=adata.X.shape[1]))).T
    adata.X = current_layer
    adata.layers["binary_residual_normalization"] = X_res
    print("binary_residual_normalization layer added")


def normalize(adata: ad.AnnData):
    adata.layers['raw'] = adata.X
    adata.layers['binarized'] = epi.pp.binarize(adata, copy=True).X
    for layer in ['raw', 'binarized']:
        log1pPF(adata, layer=layer)
        scran(adata, layer=layer)
    sctransform(adata)
    # norm.binary_residuals(adata_sim)
