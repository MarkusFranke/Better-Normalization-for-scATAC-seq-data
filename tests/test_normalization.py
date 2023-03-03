import pytest
import epinorm.data as data
import epinorm.normalization as norm
from pathlib import Path
import anndata as ad
import episcanpy as epi

@pytest.fixture
def mouse_brain_adata() -> ad.AnnData:
    epidata = data.MouseBrainDataset()
    return epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))


def test_log1pPF(mouse_brain_adata):
    mouse_brain_adata.layers['raw'] = mouse_brain_adata.X
    norm.log1pPF(mouse_brain_adata)
    assert "raw_log1pPF_normalization" in mouse_brain_adata.layers.keys()


def test_scran(mouse_brain_adata):
    mouse_brain_adata.layers['raw'] = mouse_brain_adata.X
    norm.scran(mouse_brain_adata)
    assert "raw_scran_normalization" in mouse_brain_adata.layers.keys()


def test_sctransform(mouse_brain_adata):
    norm.sctransform(mouse_brain_adata)
    assert "scTransform_normalization" in mouse_brain_adata.layers.keys()

def test_binary_residuals(mouse_brain_adata):
    mouse_brain_adata.layers['binarized'] = epi.pp.binarize(mouse_brain_adata, copy=True).X
    norm.binary_residuals(mouse_brain_adata)
    assert "binary_residual_normalization" in mouse_brain_adata.layers.keys()