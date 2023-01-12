import data
from pathlib import Path

MOUSE_BRAIN_BIN_MTX_SHAPE = (3880, 30877)


def test_get_annotation_per_barcode():
    epidata = data.MouseBrainDataset()
    anno = epidata.get_annotation_per_barcode()
    print(anno)
    assert 'cell_type' in set(anno.columns)

def test_mousebrain_window_mtx_new():
    epidata = data.MouseBrainDataset()
    adata = epidata.get_window_mtx()
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE

def test_mousebrain_window_mtx_save():
    epidata = data.MouseBrainDataset()
    adata = epidata.get_window_mtx(save=True, bin_matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE

def test_mousebrain_window_mtx_load():
    epidata = data.MouseBrainDataset()
    adata = epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE

def test_mousebrain_window_mtx_load_annot():
    vals = ['cell_type', 'protocol', 'domain']
    epidata = data.MouseBrainDataset()
    adata = epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'), annotations=vals)
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE
    assert list(adata.obs.columns) == vals
