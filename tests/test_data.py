import data
from src.data import load_mousebrain_window_mtx

MOUSE_BRAIN_BIN_MTX_SHAPE = (3880, 30877)


def test_get_annotation_per_barcode():
    epidata = data.load_mousebrain_dataset()
    anno = epidata.get_annotation_per_barcode()
    print(anno)
    assert 'cell_type' in set(anno.columns)

def test_mousebrain_window_mtx_new():
    adata = load_mousebrain_window_mtx(load=False)
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE

def test_mousebrain_window_mtx_save():
    adata = load_mousebrain_window_mtx(load=False, save=True)
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE

def test_mousebrain_window_mtx_load():
    adata = load_mousebrain_window_mtx()
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE

def test_mousebrain_window_mtx_load_annot():
    vals = ['cell_type', 'protocol', 'domain']
    adata = load_mousebrain_window_mtx(annotations=vals)
    assert adata.shape == MOUSE_BRAIN_BIN_MTX_SHAPE
    assert list(adata.obs.columns) == vals
