from dataclasses import dataclass
import pandas as pd
import anndata as ad
import episcanpy as epi
from pathlib import Path
import logging

logger = logging.getLogger()

MOUSE_BRAIN_BIN_MATRIX_PATH = Path('data/mouse_brain_5k/bin_by_cell.h5ad')


@dataclass
class EpiDataset():
    peak_file: Path
    fragments_file: Path
    gtf_file: Path
    barcode_file: Path
    annotation_file: Path

    def get_valid_barcodes(self):
        barcode_info = pd.read_csv(self.barcode_file)
        return barcode_info[barcode_info.is__cell_barcode == 1].barcode.tolist()

    def get_annotation_per_barcode(self, annotations=('cell_type', 'protocol')):
        anno = ad.read(self.annotation_file)
        if not isinstance(annotations, (list, tuple)):
            annotations = [annotations]
        return anno.obs[list(annotations)]

    def get_window_mtx(self, window_size=100000, species="human", fast=True,
                       annotations=('cell_type', 'protocol')) -> ad.AnnData:
        valid_barcodes = self.get_valid_barcodes()
        adata = epi.ct.window_mtx(str(self.fragments_file), valid_barcodes, window_size=window_size, species=species,
                                  fast=fast)
        anno = self.get_annotation_per_barcode(annotations)
        adata.obs[anno.columns] = 'No annotation'
        common_barcodes = adata.obs.index.intersection(anno.index)
        adata.obs.loc[common_barcodes, anno.columns] = anno.loc[common_barcodes]
        return adata


def load_mousebrain_dataset() -> EpiDataset:
    name = "5k_brain"
    base_dir = Path('data/mouse_brain_5k/')

    peak_file = base_dir / "{}_peaks.narrowPeak".format(name)
    fragments_file = base_dir / "atac_v1_adult_brain_fresh_5k_fragments.tsv.gz"
    gtf_file = base_dir / "gencode.vM25.basic.annotation.gtf.gz"
    barcode_file = base_dir / "atac_v1_adult_brain_fresh_5k_singlecell.csv"
    annotation_file = base_dir / "10x-ATAC-Brain5k.h5ad"
    return EpiDataset(peak_file=peak_file, fragments_file=fragments_file, gtf_file=gtf_file, barcode_file=barcode_file,
                      annotation_file=annotation_file)


def load_mousebrain_window_mtx(load=True, save=False, annotations=('cell_type', 'protocol')) -> ad.AnnData:
    epidata = load_mousebrain_dataset()
    if not load or not MOUSE_BRAIN_BIN_MATRIX_PATH.exists():
        adata = epidata.get_window_mtx(annotations=annotations)
    else:
        logger.info("Loading bin by cell matrix from " + str(MOUSE_BRAIN_BIN_MATRIX_PATH))
        adata = ad.read_h5ad(MOUSE_BRAIN_BIN_MATRIX_PATH)
        anno = epidata.get_annotation_per_barcode(annotations)
        adata.obs = pd.DataFrame(index=adata.obs.index)
        adata.obs[anno.columns] = 'No annotation'
        common_barcodes = adata.obs.index.intersection(anno.index)
        adata.obs.loc[common_barcodes, anno.columns] = anno.loc[common_barcodes]
        save = False

    if save:
        if MOUSE_BRAIN_BIN_MATRIX_PATH.exists():
            raise FileExistsError(str(MOUSE_BRAIN_BIN_MATRIX_PATH))
        logger.info("Saving bin by cell matrix as " + str(MOUSE_BRAIN_BIN_MATRIX_PATH))
        adata.write_h5ad(MOUSE_BRAIN_BIN_MATRIX_PATH)

    return adata
