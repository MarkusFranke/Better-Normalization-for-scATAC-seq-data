import pathlib
from dataclasses import dataclass
import pandas as pd
import anndata as ad
import episcanpy as epi
from pathlib import Path
import logging

logger = logging.getLogger()


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
                       annotations=('cell_type', 'protocol'), save: bool = False,
                       matrix_path: Path = None) -> ad.AnnData:
        valid_barcodes = self.get_valid_barcodes()
        adata = epi.ct.window_mtx(str(self.fragments_file), valid_barcodes, window_size=window_size, species=species,
                                  fast=fast)
        self.add_annotation(adata, annotations=annotations)
        if save:
            EpiDataset.__save_mtx(adata, matrix_path)

        return adata

    def get_peak_mtx(self, annotations=('cell_type', 'protocol'), save: bool = False,
                     matrix_path: Path = None, normalized_peak_size=None, fast=False) -> ad.AnnData:
        valid_barcodes = self.get_valid_barcodes()
        adata = epi.ct.peak_mtx(
            str(self.fragments_file),
            str(self.peak_file),
            valid_barcodes,
            normalized_peak_size=normalized_peak_size,
            fast=fast
        )
        self.add_annotation(adata, annotations=annotations)
        if save:
            EpiDataset.__save_mtx(adata, matrix_path)

        return adata

    @staticmethod
    def __save_mtx(adata: ad.AnnData, matrix_path: Path = None):
        if matrix_path is None or not isinstance(matrix_path, pathlib.Path):
            raise TypeError("bin_matrix_path is invalid")
        if matrix_path.exists():
            raise FileExistsError(str(matrix_path))
        logger.info("Saving matrix as " + str(matrix_path))
        adata.write_h5ad(matrix_path)

    def add_annotation(self, adata: ad.AnnData, annotations=('cell_type', 'protocol')):
        anno = self.get_annotation_per_barcode(annotations)
        adata.obs[anno.columns] = 'No annotation'
        common_barcodes = adata.obs.index.intersection(anno.index)
        adata.obs.loc[common_barcodes, anno.columns] = anno.loc[common_barcodes]

    def load_mtx(self, matrix_path: Path,
                 annotations=('cell_type', 'protocol')) -> ad.AnnData:
        if matrix_path is None or not isinstance(matrix_path, pathlib.Path):
            raise TypeError("matrix_path is invalid")

        if matrix_path.exists():
            logger.info("Loading matrix from " + str(matrix_path))
            adata = ad.read_h5ad(matrix_path)
            self.add_annotation(adata, annotations=annotations)
        else:
            raise FileNotFoundError('file %s not found' % matrix_path)
        return adata


class MouseBrainDataset(EpiDataset):

    def __init__(self):
        name = "5k_brain"
        base_dir = Path('data/mouse_brain_5k/')

        peak_file = base_dir / "{}_peaks.narrowPeak".format(name)
        fragments_file = base_dir / "atac_v1_adult_brain_fresh_5k_fragments.tsv.gz"
        gtf_file = base_dir / "gencode.vM25.basic.annotation.gtf.gz"
        barcode_file = base_dir / "atac_v1_adult_brain_fresh_5k_singlecell.csv"
        annotation_file = base_dir / "10x-ATAC-Brain5k.h5ad"
        super().__init__(peak_file=peak_file, fragments_file=fragments_file, gtf_file=gtf_file,
                         barcode_file=barcode_file,
                         annotation_file=annotation_file)
