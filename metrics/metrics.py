import scib
import anndata as ad
from abc import ABC, abstractmethod
from anndata import AnnData
from typing import List, Union
from collections import defaultdict

class Metric(ABC):
    def __init__(self, name: str):
        self.name = name
    
    @abstractmethod
    def __call__(self, adata: AnnData, label_key: str) -> float:
        pass

class Silhouette(Metric):
    def __init__(self):
        super().__init__("Average Silhouette Width")
    
    def __call__(self, adata: AnnData, label_key: str) -> float:
        return scib.me.silhouette(adata, label_key, "X_pca")

class ARI(Metric):
    def __init__(self):
        super().__init__("Adjusted Rand Index")
    
    def __call__(self, adata: AnnData, label_key: str) -> float:
        return scib.me.ari(adata, "cluster", label_key)

def apply_metrics(metrics: List[Metric], anndata: List[AnnData], label_key: Union[str, List[str]]):
    result = defaultdict(list)

    if isinstance(label_key, str):
        label_key = [label_key for i in range(len(anndata))]

    for adata, label in zip(anndata, label_key):
        scib.pp.reduce_data(adata, overwrite_hvg= False, umap = True)
        scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key=label)
        for metric in metrics:
            result[metric.name].append(metric(adata, label))

    return result


