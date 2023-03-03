from __future__ import annotations

import dataclasses
import pickle
from typing import List, Dict
from sklearn import metrics
import numpy as np
from pathlib import Path


@dataclasses.dataclass
class DAResults:
    p_vals: np.array
    scores: np.array
    log_fcs: np.array
    target: np.array
    run: str
    scenario: str
    method: str
    cell_group: str
    feature_group: str = None
    fdr: float = 0.05
    no_da_class: str = 'noDA'
    da_classes: List[str] = None

    def __post_init__(self):
        if self.da_classes is None:
            classes = np.unique(self.target)
        else:
            if not isinstance(self.da_classes, list):
                raise TypeError("da_classes must be of type list or None")
            classes = [self.no_da_class] + self.da_classes
        self.key = f"{self.run}_{self.scenario}_{self.cell_group}_{self.method}_{self.fdr}"
        idxs = np.where(np.isin(self.target, classes))
        self.p_vals = self.p_vals[idxs]
        self.log_fcs = self.log_fcs[idxs]
        self.target = self.target[idxs]
        self.target = ~(self.target == self.no_da_class)
        self.performance = self._performance(self.fdr)
        self.roc = self._roc_curve()

    def _performance(self, fdr=0.05):
        pred = self.p_vals <= fdr
        target = self.target

        tn, fp, fn, tp = metrics.confusion_matrix(target, pred).ravel()
        precision = metrics.precision_score(target, pred)
        accuracy = metrics.accuracy_score(target, pred)
        recall = metrics.recall_score(target, pred)
        f1 = metrics.f1_score(target, pred)
        mcc = metrics.matthews_corrcoef(target, pred)
        return Performance(tp=tp, fp=fp, fn=fn, tn=tn, precision=precision,
                           accuracy=accuracy, recall=recall, f1=f1, mcc=mcc,
                           fdr=fdr, method=self.method, feature_group=self.feature_group,
                           cell_group=self.cell_group, run=self.run, scenario=self.scenario)

    def _roc_curve(self):
        pred = 1 - self.p_vals
        target = self.target

        fpr, tpr, thresholds = metrics.roc_curve(target, pred, pos_label=True)
        auc = metrics.auc(fpr, tpr)
        return Roc(fpr, tpr, auc, method=self.method, feature_group=self.feature_group,
                   cell_group=self.cell_group, run=self.run, scenario=self.scenario)

    def to_pickle(self, path: Path):
        with open(path, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    @classmethod
    def from_pickle(cls, path: Path):
        with open(path, 'rb') as inp:
            res = pickle.load(inp)
            if not isinstance(res, cls):
                raise TypeError("the given file is not of type " + str(cls))
        return res


@dataclasses.dataclass
class Performance:
    tp: float
    fp: float
    tn: float
    fn: float
    accuracy: float
    recall: float
    precision: float
    f1: float
    mcc: float
    fdr: float
    run: str
    scenario: str
    method: str
    cell_group: str
    feature_group: str = None

    def __post_init__(self):
        self.key = f"{self.run}_{self.scenario}_{self.cell_group}_{self.method}_{self.fdr}"


@dataclasses.dataclass
class Roc:
    fpr: np.array
    tpr: np.array
    auc: float
    run: str
    scenario: str
    method: str
    cell_group: str
    feature_group: str = None

    def __post_init__(self):
        self.key = f"{self.run}_{self.scenario}_{self.cell_group}_{self.method}"
