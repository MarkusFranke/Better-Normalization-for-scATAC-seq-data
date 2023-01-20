from __future__ import annotations

import dataclasses
from typing import List, Dict

import anndata as ad
import anndata2ri
import pandas as pd
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import logging
from scipy import sparse
import numpy as np
import copy

logger = logging.getLogger()
rcb.logger.setLevel(logging.ERROR)


class Simulator:
    """
    Wrapper for the R simATAC simulator
    """

    def __init__(self):
        ro.r("library(simATAC)")
        self.__r_sim_obj_name = 'simObj'
        self.__r_original_name = 'origCSE'
        ro.globalenv[self.__r_sim_obj_name] = ro.r("newsimATACCount()")

    def __enter__(self):
        ro.pandas2ri.activate()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        logger.info("deleting R objects")
        if self.__r_sim_obj_name in ro.globalenv:
            del ro.globalenv[self.__r_sim_obj_name]
        if self.__r_original_name in ro.globalenv:
            del ro.globalenv[self.__r_original_name]
        ro.pandas2ri.deactivate()

    @property
    def parameters(self) -> Params:
        """
        returns the simulation parameters
        :return:
        """
        sim_obj = ro.globalenv[self.__r_sim_obj_name]
        slots = sim_obj.slots
        return Params(
            nbins=int(slots['nBins']),
            ncells=int(slots['nCells']),
            seed=int(slots['seed']),
            lib_mean1=float(slots['lib.mean1']),
            lib_mean2=float(slots['lib.mean2']),
            lib_sd1=float(slots['lib.sd1']),
            lib_sd2=float(slots['lib.sd2']),
            lib_prob=float(slots['lib.prob']),
            non_zero_prob=slots['non.zero.pro'],
            mean_coef0=float(slots['mean.coef0']),
            mean_coef1=float(slots['mean.coef1']),
            mean_coef2=float(slots['mean.coef2']),
            noise_mean=float(slots['noise.mean']),
            noise_sd=float(slots['noise.sd']),
        )

    @parameters.setter
    def parameters(self, params: Params):
        def _set_rparam(param_name, param_value):
            if param_value is not None:
                ro.globalenv['tmp'] = param_value
                ro.globalenv[self.__r_sim_obj_name] = ro.r(
                    'setParameters(%s, %s = tmp)' % (self.__r_sim_obj_name, param_name))

        _set_rparam('nBins', params.nbins)
        _set_rparam('nCells', params.ncells)
        _set_rparam('seed', params.seed)
        _set_rparam('lib.mean1', params.lib_mean1)
        _set_rparam('lib.mean2', params.lib_mean2)
        _set_rparam('lib.sd1', params.lib_sd1)
        _set_rparam('lib.sd2', params.lib_sd2)
        _set_rparam('lib.prob', params.lib_prob)
        if params.non_zero_prob is not None:
            ro.globalenv[self.__r_sim_obj_name] = ro.r(
                'setParameters(%s, default = FALSE)' % self.__r_sim_obj_name)
            _set_rparam('non.zero.pro', ro.FloatVector(params.non_zero_prob))
        _set_rparam('mean.coef0', params.mean_coef0)
        _set_rparam('mean.coef1', params.mean_coef1)
        _set_rparam('mean.coef2', params.mean_coef2)
        _set_rparam('noise.mean', params.noise_mean)
        _set_rparam('noise.sd', params.noise_sd)

    def estimate(self, adata: ad.AnnData) -> Params:
        """
        Estimates the simulation parameters using the given anndata
        :param adata:
        :return:
        """
        ro.globalenv[self.__r_original_name] = anndata2ri.py2rpy(adata)
        ro.r('%s <- simATACEstimate(%s)' % (self.__r_sim_obj_name, self.__r_original_name))
        return self.parameters

    def simulate_window_mtx(self) -> ad.AnnData:
        """
        Simulate bin by cell matrix.
        :return: AnnData object that contains the given simulation
        """
        tmp: ad.AnnData = anndata2ri.rpy2py(ro.r('simATACSimulate(%s)' % self.__r_sim_obj_name))
        return tmp

    def simulate_peak_mtx(self, peak_num=5000) -> ad.AnnData:
        """
        Simulate peak by cell matrix.
        :return: AnnData object that contains the given simulation
        """
        ro.r('sim <- simATACSimulate(%s)' % self.__r_sim_obj_name)
        tmp: ad.AnnData = anndata2ri.rpy2py(ro.globalenv['sim'])
        try:
            ro.r('mtx <- simATACgetPeakByCell(sim, peak.num=%s)' % peak_num)
            mtx: sparse.spmatrix = anndata2ri.rpy2py(ro.globalenv['mtx']).T
            peak_names = list(ro.r('rownames(mtx)'))
            adata = ad.AnnData(mtx, var=pd.DataFrame(index=peak_names), obs=tmp.obs.copy())
        finally:
            del ro.globalenv['sim']
            del ro.globalenv['mtx']
        return adata


@dataclasses.dataclass()
class Params:
    nbins: int = None
    ncells: int = None
    seed: int = None
    lib_mean1: float = None
    lib_mean2: float = None
    lib_sd1: float = None
    lib_sd2: float = None
    lib_prob: float = None
    non_zero_prob: np.ndarray = None
    mean_coef0: float = None
    mean_coef1: float = None
    mean_coef2: float = None
    noise_mean: float = None
    noise_sd: float = None


@dataclasses.dataclass
class DACellGroup:
    id: str
    ncells: int
    lib_mean: float
    accessibilities: List[DAAccessibility] = None


@dataclasses.dataclass
class DAFeatureGroup:
    """
    """
    id: str
    proportion: float


@dataclasses.dataclass
class DAAccessibility:
    """
    """
    feature_group_id: str
    score: float

    def __post_init__(self):
        if 0 > self.score or self.score > 1:
            raise ValueError('accessibility score has to be between 0 and 1')


@dataclasses.dataclass
class DAConfig:
    cell_groups: Dict[str, DACellGroup]
    feature_groups: Dict[str, DAFeatureGroup]

    def __post_init__(self):
        proportions = [fg.proportion for fg in self.feature_groups.values()]
        for prop in proportions:
            if 1 < prop or prop < 0:
                raise ValueError('proportions cant be negative or grater than 1')
        if sum(proportions) > 1:
            raise ValueError('the total proportions need to sum up to maximum 1')

        # validate feature group ids
        for fg_id in [acc.feature_group_id for cg in self.cell_groups.values() for acc in cg.accessibilities]:
            if fg_id not in self.feature_groups:
                raise ValueError('feature group with id %s has not been defined' % fg_id)

        for acc in self.cell_groups.values():
            if len(acc.accessibilities) != len(self.feature_groups):
                raise ValueError('each cell group should have the same amount of accesssibilities '
                                 'as the number of feature groups')


def simulate_da(params: Params, config: DAConfig, seed=42) -> ad.AnnData:
    # define feature groups idxs
    rng = np.random.default_rng(seed)

    fg_to_idxs = {}
    total_idxs = np.array(range(len(params.non_zero_prob)))
    feature_group_names = np.array(['noDA'] * len(total_idxs))
    for fg in config.feature_groups.values():
        fg_idxs = rng.choice(total_idxs, int(len(total_idxs) * fg.proportion))
        fg_to_idxs[fg.id] = fg_idxs
        feature_group_names[fg_idxs] = fg.id
        total_idxs = np.setdiff1d(total_idxs, fg_idxs)

    results = []
    for cell_group in config.cell_groups.values():
        accs = cell_group.accessibilities
        params_copy = copy.deepcopy(params)
        params_copy.ncells = cell_group.ncells
        params_copy.lib_mean1 = cell_group.lib_mean
        params_copy.lib_prob = 1
        for acc in accs:
            score = acc.score
            idxs = fg_to_idxs[acc.feature_group_id]
            params_copy.non_zero_prob[idxs] = score

        with Simulator() as simulator:
            simulator.parameters = params_copy
            tmp = simulator.simulate_window_mtx()
            tmp.obs['da_group'] = cell_group.id
            results.append(tmp)
    tmp = ad.concat(results, axis=0, index_unique='-')
    tmp.var['da_group'] = feature_group_names
    tmp.uns['simulation_params'] = params
    tmp.uns['da_config'] = config
    return tmp


def simulate(adata: ad.AnnData, params: Params, is_peak_mtx: bool = False, peak_num=5000):
    with Simulator() as simulator:
        simulator.estimate(adata)
        simulator.parameters = params
        if is_peak_mtx:
            adata_sim = simulator.simulate_peak_mtx(peak_num=peak_num)
        else:
            adata_sim = simulator.simulate_window_mtx()
        return adata_sim


def simulate_by_annotation(adata: ad.AnnData, params: Params, annotation_name: str = 'cell_type', is_peak_mtx: bool = False,
                           peak_num=5000) -> ad.AnnData:
    ncells = params.ncells
    if annotation_name not in adata.obs:
        raise ValueError('adata does not have an annotation called %s' % annotation_name)
    annotation_values = set(adata.obs[annotation_name]) - {'No Annotation'}
    simulations = []
    for annotation_value in annotation_values:
        params_copy = copy.deepcopy(params)
        logger.info(f'Simulating data for annotation {annotation_value}')
        adata_subset = adata[adata.obs[annotation_name] == annotation_value]
        with Simulator() as simulator:
            simulator.estimate(adata_subset)
            simulator.parameters = params_copy
            if is_peak_mtx:
                adata_sim = simulator.simulate_peak_mtx(peak_num=peak_num)
            else:
                adata_sim = simulator.simulate_window_mtx()
            adata_sim.obs[annotation_name] = [annotation_value] * ncells
            simulations.append(adata_sim)
    adata_sim = ad.concat(simulations, axis=0, index_unique='-')
    if not is_peak_mtx:
        adata_sim.var = adata.var.copy()
    return adata_sim
