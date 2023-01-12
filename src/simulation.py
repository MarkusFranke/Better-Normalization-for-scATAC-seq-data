import anndata as ad
import anndata2ri
import pandas as pd
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import logging
from scipy import sparse

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
    def parameters(self):
        """
        returns the simulation parameters
        :return:
        """
        sim_obj = ro.globalenv[self.__r_sim_obj_name]
        return {k: v for k, v in sim_obj.slots.items() if
                type(v) not in {ro.vectors.StrVector, ro.vectors.BoolVector}}

    def __setitem__(self, key, value):
        """
        Sets the given simulation parameter.
        :param key:
        :param value:
        :return:
        """
        ro.globalenv[self.__r_sim_obj_name] = ro.r('setParameters(%s, %s = %s)' % (self.__r_sim_obj_name, key, value))

    def __getitem__(self, item: str):
        """
        Gets the given simulation parameter.
        :param item:
        :return:
        """
        return ro.r('getParameters(%s, "%s")' % (self.__r_sim_obj_name, item))[0]

    def estimate(self, adata: ad.AnnData):
        """
        Estimates the simulation parameters using the given anndata
        :param adata:
        :return:
        """
        ro.globalenv[self.__r_original_name] = anndata2ri.py2rpy(adata)
        ro.globalenv[self.__r_sim_obj_name] = ro.r('simATACEstimate(%s)' % self.__r_original_name)

    def simulate_window_mtx(self) -> ad.AnnData:
        """
        Simulate bin by cell matrix.
        :return: AnnData object that contains the given simulation
        """
        return anndata2ri.rpy2py(ro.r('simATACSimulate(%s)' % self.__r_sim_obj_name))

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


def simulate(adata: ad.AnnData, is_peak_mtx: bool = False, peak_num=5000, **sim_parameters):
    with Simulator() as simulator:
        simulator.estimate(adata)
        for param_name, param_val in sim_parameters.items():
            simulator[param_name] = param_val
        if is_peak_mtx:
            adata_sim = simulator.simulate_peak_mtx(peak_num=peak_num)
        else:
            adata_sim = simulator.simulate_window_mtx()
        return adata_sim


def simulate_by_annotation(adata: ad.AnnData, annotation_name: str = 'cell_type', is_peak_mtx: bool = False,
                           peak_num=5000, **sim_parameters) -> ad.AnnData:
    ncells = 100
    if annotation_name not in adata.obs:
        raise ValueError('adata does not have an annotation called %s' % annotation_name)
    annotation_values = set(adata.obs[annotation_name]) - {'No Annotation'}
    simulations = []
    for annotation_value in annotation_values:
        logger.info(f'Simulating data for annotation {annotation_value}')
        adata_subset = adata[adata.obs[annotation_name] == annotation_value]
        with Simulator() as simulator:
            simulator.estimate(adata_subset)
            for param_name, param_val in sim_parameters.items():
                if param_name == 'nCells':
                    ncells = param_val
                simulator[param_name] = param_val
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
