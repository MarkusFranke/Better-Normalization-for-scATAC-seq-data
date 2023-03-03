from typing import Dict, List

import anndata
import episcanpy as epi
import anndata as ad
from . import performance, config
from pathlib import Path

from .. import simulation as sim
from .. import normalization as norm


def rank_features(adata: ad.AnnData, groupby='cell_group', reference='C_CONTROL'):
    epi.tl.rank_features(adata, groupby=groupby, omic='ATAC', use_raw=False, reference=reference,
                         n_features=None)
    ranks = adata.uns['rank_features_groups']
    cell_groups = ranks['names'].dtype.names
    stats = ['scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    for cg in cell_groups:
        for stat in stats:
            var_name = f'{stat}_{cg}'
            adata.var.loc[ranks['names'][cg], var_name] = ranks[stat][cg]
            print(f'{var_name} var added')


def get_rank_results(adata: ad.AnnData, run, scenario, fdr=0.05, method="raw") -> List[performance.DAResults]:
    cell_groups = list(set(adata.uns["da_config"].cell_groups.keys()) - {'C_CONTROL'})
    feature_groups = list(adata.uns["da_config"].feature_groups.keys())

    results = []
    for cg in cell_groups:
        p_vals = adata.var['pvals_adj_%s' % cg].array
        scores = adata.var['scores_%s' % cg].array
        fcs = adata.var['logfoldchanges_%s' % cg].array
        target = adata.var['feature_group'].array
        results.append(performance.DAResults(log_fcs=fcs, p_vals=p_vals, target=target, scores=scores,
                                             cell_group=cg, method=method, fdr=fdr, run=run, scenario=scenario))
        # for fg in feature_groups:
        #     name = f'{cg}_{fg}'
        #     results[name] = performance.DAResults(log_fcs=fcs, p_vals=p_vals, target=target, da_classes=[fg],
        #                                           name=name, level='feature_group', fdr=fdr)

    return results


def rank_normalization_methods(adata: ad.AnnData, run, scenario) -> List[
    performance.DAResults]:
    layers = adata.layers.keys()
    results = []
    for layer in layers:
        adata_cp = adata.copy()
        adata_cp.X = adata_cp.layers[layer]
        rank_features(adata_cp, groupby='cell_group', reference='C_CONTROL')
        results.extend(get_rank_results(adata_cp, method=layer, run=run, scenario=scenario))
    return results


def load_results(path: Path) -> List[performance.DAResults]:
    results = []
    for file in path.rglob('*.pkl'):
        results.append(performance.DAResults.from_pickle(file))
    results.sort(key=lambda x: f'{x.scenario}_{x.cell_group}_{x.method}_{x.fdr}_{x.run}')
    return results


def run_da(adata: anndata.AnnData, configs: Dict[str, config.DAConfig], output_dir: Path):
    if 'cell_type' not in adata.obs.keys():
        raise ValueError('cell_type must be included in obs')
    output_dir.mkdir(exist_ok=True)
    cell_types = list(set(adata.obs['cell_type']))
    cell_types.sort()

    print(cell_types)
    for cell_type in cell_types:
        print(f'Simulating cell type {cell_type}')
        cell_output_path = output_dir / cell_type.replace('/', '_').replace(' ', '')
        cell_output_path.mkdir(exist_ok=True)

        adata_cp = adata.copy()
        adata_cp = adata_cp[adata_cp.obs['cell_type'] == cell_type]
        params = sim.estimate_params(adata_cp)
        for scenario, cfg in configs.items():
            adata_sim = sim.simulate_da(params=params, config=cfg)
            norm.normalize(adata=adata_sim)
            results = rank_normalization_methods(adata_sim, scenario=scenario, run=cell_type.replace('/', '_').replace(' ', ''))
            for result in results:
                result.to_pickle(cell_output_path / f'{result.key}.pkl')
