from typing import Dict

import epinorm.diffacc as diffacc
import anndata as ad
import pytest
import epinorm.simulation as sim
import epinorm.data as data
from pathlib import Path
from epinorm import normalization as norm
from epinorm.diffacc.config import DAConfig


@pytest.fixture
def da_template() -> str:
    return "data/simulation/templates/template20220205.json"


@pytest.fixture
def mouse_brain_adata() -> ad.AnnData:
    epidata = data.MouseBrainDataset()
    return epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))


@pytest.fixture
def da_config(da_template) -> diffacc.config.DAConfig:
    cfgs = diffacc.config.DAConfig.from_template(da_template)
    return cfgs['S0']


@pytest.fixture
def estimated_params(mouse_brain_adata) -> sim.Params:
    mouse_brain_adata = mouse_brain_adata[mouse_brain_adata.obs['cell_type'] == 'L5 IT']
    with sim.Simulator() as simulator:
        return simulator.estimate(mouse_brain_adata)


def test_read_da_from_file(da_template):
    cfgs = diffacc.config.DAConfig.from_template(da_template)
    for k, cg in cfgs['S0'].cell_groups.items():
        assert cg.lib_mean == 14
        for acc in cg.accessibilities:
            if cg.id == 'C_CONTROL' and acc.feature_group_id == 'F_1':
                assert acc.score == 0.01


def test_da(da_config, estimated_params):
    adata_sim = sim.simulate_da(params=estimated_params, config=da_config)
    norm.normalize(adata_sim)
    diffacc.utils.rank_features(adata_sim, groupby='cell_group', reference='C_CONTROL')
    assert "pvals_adj_C_1" in adata_sim.var.keys()
    print(1)
