import pytest
import src.data as data
import src.simulation as sim
from pathlib import Path
import numpy as np
import anndata as ad


@pytest.fixture
def simple_da_config() -> sim.DAConfig:
    cg_ids = ['CA', 'CB', 'CC']
    ncells = [50, 50, 50]
    lib_means = [20, 20, 20]
    cgs = {}
    for idx in range(len(cg_ids)):
        cgs[cg_ids[idx]] = sim.DACellGroup(id=cg_ids[idx], ncells=ncells[idx], lib_mean=lib_means[idx])

    proportions = [0.3, 0.2, 0.1]
    fg_ids = ['FA', 'FB', 'FC']
    fgs = {}
    for idx in range(len(fg_ids)):
        fgs[fg_ids[idx]] = sim.DAFeatureGroup(id=fg_ids[idx], proportion=proportions[idx])

    scores = np.array([0.0, 0, 0.0])
    for cg in cgs.values():
        accs = []
        for idx, fg in enumerate(fgs.values()):
            score = scores[idx]
            accs.append(sim.DAAccessibility(score=score, feature_group_id=fg.id))
        cg.accessibilities = accs
        scores = np.minimum(scores + 1, np.array([0.8] * len(scores)))

    return sim.DAConfig(cell_groups=cgs, feature_groups=fgs)


@pytest.fixture
def mouse_brain_adata() -> ad.AnnData:
    epidata = data.MouseBrainDataset()
    return epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))


@pytest.fixture
def da_template() -> str:
    return "data/simulation/templates/template20220205.json"


@pytest.fixture
def estimated_params(mouse_brain_adata) -> sim.Params:
    mouse_brain_adata = mouse_brain_adata[mouse_brain_adata.obs['cell_type'] == 'L5 IT']
    with sim.Simulator() as simulator:
        return simulator.estimate(mouse_brain_adata)


def test_parameters():
    with sim.Simulator() as simulator:
        print(simulator.parameters)
        assert simulator.parameters.ncells == 500


def test_get_set_parameter():
    with sim.Simulator() as simulator:
        val = 10000
        params = sim.Params(ncells=val)
        simulator.parameters = params
        assert simulator.parameters.ncells == val


def test_estimate_simulate(mouse_brain_adata):
    with sim.Simulator() as simulator:
        simulator.estimate(mouse_brain_adata)
        key = 'nCells'
        simulator[key] = mouse_brain_adata.shape[0]
        print(simulator.parameters)
        assert mouse_brain_adata.shape[0] == simulator[key]
        adata_sim = simulator.simulate()
        assert adata_sim.shape == mouse_brain_adata.shape


def test_simulate_by_annotation(mouse_brain_adata):
    annotation_name = 'cell_type'
    ncells = 20
    simulation_params = sim.Params(ncells=ncells, noise_mean=-0.3, noise_sd=0.3)
    adata_sim = sim.simulate_by_annotation(mouse_brain_adata, simulation_params, annotation_name=annotation_name)
    num_simulations = len(set(mouse_brain_adata.obs[annotation_name]) - {'No Annotation'})
    assert adata_sim.shape == (ncells * num_simulations, mouse_brain_adata.shape[1])


def test_get_simulation_params():
    with sim.Simulator() as simulator:
        assert simulator.parameters.noise_mean == 0


def test_set_simulation_params():
    with sim.Simulator() as simulator:
        params = sim.Params(ncells=1)
        simulator.parameters = params
        assert simulator.parameters.ncells == 1


def test_set_simulation_params_estimate(mouse_brain_adata):
    with sim.Simulator() as simulator:
        params = simulator.estimate(mouse_brain_adata)
        params.non_zero_prob[1] = 0.5
        simulator.parameters = params
        assert simulator.parameters.non_zero_prob[1] == 0.5


def test_simulate_da(simple_da_config, mouse_brain_adata, estimated_params):
    adata_sim = sim.simulate_da(params=estimated_params, config=simple_da_config)


def test_read_da_from_file(da_template):
    cfgs = sim.DAConfig.from_template(da_template)
    for k, cg in cfgs['S0'].cell_groups.items():
        assert cg.lib_mean == 14
        for acc in cg.accessibilities:
            if cg.id == 'CA' and acc.feature_group_id == 'FA':
                assert acc.score == 0.01
