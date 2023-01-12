import src.data as data
import src.simulation as sim
from pathlib import Path


def test_parameters():
    with sim.Simulator() as simulator:
        print(simulator.parameters)
        assert len(simulator.parameters) == 15


def test_get_set_parameter():
    with sim.Simulator() as simulator:
        key = 'nCells'
        val = 10000
        simulator[key] = val
        assert simulator[key] == val


def test_estimate_simulate():
    with sim.Simulator() as simulator:
        epidata = data.MouseBrainDataset()
        adata = epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))
        simulator.estimate(adata)
        key = 'nCells'
        simulator[key] = adata.shape[0]
        print(simulator.parameters)
        assert adata.shape[0] == simulator[key]
        adata_sim = simulator.simulate()
        assert adata_sim.shape == adata.shape


def test_estimate_simulate_peak():
    epidata = data.MouseBrainDataset()
    adata = epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))
    simulation_params = {'nCells': adata.shape[0], 'noise.mean': -0.3, 'noise.sd': 0.3, 'is_peak_mtx': True, 'peak_num': 5000}
    adata_sim = sim.simulate(adata, **simulation_params)
    assert adata_sim.shape == (adata.shape[0], 5000)



def test_simulate_by_annotation_peak():
    annotation_name = 'cell_type'
    ncells = 20
    epidata = data.MouseBrainDataset()
    adata = epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))
    adata_sim = sim.simulate_by_annotation(adata, annotation_name=annotation_name, nCells=ncells, is_peak_mtx=True,
                                           peak_num=5000)
    num_simulations = len(set(adata.obs[annotation_name]) - {'No Annotation'})
    # assert adata_sim.shape == (ncells * num_simulations, 5000)


def test_simulate_by_annotation():
    annotation_name = 'cell_type'
    ncells = 20
    epidata = data.MouseBrainDataset()
    adata = epidata.load_mtx(matrix_path=Path('data/mouse_brain_5k/bin_by_cell.h5ad'))
    adata_sim = sim.simulate_by_annotation(adata, annotation_name=annotation_name, nCells=ncells)
    num_simulations = len(set(adata.obs[annotation_name]) - {'No Annotation'})
    assert adata_sim.shape == (ncells * num_simulations, adata.shape[1])