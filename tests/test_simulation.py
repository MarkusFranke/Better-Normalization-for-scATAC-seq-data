import src.data as data
import src.simulation as sim


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
        adata = data.load_mousebrain_window_mtx()[:1000, :].copy()
        simulator.estimate(adata)
        key = 'nCells'
        simulator[key] = adata.shape[0]
        print(simulator.parameters)
        assert adata.shape[0] == simulator[key]
        adata_sim = simulator.simulate()
        assert adata_sim.shape == adata.shape


def test_simulate_by_annotation():
    annotation_name = 'cell_type'
    ncells = 20
    adata = data.load_mousebrain_window_mtx()
    adata_sim = sim.simulate_by_annotation(adata, annotation_name=annotation_name, nCells=ncells,)
    num_simulations = len(set(adata.obs[annotation_name]) - {'No Annotation'})
    assert adata_sim.shape == (ncells * num_simulations, adata.shape[1])
