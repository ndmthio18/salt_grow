"""The 5 x 5 DuMuX comparison grid with a complete Dirichlet reservoir row."""

import numpy as np
import pytest

from saltv2.pnm import DGF_DUMUX_1D, DGF_DUMUX_5X5_VARTHROAT, PoreLayer


def test_the_minimum_y_row_is_eliminated_without_inventing_reservoir_throats():
    layer = PoreLayer.from_dgf_dirichlet_boundary(DGF_DUMUX_5X5_VARTHROAT)

    assert layer.dgf_vertices.shape == (30, 5)
    assert layer.dgf_edges.shape == (45, 5)
    assert layer.N == 25
    assert len(layer.lat) == 36
    assert layer.reservoir_file_indices.tolist() == [0, 5, 10, 15, 20]
    assert layer.reservoir_internal_file_edge_indices.tolist() == [20, 25, 30, 35]
    assert np.count_nonzero(layer.rt_res) == 5
    assert np.count_nonzero(layer.res_file_edge_index >= 0) == 5


def test_the_top_row_and_file_index_map_are_exact():
    layer = PoreLayer.from_dgf_dirichlet_boundary(DGF_DUMUX_5X5_VARTHROAT)

    assert layer.file_vertex_index[layer.surface].tolist() == [25, 26, 27, 28, 29]
    assert layer.file_vertex_index[layer.bottom].tolist() == [1, 6, 11, 16, 21]
    assert layer.n_surface == 5 and layer.n_layers == 5 and layer.n_side == 5
    assert np.allclose(layer.height[layer.surface], 3.5e-3)
    assert np.allclose(layer.height[layer.bottom], 0.7e-3)
    assert np.all(layer.rt_res[layer.bottom] > 0.0)
    assert np.all(np.delete(layer.rt_res, layer.bottom) == 0.0)
    assert sorted(layer.lat_file_edge_index.tolist()) == [
        edge for edge in range(45)
        if edge not in {0, 4, 8, 12, 16, 20, 25, 30, 35}
    ]


def test_the_exact_dgf_geometry_reaches_the_kernels():
    layer = PoreLayer.from_dgf_dirichlet_boundary(DGF_DUMUX_5X5_VARTHROAT)
    g_lat, g_res = layer.conductances()
    d_lat, d_res = layer.diffusive()

    assert g_lat.shape == d_lat.shape == (36,)
    assert g_res.shape == d_res.shape == (25,)
    assert np.all(g_lat > 0.0) and np.all(d_lat > 0.0)
    assert np.all(g_res[layer.bottom] > 0.0) and np.all(d_res[layer.bottom] > 0.0)
    assert np.all(np.delete(g_res, layer.bottom) == 0.0)
    assert np.all(np.delete(d_res, layer.bottom) == 0.0)


def test_invalid_boundary_selection_is_rejected():
    with pytest.raises(ValueError, match="boundary_axis"):
        PoreLayer.from_dgf_dirichlet_boundary(DGF_DUMUX_5X5_VARTHROAT, boundary_axis=2)
    with pytest.raises(ValueError, match="boundary_side"):
        PoreLayer.from_dgf_dirichlet_boundary(DGF_DUMUX_5X5_VARTHROAT, boundary_side="left")


def test_the_legacy_one_reservoir_column_is_unchanged():
    layer = PoreLayer.from_dgf(DGF_DUMUX_1D)
    assert layer.N == 5 and layer.n_surface == 1
    assert np.count_nonzero(layer.rt_res) == 1

