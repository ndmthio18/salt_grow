from salt_phase2 import build_geometry_benchmarks, evaluate_geometry_benchmark


def test_synthetic_geometry_benchmarks_match_reference_answers():
    cases = build_geometry_benchmarks()
    assert [case.name for case in cases] == [
        "solid_pillar_exterior_bypass",
        "porous_slab_through_channel",
        "two_grains_connected_intergranular_gap",
        "closed_internal_cavity",
        "externally_connected_blind_pore",
    ]
    for case in cases:
        result = evaluate_geometry_benchmark(case)
        assert (
            result["window_bottom_to_top"]
            == result["expected_window_bottom_to_top"]
        )
        assert (
            result["material_bottom_to_top"]
            == result["expected_material_bottom_to_top"]
        )
        assert (
            result["material_external_void_voxels"]
            == result["expected_material_external_void_voxels"]
        )
        assert (
            result["material_enclosed_void_voxels"]
            == result["expected_material_enclosed_void_voxels"]
        )
        assert (
            result["accessible_intergranular_void_voxels"]
            == result["expected_accessible_intergranular_void_voxels"]
        )


def test_exterior_bypass_is_not_misclassified_as_pillar_porosity():
    pillar = build_geometry_benchmarks()[0]
    result = evaluate_geometry_benchmark(pillar)
    assert result["window_bottom_to_top"]
    assert not result["material_bottom_to_top"]
    assert result["exterior_bypass_detected"]


def test_blind_pore_is_external_but_not_through_connected():
    blind = build_geometry_benchmarks()[-1]
    result = evaluate_geometry_benchmark(blind)
    assert result["material_external_void_voxels"] == 3
    assert not result["material_bottom_to_top"]
    assert result["material_enclosed_void_voxels"] == 0
