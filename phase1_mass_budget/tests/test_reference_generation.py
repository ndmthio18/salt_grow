import hashlib
from pathlib import Path
import sys


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))

from run_reference_case import build_time_grid, generate  # noqa: E402


def test_time_grid_uses_configured_duration_and_timestep():
    assert build_time_grid(120.0, 30.0) == [0.0, 30.0, 60.0, 90.0, 120.0]
    assert build_time_grid(125.0, 60.0) == [0.0, 60.0, 120.0, 125.0]


def test_reference_generation_is_byte_reproducible_in_same_environment(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    first_summary = generate(first, total_duration_s=120.0, timestep_s=30.0)
    second_summary = generate(second, total_duration_s=120.0, timestep_s=30.0)
    assert first_summary["time_grid_point_count"] == 5
    assert first_summary == second_summary
    for path in sorted(first.iterdir()):
        peer = second / path.name
        assert peer.exists()
        assert hashlib.sha256(path.read_bytes()).digest() == hashlib.sha256(
            peer.read_bytes()
        ).digest()
