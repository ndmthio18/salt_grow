"""Test-path setup for Phase 2A and its verified Phase 1 dependency."""

from pathlib import Path
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "phase2_ca"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase2_ca" / "src"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))
