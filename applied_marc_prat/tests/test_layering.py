"""The package's dependency rules, as tests (CLAUDE.md, rules R1 and R5).

    materials  <-  pnm          materials  <-  growth          pnm and growth: strangers
    {materials, pnm, growth}  <-  coupling  <-  cases  <-  runs/, tests/

Read from the source with ``ast``, so a lazy import inside a function counts as well.
"""
import ast
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parents[1] / "saltv2"
LAYERS = ("materials", "pnm", "growth", "coupling", "cases")
ALLOWED = {
    "materials": set(),
    "pnm": {"materials"},
    "growth": {"materials"},
    "coupling": {"materials", "pnm", "growth"},
    "cases": {"materials", "pnm", "growth", "coupling"},
}


def _imports(path):
    """Absolute dotted names of everything ``path`` imports, relative imports resolved."""
    here = ("saltv2",) + path.relative_to(PACKAGE).parent.parts
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom):
            base = here[:len(here) - (node.level - 1)] if node.level else ()
            module = ".".join(base + ((node.module,) if node.module else ()))
            for alias in node.names:                     # ``from . import brine`` names a module
                yield f"{module}.{alias.name}"


@pytest.mark.parametrize("layer", LAYERS)
def test_a_layer_imports_only_the_layers_below_it(layer):
    offenders = []
    for path in sorted((PACKAGE / layer).rglob("*.py")):
        for name in _imports(path):
            parts = name.split(".")
            if parts[0] == "saltv2" and len(parts) > 1 and parts[1] in LAYERS:
                if parts[1] != layer and parts[1] not in ALLOWED[layer]:
                    offenders.append(f"{path.relative_to(PACKAGE)} imports {name}")
    assert not offenders, "\n".join(offenders)


def test_the_package_never_imports_the_run_scripts():
    """R1.  The old flat module names are listed so that a revived copy is caught too."""
    forbidden = {"runs", "pre_experiment", "pre_experiment_vapor", "feedback2d"}
    offenders = [f"{path.relative_to(PACKAGE)} imports {name}"
                 for path in sorted(PACKAGE.rglob("*.py")) for name in _imports(path)
                 if name.split(".")[0] in forbidden]
    assert not offenders, "\n".join(offenders)


def test_the_kernels_do_not_cite_the_run_directories():
    """R1, second half: a kernel's docstring states its own evidence.  ``cases`` is exempt -- a
    case value's provenance IS the execution record under runs/."""
    offenders = []
    for layer in ("materials", "pnm", "growth", "coupling"):
        for path in sorted((PACKAGE / layer).rglob("*.py")):
            for n, line in enumerate(path.read_text().splitlines(), 1):
                if "runs/" in line or "pre_experiment" in line or "feedback2d" in line:
                    offenders.append(f"{path.relative_to(PACKAGE)}:{n}: {line.strip()}")
    assert not offenders, "\n".join(offenders)
