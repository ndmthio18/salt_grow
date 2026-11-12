"""The Prat (2026) growth law, as implemented and validated in ``prat2026/core.py``.

``prat2026/`` is a verbatim copy of a flat, uninstalled package (its modules import each other as
``core`` / ``params``), so it has to be put on ``sys.path`` to be used.  This module is the ONE
place that does it; everything in ``saltv2`` writes ``from saltv2.growth.prat import prat``.
"""
import sys
from pathlib import Path

PRAT2026 = Path(__file__).resolve().parents[2] / "prat2026"
if not (PRAT2026 / "core.py").is_file():
    raise ImportError(f"the Prat (2026) growth law is expected at {PRAT2026 / 'core.py'}")
if str(PRAT2026) not in sys.path:
    sys.path.insert(0, str(PRAT2026))
import core as prat                                                            # noqa: E402

__all__ = ["prat"]
