"""Put the prat2026 package directory on sys.path so tests can import core."""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
