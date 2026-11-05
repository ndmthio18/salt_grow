"""Fail-closed cache identity for configuration, model source and data bytes."""
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
import json


def source_hash() -> str:
    digest = sha256()
    for path in sorted(Path(__file__).parent.glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def identity(p, n, closure, initial_profile="equilibrium") -> dict:
    config = dict(physical=asdict(p), numerical=asdict(n), closure=closure,
                  initial_profile=initial_profile)
    return dict(config_sha256=sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
                code_sha256=source_hash())


def cache_matches(metadata_path, data_path, expected) -> bool:
    if not metadata_path.exists() or not data_path.exists():
        return False
    try:
        record = json.loads(metadata_path.read_text(encoding="utf-8"))
        return (record.get("identity") == expected and
                record.get("data_sha256") == sha256(data_path.read_bytes()).hexdigest())
    except (OSError, ValueError):
        return False
