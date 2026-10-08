"""YAML configuration with defaults. Configs are plain nested dicts, deep-merged over DEFAULTS."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

DEFAULTS: dict[str, Any] = {
    "oligo_len": 200,
    "header_copies": 5,
    "global_seed": 0,
    "primers": "primers.yaml",  # relative to the config file's directory
    "container": {"zstd_level": 19, "max_ratio": 0.98},
    "ecc": {"parity_permille": 150, "k_max": None},
    "codec": {"name": "goldman", "params": {}},
    "screening": {},  # see dnastore.screening.ScreeningConfig
    "compute_mfe": False,
    "workers": None,
}


def deep_merge(base: dict[str, Any], over: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict) and k != "params":
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_config(path: str | Path | None = None, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    """Load a YAML config, merge over defaults, and resolve the primer path."""
    cfg = copy.deepcopy(DEFAULTS)
    base_dir = Path(__file__).resolve().parent / "resources"
    if path is not None:
        path = Path(path)
        base_dir = path.parent
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}
        if not isinstance(raw, dict):
            raise ValueError("configuration must be a YAML mapping")
        unknown = set(raw) - (set(DEFAULTS) | {"seed_trits"})
        if unknown:
            raise ValueError(f"unknown configuration keys: {sorted(unknown)}")
        cfg = deep_merge(cfg, raw)
    cfg = deep_merge(cfg, overrides or {})
    p = Path(cfg["primers"])
    cfg["primers"] = str(p if p.is_absolute() else base_dir / p)
    for key in ("container", "ecc", "codec", "screening"):
        if not isinstance(cfg[key], dict):
            raise ValueError(f"{key} must be a mapping")
    if cfg["workers"] is not None and (not isinstance(cfg["workers"], int) or cfg["workers"] < 1):
        raise ValueError("workers must be null or a positive integer")
    return cfg
