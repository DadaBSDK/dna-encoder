from __future__ import annotations

from pathlib import Path

import pytest

from dnastore.config import load_config
from dnastore.primers import load_primers

ROOT = Path(__file__).resolve().parent.parent
CODEC_CONFIGS = [
    ("naive2bit", {}),
    ("naive2bit", {"whiten": True}),
    ("goldman", {}),
]


@pytest.fixture(scope="session")
def primers():
    return load_primers(ROOT / "configs" / "primers.yaml")


def make_cfg(name: str, params: dict, **over) -> dict:
    return load_config(ROOT / "configs" / "default.yaml", {"codec": {"name": name, "params": params}, **over})
