#!/usr/bin/env python
"""Smoke test: every cnvlib submodule imports cleanly.

The per-command unit tests now live in test_<command>.py files; this module
keeps the cheap import-time smoke check that previously rode along with them.
"""

import importlib

import pytest

CNVLIB_SUBMODULES = [
    "access",
    "antitarget",
    "autobin",
    "batch",
    "bintest",
    "call",
    "cluster",
    "cmdutil",
    "cnary",
    "commands",
    "core",
    "coverage",
    "diagram",
    "export",
    "fix",
    "heatmap",
    "import_rna",
    "importers",
    "metrics",
    "parallel",
    "params",
    "plots",
    "purity",
    "reference",
    "reports",
    "samutil",
    "scatter",
    "segfilters",
    "segmentation",
    "segmetrics",
    "smoothing",
    "vary",
]


@pytest.mark.parametrize("submodule", CNVLIB_SUBMODULES)
def test_submodule_imports(submodule):
    importlib.import_module(f"cnvlib.{submodule}")
