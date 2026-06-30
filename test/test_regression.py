#!/usr/bin/env python
"""End-to-end numerical regression test on checked-in fixtures (#1039).

Runs the canonical CNVkit pipeline -- ``import-picard`` -> ``reference`` ->
``fix`` -> ``segment`` -> ``call`` -- on the small Picard-derived fixtures under
``test/picard/`` and compares the ``.cnr``/``.cns``/``.call.cns`` outputs against
frozen reference files checked into ``test/formats/regression/``.

The pipeline uses the *haar* segmenter (pure-Python, deterministic, no R
dependency) so the result is reproducible in CI without the stochasticity of CBS.
This guards the whole numerical stack at once: a change in import, pooled-
reference aggregation, bias correction, segmentation, or copy-number calling that
moves any output value beyond tolerance makes this test fail with a pointed diff.

Comparison policy
-----------------
- Float columns (any column read back as floating-point -- ``log2``, ``depth``,
  ``weight``): ``np.isclose`` with ``rtol=1e-5, atol=1e-6``. The frozen fixtures
  are stored at tabio's ``%.6g`` precision, whose quantization alone reaches
  ~5e-6 relative error (measured); ``rtol=1e-5`` leaves ~2x headroom while still
  catching genuine numerical drift, and absorbs cross-BLAS/cross-OS
  floating-point differences. (Verified empirically: identical structure and
  sub-tolerance floats on Linux x86_64 / numpy 2.5 vs macOS arm64 / numpy 2.4.)
  A corollary is the sensitivity floor: a uniform drift smaller than ~1e-5
  relative is below the fixture's own resolution and cannot be detected here; for
  finer sensitivity, store fixtures at higher precision and tighten the tolerance.
- Every other column (``chromosome``, ``start``, ``end``, ``gene`` and the
  integer ``cn``, ``probes``): **exact** match. Coordinates are copied through
  from the inputs and the integers are discrete decisions; any change is a real
  regression, not floating-point noise. Haar segment endpoints are deterministic
  and were verified bit-identical across OS/arch.

Regenerating the fixtures
-------------------------
Deliberate, never on a routine test run::

    cd test && python test_regression.py --regenerate
    # or, equivalently:
    make -C test regression-fixtures

A PR that regenerates these fixtures MUST explain the intended numerical change
in its description: an unexplained fixture diff is a red flag to scrutinize, not
a rubber stamp.
"""

import logging
import os
import sys
import tempfile
import unittest
import warnings

import numpy as np
import pytest

import cnvlib
from cnvlib import commands, importers
from skgenome import tabio

logging.basicConfig(level=logging.ERROR, format="%(message)s")
warnings.filterwarnings("ignore", category=DeprecationWarning)

# The full pipeline run takes several seconds (see the >5s `slow` marker in
# pyproject.toml), so it is skipped by the `pytest -m 'not slow'` fast path.
pytestmark = pytest.mark.slow

# Resolve fixture locations relative to this file so the module works both under
# pytest (conftest chdirs into test/) and when run standalone for regeneration.
HERE = os.path.dirname(os.path.abspath(__file__))
PICARD_DIR = os.path.join(HERE, "picard")
FROZEN_DIR = os.path.join(HERE, "formats", "regression")

# Normal samples (the "_5" replicates) pooled into the reference, mirroring
# test/Makefile's `reference build/p2-*_5.*targetcoverage.cnn`.
NORMAL_IDS = ("p2-5_5", "p2-9_5", "p2-20_5")
# Tumor sample carrying real copy-number events across all 24 chromosomes.
SAMPLE_ID = "p2-9_2"

# Tolerance for float columns; set by %.6g fixture quantization (~5e-6 rel),
# not raw FP drift. Float vs exact comparison is chosen per-column by dtype.
RTOL = 1e-5
ATOL = 1e-6


def _import_picard(sample_id, kind):
    """Import one Picard CalculateHsMetrics CSV into a GenomicArray."""
    path = os.path.join(PICARD_DIR, f"{sample_id}.{kind}coverage.csv")
    return importers.do_import_picard(path)


def run_pipeline():
    """Run import -> reference -> fix -> segment -> call; return the three stages.

    The imported coverages are round-tripped through ``.cnn`` files because
    ``do_reference`` accepts only file paths (it reads the normals itself), and
    because ``do_import_picard`` yields a ``GenomicArray`` whereas ``do_fix``
    needs the ``CopyNumArray`` that ``cnvlib.read`` produces. This also mirrors
    the CLI flow exactly (``import-picard`` writes ``.cnn``; ``reference``/``fix``
    read them back).

    Returns a dict mapping fixture basename ('p2-9_2.cnr', '.cns', '.call.cns')
    to the freshly computed CopyNumArray.
    """
    with tempfile.TemporaryDirectory() as tmp:

        def import_to_cnn(sample_id, kind):
            path = os.path.join(tmp, f"{sample_id}.{kind}coverage.cnn")
            tabio.write(_import_picard(sample_id, kind), path)
            return path

        tgt_fnames = [import_to_cnn(s, "target") for s in NORMAL_IDS]
        anti_fnames = [import_to_cnn(s, "antitarget") for s in NORMAL_IDS]
        reference = commands.do_reference(
            tgt_fnames, anti_fnames, is_haploid_x_reference=True
        )
        sample_target = cnvlib.read(import_to_cnn(SAMPLE_ID, "target"))
        sample_antitarget = cnvlib.read(import_to_cnn(SAMPLE_ID, "antitarget"))

    cnr = commands.do_fix(sample_target, sample_antitarget, reference)
    cns = commands.do_segmentation(cnr, "haar")
    call_cns = commands.do_call(cns, method="threshold", is_haploid_x_reference=True)
    return {
        f"{SAMPLE_ID}.cnr": cnr,
        f"{SAMPLE_ID}.cns": cns,
        f"{SAMPLE_ID}.call.cns": call_cns,
    }


def regenerate():
    """Write the frozen reference fixtures. Invoked deliberately, never in CI."""
    os.makedirs(FROZEN_DIR, exist_ok=True)
    for basename, cna in run_pipeline().items():
        path = os.path.join(FROZEN_DIR, basename)
        tabio.write(cna, path)
        # Surface any NaNs so a maintainer consciously confirms they are intended
        # before checking the fixture in (a NaN frozen as canonical would be
        # accepted on every later run; see test docstring comparison policy).
        nan_cols = {
            col: int(cna.data[col].isna().sum())
            for col in cna.data.columns
            if cna.data[col].dtype.kind == "f" and cna.data[col].isna().any()
        }
        note = f"  WARNING NaN present: {nan_cols}" if nan_cols else ""
        print(f"wrote {path} ({len(cna)} rows){note}")


class RegressionTests(unittest.TestCase):
    """Pipeline outputs must match the frozen fixtures within tolerance."""

    def _compare(self, basename, fresh):
        frozen_path = os.path.join(FROZEN_DIR, basename)
        self.assertTrue(
            os.path.exists(frozen_path),
            f"Missing frozen fixture {frozen_path}; regenerate with "
            "`python test_regression.py --regenerate`.",
        )
        frozen = cnvlib.read(frozen_path)
        self.assertGreater(
            len(frozen), 0, f"{basename}: frozen fixture is empty -- regenerate it."
        )
        self.assertEqual(
            len(fresh),
            len(frozen),
            f"{basename}: row count changed {len(frozen)} -> {len(fresh)} "
            "(segmentation/calling regression). If intended, regenerate the "
            "fixtures and explain the change in the PR.",
        )
        self.assertEqual(
            list(fresh.data.columns),
            list(frozen.data.columns),
            f"{basename}: column set changed.",
        )

        for col in frozen.data.columns:
            fresh_col = fresh.data[col].to_numpy()
            frozen_col = frozen.data[col].to_numpy()
            # Float columns get tolerance; coordinates, gene labels, and integer
            # calls are compared exactly. Detect floats by dtype so a future
            # output column cannot silently fall into the wrong path.
            if fresh_col.dtype.kind == "f" or frozen_col.dtype.kind == "f":
                bad = ~np.isclose(
                    fresh_col, frozen_col, rtol=RTOL, atol=ATOL, equal_nan=True
                )
                detail = f"drifted beyond tolerance (rtol={RTOL}, atol={ATOL})"
            else:
                bad = fresh_col != frozen_col
                detail = "changed"
            if bad.any():
                i = int(np.argmax(bad))
                self.fail(
                    f"{basename}: column '{col}' {detail} at row {i}: "
                    f"frozen={frozen_col[i]!r} fresh={fresh_col[i]!r} "
                    f"({int(bad.sum())} of {len(bad)} rows differ). If intended, "
                    "regenerate the fixtures and explain the change in the PR."
                )

    def test_pipeline_matches_frozen(self):
        """fix -> segment -> call outputs match the frozen references.

        One pipeline run, three fixtures compared via subtests so a failure in
        any stage is reported distinctly while the others still run.
        """
        outputs = run_pipeline()
        for basename, fresh in outputs.items():
            with self.subTest(fixture=basename):
                self._compare(basename, fresh)


if __name__ == "__main__":
    if "--regenerate" in sys.argv:
        regenerate()
    else:
        unittest.main()
