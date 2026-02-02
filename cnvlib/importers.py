"""Import from other formats to the CNVkit format."""

from __future__ import annotations
import logging
from typing import TYPE_CHECKING, Union

import numpy as np
from skgenome import tabio

from . import params

if TYPE_CHECKING:
    from collections.abc import Iterator
    from cnvlib.cnary import CopyNumArray


# __________________________________________________________________________
# import-picard


def do_import_picard(fname, too_many_no_coverage=100):
    """Import coverage data from Picard CalculateHsMetrics output.

    Converts Picard CalculateHsMetrics output to CNVkit format by reading the
    interval-level coverage data, cleaning up gene names, and calculating log2
    ratios.

    Parameters
    ----------
    fname : str
        Path to Picard CalculateHsMetrics output file (typically ends in
        .per_target_coverage or .interval_summary).
    too_many_no_coverage : int, optional
        Threshold for warning about excessive zero-coverage bins.
        Default: 100.

    Returns
    -------
    CopyNumArray
        CNVkit-compatible genomic array with columns including 'gene', 'log2',
        and 'ratio' (coverage relative to average).

    Notes
    -----
    Bins with zero coverage are assigned the NULL_LOG2_COVERAGE value to avoid
    math domain errors in log2 transformation.

    Gene names from overlapping intervals are deduplicated (see `unpipe_name`).

    See Also
    --------
    unpipe_name : Cleans up duplicated gene names from Picard output
    """
    garr = tabio.read(fname, "picardhs")
    garr["gene"] = garr["gene"].apply(unpipe_name)
    # Create log2 column from coverages, avoiding math domain error
    coverages = garr["ratio"].copy()
    no_cvg_idx = coverages == 0
    if no_cvg_idx.sum() > too_many_no_coverage:
        logging.warning(
            "WARNING: Sample %s has >%d bins with no coverage",
            garr.sample_id,
            too_many_no_coverage,
        )
    coverages[no_cvg_idx] = 2**params.NULL_LOG2_COVERAGE
    garr["log2"] = np.log2(coverages)
    return garr


def unpipe_name(name):
    """Fix the duplicated gene names Picard spits out.

    Return a string containing the single gene name, sans duplications and pipe
    characters.

    Picard CalculateHsMetrics combines the labels of overlapping intervals
    by joining all labels with '|', e.g. 'BRAF|BRAF' -- no two distinct
    targeted genes actually overlap, though, so these dupes are redundant.
    Meaningless target names are dropped, e.g. 'CGH|FOO|-' resolves as 'FOO'.
    In case of ambiguity, the longest name is taken, e.g. "TERT|TERT Promoter"
    resolves as "TERT Promoter".
    """
    if "|" not in name:
        return name
    gene_names = set(name.split("|"))
    if len(gene_names) == 1:
        return gene_names.pop()
    cleaned_names = gene_names.difference(params.IGNORE_GENE_NAMES)
    if cleaned_names:
        gene_names = cleaned_names
    new_name = sorted(gene_names, key=len, reverse=True)[0]
    if len(gene_names) > 1:
        logging.warning("WARNING: Ambiguous gene name %r; using %r", name, new_name)
    return new_name


# __________________________________________________________________________
# import-theta


def do_import_theta(
    segarr: CopyNumArray, theta_results_fname: str, ploidy: int = 2
) -> Iterator[CopyNumArray]:
    """Import absolute copy number calls from THetA results.

    THetA (Tumor Heterogeneity Analysis) infers tumor purity and clonal/subclonal
    copy number profiles. This function converts THetA's output into CNVkit
    segment format with absolute copy numbers.

    Parameters
    ----------
    segarr : CopyNumArray
        Input segmented copy number data (.cns file). Typically the segmentation
        results before calling absolute copy numbers.
    theta_results_fname : str
        Path to THetA results file (.results or .withBounds file).
    ploidy : int, optional
        Expected baseline ploidy (usually 2 for diploid). Used to calculate
        log2 ratios from absolute copy numbers. Default: 2.

    Yields
    ------
    CopyNumArray
        One or more segment arrays with 'cn' (absolute copy number) and
        recalculated 'log2' values. THetA may output multiple solutions,
        each yielded separately.

    Notes
    -----
    - Only autosomal segments are processed; sex chromosomes are excluded
      because THetA doesn't handle them well.
    - Segments with unknown copy number (marked as None/X in THetA output)
      are dropped.
    - Copy number 0 is treated as 0.5 for log2 calculation to avoid -inf.

    See Also
    --------
    parse_theta_results : Parses the THetA results file format
    """
    theta = parse_theta_results(theta_results_fname)
    # THetA doesn't handle sex chromosomes well
    segarr = segarr.autosomes()
    for copies in theta["C"]:
        if len(copies) != len(segarr):
            copies = copies[: len(segarr)]
        # Drop any segments where the C value is None
        mask_drop = np.array([c is None for c in copies], dtype="bool")
        segarr = segarr[~mask_drop].copy()
        ok_copies = np.asarray([c for c in copies if c is not None], dtype=float)
        # Replace remaining segment values with these integers
        segarr["cn"] = ok_copies.astype("int")
        ok_copies[ok_copies == 0] = 0.5
        segarr["log2"] = np.log2(ok_copies / ploidy)
        segarr.sort_columns()
        yield segarr


def parse_theta_results(
    fname: str,
) -> dict[str, Union[float, list[float], list[list[int]], list[list[float]]]]:
    """Parse THetA results into a data structure.

    Columns: NLL, mu, C, p*
    """
    with open(fname) as handle:
        header = next(handle).rstrip().split("\t")
        body = next(handle).rstrip().split("\t")
        assert len(body) == len(header) == 4

        # NLL
        nll = float(body[0])

        # mu
        mu = body[1].split(",")
        mu_normal = float(mu[0])
        mu_tumors = list(map(float, mu[1:]))

        # C
        copies = body[2].split(":")
        if len(mu_tumors) == 1:
            # 1D array of integers
            # Replace X with None for "missing"
            copies = [[int(c) if c.isdigit() else None for c in copies]]
        else:
            # List of lists of integer-or-None (usu. 2 x #segments)
            copies = [
                [int(c) if c.isdigit() else None for c in subcop]
                for subcop in zip(*[c.split(",") for c in copies], strict=False)
            ]

        # p*
        probs = body[3].split(",")
        if len(mu_tumors) == 1:
            # 1D array of floats, or None for "X" (missing/unknown)
            probs = [float(p) if not p.isalpha() else None for p in probs]
        else:
            probs = [
                [float(p) if not p.isalpha() else None for p in subprob]
                for subprob in zip(*[p.split(",") for p in probs], strict=False)
            ]
    return {
        "NLL": nll,
        "mu_normal": mu_normal,
        "mu_tumors": mu_tumors,
        "C": copies,
        "p*": probs,
    }
