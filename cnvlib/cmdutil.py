"""Functions reused within command-line implementations."""

from __future__ import annotations
import logging
import sys

import numpy as np

from skgenome import tabio

from .cnary import CopyNumArray as CNA
from skgenome import GenomicArray as GA
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

    import pandas as pd

    from cnvlib.cnary import CopyNumArray
    from cnvlib.vary import VariantArray
    from skgenome.gary import GenomicArray

# Median |alt_freq - 0.5| above which the variant set looks unfit for BAF
# analysis (heterozygous germline SNPs cluster tightly around 0.5).
_BAF_INPUT_MEDIAN_TOL = 0.2
# Below this many heterozygous variants, distribution checks are unreliable.
_BAF_INPUT_MIN_HET = 50


def read_cna(
    infile: str, sample_id: str | None = None, meta: dict[str, str] | None = None
) -> CopyNumArray:
    """Read a CNVkit file (.cnn, .cnr, .cns) to create a CopyNumArray object."""
    return tabio.read(infile, into=CNA, sample_id=sample_id, meta=meta)  # type: ignore[return-value]


def read_ga(
    infile: str, sample_id: str | None = None, meta: dict[str, str] | None = None
) -> GenomicArray:
    """Read a CNVkit file (.cnn, .cnr, .cns) to create a GenomicArray (!) object."""
    return tabio.read(infile, into=GA, sample_id=sample_id, meta=meta)  # type: ignore[return-value]


def load_het_snps(
    vcf_fname: str,
    sample_id: str | None = None,
    normal_id: str | None = None,
    min_variant_depth: int = 20,
    zygosity_freq: float | None = None,
    tumor_boost: bool = False,
) -> VariantArray:
    if vcf_fname is None:
        return None  # type: ignore[return-value,unreachable]
    varr = tabio.read(
        vcf_fname,
        "vcf",
        sample_id=sample_id,
        normal_id=normal_id,
        min_depth=min_variant_depth,
        skip_somatic=True,
    )
    if zygosity_freq is None and "n_zygosity" in varr and not varr["n_zygosity"].any():
        # Mutect2 sets all normal genotypes to 0/0 -- work around it
        logging.warning(
            "VCF normal sample's genotypes are all 0/0 or missing; "
            "inferring genotypes from allele frequency instead"
        )
        zygosity_freq = 0.25
    if zygosity_freq is not None:
        varr = varr.zygosity_from_freq(zygosity_freq, 1 - zygosity_freq)  # type: ignore[attr-defined]
    if "n_zygosity" in varr:
        # Infer & drop (more) somatic loci based on genotype
        somatic_idx = (varr["zygosity"] != 0.0) & (varr["n_zygosity"] == 0.0)
        if somatic_idx.any() and not somatic_idx.all():
            logging.info(
                f"Skipping {somatic_idx.sum()} additional somatic records "
                + "based on T/N genotypes",
            )
        varr = varr[~somatic_idx]
    orig_len = len(varr)
    varr = varr.heterozygous()  # type: ignore[attr-defined]
    logging.info("Kept %d heterozygous of %d VCF records", len(varr), orig_len)
    _warn_if_baf_input_suspicious(
        varr["alt_freq"] if "alt_freq" in varr else None
    )
    # TODO use/explore tumor_boost option
    if tumor_boost:
        varr["alt_freq"] = varr.tumor_boost()
    return varr  # type: ignore[no-any-return]


def _warn_if_baf_input_suspicious(alt_freqs: pd.Series | None) -> None:
    """Log warnings when the heterozygous-SNP set looks unfit for BAF analysis.

    Triggered for two scenarios that commonly produce nonsensical downstream
    BAF and allele-specific copy number values:

    - No heterozygous variants survive filtering (somatic-only VCF, wrong
      sample IDs, or empty input).
    - The median allele frequency is far from 0.5, indicating the variants
      are unlikely to be heterozygous germline SNPs (e.g. somatic variants
      lacking a SOMATIC INFO tag, homozygous-only VCF).
    """
    if alt_freqs is None or len(alt_freqs) == 0:
        logging.warning(
            "No heterozygous variants remain after filtering. BAF and "
            "allele-specific copy number cannot be computed. Check that "
            "the VCF includes germline (non-somatic) heterozygous SNPs and "
            "that the correct sample IDs are given (-i/--sample-id, "
            "-n/--normal-id)."
        )
        return
    if len(alt_freqs) >= _BAF_INPUT_MIN_HET:
        median_af = float(np.nanmedian(alt_freqs))
        if abs(median_af - 0.5) > _BAF_INPUT_MEDIAN_TOL:
            logging.warning(
                "Median allele frequency %.2f is far from the 0.5 expected "
                "for heterozygous germline SNPs. The VCF may contain "
                "non-germline variants, or the sample IDs (-i/--sample-id, "
                "-n/--normal-id) may be incorrect. BAF and allele-specific "
                "copy number output may be unreliable.",
                median_af,
            )


def verify_sample_sex(
    cnarr: CopyNumArray,
    sex_arg: str | None,
    is_haploid_x_reference: bool,
    diploid_parx_genome: str | None,
) -> bool:
    is_sample_female: bool = cnarr.guess_xx(  # type: ignore[assignment]
        is_haploid_x_reference, diploid_parx_genome, verbose=False
    )
    if sex_arg:
        is_sample_female_given = sex_arg.lower() not in ["y", "m", "male"]
        if is_sample_female != is_sample_female_given:
            logging.warning(
                "Sample sex specified as %s but chromosomal X/Y ploidy looks like %s",
                "female" if is_sample_female_given else "male",
                "female" if is_sample_female else "male",
            )
            is_sample_female = is_sample_female_given
    logging.info(
        "Treating sample %s as %s",
        cnarr.sample_id or "",
        "female" if is_sample_female else "male",
    )
    return is_sample_female


def write_tsv(
    outfname: str | None,
    rows: Iterable[Iterable[object]],
    colnames: Iterable[str] | None = None,
) -> None:
    """Write rows, with optional column header, to tabular file."""
    with tabio.safe_write(outfname or sys.stdout) as handle:  # type: ignore[arg-type]
        if colnames:
            header = "\t".join(colnames) + "\n"
            handle.write(header)
        handle.writelines("\t".join(map(str, row)) + "\n" for row in rows)


def write_text(outfname: str | None, text: str, *more_texts: str) -> None:
    """Write one or more strings (blocks of text) to a file."""
    with tabio.safe_write(outfname or sys.stdout) as handle:  # type: ignore[arg-type]
        handle.write(text)
        if more_texts:
            for mtext in more_texts:
                handle.write(mtext)


def write_dataframe(
    outfname: str | None, dframe: pd.DataFrame, header: bool = True
) -> None:
    """Write a pandas.DataFrame to a tabular file."""
    with tabio.safe_write(outfname or sys.stdout) as handle:  # type: ignore[arg-type]
        dframe.to_csv(handle, header=header, index=False, sep="\t", float_format="%.6g")
