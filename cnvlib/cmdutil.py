"""Functions reused within command-line implementations."""

from __future__ import annotations

import logging
import sys
from typing import TYPE_CHECKING

import numpy as np

from skgenome import GenomicArray as GA
from skgenome import tabio
from skgenome.chromnames import infer_sex_chrom_labels

from .cnary import CopyNumArray as CNA
from .cnary import is_female_default
from .vary import chrx_het_density_rejects_haploid

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


def _resolve_zygosity_freq(
    germline_varr: VariantArray, zygosity_freq: float | None, *, warn: bool
) -> float | None:
    """Apply the Mutect2 all-zero-normal workaround heuristic.

    When the VCF's normal-sample genotypes are all 0/0 or missing (Mutect2
    quirk), infer zygosity from allele frequency by setting ``zygosity_freq``
    to 0.25. Caller must pass the germline subset only (e.g. ``varr`` from a
    ``skip_somatic=True`` read), so that the all-zero check reflects the
    germline rows and is not masked by interleaved somatic calls.
    """
    if (
        zygosity_freq is None
        and "n_zygosity" in germline_varr
        and not germline_varr["n_zygosity"].any()
    ):
        if warn:
            logging.warning(
                "VCF normal sample's genotypes are all 0/0 or missing; "
                "inferring genotypes from allele frequency instead"
            )
        return 0.25
    return zygosity_freq


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
    zygosity_freq = _resolve_zygosity_freq(varr, zygosity_freq, warn=True)  # type: ignore[arg-type]
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
    # Count chrX SNPs (and chrX hets) BEFORE the global het filter, so the
    # downstream chrX-het-density confirmer (verify_sample_sex with VCF, #341)
    # has both numerator and denominator. Detect the chrX label via
    # skgenome's context-aware classifier rather than an inline regex.
    chrx_label, _ = infer_sex_chrom_labels(set(varr.chromosome.unique()))
    if chrx_label is not None:
        on_chrx = varr.chromosome == chrx_label
        n_chrx_total = int(on_chrx.sum())
        # Match VariantArray.heterozygous()'s column-selection: when a matched
        # normal VCF is loaded both ``zygosity`` (tumor) and ``n_zygosity``
        # (normal) are present, and the germline-het signal lives in the
        # normal. Use the same precedence so the binomial test's numerator
        # comes from the same population the het-filter would have selected.
        if "n_zygosity" in varr:
            zyg = varr["n_zygosity"]
        elif "zygosity" in varr:
            zyg = varr["zygosity"]
        else:
            zyg = None
        if zyg is not None:
            n_chrx_het = int((on_chrx & (zyg != 0.0) & (zyg != 1.0)).sum())
        else:
            n_chrx_het = 0
    else:
        n_chrx_total = 0
        n_chrx_het = 0

    orig_len = len(varr)
    varr = varr.heterozygous()  # type: ignore[attr-defined]
    logging.info("Kept %d heterozygous of %d VCF records", len(varr), orig_len)
    _warn_if_baf_input_suspicious(varr["alt_freq"] if "alt_freq" in varr else None)
    # TODO use/explore tumor_boost option
    if tumor_boost:
        varr["alt_freq"] = varr.tumor_boost()
    # Stash the pre-filter chrX counts so verify_sample_sex can run the
    # chrX-het-density confirmer below.
    varr.meta["chrx_snp_total"] = n_chrx_total
    varr.meta["chrx_het_count"] = n_chrx_het
    return varr  # type: ignore[no-any-return]


def _partition_snp_extras(
    varr: VariantArray,
    *,
    include_loh: bool,
    include_somatic: bool,
) -> tuple[VariantArray | None, VariantArray | None]:
    """Split a VariantArray into LOH-evidence and somatic-SNV subsets.

    ``varr`` should be the full record set as returned by ``tabio.read`` with
    ``skip_somatic=False`` (post-zygosity_freq inference if applicable). The
    function does not modify ``varr``.

    LOH evidence = tumor-homozygous loci excluding any somatic call. When a
    matched normal column (``n_zygosity``) is present, additionally restrict
    to normal-heterozygous loci (the stricter, T/N-aware definition). In
    tumor-only flows the column is absent and any tumor-homozygous locus
    qualifies; this is the same convention as the 2017 #290 request, which
    treated tumor-homozygous SNVs as LOH evidence in tumor-only mode.

    Somatic = the union of (a) the VCF ``SOMATIC`` flag and (b) T/N-inferred
    somatic loci (tumor-het with normal-homref), matching the loci that
    ``load_het_snps`` excludes when building the heterozygous subset.
    """
    if not (include_loh or include_somatic):
        return None, None
    vcf_somatic = varr["somatic"]
    if "n_zygosity" in varr:
        tn_somatic = (varr["zygosity"] != 0.0) & (varr["n_zygosity"] == 0.0)
        all_somatic = vcf_somatic | tn_somatic
    else:
        all_somatic = vcf_somatic
    loh: VariantArray | None = None
    if include_loh:
        non_somatic = varr[~all_somatic]
        tumor_hom = (non_somatic["zygosity"] == 0.0) | (non_somatic["zygosity"] == 1.0)
        if "n_zygosity" in non_somatic:
            loh_mask = tumor_hom & (non_somatic["n_zygosity"] == 0.5)
        else:
            loh_mask = tumor_hom
        loh = non_somatic[loh_mask]
    somatic: VariantArray | None = None
    if include_somatic:
        somatic = varr[all_somatic]
    return loh, somatic


def load_snp_subsets(
    vcf_fname: str,
    sample_id: str | None = None,
    normal_id: str | None = None,
    min_variant_depth: int = 20,
    zygosity_freq: float | None = None,
    *,
    include_loh: bool = False,
    include_somatic: bool = False,
) -> tuple[VariantArray, VariantArray | None, VariantArray | None]:
    """Load heterozygous SNPs plus optional LOH-evidence and somatic subsets.

    The ``het`` subset is identical to ``load_het_snps``'s return -- same BAF
    distribution, same chrX-counting in ``.meta``, same warnings. The extras
    are surfaced for display only (e.g. scatter plot ``--show-snvs`` modes,
    #290) and are excluded from segment-overlay calculations to preserve the
    BAF math that the het subset drives.

    When both ``include_*`` flags are False, only ``het`` is loaded -- the
    second VCF read is skipped. The extras are ``None`` in that case.
    """
    het = load_het_snps(
        vcf_fname, sample_id, normal_id, min_variant_depth, zygosity_freq
    )
    if not (include_loh or include_somatic):
        return het, None, None
    # Re-read the VCF without skip_somatic to recover the SOMATIC-flagged
    # loci. The chrX-counting and BAF-distribution warnings already fired
    # during the load_het_snps call, so do not re-emit them here.
    varr_all = tabio.read(
        vcf_fname,
        "vcf",
        sample_id=sample_id,
        normal_id=normal_id,
        min_depth=min_variant_depth,
        skip_somatic=False,
    )
    # Resolve the Mutect2 workaround against the germline subset only --
    # interleaved somatic rows can mask the all-zero-normal heuristic that
    # load_het_snps sees, which would leave germline rows with raw all-zero
    # zygosity here while load_het_snps applied the workaround. Restricting
    # to ~somatic mirrors load_het_snps's skip_somatic=True input.
    germline = varr_all[~varr_all["somatic"]]  # type: ignore[index]
    zygosity_freq = _resolve_zygosity_freq(germline, zygosity_freq, warn=False)  # type: ignore[arg-type]
    if zygosity_freq is not None:
        varr_all = varr_all.zygosity_from_freq(zygosity_freq, 1 - zygosity_freq)  # type: ignore[attr-defined]
    loh, somatic = _partition_snp_extras(
        varr_all,  # type: ignore[arg-type]
        include_loh=include_loh,
        include_somatic=include_somatic,
    )
    return het, loh, somatic


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
    variants: VariantArray | None = None,
) -> bool:
    guess = cnarr.guess_xx(is_haploid_x_reference, diploid_parx_genome, verbose=False)
    # None means "couldn't tell" (no chrX / degenerate); default to female so an
    # indeterminate sample is never silently treated as male (#360 family).
    is_sample_female = is_female_default(guess)
    if sex_arg:
        is_sample_female_given = sex_arg.lower() not in ["y", "m", "male"]
        # Only warn on a real disagreement -- not when we merely defaulted.
        if guess is not None and is_sample_female != is_sample_female_given:
            logging.warning(
                "Sample sex specified as %s but chromosomal X/Y ploidy looks like %s",
                "female" if is_sample_female_given else "male",
                "female" if is_sample_female else "male",
            )
        is_sample_female = is_sample_female_given
    elif (
        variants is not None and not is_sample_female and cnarr.chr_x_label is not None
    ):
        # No user override + coverage inferred male: ask the VCF for an
        # independent chrX-heterozygous-SNP confirmer. True haploid X has
        # ~no het SNPs, so observing many of them rejects the haploid-X
        # null and we override to female -- the chrX expected ploidy used
        # downstream (purity rescaling, diagram shift_xx, etc.) is then
        # diploid, which is the correct baseline for calls and LOH (#341).
        # Coverage-based call stays untouched when the test isn't decisive,
        # so this is a one-way upgrade toward female (the safe direction
        # for #360-family indeterminacy).
        n_total = variants.meta.get("chrx_snp_total", 0)
        n_het = variants.meta.get("chrx_het_count", 0)
        rejected, p_value = chrx_het_density_rejects_haploid(n_total, n_het)
        if rejected:
            logging.warning(
                "Coverage inferred male, but chrX heterozygous-SNP density "
                "(%d het of %d total chrX SNPs in the VCF) rejects the "
                "haploid-X null at binomial p=%.2g; treating sample as "
                "female so chrX copy-number calls use a diploid-X "
                "baseline. Pass --sample-sex male to suppress this "
                "override.",
                n_het,
                n_total,
                p_value,
            )
            is_sample_female = True
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
