"""Variant Call Format (VCF) for SNV loci."""

from __future__ import annotations
import collections
import logging
import os

from itertools import chain

import numpy as np
import pandas as pd
from typing import TYPE_CHECKING, Any

from skgenome._pysam import PYSAM_INSTALL_MSG

if TYPE_CHECKING:
    import pysam
    from collections.abc import Iterator
    from pysam.libcbcf import (
        VariantFile,
        VariantRecord,
        VariantRecordInfo,
        VariantRecordSample,
    )


# Strelka reports per-base counts (tier1, tier2) instead of GT/AD (#943)
_STRELKA_TIER1_FIELD = {"A": "AU", "C": "CU", "G": "GU", "T": "TU"}


def read_vcf(
    infile: str,
    sample_id: str | None = None,
    normal_id: str | None = None,
    min_depth: int | None = None,
    skip_reject: bool = False,
    skip_somatic: bool = False,
) -> pd.DataFrame:
    """Read one tumor-normal pair or unmatched sample from a VCF file.

    By default, return the first tumor-normal pair or unmatched sample in the
    file.  If `sample_id` is a string identifier, return the (paired or single)
    sample  matching that ID.  If `sample_id` is a positive integer, return the
    sample or pair at that index position, counting from 0.

    Samples without a ``GT`` genotype field (e.g. Strelka somatic calls) are
    tolerated: zygosity is then inferred from the alternate-allele frequency.
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            f"pysam is required for reading VCF files. {PYSAM_INSTALL_MSG}"
        ) from None
    if not isinstance(infile, (str, os.PathLike)):
        # pysam requires a path on disk; an open file handle can't be parsed
        raise ValueError(
            f"Must give a VCF filename, not an open file handle: {infile!r}"
        )
    try:
        vcf_reader = pysam.VariantFile(infile)
    except (OSError, ValueError) as exc:
        # e.g. a malformed header (FORMAT column declared but no samples, #680)
        raise ValueError(f"Could not read VCF file {infile!r}: {exc}") from exc
    if vcf_reader.header.samples:
        sid, nid = _choose_samples(vcf_reader, sample_id, normal_id)
        logging.info(
            "Selected test sample %s and control sample %s",
            sid,
            nid if nid else "",
        )
        # NB: in-place
        vcf_reader.subset_samples(list(filter(None, (sid, nid))))
    else:
        logging.warning("VCF file %s has no sample genotypes", infile)
        sid = sample_id  # type: ignore[assignment]
        nid = None

    columns = [
        "chromosome",
        "start",
        "end",
        "ref",
        "alt",
        "somatic",
        "zygosity",
        "depth",
        "alt_count",
    ]
    if nid:
        columns.extend(["n_zygosity", "n_depth", "n_alt_count"])

    rows = _parse_records(vcf_reader, sid, nid, skip_reject)  # type: ignore[arg-type]
    table = pd.DataFrame.from_records(rows, columns=columns)
    table["alt_freq"] = table["alt_count"] / table["depth"]
    if nid:
        table["n_alt_freq"] = table["n_alt_count"] / table["n_depth"]
    # Fill missing genotype/depth/count fields with 0, but leave allele
    # *frequencies* as NaN: a het call with no allele counts has an unknown
    # (not zero) frequency, and coercing it to 0 pins mirrored BAF to 0 over
    # the whole region (#407). NaN is excluded by the np.nanmedian aggregation.
    freq_cols = {"alt_freq", "n_alt_freq"}
    fill_cols = [c for c in table.columns[6:] if c not in freq_cols]
    table = table.fillna({col: 0.0 for col in fill_cols})
    # Filter out records as requested
    cnt_depth = cnt_som = 0
    if min_depth:
        if table["depth"].any():
            dkey = "n_depth" if "n_depth" in table.columns else "depth"
            idx_depth = table[dkey] >= min_depth
            cnt_depth = (~idx_depth).sum()
            table = table[idx_depth]
        else:
            logging.warning("Depth info not available for filtering")
    if skip_somatic:
        idx_som = table["somatic"]
        cnt_som = idx_som.sum()
        table = table[~idx_som]
    logging.info(
        "Loaded %d records; skipped: %d somatic, %d depth",
        len(table),
        cnt_som,
        cnt_depth,
    )
    # return sid, nid, table
    return table


def _choose_samples(
    vcf_reader: VariantFile, sample_id: str | None, normal_id: str | None
) -> tuple[str, str]:
    """Emit the sample IDs of all samples or tumor-normal pairs in the VCF.

    Determine tumor-normal pairs from the PEDIGREE tag(s). If no PEDIGREE tag
    is present, use the specified sample_id and normal_id as the pair, or if
    unspecified, emit all samples as unpaired tumors.
    """
    vcf_samples = list(vcf_reader.header.samples)
    if isinstance(sample_id, int):  # type: ignore[unreachable]
        sample_id = vcf_samples[sample_id]  # type: ignore[unreachable]
    if isinstance(normal_id, int):  # type: ignore[unreachable]
        normal_id = vcf_samples[normal_id]  # type: ignore[unreachable]
    for sid in (sample_id, normal_id):
        if sid and sid not in vcf_samples:
            raise IndexError(f"Specified sample {sid} not in VCF file")
    pairs = None
    peds = list(_parse_pedigrees(vcf_reader))
    if peds:
        # Trust the PEDIGREE tag
        pairs = peds
    elif normal_id:
        # All/any other samples are tumors paired with this normal
        other_ids = [s for s in vcf_samples if s != normal_id]
        if not other_ids:
            raise IndexError(
                f"No other sample in VCF besides the specified normal {normal_id}; "
                + "did you mean to use this as the sample_id instead?"
            )
        pairs = [(oid, normal_id) for oid in other_ids]
    else:
        # All samples are unpaired tumors
        pairs = [(sid, None) for sid in vcf_samples]  # type: ignore[misc]
    if sample_id:
        # Keep only the specified tumor/test sample
        pairs = [(s, n) for s, n in pairs if s == sample_id]
    if not pairs:
        # sample_id refers to a normal/control sample -- salvage it
        pairs = [(sample_id, None)]  # type: ignore[list-item]
    for sid in set(chain(*pairs)) - {None}:
        _confirm_unique(sid, vcf_samples)

    sid, nid = pairs[0]
    if len(pairs) > 1:
        if nid:
            logging.warning(
                "WARNING: VCF file contains multiple tumor-normal "
                "pairs; returning the first pair '%s' / '%s'",
                sid,
                nid,
            )
        else:
            logging.warning(
                "WARNING: VCF file contains multiple samples; "
                "returning the first sample '%s'",
                sid,
            )

    return sid, nid


def _parse_pedigrees(vcf_reader: VariantFile) -> Iterator[tuple[str, str]]:
    """Extract tumor/normal pair sample IDs from the VCF header.

    Return an iterable of (tumor sample ID, normal sample ID).
    """
    meta = collections.defaultdict(list)
    for hr in vcf_reader.header.records:
        if hr.key and hr.key not in ("ALT", "FILTER", "FORMAT", "INFO", "contig"):
            meta[hr.key].append(dict(hr.items()))
    # Prefer the standard tag
    if "PEDIGREE" in meta:
        for tag in meta["PEDIGREE"]:
            if "Derived" in tag:
                sample_id = tag["Derived"]
                normal_id = tag["Original"]
                logging.debug(
                    "Found tumor sample %s and normal sample %s "
                    "in the VCF header PEDIGREE tag",
                    sample_id,
                    normal_id,
                )
                yield sample_id, normal_id
    # GATK Mutect and Mutect2 imply paired tumor & normal IDs
    elif "GATKCommandLine" in meta:
        # GATK 3.0(?) and earlier
        for tag in meta["GATKCommandLine"]:
            if tag.get("ID") == "MuTect":  # any others OK?
                options = dict(
                    kv.split("=", 1)
                    for kv in (tag["CommandLineOptions"].strip('"').split())
                    if "=" in kv
                )
                sample_id = options.get("tumor_sample_name")  # type: ignore[assignment]
                normal_id = options["normal_sample_name"]
                logging.debug(
                    "Found tumor sample %s and normal sample "
                    "%s in the MuTect VCF header",
                    sample_id,
                    normal_id,
                )
                yield sample_id, normal_id  # type: ignore[misc]
    elif "GATKCommandLine.MuTect2" in meta:
        # GATK 3+ metadata is suboptimal.
        # Apparent T/N convention: The samples are just renamed TUMOR and
        # NORMAL, listed in arbitrary order. Metadata is data! We have always
        # been at war with metadata! (#195)
        # Mutect2 can also run in tumor-only mode (safe fallback)
        if len(vcf_reader.header.samples) == 2:
            sample_ids = tuple(vcf_reader.header.samples)
            if sample_ids == ("NORMAL", "TUMOR"):
                yield ("TUMOR", "NORMAL")
            else:
                yield sample_ids  # type: ignore[misc]


def _confirm_unique(sample_id: str, samples: list[str]) -> None:
    occurrences = [s for s in samples if s == sample_id]
    if len(occurrences) != 1:
        raise IndexError(f"Did not find a single sample ID '{sample_id}' in: {samples}")


def _parse_records(
    records: VariantFile, sample_id: str, normal_id: str, skip_reject: bool
) -> Iterator[Any]:
    """Parse VCF records into DataFrame rows.

    Apply filters to skip records with low depth, homozygosity, the REJECT
    flag, or the SOMATIC info field.
    """
    cnt_reject = 0  # For logging
    for record in records:
        if (
            skip_reject
            and record.filter
            and len(record.filter) > 0
            and len(set(record.filter) - {".", "PASS", "KEEP"})
        ):
            cnt_reject += 1
            continue

        if record.samples:
            sample = record.samples[sample_id]
            try:
                depth, zygosity, alt_count = _extract_genotype(sample, record)
                if normal_id:
                    normal = record.samples[normal_id]
                    n_depth, n_zygosity, n_alt_count = _extract_genotype(normal, record)
            except Exception as exc:
                logging.error(
                    "Skipping %s:%d %s @ %s; %s",
                    record.chrom,
                    record.pos,
                    record.ref,
                    sample_id,
                    exc,
                )
                raise
        else:
            # Assume unpaired tumor; take DP, AF from INFO (e.g. LoFreq)
            depth = record.info.get("DP", 0.0)  # type: ignore[assignment,arg-type]
            if "AF" in record.info:
                alt_freq = record.info["AF"]
                alt_count = round(alt_freq * depth)
                # NB: No genotype, so crudely guess from allele frequency
                zygosity = _zygosity_from_freq(alt_freq)
            else:
                alt_count = 0
                zygosity = 0.0

        is_som = "SOMATIC" in record.info and bool(record.info.get("SOMATIC"))
        # Split multiallelics?
        # XXX Ensure sample genotypes are handled properly
        start = record.start
        if record.alts:
            for alt in record.alts:
                if alt == "<NON_REF>":
                    # gVCF placeholder -- not a real allele
                    continue
                end = _get_end(start, alt, record.info)
                row = (
                    record.chrom,
                    start,
                    end,
                    record.ref,
                    alt,
                    is_som,
                    zygosity,
                    depth,
                    alt_count,
                )
                if normal_id:
                    row += (n_zygosity, n_depth, n_alt_count)  # type: ignore[assignment]
                yield row

    if cnt_reject:
        logging.info("Filtered out %d records", cnt_reject)


def _extract_genotype(
    sample: VariantRecordSample, record: VariantRecord
) -> tuple[int | float, float, int | float]:
    depth = _get_depth(sample, record)
    alt_count = _get_alt_count(sample, record)
    zygosity = _get_zygosity(sample, depth, alt_count)
    return depth, zygosity, alt_count


def _get_depth(sample: VariantRecordSample, record: VariantRecord) -> int | float:
    """Get the total read depth at a sample's locus."""
    if "DP" in sample:
        return sample["DP"]  # type: ignore[no-any-return]
    if "AD" in sample and isinstance(sample["AD"], tuple):
        return _safesum(sample["AD"])
    if "DP" in record.info:
        return record.info["DP"]  # type: ignore[no-any-return]
    # SV or not called, probably
    return np.nan


def _get_zygosity(
    sample: VariantRecordSample, depth: int | float, alt_count: int | float
) -> float:
    """Get the zygosity (0=hom-ref, 0.5=het, 1=hom-alt) of a sample's genotype.

    Use the explicit GT call when present and complete; otherwise -- when GT is
    absent (e.g. Strelka, #943) or a no-call (``./.``) -- infer it from the
    alternate-allele frequency.
    """
    if "GT" in sample:
        gt = tuple(sample["GT"])
        # pysam reports a no-call allele as None; it is a missing-allele
        # placeholder, not a distinct allele, so don't count it as one
        # (gh-9vv). Classify by the alleles that were actually called.
        called = set(gt) - {None}
        has_nocall = None in gt
        # A complete no-call ('./.', i.e. nothing called) carries no genotype,
        # so fall through to frequency inference rather than guessing.
        if called:
            if 0 in called:
                # A REF allele was called -- het only if an ALT was too
                # ('0/1'); '0/0' and the partial '0/.' (no ALT evidence)
                # are hom-ref.
                return 0.5 if len(called) > 1 else 0.0
            # Only ALT allele(s) called. '1/1' (and '1/2') are unambiguous,
            # but a partial '1/.' leaves the other allele unknown; real
            # ensemble VCFs show these lean heterozygous, so keep them het.
            if has_nocall:
                return 0.5
            return 0.5 if len(called) > 1 else 1.0
    # No (complete) genotype call -- guess from allele frequency, if we have it
    if (
        depth
        and not np.isnan(depth)
        and alt_count is not None
        and not np.isnan(alt_count)
    ):
        return _zygosity_from_freq(alt_count / depth)
    return 0.0


def _zygosity_from_freq(alt_freq: float) -> float:
    """Crudely guess zygosity from an alternate-allele frequency."""
    if alt_freq < 0.25:
        return 0.0
    if alt_freq < 0.75:
        return 0.5
    return 1.0


def _get_alt_count(sample: VariantRecordSample, record: VariantRecord) -> int | float:
    """Get the alternative allele count from a sample in a VCF record."""
    if sample.get("AD") not in (None, (None,)):
        # GATK and other callers: (ref depth, alt depth)
        if isinstance(sample["AD"], tuple):
            # Ensure we have alternative alleles and thus two AD values
            # ref only calls in GATK can be missing these
            # 1       49515   .       G       .       50.80   .       AN=2;DP=34;MQ=40.01  GT:AD:DP:MMQ  0/0:34:34:.
            alt_count = sample["AD"][1] if len(sample["AD"]) > 1 else 0.0
        # VarScan
        else:
            alt_count = sample["AD"]
    elif sample.get("CLCAD2") not in (None, (None,)):
        # Qiagen CLC Genomics Server -- similar to GATK's AD
        alt_count = sample["CLCAD2"][1]
    elif "AO" in sample:
        if isinstance(sample["AO"], tuple):
            alt_count = _safesum(sample["AO"])
        elif sample["AO"]:
            alt_count = sample["AO"]
        else:
            alt_count = 0.0
    else:
        alt_count = _strelka_alt_count(sample, record)
    return alt_count


def _strelka_alt_count(
    sample: VariantRecordSample, record: VariantRecord
) -> int | float:
    """Get the alt-allele count from Strelka's per-base tier-1 counts.

    Strelka SNV VCFs report AU/CU/GU/TU (counts of A/C/G/T alleles in tiers
    1,2) rather than AD. Return the tier-1 count for the (first) ALT base, or
    NaN if unavailable (e.g. indels, which use a different layout).
    """
    if not record.alts:
        return np.nan
    field = _STRELKA_TIER1_FIELD.get(record.alts[0].upper())
    if field is not None and field in sample:
        value = sample[field]
        # tier-1 count is the first of the (tier1, tier2) pair
        return value[0] if isinstance(value, tuple) else value  # type: ignore[no-any-return]
    return np.nan


def _safesum(tup: tuple[None] | tuple[int, int] | tuple[int]) -> int:
    return sum(filter(None, tup))


def _get_end(posn: int, alt: str, info: VariantRecordInfo) -> int:
    """Get record end position."""
    if "END" in info:
        # Structural variant
        return info["END"]  # type: ignore[no-any-return]
    return posn + len(alt)


# _____________________________________________________________________


def write_vcf(dframe):
    """Variant Call Format (VCF) for SV loci."""
    return NotImplemented
    # See export.export_vcf()
