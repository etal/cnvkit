#!/usr/bin/env python
"""Genotype a sample at known SNP sites, for allele-frequency (BAF) analysis.

Pile up the reads at each site listed in a VCF -- for example, common SNPs from
dbSNP or the 1000 Genomes Project, or a normal sample's germline calls -- and
call genotypes there with bcftools. The output VCF carries each sample's
genotype (GT) and allele depths (AD, DP), which is what CNVkit reads to find
heterozygous germline SNPs and their B-allele frequencies.

This is useful when the available VCF holds only somatic variants, as from
GATK Mutect2 run without germline output, and so has no heterozygous SNPs.

Give the matched normal BAM, if there is one, with -n. Both samples are then
called jointly, and the output's PEDIGREE header pairs them so that CNVkit
judges zygosity from the normal sample, as it does for a tumor-normal VCF from
a variant caller.

Sample names in the output are the BAM filenames without directory or
extension. Chromosome names must match between the BAM, the FASTA and the
sites VCF; sites on chromosomes missing from the BAM or the FASTA are skipped.
"""

import argparse
import logging
import os
import sys
import tempfile
from collections import Counter
from collections.abc import Iterable, Iterator

import pysam
import pysam.bcftools

from skgenome import tabio

from ..core import fbase
from ..samutil import ensure_bam_index, get_bam_chroms

# Sites per bcftools run. pysam returns each run's output as one string, so
# a whole-genome site list is processed in blocks to bound memory use.
# Splitting the list changes no allele depth (AD, DP). At a site deeper than
# 255 reads, though, bcftools computes the genotype likelihoods from a random
# subsample of the bases, drawn from a generator that runs on from site to
# site. PL and QUAL there -- and so the genotype, at a borderline site --
# therefore depend on which sites preceded it in the same run, and so on where
# the blocks begin.
CHUNK_SIZE = 50_000
# Per-BAM read depth cap. bcftools' default of 250 would subsample the reads
# at a deeply sequenced site and so coarsen its allele frequency; this matches
# the default of samtools mpileup.
MAX_DEPTH = 8000


def argument_parsing() -> argparse.Namespace:
    AP = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    AP.add_argument("sample_bam", help="Sample (e.g. tumor) BAM or CRAM file.")
    AP.add_argument(
        "sites",
        help="""VCF of the SNP sites to genotype. Only positions are used, so
                a sites-only VCF without sample columns is fine.""",
    )
    AP.add_argument(
        "-f",
        "--fasta",
        metavar="FILENAME",
        required=True,
        help="""Reference genome, FASTA format: the one the reads were aligned
                to. An index (.fai) is created beside it if absent.""",
    )
    AP.add_argument(
        "-n",
        "--normal",
        metavar="BAM",
        help="Matched normal BAM or CRAM file, to be genotyped jointly.",
    )
    AP.add_argument(
        "-o", "--output", metavar="FILENAME", help="Output VCF. [Default: stdout]"
    )
    return AP.parse_args()


def snpfilter(
    sample_bam: str,
    sites_fname: str,
    fasta: str,
    normal_bam: str | None = None,
    chunk_size: int = CHUNK_SIZE,
) -> Iterator[str]:
    """Genotype a sample, and optionally its matched normal, at the given sites.

    Yield the lines of the output VCF.
    """
    bams = [sample_bam] if normal_bam is None else [sample_bam, normal_bam]
    names = [fbase(bam) for bam in bams]
    if len(set(names)) < len(names):
        raise ValueError(
            f"The sample and normal BAMs would both get the same sample name "
            f"'{names[0]}'; rename one of the files"
        )
    for bam in bams:
        ensure_bam_index(bam)
    bam_contigs = set.intersection(*(get_bam_chroms(bam, fasta) for bam in bams))
    with pysam.FastaFile(fasta) as fa:
        fasta_contigs = set(fa.references)
    contigs = bam_contigs & fasta_contigs
    if not contigs:
        raise ValueError(
            "The BAM and the FASTA have no chromosome name in common "
            f"(BAM: {', '.join(_abbreviate(sorted(bam_contigs)))}; "
            f"FASTA: {', '.join(_abbreviate(sorted(fasta_contigs)))})"
        )
    n_sites = n_records = 0
    with tempfile.TemporaryDirectory() as tmpdir:
        regions_fname = os.path.join(tmpdir, "sites.tsv")
        pileup_fname = os.path.join(tmpdir, "pileup.vcf")
        positions = _read_positions(sites_fname, contigs)
        for i, chunk in enumerate(_chunks(positions, chunk_size)):
            n_sites += len(chunk)
            with open(regions_fname, "w") as handle:
                handle.writelines(f"{chrom}\t{pos}\n" for chrom, pos in chunk)
            # Without split_lines, pysam's dispatcher returns standard output
            # as one string
            pileup: str = pysam.bcftools.mpileup(  # type: ignore[assignment]
                "--fasta-ref", fasta,
                "--regions-file", regions_fname,
                "--annotate", "FORMAT/AD,FORMAT/DP",
                "--max-depth", str(MAX_DEPTH),
                "--skip-indels",
                "--count-orphans",
                "--ignore-RG",
                "--output-type", "v",
                *bams,
            )  # fmt: skip
            with open(pileup_fname, "w") as handle:
                handle.write(pileup)
            # Keep sites where no sample carries an ALT allele (no --variants-
            # only), so that the output covers every site that had reads.
            calls: str = pysam.bcftools.call(  # type: ignore[assignment]
                "--multiallelic-caller", "--output-type", "v", pileup_fname
            )
            for line in calls.splitlines(keepends=True):
                if not line.startswith("#"):
                    n_records += 1
                    yield line
                elif i == 0:
                    if line.startswith("#CHROM"):
                        if normal_bam is not None:
                            yield f"##PEDIGREE=<Derived={names[0]},Original={names[1]}>\n"
                        fields = line.rstrip("\n").split("\t")
                        line = "\t".join(fields[:9] + names) + "\n"
                    yield line
    if n_records:
        logging.info("%d of %d sites had reads to genotype", n_records, n_sites)
    else:
        logging.warning(
            "WARNING: None of the %d sites had any reads; the VCF has no records",
            n_sites,
        )


def _read_positions(sites_fname: str, contigs: set[str]) -> Iterator[tuple[str, int]]:
    """Yield the distinct (chromosome, 1-based position) pairs of the sites.

    Sites on chromosomes outside ``contigs`` are skipped: bcftools fails on a
    region it cannot find in the BAM, and where the FASTA lacks the chromosome
    it piles up against a reference base of N, which the caller then drops
    without notice.
    Consecutive records at one position, such as the split alleles of a
    multiallelic site, give a single position, since bcftools would otherwise
    genotype the site once for each.
    """
    skipped: Counter[str] = Counter()
    n_kept = 0
    last = None
    with pysam.VariantFile(sites_fname, drop_samples=True) as sites:
        for record in sites:
            if record.chrom not in contigs:
                skipped[record.chrom] += 1
                continue
            position = (record.chrom, record.pos)
            if position != last:
                n_kept += 1
                yield position
                last = position
    if skipped:
        logging.info(
            "Skipped %d sites on %d chromosomes not in both the BAM and the FASTA: %s",
            skipped.total(),
            len(skipped),
            ", ".join(_abbreviate(skipped)),
        )
    if not n_kept:
        if not skipped:
            raise ValueError(f"No sites found in {sites_fname}")
        raise ValueError(
            f"None of the sites in {sites_fname} is on a chromosome in both the "
            "BAM and the FASTA. Check that the chromosome names match "
            f"(sites VCF: {', '.join(_abbreviate(skipped))}; "
            f"BAM and FASTA: {', '.join(_abbreviate(sorted(contigs)))})"
        )


def _abbreviate(names: Iterable[str], limit: int = 5) -> list[str]:
    names = list(names)
    return names if len(names) <= limit else [*names[:limit], "..."]


def _chunks(
    positions: Iterable[tuple[str, int]], size: int
) -> Iterator[list[tuple[str, int]]]:
    chunk: list[tuple[str, int]] = []
    for position in positions:
        chunk.append(position)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = argument_parsing()
    if args.output and os.path.exists(args.output):
        # The inputs are still being read while the output is written. An input
        # that does not exist cannot be compared, and pysam reports it instead.
        for input_fname in (args.sites, args.sample_bam, args.normal, args.fasta):
            if (
                input_fname
                and os.path.exists(input_fname)
                and os.path.samefile(args.output, input_fname)
            ):
                raise ValueError(
                    f"Output file {args.output} is also an input file ({input_fname})"
                )
    lines = snpfilter(args.sample_bam, args.sites, args.fasta, args.normal)
    # snpfilter() is a generator: its input checks and its first block of sites
    # run at this next(), before the output file is opened. A failure in a
    # later block truncates the output, as for other CNVkit commands.
    first_line = next(lines)
    with tabio.safe_write(args.output or sys.stdout) as handle:  # type: ignore[arg-type]
        handle.write(first_line)
        handle.writelines(lines)


if __name__ == "__main__":
    main()
