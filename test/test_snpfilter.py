#!/usr/bin/env python
"""Tests for snpfilter.py, which genotypes a sample at known SNP sites.

Each test builds a small synthetic genome in ``tmp_path``: a random reference
sequence, a sites VCF, and tumor and normal BAMs whose allele counts at every
site are set exactly. The pipeline runs on the bcftools bundled with pysam, so
no external binary is needed. The output is then read back through
``cmdutil.load_het_snps``, the loader that the BAF-consuming commands use.
"""

import logging
import random
import sys

import numpy as np
import pysam
import pytest

from cnvlib import cmdutil
from cnvlib.cli import snpfilter

READ_LEN = 100
# Above the 250 reads per file to which bcftools mpileup subsamples by default
DEPTH = 300
CONTIGS = {"chr1": 8000, "chr2": 6000}
# (normal ALT fraction, tumor ALT fraction, number of sites)
DESIGN = {
    "het": (0.5, 0.5, 50),
    "imbalanced": (0.5, 0.8, 6),
    "loh": (0.5, 1.0, 6),
    "hom_ref": (0.0, 0.0, 3),
    "hom_alt": (1.0, 1.0, 3),
}


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    """Write the synthetic reference, sites VCF and BAMs; return their paths."""
    out = tmp_path_factory.mktemp("snpfilter")
    rng = random.Random(0)
    seqs = {
        c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CONTIGS.items()
    }
    fasta = out / "ref.fa"
    with open(fasta, "w") as handle:
        for chrom, seq in seqs.items():
            handle.write(f">{chrom}\n{seq}\n")
    pysam.faidx(str(fasta))

    # Sites 200 bp apart, so that no read reaches a neighboring site
    slots = [
        (chrom, pos)
        for chrom, size in CONTIGS.items()
        for pos in range(150, size - 150, 200)
    ]
    kinds = [kind for kind, (_, _, n) in DESIGN.items() for _ in range(n)]
    rng.shuffle(kinds)
    sites = []  # (chrom, 0-based position, ref, alt, kind)
    for (chrom, pos), kind in zip(slots, kinds, strict=True):
        ref = seqs[chrom][pos]
        alt = rng.choice([b for b in "ACGT" if b != ref])
        sites.append((chrom, pos, ref, alt, kind))

    sites_vcf = out / "sites.vcf"
    write_sites(sites_vcf, [(c, p + 1, r, a) for c, p, r, a, _ in sites])
    bams = {
        name: write_bam(out / f"{name}.bam", seqs, sites, column)
        for column, name in enumerate(("normal", "tumor"))
    }
    return {
        "fasta": str(fasta),
        "sites": str(sites_vcf),
        "tumor": bams["tumor"],
        "normal": bams["normal"],
        "dir": out,
        "kinds": {(c, p + 1): kind for c, p, _, _, kind in sites},
    }


def write_sites(path, records):
    """Write a sites-only VCF of (chrom, 1-based position, ref, alt) records."""
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom, size in CONTIGS.items():
            handle.write(f"##contig=<ID={chrom},length={size}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in records:
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")


def write_bam(path, seqs, sites, column):
    """Write a sorted, indexed BAM with exactly the designed allele counts."""
    # Both BAMs name the same sample in their read groups, as a pipeline might
    # for one patient. Without --ignore-RG, bcftools would merge the two files
    # into a single sample column.
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": n} for c, n in CONTIGS.items()],
        "RG": [{"ID": "rg1", "SM": "patient"}],
    }
    tids = {c: i for i, c in enumerate(CONTIGS)}
    reads = []
    for chrom, pos, _ref, alt, kind in sites:
        n_alt = round(DESIGN[kind][column] * DEPTH)
        for i in range(DEPTH):
            start = pos - (i * 37 % READ_LEN)
            bases = list(seqs[chrom][start : start + READ_LEN])
            if i < n_alt:
                bases[pos - start] = alt
            reads.append((tids[chrom], start, "".join(bases)))
    reads.sort()
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, (tid, start, seq) in enumerate(reads):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"r{i}"
            read.reference_id = tid
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = f"{READ_LEN}M"
            read.query_sequence = seq
            read.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            read.set_tag("RG", "rg1")
            bam.write(read)
    pysam.index(str(path))
    return str(path)


def run(path, *args, **kwargs):
    """Run the helper and write its output VCF to ``path``."""
    with open(path, "w") as handle:
        handle.writelines(snpfilter.snpfilter(*args, **kwargs))
    return str(path)


def records(vcf_fname):
    with open(vcf_fname) as handle:
        return [line for line in handle if not line.startswith("#")]


def sample_fields(record, *keys):
    """The named FORMAT fields of each sample in a VCF record line."""
    fields = record.rstrip("\n").split("\t")
    names = fields[8].split(":")
    return [
        tuple(dict(zip(names, sample.split(":"), strict=True))[key] for key in keys)
        for sample in fields[9:]
    ]


def het_kinds(varr, kinds):
    return sorted(
        kinds[(c, s + 1)] for c, s in zip(varr.chromosome, varr.start, strict=True)
    )


def test_tumor_only_round_trip(genome, tmp_path, caplog):
    """Germline hets come back as hets with BAF near 0.5."""
    vcf = run(tmp_path / "calls.vcf", genome["tumor"], genome["sites"], genome["fasta"])
    with open(vcf) as handle:
        sample_line = next(line for line in handle if line.startswith("#CHROM"))
    assert sample_line.rstrip("\n").split("\t")[9:] == ["tumor"]

    # Only load_het_snps's own warnings count, not snpfilter's log records
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        varr = cmdutil.load_het_snps(vcf)
    # Enough hets for load_het_snps to check the median, and it did not object
    assert len(varr) >= cmdutil._BAF_INPUT_MIN_HET
    assert not caplog.records
    assert np.isfinite(varr["alt_freq"]).all()
    assert abs(np.median(varr["alt_freq"]) - 0.5) < 0.05
    # Genotyped alone, the tumor's imbalanced and LOH sites are called
    # homozygous at this depth. Classifying by allele frequency instead
    # recovers the imbalanced ones, as doc/baf.rst advises.
    assert het_kinds(varr, genome["kinds"]) == ["het"] * 50
    varr = cmdutil.load_het_snps(vcf, zygosity_freq=0.1)
    assert het_kinds(varr, genome["kinds"]) == ["het"] * 50 + ["imbalanced"] * 6


def test_paired_round_trip(genome, tmp_path):
    """A matched normal is paired automatically and recovers the LOH sites."""
    vcf = run(
        tmp_path / "calls.vcf",
        genome["tumor"],
        genome["sites"],
        genome["fasta"],
        normal_bam=genome["normal"],
    )
    # Each sample's own depth, not the joint INFO/DP that the reader would
    # otherwise fall back to, and not subsampled
    for record in records(vcf):
        assert sample_fields(record, "DP") == [(str(DEPTH),)] * 2
        assert all(gt != "./." for (gt,) in sample_fields(record, "GT"))
    varr = cmdutil.load_het_snps(vcf)
    assert "n_zygosity" in varr
    assert (varr["n_zygosity"] == 0.5).all()
    assert het_kinds(varr, genome["kinds"]) == (
        ["het"] * 50 + ["imbalanced"] * 6 + ["loh"] * 6
    )
    loh = varr[varr["alt_freq"] > 0.9]
    assert len(loh) == 6


def test_chunking_does_not_change_allele_depths(genome, tmp_path):
    """Splitting the site list into blocks changes no site, genotype or depth.

    Above 255 reads, bcftools subsamples the bases for its genotype likelihoods
    with a generator whose state carries over from site to site, so PL and QUAL
    can differ between the runs; at these unambiguous sites, GT cannot.
    """
    args = (genome["tumor"], genome["sites"], genome["fasta"])
    whole = run(tmp_path / "whole.vcf", *args, normal_bam=genome["normal"])
    split = run(
        tmp_path / "split.vcf", *args, normal_bam=genome["normal"], chunk_size=7
    )

    def contents(vcf_fname):
        with open(vcf_fname) as handle:
            lines = handle.readlines()
        # bcftools' command lines name temporary files
        header = [line for line in lines if line.startswith("#")]
        header = [line for line in header if not line.startswith("##bcftools")]
        calls = [
            (line.split("\t")[:5], sample_fields(line, "GT", "AD", "DP"))
            for line in lines
            if not line.startswith("#")
        ]
        return header, calls

    assert contents(split) == contents(whole)
    assert len(records(whole)) == sum(n for _, _, n in DESIGN.values())


def test_sites_off_the_bam_contigs(genome, tmp_path):
    """Sites on contigs the BAM lacks are skipped; repeated positions collapse."""
    lines = (genome["dir"] / "sites.vcf").read_text().splitlines(keepends=True)
    meta = [line for line in lines if line.startswith("##")]
    header = [line for line in lines if line.startswith("#CHROM")]
    body = [line for line in lines if not line.startswith("#")]
    # A second allele at the first site, as in a split multiallelic record
    first = body[0].split("\t")
    other_alt = next(b for b in "ACGT" if b not in first[3:5])
    second_allele = "\t".join([*first[:4], other_alt, *first[5:]])
    extra = tmp_path / "extra.vcf"
    extra.write_text(
        "".join(
            [
                *meta,
                "##contig=<ID=chrUn,length=1000>\n",
                *header,
                body[0],
                second_allele,
                *body[1:],
                "chrUn\t500\t.\tA\tG\t.\t.\t.\n",
            ]
        )
    )
    plain = run(
        tmp_path / "plain.vcf", genome["tumor"], genome["sites"], genome["fasta"]
    )
    padded = run(tmp_path / "padded.vcf", genome["tumor"], str(extra), genome["fasta"])
    assert records(padded) == records(plain)


def test_no_sites_on_the_bam_contigs(genome, tmp_path):
    """Mismatched chromosome naming is reported, not written as an empty VCF."""
    with open(genome["sites"]) as handle:
        renamed = handle.read().replace("chr", "")
    sites = tmp_path / "nochr.vcf"
    sites.write_text(renamed)
    with pytest.raises(ValueError, match="chromosome names"):
        list(snpfilter.snpfilter(genome["tumor"], str(sites), genome["fasta"]))

    empty = tmp_path / "empty.vcf"
    write_sites(empty, [])
    with pytest.raises(ValueError, match="No sites"):
        list(snpfilter.snpfilter(genome["tumor"], str(empty), genome["fasta"]))


def test_fasta_contigs(genome, tmp_path, caplog):
    """Chromosomes the FASTA lacks are skipped and reported."""
    with pysam.FastaFile(genome["fasta"]) as fa:
        chr1 = fa.fetch("chr1")
    partial = tmp_path / "chr1.fa"
    partial.write_text(f">chr1\n{chr1}\n")
    with caplog.at_level(logging.INFO):
        vcf = run(tmp_path / "chr1.vcf", genome["tumor"], genome["sites"], str(partial))
    n_chr1 = sum(1 for chrom, _pos in genome["kinds"] if chrom == "chr1")
    assert [rec.split("\t")[0] for rec in records(vcf)] == ["chr1"] * n_chr1
    n_chr2 = len(genome["kinds"]) - n_chr1
    assert f"Skipped {n_chr2} sites on 1 chromosomes" in caplog.text

    renamed = tmp_path / "nochr.fa"
    renamed.write_text(f">1\n{chr1}\n")
    with pytest.raises(ValueError, match="no chromosome name in common"):
        list(snpfilter.snpfilter(genome["tumor"], genome["sites"], str(renamed)))


def test_refused_run_leaves_files_alone(genome, tmp_path, monkeypatch):
    """With -o, an input error or an input named as the output changes nothing."""

    # A copy, since a failed guard would overwrite it, and the genome fixture is
    # shared with the other tests in this module
    reference = (genome["dir"] / "ref.fa").read_bytes()
    fasta = tmp_path / "ref.fa"
    fasta.write_bytes(reference)

    def run_main(sites, output):
        argv = ["snpfilter.py", genome["tumor"], str(sites), "-f", str(fasta)]
        monkeypatch.setattr(sys, "argv", [*argv, "-o", str(output)])
        snpfilter.main()

    output = tmp_path / "calls.vcf"
    run_main(genome["sites"], output)
    complete = output.read_text()
    assert len(records(output)) == len(genome["kinds"])

    nochr = tmp_path / "nochr.vcf"
    nochr.write_text((genome["dir"] / "sites.vcf").read_text().replace("chr", ""))
    with pytest.raises(ValueError, match="chromosome names"):
        run_main(nochr, output)
    # The earlier output, spelled differently, is a valid sites VCF that would be
    # truncated while still being read. The ".." keeps the two paths textually
    # distinct, so that os.path.samefile must decide; pathlib removes a ".".
    with pytest.raises(ValueError, match="also an input"):
        run_main(tmp_path / ".." / tmp_path.name / "calls.vcf", output)
    assert output.read_text() == complete

    alias = tmp_path / "alias.fa"
    alias.symlink_to(fasta)
    with pytest.raises(ValueError, match="also an input"):
        run_main(genome["sites"], alias)
    assert fasta.read_bytes() == reference


def test_same_sample_name_refused(genome, tmp_path):
    """Two BAMs that would get the same sample name are refused."""
    normal = tmp_path / "tumor.bam"
    for suffix in ("", ".bai"):
        (tmp_path / f"tumor.bam{suffix}").write_bytes(
            (genome["dir"] / f"normal.bam{suffix}").read_bytes()
        )
    with pytest.raises(ValueError, match="same sample name"):
        list(
            snpfilter.snpfilter(
                genome["tumor"],
                genome["sites"],
                genome["fasta"],
                normal_bam=str(normal),
            )
        )
