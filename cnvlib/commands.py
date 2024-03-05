"""Command-line interface and corresponding API for CNVkit."""
# NB: argparse CLI definitions and API functions are interwoven:
#   "_cmd_*" handles I/O and arguments processing for the command
#   "do_*" runs the command's functionality as an API
import argparse
import logging
import multiprocessing
import os
import sys
import warnings
from typing import List

# Filter spurious Cython warnings re: numpy
# https://github.com/numpy/numpy/pull/432
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# Choose a safe plotting backend supported by the current platform
import matplotlib

if not os.environ.get("DISPLAY"):
    # If running headless, use a GUI-less backend
    matplotlib.use("Agg")
elif sys.platform == "darwin":
    # Prevent crash on OS X
    # https://github.com/MTG/sms-tools/issues/36
    matplotlib.use("TkAgg")

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

pyplot.ioff()

import pandas as pd
from skgenome import tabio, GenomicArray as _GA
from skgenome.rangelabel import to_label

from . import (
    access,
    antitarget,
    autobin,
    batch,
    bintest,
    call,
    core,
    coverage,
    diagram,
    export,
    fix,
    heatmap,
    import_rna,
    importers,
    metrics,
    parallel,
    reference,
    reports,
    scatter,
    segmentation,
    segmetrics,
    target,
)
from .cmdutil import (
    load_het_snps,
    read_cna,
    verify_sample_sex,
    write_tsv,
    write_text,
    write_dataframe,
)

from ._version import __version__


__all__ = []  # type: List[str]


def public(fn):
    __all__.append(fn.__name__)
    return fn


AP = argparse.ArgumentParser(
    description="CNVkit, a command-line toolkit for copy number analysis.",
    epilog="See the online manual for details: https://cnvkit.readthedocs.io",
)
# Print help and exit if run without arguments
AP.set_defaults(func=lambda args: AP.print_help())
AP_subparsers = AP.add_subparsers(help="Sub-commands (use with -h for more info)")


# _____________________________________________________________________________
# Shared parameters
def add_diploid_parx_genome(P):
    P.add_argument(
        "--diploid-parx-genome",
        type=str,
        help="Considers the given human genome's PAR of chromosome X as autosomal. Example: 'grch38'",
    )

# _____________________________________________________________________________
# Core pipeline

# batch -----------------------------------------------------------------------


def _cmd_batch(args):
    """Run the complete CNVkit pipeline on one or more BAM files."""
    logging.info("CNVkit %s", __version__)
    # Validate/restrict options, beyond what argparse mutual exclusion can do
    bad_args_msg = ""
    if args.reference:
        bad_flags = [
            flag
            for is_used, flag in (
                (args.normal is not None, "-n/--normal"),
                (args.fasta, "-f/--fasta"),
                (args.targets, "-t/--targets"),
                (args.antitargets, "-a/--antitargets"),
                (args.access, "-g/--access"),
                (args.annotate, "--annotate"),
                (args.short_names, "--short-names"),
                (args.target_avg_size, "--target-avg-size"),
                (args.antitarget_avg_size, "--antitarget-avg-size"),
                (args.antitarget_min_size, "--antitarget-min-size"),
            )
            if is_used
        ]
        if bad_flags:
            bad_args_msg = (
                "If -r/--reference is given, options to construct "
                "a new reference (%s) should not be used." % ", ".join(bad_flags)
            )
    elif args.normal is None:
        bad_args_msg = (
            "Option -n/--normal must be given to build a new "
            "reference if -r/--reference is not used."
        )
    elif args.seq_method in ("hybrid", "amplicon") and not args.targets:
        bad_args_msg = (
            "For the '%r' sequencing method, option -t/--targets "
            "(at least) must be given to build a new reference if "
            "-r/--reference is not used." % args.seq_method
        )
    if bad_args_msg:
        sys.exit(bad_args_msg + "\n(See: cnvkit.py batch -h)")

    # Ensure sample IDs are unique to avoid overwriting outputs
    seen_sids = {}
    for fname in (args.bam_files or []) + (args.normal or []):
        sid = core.fbase(fname)
        if sid in seen_sids:
            sys.exit(f"Duplicate sample ID {sid!r} (from {fname} and {seen_sids[sid]})")
        seen_sids[sid] = fname

    if args.processes < 1:
        args.processes = multiprocessing.cpu_count()

    if not args.reference:
        # Build a copy number reference; update (anti)targets upon request
        args.reference, args.targets, args.antitargets = batch.batch_make_reference(
            args.normal,
            args.targets,
            args.antitargets,
            args.male_reference,
            args.diploid_parx_genome,
            args.fasta,
            args.annotate,
            args.short_names,
            args.target_avg_size,
            args.access,
            args.antitarget_avg_size,
            args.antitarget_min_size,
            args.output_reference,
            args.output_dir,
            args.processes,
            args.count_reads,
            args.seq_method,
            args.cluster,
        )
    elif args.targets is None and args.antitargets is None:
        # Extract (anti)target BEDs from the given, existing CN reference
        ref_arr = read_cna(args.reference)
        targets, antitargets = reference.reference2regions(ref_arr)
        ref_pfx = os.path.join(args.output_dir, core.fbase(args.reference))
        args.targets = ref_pfx + ".target-tmp.bed"
        args.antitargets = ref_pfx + ".antitarget-tmp.bed"
        tabio.write(targets, args.targets, "bed4")
        tabio.write(antitargets, args.antitargets, "bed4")

    if args.bam_files:
        if args.processes == 1:
            procs_per_bam = 1
            logging.info("Running %d samples in serial", len(args.bam_files))
        else:
            procs_per_bam = max(1, args.processes // len(args.bam_files))
            logging.info(
                "Running %d samples in %d processes (that's %d processes per bam)",
                len(args.bam_files),
                args.processes,
                procs_per_bam,
            )

        with parallel.pick_pool(args.processes) as pool:
            for bam in args.bam_files:
                pool.submit(
                    batch.batch_run_sample,
                    bam,
                    args.targets,
                    args.antitargets,
                    args.reference,
                    args.output_dir,
                    args.male_reference,
                    args.diploid_parx_genome,
                    args.scatter,
                    args.diagram,
                    args.rscript_path,
                    args.count_reads,
                    args.drop_low_coverage,
                    args.seq_method,
                    args.segment_method,
                    procs_per_bam,
                    args.cluster,
                    args.fasta,
                    args.diploid_parx_genome,
                )
    else:
        logging.info(
            "No tumor/test samples (but %d normal/control samples) "
            "specified on the command line.",
            len(args.normal),
        )


P_batch = AP_subparsers.add_parser("batch", help=_cmd_batch.__doc__)
P_batch.add_argument("bam_files", nargs="*", help="Mapped sequence reads (.bam)")
P_batch.add_argument(
    "-m",
    "--seq-method",
    "--method",
    choices=("hybrid", "amplicon", "wgs"),
    default="hybrid",
    help="""Sequencing assay type: hybridization capture ('hybrid'), targeted amplicon
            sequencing ('amplicon'), or whole genome sequencing ('wgs'). Determines
            whether and how to use antitarget bins. [Default: %(default)s]""",
)
P_batch.add_argument(
    "--segment-method",
    choices=segmentation.SEGMENT_METHODS,
    default="cbs",
    help="""Method used in the 'segment' step. [Default: %(default)s]""",
)
P_batch.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Use or assume a male reference (i.e. female samples will have +1 log-CNR of
            chrX; otherwise male samples would have -1 chrX).""",
)
P_batch.add_argument(
    "-c",
    "--count-reads",
    action="store_true",
    help="""Get read depths by counting read midpoints within each bin. (An alternative
            algorithm).""",
)
P_batch.add_argument(
    "--drop-low-coverage",
    action="store_true",
    help="""Drop very-low-coverage bins before segmentation to avoid false-positive
            deletions in poor-quality tumor samples.""",
)
P_batch.add_argument(
    "-p",
    "--processes",
    nargs="?",
    type=int,
    const=0,
    default=1,
    help="""Number of subprocesses used to running each of the BAM files in parallel.
            Without an argument, use the maximum number of available CPUs. [Default:
            process each BAM in serial]""",
)
P_batch.add_argument(
    "--rscript-path",
    metavar="PATH",
    default="Rscript",
    help="""Path to the Rscript executable to use for running R code. Use this option
            to specify a non-default R installation. [Default: %(default)s]""",
)
add_diploid_parx_genome(P_batch)


# Reference-building options
P_batch_newref = P_batch.add_argument_group("To construct a new copy number reference")
P_batch_newref.add_argument(
    "-n",
    "--normal",
    nargs="*",
    metavar="FILES",
    help="""Normal samples (.bam) used to construct the pooled, paired, or flat
            reference. If this option is used but no filenames are given, a "flat"
            reference will be built. Otherwise, all filenames following this option will
            be used.""",
)
P_batch_newref.add_argument(
    "-f",
    "--fasta",
    metavar="FILENAME",
    help="Reference genome, FASTA format (e.g. UCSC hg19.fa)",
)
P_batch_newref.add_argument(
    "-t", "--targets", metavar="FILENAME", help="Target intervals (.bed or .list)"
)
P_batch_newref.add_argument(
    "-a",
    "--antitargets",
    metavar="FILENAME",
    help="Antitarget intervals (.bed or .list)",
)
# For pre-processing targets
P_batch_newref.add_argument(
    "--annotate",
    metavar="FILENAME",
    help="""Use gene models from this file to assign names to the target regions.
            Format: UCSC refFlat.txt or ensFlat.txt file (preferred), or BED, interval
            list, GFF, or similar.""",
)
P_batch_newref.add_argument(
    "--short-names",
    action="store_true",
    help="Reduce multi-accession bait labels to be short and consistent.",
)
P_batch_newref.add_argument(
    "--target-avg-size",
    type=int,
    help="Average size of split target bins (results are approximate).",
)
# For antitargets:
P_batch_newref.add_argument(
    "-g",
    "--access",
    metavar="FILENAME",
    help="""Regions of accessible sequence on chromosomes (.bed), as output by the
            'access' command.""",
)
P_batch_newref.add_argument(
    "--antitarget-avg-size",
    type=int,
    help="Average size of antitarget bins (results are approximate).",
)
P_batch_newref.add_argument(
    "--antitarget-min-size",
    type=int,
    help="Minimum size of antitarget bins (smaller regions are dropped).",
)
P_batch_newref.add_argument(
    "--output-reference",
    metavar="FILENAME",
    help="""Output filename/path for the new reference file being created. (If given,
            ignores the -o/--output-dir option and will write the file to the given
            path. Otherwise, \"reference.cnn\" will be created in the current directory
            or specified output directory.)""",
)
P_batch_newref.add_argument(
    "--cluster",
    action="store_true",
    help="""Calculate and use cluster-specific summary stats in the reference pool to
            normalize samples.""",
)

P_batch_oldref = P_batch.add_argument_group("To reuse an existing reference")
P_batch_oldref.add_argument(
    "-r", "--reference", help="Copy number reference file (.cnn)."  # required=True,
)

# Reporting options
P_batch_report = P_batch.add_argument_group("Output options")
P_batch_report.add_argument(
    "-d", "--output-dir", metavar="DIRECTORY", default=".", help="Output directory."
)
P_batch_report.add_argument(
    "--scatter",
    action="store_true",
    help="Create a whole-genome copy ratio profile as a PDF scatter plot.",
)
P_batch_report.add_argument(
    "--diagram",
    action="store_true",
    help="Create an ideogram of copy ratios on chromosomes as a PDF.",
)

P_batch.set_defaults(func=_cmd_batch)


# target ----------------------------------------------------------------------

do_target = public(target.do_target)


def _cmd_target(args):
    """Transform bait intervals into targets more suitable for CNVkit."""
    regions = tabio.read_auto(args.interval)
    regions = target.do_target(
        regions, args.annotate, args.short_names, args.split, args.avg_size
    )
    tabio.write(regions, args.output, "bed4")


P_target = AP_subparsers.add_parser("target", help=_cmd_target.__doc__)
P_target.add_argument(
    "interval", help="""BED or interval file listing the targeted regions."""
)
P_target.add_argument(
    "--annotate",
    help="""Use gene models from this file to assign names to the target regions.
            Format: UCSC refFlat.txt or ensFlat.txt file (preferred), or BED, interval
            list, GFF, or similar.""",
)
P_target.add_argument(
    "--short-names",
    action="store_true",
    help="Reduce multi-accession bait labels to be short and consistent.",
)
P_target.add_argument(
    "--split",
    action="store_true",
    help="Split large tiled intervals into smaller, consecutive targets.",
)
# Exons: [114--188==203==292--21750], mean=353 -> outlier=359, extreme=515
#   NV2:  [65--181==190==239--12630], mean=264 -> outlier=277, extreme=364
# Default avg_size chosen s.t. minimum bin size after split is ~= median
P_target.add_argument(
    "-a",
    "--avg-size",
    type=int,
    default=200 / 0.75,
    help="""Average size of split target bins (results are approximate).
            [Default: %(default)s]""",
)
P_target.add_argument(
    "-o", "--output", metavar="FILENAME", help="""Output file name."""
)
P_target.set_defaults(func=_cmd_target)


# access ----------------------------------------------------------------------

do_access = public(access.do_access)


def _cmd_access(args):
    """List the locations of accessible sequence regions in a FASTA file."""
    access_arr = access.do_access(args.fa_fname, args.exclude, args.min_gap_size)
    tabio.write(access_arr, args.output, "bed3")


P_access = AP_subparsers.add_parser("access", help=_cmd_access.__doc__)
P_access.add_argument("fa_fname", help="Genome FASTA file name")
P_access.add_argument(
    "-s",
    "--min-gap-size",
    type=int,
    default=5000,
    help="""Minimum gap size between accessible sequence regions. Regions separated by
            less than this distance will be joined together. [Default: %(default)s]""",
)
P_access.add_argument(
    "-x",
    "--exclude",
    action="append",
    default=[],
    help="""Additional regions to exclude, in BED format. Can be used multiple
            times.""",
)
P_access.add_argument(
    "-o",
    "--output",
    metavar="FILENAME",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="Output file name",
)
P_access.set_defaults(func=_cmd_access)


# antitarget ------------------------------------------------------------------

do_antitarget = public(antitarget.do_antitarget)


def _cmd_antitarget(args):
    """Derive off-target ("antitarget") bins from target regions."""
    targets = tabio.read_auto(args.targets)
    access = tabio.read_auto(args.access) if args.access else None
    out_arr = antitarget.do_antitarget(targets, access, args.avg_size, args.min_size)
    if not args.output:
        base, ext = args.interval.rsplit(".", 1)
        args.output = base + ".antitarget." + ext
    tabio.write(out_arr, args.output, "bed4")


P_anti = AP_subparsers.add_parser("antitarget", help=_cmd_antitarget.__doc__)
P_anti.add_argument(
    "targets", help="""BED or interval file listing the targeted regions."""
)
P_anti.add_argument(
    "-g",
    "--access",
    metavar="FILENAME",
    help="""Regions of accessible sequence on chromosomes (.bed), as output by
            genome2access.py.""",
)
P_anti.add_argument(
    "-a",
    "--avg-size",
    type=int,
    default=150000,
    help="""Average size of antitarget bins (results are approximate).
            [Default: %(default)s]""",
)
P_anti.add_argument(
    "-m",
    "--min-size",
    type=int,
    help="""Minimum size of antitarget bins (smaller regions are dropped).
            [Default: 1/16 avg size, calculated]""",
)
P_anti.add_argument("-o", "--output", metavar="FILENAME", help="""Output file name.""")
P_anti.set_defaults(func=_cmd_antitarget)


# autobin ---------------------------------------------------------------------

do_autobin = public(autobin.do_autobin)


def _cmd_autobin(args):
    """Quickly calculate reasonable bin sizes from BAM read counts."""
    if args.method in ("hybrid", "amplicon") and not args.targets:
        raise RuntimeError(f"Sequencing method {args.method!r} requires targets (-t)")
    if args.method == "wgs":
        if not args.access:
            raise RuntimeError(
                "Sequencing method 'wgs' requires accessible regions (-g)"
            )
        if args.targets:
            logging.warning("Targets will be ignored: %s", args.targets)
    if args.method == "amplicon" and args.access:
        logging.warning(
            "Sequencing-accessible regions will be ignored: %s", args.access
        )

    def read_regions(bed_fname):
        if bed_fname:
            regions = tabio.read_auto(bed_fname)
            if len(regions):
                return regions
            logging.warning(
                "No regions to estimate depth from %s",
                regions.meta.get("filename", ""),
            )
        return None

    tgt_arr = read_regions(args.targets)
    access_arr = read_regions(args.access)
    bam_fname = autobin.midsize_file(args.bams)
    fields = autobin.do_autobin(
        bam_fname,
        args.method,
        tgt_arr,
        access_arr,
        args.bp_per_bin,
        args.target_min_size,
        args.target_max_size,
        args.antitarget_min_size,
        args.antitarget_max_size,
        args.fasta,
    )
    (_tgt_depth, tgt_bin_size), (_anti_depth, anti_bin_size) = fields

    # Create & write BED files
    target_out_arr = target.do_target(
        access_arr if args.method == "wgs" else tgt_arr,
        args.annotate,
        args.short_names,
        do_split=True,
        avg_size=tgt_bin_size,
    )
    tgt_name_base = tgt_arr.sample_id if tgt_arr else core.fbase(bam_fname)
    target_bed = args.target_output_bed or tgt_name_base + ".target.bed"
    tabio.write(target_out_arr, target_bed, "bed4")
    if args.method == "hybrid" and anti_bin_size:
        # Build antitarget BED from the given targets
        anti_arr = antitarget.do_antitarget(
            target_out_arr,
            access=access_arr,
            avg_bin_size=anti_bin_size,
            min_bin_size=args.antitarget_min_size,
        )
    else:
        # No antitargets for wgs, amplicon
        anti_arr = _GA([])
    antitarget_bed = args.antitarget_output_bed or tgt_name_base + ".antitarget.bed"
    tabio.write(anti_arr, antitarget_bed, "bed4")

    # Print depths & bin sizes as a table on stdout
    labels = ("Target", "Antitarget")
    width = max(map(len, labels)) + 1
    print(" " * width, "Depth", "Bin size", sep="\t")
    for label, (depth, binsize) in zip(labels, fields):
        if depth is not None:
            print((label + ":").ljust(width), format(depth, ".3f"), binsize, sep="\t")


P_autobin = AP_subparsers.add_parser("autobin", help=_cmd_autobin.__doc__)
P_autobin.add_argument(
    "bams", nargs="+", help="""Sample BAM file(s) to test for target coverage"""
)
P_autobin.add_argument(
    "-f",
    "--fasta",
    metavar="FILENAME",
    help="Reference genome, FASTA format (e.g. UCSC hg19.fa)",
)
P_autobin.add_argument(
    "-m",
    "--method",
    choices=("hybrid", "amplicon", "wgs"),
    default="hybrid",
    help="""Sequencing protocol: hybridization capture ('hybrid'), targeted amplicon
            sequencing ('amplicon'), or whole genome sequencing ('wgs'). Determines
            whether and how to use antitarget bins. [Default: %(default)s]""",
)
P_autobin.add_argument(
    "-g",
    "--access",
    metavar="FILENAME",
    help="""Sequencing-accessible genomic regions, or exons to use as possible targets
            (e.g. output of refFlat2bed.py)""",
)
P_autobin.add_argument(
    "-t",
    "--targets",
    help="""Potentially targeted genomic regions, e.g. all possible exons for the
            reference genome. Format: BED, interval list, etc.""",
)
P_autobin.add_argument(
    "-b",
    "--bp-per-bin",
    type=float,
    default=100000.0,
    help="""Desired average number of sequencing read bases mapped to each bin.
            [Default: %(default)s]""",
)

P_autobin.add_argument(
    "--target-max-size",
    metavar="BASES",
    type=int,
    default=20000,
    help="Maximum size of target bins. [Default: %(default)s]",
)
P_autobin.add_argument(
    "--target-min-size",
    metavar="BASES",
    type=int,
    default=20,
    help="Minimum size of target bins. [Default: %(default)s]",
)
P_autobin.add_argument(
    "--antitarget-max-size",
    metavar="BASES",
    type=int,
    default=500000,
    help="Maximum size of antitarget bins. [Default: %(default)s]",
)
P_autobin.add_argument(
    "--antitarget-min-size",
    metavar="BASES",
    type=int,
    default=500,
    help="Minimum size of antitarget bins. [Default: %(default)s]",
)

P_autobin.add_argument(
    "--annotate",
    metavar="FILENAME",
    help="""Use gene models from this file to assign names to the target regions.
            Format: UCSC refFlat.txt or ensFlat.txt file (preferred), or BED, interval
            list, GFF, or similar.""",
)
P_autobin.add_argument(
    "--short-names",
    action="store_true",
    help="Reduce multi-accession bait labels to be short and consistent.",
)
P_autobin.add_argument(
    "--target-output-bed",
    metavar="FILENAME",
    help="Filename for target BED output. If not specified, constructed from the input file basename.",
)
P_autobin.add_argument(
    "--antitarget-output-bed",
    metavar="FILENAME",
    help="Filename for antitarget BED output. If not specified, constructed from the input file basename.",
)
# Option: --dry-run to not write BED files?
P_autobin.set_defaults(func=_cmd_autobin)


# coverage --------------------------------------------------------------------

do_coverage = public(coverage.do_coverage)


def _cmd_coverage(args):
    """Calculate coverage in the given regions from BAM read depths."""
    pset = coverage.do_coverage(
        args.interval,
        args.bam_file,
        args.count,
        args.min_mapq,
        args.processes,
        args.fasta,
    )
    if not args.output:
        # Create an informative but unique name for the coverage output file
        bambase = core.fbase(args.bam_file)
        bedbase = core.fbase(args.interval)
        tgtbase = (
            "antitargetcoverage" if "anti" in bedbase.lower() else "targetcoverage"
        )
        args.output = f"{bambase}.{tgtbase}.cnn"
        if os.path.exists(args.output):
            args.output = f"{bambase}.{bedbase}.cnn"
    core.ensure_path(args.output)
    tabio.write(pset, args.output)


P_coverage = AP_subparsers.add_parser("coverage", help=_cmd_coverage.__doc__)
P_coverage.add_argument("bam_file", help="Mapped sequence reads (.bam)")
P_coverage.add_argument("interval", help="Intervals (.bed or .list)")
P_coverage.add_argument(
    "-f",
    "--fasta",
    metavar="FILENAME",
    help="Reference genome, FASTA format (e.g. UCSC hg19.fa)",
)
P_coverage.add_argument(
    "-c",
    "--count",
    action="store_true",
    help="""Get read depths by counting read midpoints within each bin.
            (An alternative algorithm).""",
)
P_coverage.add_argument(
    "-q",
    "--min-mapq",
    type=int,
    default=0,
    help="""Minimum mapping quality score (phred scale 0-60) to count a read for
            coverage depth. [Default: %(default)s]""",
)
P_coverage.add_argument(
    "-o", "--output", metavar="FILENAME", help="""Output file name."""
)
P_coverage.add_argument(
    "-p",
    "--processes",
    nargs="?",
    type=int,
    const=0,
    default=1,
    help="""Number of subprocesses to calculate coverage in parallel.
            Without an argument, use the maximum number of available CPUs.
            [Default: use 1 process]""",
)
P_coverage.set_defaults(func=_cmd_coverage)


# reference -------------------------------------------------------------------

do_reference = public(reference.do_reference)
do_reference_flat = public(reference.do_reference_flat)


def _cmd_reference(args):
    """Compile a coverage reference from the given files (normal samples)."""
    usage_err_msg = "Give .cnn samples OR targets and (optionally) antitargets."
    if args.targets:
        # Flat reference
        assert not args.references, usage_err_msg
        ref_probes = reference.do_reference_flat(
            args.targets, args.antitargets, args.fasta, args.male_reference
        )
    elif args.references:
        # Pooled reference
        assert not args.targets and not args.antitargets, usage_err_msg
        filenames = []
        for path in args.references:
            if os.path.isdir(path):
                filenames.extend(
                    os.path.join(path, f)
                    for f in os.listdir(path)
                    if f.endswith("targetcoverage.cnn")
                )
            else:
                filenames.append(path)
        targets = [f for f in filenames if "antitarget" not in f]
        antitargets = [f for f in filenames if "antitarget" in f]
        logging.info(
            "Number of target and antitarget files: %d, %d",
            len(targets),
            len(antitargets),
        )
        female_samples = (
            (args.sample_sex.lower() not in ["y", "m", "male"])
            if args.sample_sex
            else None
        )
        ref_probes = reference.do_reference(
            targets,
            antitargets,
            args.fasta,
            args.male_reference,
            args.diploid_parx_genome,
            female_samples,
            args.do_gc,
            args.do_edge,
            args.do_rmask,
            args.cluster,
            args.min_cluster_size,
        )
    else:
        raise ValueError(usage_err_msg)

    ref_fname = args.output or "cnv_reference.cnn"
    core.ensure_path(ref_fname)
    tabio.write(ref_probes, ref_fname)


P_reference = AP_subparsers.add_parser("reference", help=_cmd_reference.__doc__)
P_reference.add_argument(
    "references",
    nargs="*",
    help="""Normal-sample target or antitarget .cnn files, or the directory that
            contains them.""",
)
P_reference.add_argument(
    "-f", "--fasta", help="Reference genome, FASTA format (e.g. UCSC hg19.fa)"
)
P_reference.add_argument("-o", "--output", metavar="FILENAME", help="Output file name.")
P_reference.add_argument(
    "-c",
    "--cluster",
    action="store_true",
    help="""Calculate and store summary stats for clustered subsets of the normal
            samples with similar coverage profiles.""",
)
P_reference.add_argument(
    "--min-cluster-size",
    metavar="NUM",
    type=int,
    default=4,
    help="""Minimum cluster size to keep in reference profiles. [Default: %(default)s]""",
)
P_reference.add_argument(
    "-x",
    "--sample-sex",
    "-g",
    "--gender",
    dest="sample_sex",
    choices=("m", "y", "male", "Male", "f", "x", "female", "Female"),
    help="""Specify the chromosomal sex of all given samples as male or female.
            (Default: guess each sample from coverage of X and Y chromosomes).""",
)
P_reference.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Create a male reference: shift female samples' chrX log-coverage by -1, so
            the reference chrX average is -1. Otherwise, shift male samples' chrX by +1,
            so the reference chrX average is 0.""",
)
add_diploid_parx_genome(P_reference)

P_reference_flat = P_reference.add_argument_group(
    'To construct a generic, "flat" copy number reference with neutral '
    "expected coverage"
)
P_reference_flat.add_argument(
    "-t", "--targets", help="Target intervals (.bed or .list)"
)
P_reference_flat.add_argument(
    "-a", "--antitargets", help="Antitarget intervals (.bed or .list)"
)

P_reference_bias = P_reference.add_argument_group(
    "To disable specific automatic bias corrections"
)
P_reference_bias.add_argument(
    "--no-gc", dest="do_gc", action="store_false", help="Skip GC correction."
)
P_reference_bias.add_argument(
    "--no-edge",
    dest="do_edge",
    action="store_false",
    help="Skip edge-effect correction.",
)
P_reference_bias.add_argument(
    "--no-rmask",
    dest="do_rmask",
    action="store_false",
    help="Skip RepeatMasker correction.",
)
P_reference.set_defaults(func=_cmd_reference)


# fix -------------------------------------------------------------------------

do_fix = public(fix.do_fix)


def _cmd_fix(args):
    """Combine target and antitarget coverages and correct for biases.

    Adjust raw coverage data according to the given reference, correct potential
    biases and re-center.
    """
    # Verify that target and antitarget are from the same sample
    tgt_raw = read_cna(args.target, sample_id=args.sample_id)
    anti_raw = read_cna(args.antitarget, sample_id=args.sample_id)
    if len(anti_raw) and tgt_raw.sample_id != anti_raw.sample_id:
        raise ValueError(
            "Sample IDs do not match: "
            + f"'{tgt_raw.sample_id}' (target) vs. '{anti_raw.sample_id}' (antitarget)"
        )
    target_table = fix.do_fix(
        tgt_raw,
        anti_raw,
        read_cna(args.reference),
        args.diploid_parx_genome,
        args.do_gc,
        args.do_edge,
        args.do_rmask,
        args.cluster,
        args.smoothing_window_fraction,
    )
    tabio.write(target_table, args.output or tgt_raw.sample_id + ".cnr")


P_fix = AP_subparsers.add_parser("fix", help=_cmd_fix.__doc__)
P_fix.add_argument("target", help="Target coverage file (.targetcoverage.cnn).")
P_fix.add_argument(
    "antitarget", help="Antitarget coverage file (.antitargetcoverage.cnn)."
)
P_fix.add_argument("reference", help="Reference coverage (.cnn).")
P_fix.add_argument(
    "-c",
    "--cluster",
    action="store_true",
    help="""Compare and use cluster-specific values present in the reference profile.
            (Requires that the reference profile was built with the --cluster
            option.)""",
)
P_fix.add_argument(
    "-i",
    "--sample-id",
    help="Sample ID for target/antitarget files. Otherwise inferred from file names.",
)
# P_fix.add_argument('--do-gc', action='store_true', default=True,
#         help="Do GC correction.")
# P_fix.add_argument('--do-edge', action='store_true',
#         help="Do edge-effect correction.")
# P_fix.add_argument('--do-size', action='store_true',
#         help="Do interval-size correction.")
P_fix.add_argument(
    "--no-gc", dest="do_gc", action="store_false", help="Skip GC correction."
)
P_fix.add_argument(
    "--no-edge",
    dest="do_edge",
    action="store_false",
    help="Skip edge-effect correction.",
)
P_fix.add_argument(
    "--no-rmask",
    dest="do_rmask",
    action="store_false",
    help="Skip RepeatMasker correction.",
)
P_fix.add_argument("-o", "--output", metavar="FILENAME", help="Output file name.")
add_diploid_parx_genome(P_fix)
P_fix.add_argument(
    "--smoothing-window-fraction",
    type=float,
    help=(
        "If specified, sets the smoothing window fraction for rolling median bias smoothing"
        " based on traits. Otherwise, defaults to 1/sqrt(len(data))."
    ),
    default=None
)
P_fix.set_defaults(func=_cmd_fix)


# segment ---------------------------------------------------------------------

do_segmentation = public(segmentation.do_segmentation)


def _cmd_segment(args):
    """Infer copy number segments from the given coverage table."""
    cnarr = read_cna(args.filename)
    variants = load_het_snps(
        args.vcf,
        args.sample_id,
        args.normal_id,
        args.min_variant_depth,
        args.zygosity_freq,
    )
    results = segmentation.do_segmentation(
        cnarr,
        args.method,
        args.diploid_parx_genome,
        args.threshold,
        variants=variants,
        skip_low=args.drop_low_coverage,
        skip_outliers=args.drop_outliers,
        save_dataframe=bool(args.dataframe),
        rscript_path=args.rscript_path,
        processes=args.processes,
        smooth_cbs=args.smooth_cbs,
    )

    if args.dataframe:
        segments, dframe = results
        with open(args.dataframe, "w") as handle:
            handle.write(dframe)
        logging.info("Wrote %s", args.dataframe)
    else:
        segments = results
    tabio.write(segments, args.output or segments.sample_id + ".cns")


P_segment = AP_subparsers.add_parser("segment", help=_cmd_segment.__doc__)
P_segment.add_argument(
    "filename", help="Bin-level log2 ratios (.cnr file), as produced by 'fix'."
)
P_segment.add_argument(
    "-o",
    "--output",
    metavar="FILENAME",
    help="Output table file name (CNR-like table of segments, .cns).",
)
P_segment.add_argument(
    "-d",
    "--dataframe",
    help="""File name to save the raw R dataframe emitted by CBS or Fused Lasso. (Useful
            for debugging.)""",
)
P_segment.add_argument(
    "-m",
    "--method",
    choices=segmentation.SEGMENT_METHODS,
    default="cbs",
    help="""Segmentation method (see docs), or 'none' for chromosome arm-level averages
            as segments. [Default: %(default)s]""",
)
P_segment.add_argument(
    "-t",
    "--threshold",
    type=float,
    help="""Significance threshold (p-value or FDR, depending on method) to accept
            breakpoints during segmentation. For HMM methods, this is the smoothing
            window size.""",
)
P_segment.add_argument(
    "--drop-low-coverage",
    action="store_true",
    help="""Drop very-low-coverage bins before segmentation to avoid false-positive
            deletions in poor-quality tumor samples.""",
)
P_segment.add_argument(
    "--drop-outliers",
    metavar="FACTOR",
    type=float,
    default=10,
    help="""Drop outlier bins more than this many multiples of the 95th quantile away
            from the average within a rolling window. Set to 0 for no outlier filtering.
            [Default: %(default)g]""",
)
P_segment.add_argument(
    "--rscript-path",
    metavar="PATH",
    default="Rscript",
    help="""Path to the Rscript executable to use for running R code. Use this option
            to specify a non-default R installation. [Default: %(default)s]""",
)
P_segment.add_argument(
    "-p",
    "--processes",
    nargs="?",
    type=int,
    const=0,
    default=1,
    help="""Number of subprocesses to segment in parallel. Give 0 or a negative value to
            use the maximum number of available CPUs. [Default: use 1 process]""",
)
P_segment.add_argument(
    "--smooth-cbs",
    action="store_true",
    help="""Perform an additional smoothing before CBS segmentation, which in some cases
            may increase the sensitivity. Used only for CBS method.""",
)
add_diploid_parx_genome(P_segment)

P_segment_vcf = P_segment.add_argument_group(
    "To additionally segment SNP b-allele frequencies"
)
P_segment_vcf.add_argument(
    "-v",
    "--vcf",
    metavar="FILENAME",
    help="""VCF file name containing variants for segmentation by allele
            frequencies.""",
)
P_segment_vcf.add_argument(
    "-i",
    "--sample-id",
    help="""Specify the name of the sample in the VCF (-v/--vcf) to use for b-allele
            frequency extraction and as the default plot title.""",
)
P_segment_vcf.add_argument(
    "-n",
    "--normal-id",
    help="""Corresponding normal sample ID in the input VCF (-v/--vcf). This sample is
            used to select only germline SNVs to plot b-allele frequencies.""",
)
P_segment_vcf.add_argument(
    "--min-variant-depth",
    type=int,
    default=20,
    help="""Minimum read depth for a SNV to be displayed in the b-allele frequency plot.
            [Default: %(default)s]""",
)
P_segment_vcf.add_argument(
    "-z",
    "--zygosity-freq",
    metavar="ALT_FREQ",
    nargs="?",
    type=float,
    const=0.25,
    help="""Ignore VCF's genotypes (GT field) and instead infer zygosity from allele
            frequencies. [Default if used without a number: %(const)s]""",
)

P_segment.set_defaults(func=_cmd_segment)


# call ------------------------------------------------------------------------

do_call = public(call.do_call)


def _cmd_call(args):
    """Call copy number variants from segmented log2 ratios."""
    if args.purity and not 0.0 < args.purity <= 1.0:
        raise RuntimeError("Purity must be between 0 and 1.")

    cnarr = read_cna(args.filename)
    if args.center_at:
        logging.info("Shifting log2 ratios by %f", -args.center_at)
        cnarr["log2"] -= args.center_at
    elif args.center:
        cnarr.center_all(args.center, skip_low=args.drop_low_coverage, verbose=True, diploid_parx_genome=args.diploid_parx_genome)

    varr = load_het_snps(
        args.vcf,
        args.sample_id,
        args.normal_id,
        args.min_variant_depth,
        args.zygosity_freq,
    )
    is_sample_female = (
        verify_sample_sex(cnarr, args.sample_sex, args.male_reference, args.diploid_parx_genome)
        if args.purity and args.purity < 1.0
        else None
    )
    cnarr = call.do_call(
        cnarr,
        varr,
        args.method,
        args.ploidy,
        args.purity,
        args.male_reference,
        args.diploid_parx_genome,
        is_sample_female,
        args.filters,
        args.thresholds,
    )
    tabio.write(cnarr, args.output or cnarr.sample_id + ".call.cns")


def csvstring(text):
    return tuple(map(float, text.split(",")))


P_call = AP_subparsers.add_parser("call", help=_cmd_call.__doc__)
P_call.add_argument("filename", help="Copy ratios (.cnr or .cns).")
P_call.add_argument(
    "--center",
    nargs="?",
    const="median",
    choices=("mean", "median", "mode", "biweight"),
    help="""Re-center the log2 ratio values using this estimator of the center or
            average value. ('median' if no argument given.)""",
)
P_call.add_argument(
    "--center-at",
    type=float,
    help="""Subtract a constant number from all log2 ratios. For "manual" re-centering,
            in case the --center option gives unsatisfactory results.)""",
)
P_call.add_argument(
    "--filter",
    action="append",
    default=[],
    dest="filters",
    choices=(
        "ampdel",
        "cn",
        "ci",
        "sem",  # 'bic'
    ),
    help="""Merge segments flagged by the specified filter(s) with the adjacent
            segment(s).""",
)
P_call.add_argument(
    "-m",
    "--method",
    choices=("threshold", "clonal", "none"),
    default="threshold",
    help="""Calling method. [Default: %(default)s]""",
)
P_call.add_argument(
    "-t",
    "--thresholds",
    type=csvstring,
    default="-1.1,-0.25,0.2,0.7",
    help="""Hard thresholds for calling each integer copy number, separated by commas.
            Use the '=' sign on the command line, e.g.: -t=-1,0,1
            [Default: %(default)s]""",
)
P_call.add_argument(
    "--ploidy",
    type=int,
    default=2,
    help="Ploidy of the sample cells. [Default: %(default)d]",
)
P_call.add_argument(
    "--purity",
    type=float,
    help="Estimated tumor cell fraction, a.k.a. purity or cellularity.",
)
P_call.add_argument(
    "--drop-low-coverage",
    action="store_true",
    help="""Drop very-low-coverage bins before segmentation to avoid false-positive
            deletions in poor-quality tumor samples.""",
)
P_call.add_argument(
    "-x",
    "--sample-sex",
    "-g",
    "--gender",
    dest="sample_sex",
    choices=("m", "y", "male", "Male", "f", "x", "female", "Female"),
    help="""Specify the sample's chromosomal sex as male or female. (Otherwise guessed
            from X and Y coverage).""",
)
P_call.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Was a male reference used?  If so, expect half ploidy on chrX and chrY;
            otherwise, only chrY has half ploidy. In CNVkit, if a male reference was
            used, the "neutral" copy number (ploidy) of chrX is 1; chrY is haploid for
            either reference sex.""",
)
P_call.add_argument(
    "-o",
    "--output",
    metavar="FILENAME",
    help="Output table file name (CNR-like table of segments, .cns).",
)

P_call_vcf = P_call.add_argument_group(
    "To additionally process SNP b-allele frequencies for allelic copy number"
)
P_call_vcf.add_argument(
    "-v",
    "--vcf",
    metavar="FILENAME",
    help="""VCF file name containing variants for calculation of b-allele
            frequencies.""",
)
P_call_vcf.add_argument(
    "-i",
    "--sample-id",
    help="""Name of the sample in the VCF (-v/--vcf) to use for b-allele
            frequency extraction.""",
)
P_call_vcf.add_argument(
    "-n",
    "--normal-id",
    help="""Corresponding normal sample ID in the input VCF (-v/--vcf). This sample is
            used to select only germline SNVs to calculate b-allele frequencies.""",
)
P_call_vcf.add_argument(
    "--min-variant-depth",
    type=int,
    default=20,
    help="""Minimum read depth for a SNV to be used in the b-allele frequency
            calculation. [Default: %(default)s]""",
)
P_call_vcf.add_argument(
    "-z",
    "--zygosity-freq",
    metavar="ALT_FREQ",
    nargs="?",
    type=float,
    const=0.25,
    help="""Ignore VCF's genotypes (GT field) and instead infer zygosity from allele
            frequencies. [Default if used without a number: %(const)s]""",
)
add_diploid_parx_genome(P_call)
P_call.set_defaults(func=_cmd_call)


# _____________________________________________________________________________
# Plots and graphics

# diagram ---------------------------------------------------------------------


def _cmd_diagram(args):
    """Draw copy number (log2 coverages, segments) on chromosomes as a diagram.

    If both the raw probes and segments are given, show them side-by-side on
    each chromosome (segments on the left side, probes on the right side).
    """
    if not args.filename and not args.segment:
        raise ValueError(
            "Must specify a filename as an argument or with "
            "the '-s' option, or both. You did neither."
        )

    cnarr = read_cna(args.filename) if args.filename else None
    segarr = read_cna(args.segment) if args.segment else None
    if args.adjust_xy:
        is_sample_female = verify_sample_sex(
            cnarr or segarr, args.sample_sex, args.male_reference, args.diploid_parx_genome
        )
        if cnarr:
            cnarr = cnarr.shift_xx(args.male_reference, is_sample_female)
        if segarr:
            segarr = segarr.shift_xx(args.male_reference, is_sample_female)
    outfname = diagram.create_diagram(
        cnarr,
        segarr,
        args.threshold,
        args.min_probes,
        args.output,
        args.chromosome,
        args.title,
        args.show_labels,
    )
    logging.info("Wrote %s", outfname)


P_diagram = AP_subparsers.add_parser("diagram", help=_cmd_diagram.__doc__)
P_diagram.add_argument(
    "filename",
    nargs="?",
    help="""Processed coverage data file (*.cnr), the output of the 'fix'
            sub-command.""",
)
P_diagram.add_argument(
    "-s",
    "--segment",
    help="Segmentation calls (.cns), the output of the 'segment' command.",
)
P_diagram.add_argument(
    "-c",
    "--chromosome",
    help="""Chromosome to display, e.g. 'chr1' (no chromosomal range allowed)""",
)
P_diagram.add_argument(
    "-t",
    "--threshold",
    type=float,
    default=0.5,
    help="""Copy number change threshold to label genes. [Default: %(default)s]""",
)
P_diagram.add_argument(
    "-m",
    "--min-probes",
    type=int,
    default=3,
    help="""Minimum number of covered probes to label a gene. [Default: %(default)d]""",
)
P_diagram.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Assume inputs were normalized to a male reference (i.e. female samples will
            have +1 log-CNR of chrX; otherwise male samples would have -1 chrX).""",
)
P_diagram.add_argument(
    "-x",
    "--sample-sex",
    "-g",
    "--gender",
    dest="sample_sex",
    choices=("m", "y", "male", "Male", "f", "x", "female", "Female"),
    help="""Specify the sample's chromosomal sex as male or female. (Otherwise guessed
            from X and Y coverage).""",
)
P_diagram.add_argument(
    "--no-shift-xy",
    dest="adjust_xy",
    action="store_false",
    help="Don't adjust the X and Y chromosomes according to sample sex.",
)
P_diagram.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output PDF file name."
)
P_diagram_aes = P_diagram.add_argument_group("Plot aesthetics")
P_diagram_aes.add_argument(
    "--title", help="Plot title. [Default: sample ID, from filename or -i]"
)
P_diagram_aes.add_argument(
    "--no-gene-labels",
    dest="show_labels",
    action="store_false",
    help="""Disable gene_name labels on plot (useful when a lot of CNV were called).""",
)
P_diagram.set_defaults(func=_cmd_diagram)


# scatter ---------------------------------------------------------------------

do_scatter = public(scatter.do_scatter)


def _cmd_scatter(args):
    """Plot probe log2 coverages and segmentation calls together."""
    cnarr = read_cna(args.filename, sample_id=args.sample_id) if args.filename else None
    segarr = read_cna(args.segment, sample_id=args.sample_id) if args.segment else None
    varr = load_het_snps(
        args.vcf,
        args.sample_id,
        args.normal_id,
        args.min_variant_depth,
        args.zygosity_freq,
    )
    scatter_opts = {
        k: v
        for k, v in (
            ("do_trend", args.trend),
            ("by_bin", args.by_bin),
            ("window_width", args.width),
            ("y_min", args.y_min),
            ("y_max", args.y_max),
            ("fig_size", args.fig_size),
            ("antitarget_marker", args.antitarget_marker),
            ("segment_color", args.segment_color),
        )
        if v is not None
    }

    if args.range_list:
        with PdfPages(args.output) as pdf_out:
            for region in tabio.read_auto(args.range_list).coords():
                try:
                    if args.title is not None:
                        scatter_opts["title"] = f"{args.title} {region.chromosome}"
                    scatter.do_scatter(
                        cnarr, segarr, varr, show_range=region, **scatter_opts
                    )
                except ValueError as exc:
                    # Probably no bins in the selected region
                    logging.warning("Not plotting region %r: %s", to_label(region), exc)
                pdf_out.savefig()
                pyplot.close()
    else:
        if args.title is not None:
            scatter_opts["title"] = args.title
        scatter.do_scatter(
            cnarr, segarr, varr, args.chromosome, args.gene, **scatter_opts
        )
        if args.output:
            oformat = os.path.splitext(args.output)[-1].replace(".", "")
            pyplot.savefig(args.output, format=oformat, bbox_inches="tight")
            logging.info("Wrote %s", args.output)
        else:
            pyplot.show()


P_scatter = AP_subparsers.add_parser("scatter", help=_cmd_scatter.__doc__)
P_scatter.add_argument(
    "filename",
    nargs="?",
    help="""Processed bin-level copy ratios (*.cnr), the output of the 'fix'
            sub-command.""",
)
P_scatter.add_argument(
    "-s",
    "--segment",
    metavar="FILENAME",
    help="Segmentation calls (.cns), the output of the 'segment' command.",
)
P_scatter.add_argument(
    "-c",
    "--chromosome",
    metavar="RANGE",
    help="""Chromosome or chromosomal range, e.g. 'chr1' or 'chr1:2333000-2444000', to
            display. If a range is given, all targeted genes in this range will be
            shown, unless -g/--gene is also given.""",
)
P_scatter.add_argument(
    "-g", "--gene", help="Name of gene or genes (comma-separated) to display."
)
P_scatter.add_argument(
    "-l",
    "--range-list",
    help="""File listing the chromosomal ranges to display, as BED, interval list or
            'chr:start-end' text. Creates focal plots similar to -c/--chromosome for
            each listed region, combined into a multi-page PDF.  The output filename
            must also be specified (-o/--output).""",
)
P_scatter.add_argument(
    "-w",
    "--width",
    type=float,
    default=1e6,
    help="""Width of margin to show around the selected gene(s) (-g/--gene) or small
            chromosomal region (-c/--chromosome). [Default: %(default)d]""",
)
P_scatter.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output PDF file name."
)

P_scatter_aes = P_scatter.add_argument_group("Plot aesthetics")
P_scatter_aes.add_argument(
    "-a",
    "--antitarget-marker",
    metavar="CHARACTER",
    dest="antitarget_marker",
    default=None,
    help="""Plot antitargets using this symbol when plotting in a selected chromosomal
            region (-g/--gene or -c/--chromosome). [Default: same as targets]""",
)
P_scatter_aes.add_argument(
    "--by-bin",
    action="store_true",
    help="""Plot data x-coordinates by bin indices instead of genomic coordinates. All
            bins will be shown with equal width, no blank regions will be shown, and
            x-axis values indicate bin number (within chromosome) instead of genomic
            position.""",
)
P_scatter_aes.add_argument(
    "--segment-color",
    default=scatter.SEG_COLOR,
    help="""Plot segment lines in this color. Value can be any string accepted by
            matplotlib, e.g. 'red' or '#CC0000'.""",
)
P_scatter_aes.add_argument(
    "--title", help="Plot title. [Default: sample ID, from filename or -i]"
)
P_scatter_aes.add_argument(
    "-t",
    "--trend",
    action="store_true",
    help="Draw a smoothed local trendline on the scatter plot.",
)
P_scatter_aes.add_argument("--y-max", type=float, help="y-axis upper limit.")
P_scatter_aes.add_argument("--y-min", type=float, help="y-axis lower limit.")
P_scatter_aes.add_argument(
    "--fig-size",
    nargs=2,
    metavar=("WIDTH", "HEIGHT"),
    type=float,
    help="""Width and height of the plot in inches. [Default: Pre-defined in Matplotlib
            'rcParams' variable (most of the time: '6.4 4.8')]""",
)

P_scatter_vcf = P_scatter.add_argument_group("To plot SNP b-allele frequencies")
P_scatter_vcf.add_argument(
    "-v",
    "--vcf",
    metavar="FILENAME",
    help="""VCF file name containing variants to plot for SNV b-allele frequencies.""",
)
P_scatter_vcf.add_argument(
    "-i",
    "--sample-id",
    help="""Name of the sample in the VCF to use for b-allele frequency extraction and
            as the default plot title.""",
)
P_scatter_vcf.add_argument(
    "-n",
    "--normal-id",
    help="""Corresponding normal sample ID in the input VCF. This sample is used to
            select only germline SNVs to plot.""",
)
P_scatter_vcf.add_argument(
    "-m",
    "--min-variant-depth",
    type=int,
    default=20,
    help="""Minimum read depth for a SNV to be used in the b-allele frequency
            calculation. [Default: %(default)s]""",
)
P_scatter_vcf.add_argument(
    "-z",
    "--zygosity-freq",
    metavar="ALT_FREQ",
    nargs="?",
    type=float,
    const=0.25,
    help="""Ignore VCF's genotypes (GT field) and instead infer zygosity from allele
            frequencies. [Default if used without a number: %(const)s]""",
)

P_scatter.set_defaults(func=_cmd_scatter)


# heatmap ---------------------------------------------------------------------

do_heatmap = public(heatmap.do_heatmap)


def _cmd_heatmap(args):
    """Plot copy number for multiple samples as a heatmap."""
    cnarrs = []
    for fname in args.filenames:
        cnarr = read_cna(fname)
        if args.adjust_xy:
            is_sample_female = verify_sample_sex(
                cnarr, args.sample_sex, args.male_reference, args.diploid_parx_genome
            )
            cnarr = cnarr.shift_xx(args.male_reference, is_sample_female, args.diploid_parx_genome)
        cnarrs.append(cnarr)
    heatmap.do_heatmap(
        cnarrs,
        args.chromosome,
        args.desaturate,
        args.by_bin,
        args.delim_sampl,
        args.vertical,
        args.title,
    )
    if args.output:
        oformat = os.path.splitext(args.output)[-1].replace(".", "")
        pyplot.savefig(args.output, format=oformat, bbox_inches="tight")
        logging.info("Wrote %s", args.output)
    else:
        pyplot.show()


P_heatmap = AP_subparsers.add_parser("heatmap", help=_cmd_heatmap.__doc__)
P_heatmap.add_argument(
    "filenames",
    nargs="+",
    help="Sample coverages as raw probes (.cnr) or segments (.cns).",
)
P_heatmap.add_argument(
    "-c",
    "--chromosome",
    help="""Chromosome (e.g. 'chr1') or chromosomal range (e.g. 'chr1:2333000-2444000')
            to display. If a range is given, all targeted genes in this range will be
            shown, unless '--gene'/'-g' is already given.""",
)
P_heatmap.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Assume inputs were normalized to a male reference (i.e. female samples will
            have +1 log-CNR of chrX; otherwise male samples would have -1 chrX).""",
)
P_heatmap.add_argument(
    "-x",
    "--sample-sex",
    "-g",
    "--gender",
    dest="sample_sex",
    choices=("m", "y", "male", "Male", "f", "x", "female", "Female"),
    help="""Specify the chromosomal sex of all given samples as male or female.
            [Default: guess each sample from coverage of X and Y chromosomes].""",
)
P_heatmap.add_argument(
    "--no-shift-xy",
    dest="adjust_xy",
    action="store_false",
    help="Don't adjust the X and Y chromosomes according to sample sex.",
)
P_heatmap.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output PDF file name."
)
# P_heatmap.add_argument('-g', '--gene',
#         help="Name of gene to display.")

P_heatmap_aes = P_heatmap.add_argument_group("Plot aesthetics")
P_heatmap_aes.add_argument(
    "-b",
    "--by-bin",
    action="store_true",
    help="""Plot data x-coordinates by bin indices instead of genomic coordinates. All
            bins will be shown with equal width, no blank regions will be shown, and
            x-axis values indicate bin number (within chromosome) instead of genomic
            position.""",
)
P_heatmap_aes.add_argument(
    "-d",
    "--desaturate",
    action="store_true",
    help="Tweak color saturation to focus on significant changes.",
)
P_heatmap_aes.add_argument(
    "-v",
    "--vertical",
    action="store_true",
    help="Plot heatmap with samples as X-axis (instead of Y-axis).",
)
P_heatmap_aes.add_argument(
    "--delimit-samples",
    action="store_true",
    dest="delim_sampl",
    help="Add an horizontal delimitation line between each sample.",
)
P_heatmap_aes.add_argument(
    "-t", "--title", help="Plot title. [Default: Range if provided, otherwise none]"
)
P_heatmap.set_defaults(func=_cmd_heatmap)


# _____________________________________________________________________________
# Tabular outputs

# breaks ----------------------------------------------------------------------

do_breaks = public(reports.do_breaks)


def _cmd_breaks(args):
    """List the targeted genes in which a copy number breakpoint occurs."""
    cnarr = read_cna(args.filename)
    segarr = read_cna(args.segment)
    bpoints = do_breaks(cnarr, segarr, args.min_probes)
    logging.info("Found %d gene breakpoints", len(bpoints))
    write_dataframe(args.output, bpoints)


P_breaks = AP_subparsers.add_parser("breaks", help=_cmd_breaks.__doc__)
P_breaks.add_argument(
    "filename",
    help="""Processed sample coverage data file (*.cnr), the output of the 'fix'
            sub-command.""",
)
P_breaks.add_argument(
    "segment", help="Segmentation calls (.cns), the output of the 'segment' command)."
)
P_breaks.add_argument(
    "-m",
    "--min-probes",
    type=int,
    default=1,
    help="""Minimum number of within-gene probes on both sides of a breakpoint to report
            it. [Default: %(default)d]""",
)
P_breaks.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output table file name."
)
P_breaks.set_defaults(func=_cmd_breaks)


# genemetrics/gainloss --------------------------------------------------------

do_genemetrics = public(reports.do_genemetrics)


def _cmd_genemetrics(args):
    """Identify targeted genes with copy number gain or loss."""
    cnarr = read_cna(args.filename)
    segarr = read_cna(args.segment) if args.segment else None
    is_sample_female = verify_sample_sex(cnarr, args.sample_sex, args.male_reference, args.diploid_parx_genome)
    # TODO use the stats args
    table = do_genemetrics(
        cnarr,
        segarr,
        args.threshold,
        args.min_probes,
        args.drop_low_coverage,
        args.male_reference,
        is_sample_female,
        args.diploid_parx_genome,
    )
    logging.info("Found %d gene-level gains and losses", len(table))
    write_dataframe(args.output, table)


P_genemetrics = AP_subparsers.add_parser("genemetrics", help=_cmd_genemetrics.__doc__)
P_genemetrics.add_argument(
    "filename",
    help="""Processed sample coverage data file (*.cnr), the output of the 'fix'
            sub-command.""",
)
P_genemetrics.add_argument(
    "-s",
    "--segment",
    help="Segmentation calls (.cns), the output of the 'segment' command).",
)
P_genemetrics.add_argument(
    "-t",
    "--threshold",
    type=float,
    default=0.2,
    help="""Copy number change threshold to report a gene gain/loss.
            [Default: %(default)s]""",
)
P_genemetrics.add_argument(
    "-m",
    "--min-probes",
    type=int,
    default=3,
    help="""Minimum number of covered probes to report a gain/loss.
            [Default: %(default)d]""",
)
P_genemetrics.add_argument(
    "--drop-low-coverage",
    action="store_true",
    help="""Drop very-low-coverage bins before segmentation to avoid false-positive
            deletions in poor-quality tumor samples.""",
)
P_genemetrics.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Assume inputs were normalized to a male reference (i.e. female samples will
            have +1 log-coverage of chrX; otherwise male samples would have -1
            chrX).""",
)
P_genemetrics.add_argument(
    "-x",
    "--sample-sex",
    "-g",
    "--gender",
    dest="sample_sex",
    choices=("m", "y", "male", "Male", "f", "x", "female", "Female"),
    help="""Specify the sample's chromosomal sex as male or female. (Otherwise guessed
            from X and Y coverage).""",
)
P_genemetrics.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output table file name."
)
add_diploid_parx_genome(P_genemetrics)

P_genemetrics_stats = P_genemetrics.add_argument_group("Statistics available")
# Location statistics
P_genemetrics_stats.add_argument(
    "--mean",
    action="append_const",
    dest="location_stats",
    const="mean",
    help="Mean log2-ratio (unweighted).",
)
P_genemetrics_stats.add_argument(
    "--median",
    action="append_const",
    dest="location_stats",
    const="median",
    help="Median.",
)
P_genemetrics_stats.add_argument(
    "--mode",
    action="append_const",
    dest="location_stats",
    const="mode",
    help="Mode (i.e. peak density of log2 ratios).",
)
P_genemetrics_stats.add_argument(
    "--ttest",
    action="append_const",
    dest="location_stats",
    const="p_ttest",
    help="One-sample t-test of bin log2 ratios versus 0.0.",
)
# Dispersion statistics
P_genemetrics_stats.add_argument(
    "--stdev",
    action="append_const",
    dest="spread_stats",
    const="stdev",
    help="Standard deviation.",
)
P_genemetrics_stats.add_argument(
    "--sem",
    action="append_const",
    dest="spread_stats",
    const="sem",
    help="Standard error of the mean.",
)
P_genemetrics_stats.add_argument(
    "--mad",
    action="append_const",
    dest="spread_stats",
    const="mad",
    help="Median absolute deviation (standardized).",
)
P_genemetrics_stats.add_argument(
    "--mse",
    action="append_const",
    dest="spread_stats",
    const="mse",
    help="Mean squared error.",
)
P_genemetrics_stats.add_argument(
    "--iqr",
    action="append_const",
    dest="spread_stats",
    const="iqr",
    help="Inter-quartile range.",
)
P_genemetrics_stats.add_argument(
    "--bivar",
    action="append_const",
    dest="spread_stats",
    const="bivar",
    help="Tukey's biweight midvariance.",
)
# Interval statistics
P_genemetrics_stats.add_argument(
    "--ci",
    action="append_const",
    dest="interval_stats",
    const="ci",
    help="Confidence interval (by bootstrap).",
)
P_genemetrics_stats.add_argument(
    "--pi",
    action="append_const",
    dest="interval_stats",
    const="pi",
    help="Prediction interval.",
)
P_genemetrics_stats.add_argument(
    "-a",
    "--alpha",
    type=float,
    default=0.05,
    help="""Level to estimate confidence and prediction intervals; use with --ci and
            --pi. [Default: %(default)s]""",
)
P_genemetrics_stats.add_argument(
    "-b",
    "--bootstrap",
    type=int,
    default=100,
    help="""Number of bootstrap iterations to estimate confidence interval; use with
        --ci. [Default: %(default)d]""",
)
P_genemetrics_stats.set_defaults(location_stats=[], spread_stats=[], interval_stats=[])
P_genemetrics.set_defaults(func=_cmd_genemetrics)

# Shims
AP_subparsers._name_parser_map["gainloss"] = P_genemetrics
do_gainloss = public(do_genemetrics)


# sex/gender ------------------------------------------------------------------


def _cmd_sex(args):
    """Guess samples' sex from the relative coverage of chromosomes X and Y."""
    cnarrs = map(read_cna, args.filenames)
    table = do_sex(cnarrs, args.male_reference, args.diploid_parx_genome)
    write_dataframe(args.output, table, header=True)


@public
def do_sex(cnarrs, is_haploid_x_reference, diploid_parx_genome):
    """Guess samples' sex from the relative coverage of chromosomes X and Y."""

    def strsign(num):
        if num > 0:
            return "+%.3g" % num
        return "%.3g" % num

    def guess_and_format(cna):
        is_xy, stats = cna.compare_sex_chromosomes(is_haploid_x_reference, diploid_parx_genome)
        return (
            cna.meta["filename"] or cna.sample_id,
            "Male" if is_xy else "Female",
            strsign(stats["chrx_ratio"]) if stats else "NA",
            strsign(stats["chry_ratio"]) if stats else "NA",
        )

    rows = (guess_and_format(cna) for cna in cnarrs)
    columns = ["sample", "sex", "X_logratio", "Y_logratio"]
    return pd.DataFrame.from_records(rows, columns=columns)


P_sex = AP_subparsers.add_parser("sex", help=_cmd_sex.__doc__)
P_sex.add_argument(
    "filenames", nargs="+", help="Copy number or copy ratio files (*.cnn, *.cnr)."
)
P_sex.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Assume inputs were normalized to a male reference (i.e. female samples will
            have +1 log-coverage of chrX; otherwise male samples would have -1 chrX).""",
)
P_sex.add_argument("-o", "--output", metavar="FILENAME", help="Output table file name.")
add_diploid_parx_genome(P_sex)
P_sex.set_defaults(func=_cmd_sex)

# Shims
AP_subparsers._name_parser_map["gender"] = P_sex
do_gender = public(do_sex)


# metrics ---------------------------------------------------------------------

do_metrics = public(metrics.do_metrics)


def _cmd_metrics(args):
    """Compute coverage deviations and other metrics for self-evaluation."""
    if (
        len(args.cnarrays) > 1
        and args.segments
        and len(args.segments) > 1
        and len(args.cnarrays) != len(args.segments)
    ):
        raise ValueError(
            "Number of coverage/segment filenames given must be "
            "equal, if more than 1 segment file is given."
        )

    cnarrs = map(read_cna, args.cnarrays)
    if args.segments:
        args.segments = map(read_cna, args.segments)
    table = metrics.do_metrics(cnarrs, args.segments, args.drop_low_coverage)
    write_dataframe(args.output, table)


P_metrics = AP_subparsers.add_parser("metrics", help=_cmd_metrics.__doc__)
P_metrics.add_argument(
    "cnarrays",
    nargs="+",
    help="""One or more bin-level coverage data files (*.cnn, *.cnr).""",
)
P_metrics.add_argument(
    "-s",
    "--segments",
    nargs="+",
    help="""One or more segmentation data files (*.cns, output of the 'segment'
            command). If more than one file is given, the number must match the coverage
            data files, in which case the input files will be paired together in the
            given order. Otherwise, the same segments will be used for all coverage
            files.""",
)
P_metrics.add_argument(
    "--drop-low-coverage",
    action="store_true",
    help="""Drop very-low-coverage bins before calculations to reduce negative "fat
            tail" of bin log2 values in poor-quality tumor samples.""",
)
P_metrics.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output table file name."
)
P_metrics.set_defaults(func=_cmd_metrics)


# segmetrics ------------------------------------------------------------------

do_segmetrics = public(segmetrics.do_segmetrics)


def _cmd_segmetrics(args):
    """Compute segment-level metrics from bin-level log2 ratios."""
    if not 0.0 < args.alpha <= 1.0:
        raise RuntimeError("alpha must be between 0 and 1.")

    if not any((args.location_stats, args.spread_stats, args.interval_stats)):
        logging.info("No stats specified")
        return

    # Calculate all metrics
    cnarr = read_cna(args.cnarray)
    segarr = read_cna(args.segments)
    segarr = do_segmetrics(
        cnarr,
        segarr,
        args.location_stats,
        args.spread_stats,
        args.interval_stats,
        args.alpha,
        args.bootstrap,
        args.smooth_bootstrap,
        skip_low=args.drop_low_coverage,
    )
    tabio.write(segarr, args.output or segarr.sample_id + ".segmetrics.cns")


P_segmetrics = AP_subparsers.add_parser("segmetrics", help=_cmd_segmetrics.__doc__)
P_segmetrics.add_argument(
    "cnarray", help="""Bin-level copy ratio data file (*.cnn, *.cnr)."""
)
P_segmetrics.add_argument(
    "-s",
    "--segments",
    required=True,
    help="Segmentation data file (*.cns, output of the 'segment' command).",
)
P_segmetrics.add_argument(
    "--drop-low-coverage",
    action="store_true",
    help="""Drop very-low-coverage bins before calculations to avoid negative bias in
            poor-quality tumor samples.""",
)
P_segmetrics.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output table file name."
)

P_segmetrics_stats = P_segmetrics.add_argument_group("Statistics available")
# Location statistics
P_segmetrics_stats.add_argument(
    "--mean",
    action="append_const",
    dest="location_stats",
    const="mean",
    help="Mean log2 ratio (unweighted).",
)
P_segmetrics_stats.add_argument(
    "--median",
    action="append_const",
    dest="location_stats",
    const="median",
    help="Median.",
)
P_segmetrics_stats.add_argument(
    "--mode",
    action="append_const",
    dest="location_stats",
    const="mode",
    help="Mode (i.e. peak density of bin log2 ratios).",
)
P_segmetrics_stats.add_argument(
    "--t-test",
    action="append_const",
    dest="location_stats",
    const="p_ttest",
    help="One-sample t-test of bin log2 ratios versus 0.0.",
)
# Dispersion statistics
P_segmetrics_stats.add_argument(
    "--stdev",
    action="append_const",
    dest="spread_stats",
    const="stdev",
    help="Standard deviation.",
)
P_segmetrics_stats.add_argument(
    "--sem",
    action="append_const",
    dest="spread_stats",
    const="sem",
    help="Standard error of the mean.",
)
P_segmetrics_stats.add_argument(
    "--mad",
    action="append_const",
    dest="spread_stats",
    const="mad",
    help="Median absolute deviation (standardized).",
)
P_segmetrics_stats.add_argument(
    "--mse",
    action="append_const",
    dest="spread_stats",
    const="mse",
    help="Mean squared error.",
)
P_segmetrics_stats.add_argument(
    "--iqr",
    action="append_const",
    dest="spread_stats",
    const="iqr",
    help="Inter-quartile range.",
)
P_segmetrics_stats.add_argument(
    "--bivar",
    action="append_const",
    dest="spread_stats",
    const="bivar",
    help="Tukey's biweight midvariance.",
)
# Interval statistics
P_segmetrics_stats.add_argument(
    "--ci",
    action="append_const",
    dest="interval_stats",
    const="ci",
    help="Confidence interval (by bootstrap).",
)
P_segmetrics_stats.add_argument(
    "--pi",
    action="append_const",
    dest="interval_stats",
    const="pi",
    help="Prediction interval.",
)
P_segmetrics_stats.add_argument(
    "-a",
    "--alpha",
    type=float,
    default=0.05,
    help="""Level to estimate confidence and prediction intervals; use with --ci and
            --pi. [Default: %(default)s]""",
)
P_segmetrics_stats.add_argument(
    "-b",
    "--bootstrap",
    type=int,
    default=100,
    help="""Number of bootstrap iterations to estimate confidence interval; use with
            --ci. [Default: %(default)d]""",
)
P_segmetrics_stats.add_argument(
    "--smooth-bootstrap",
    action="store_true",
    help="""Apply Gaussian noise to bootstrap samples, a.k.a. smoothed bootstrap, to
            estimate confidence interval; use with --ci.""",
)

P_segmetrics_stats.set_defaults(location_stats=[], spread_stats=[], interval_stats=[])
P_segmetrics.set_defaults(func=_cmd_segmetrics)


# bintest -----------------------------------------------------------------------

do_bintest = public(bintest.do_bintest)


def _cmd_bintest(args):
    """Test for single-bin copy number alterations."""
    cnarr = read_cna(args.cnarray)
    segments = read_cna(args.segment) if args.segment else None
    sig = do_bintest(cnarr, segments, args.alpha, args.target)
    tabio.write(sig, args.output or sys.stdout)


P_bintest = AP_subparsers.add_parser("bintest", help=_cmd_bintest.__doc__)
P_bintest.add_argument(
    "cnarray", help="Bin-level log2 ratios (.cnr file), as produced by 'fix'."
)
P_bintest.add_argument(
    "-s",
    "--segment",
    metavar="FILENAME",
    help="""Segmentation calls (.cns), the output of the 'segment' command).""",
)
P_bintest.add_argument(
    "-a",
    "--alpha",
    type=float,
    default=0.005,
    help="Significance threhold. [Default: %(default)s]",
)
P_bintest.add_argument(
    "-t",
    "--target",
    action="store_true",
    help="Test target bins only; ignore off-target bins.",
)
P_bintest.add_argument("-o", "--output", help="Output filename.")
P_bintest.set_defaults(func=_cmd_bintest)


# _____________________________________________________________________________
# Other I/O and compatibility

# import-picard ---------------------------------------------------------------


def _cmd_import_picard(args):
    """Convert Picard CalculateHsMetrics tabular output to CNVkit .cnn files.

    The input file is generated by the PER_TARGET_COVERAGE option in the
    CalculateHsMetrics script in Picard tools.

    If 'antitarget' is in the input filename, the generated output filename will
    have the suffix '.antitargetcoverage.cnn', otherwise '.targetcoverage.cnn'.
    """
    for fname in args.targets:
        if not os.path.isfile(fname):
            # Legacy usage: previously accepted directory as an argument
            raise ValueError("Not a file: %s" % fname)
        garr = importers.do_import_picard(fname)
        outfname = "{}.{}targetcoverage.cnn".format(
            garr.sample_id, "anti" if "antitarget" in fname else ""
        )
        if args.output_dir:
            if not os.path.isdir(args.output_dir):
                os.mkdir(args.output_dir)
                logging.info("Created directory %s", args.output_dir)
            outfname = os.path.join(args.output_dir, outfname)
        tabio.write(garr, outfname)


P_import_picard = AP_subparsers.add_parser(
    "import-picard", help=_cmd_import_picard.__doc__
)
P_import_picard.add_argument(
    "targets", nargs="+", help="""Sample coverage .csv files (target and antitarget)."""
)
P_import_picard.add_argument(
    "-d",
    "--output-dir",
    metavar="DIRECTORY",
    default=".",
    help="Output directory name.",
)
P_import_picard.set_defaults(func=_cmd_import_picard)


# import-seg ------------------------------------------------------------------


def _cmd_import_seg(args):
    """Convert a SEG file to CNVkit .cns files."""
    from .cnary import CopyNumArray as _CNA

    if args.chromosomes:
        if args.chromosomes == "human":
            chrom_names = {"23": "X", "24": "Y", "25": "M"}
        else:
            chrom_names = dict(kv.split(":") for kv in args.chromosomes.split(","))
    else:
        chrom_names = args.chromosomes
    for sid, segtable in tabio.seg.parse_seg(
        args.segfile, chrom_names, args.prefix, args.from_log10
    ):
        segarr = _CNA(segtable, {"sample_id": sid})
        tabio.write(segarr, os.path.join(args.output_dir, sid + ".cns"))


P_import_seg = AP_subparsers.add_parser("import-seg", help=_cmd_import_seg.__doc__)
P_import_seg.add_argument(
    "segfile", help="""Input file in SEG format. May contain multiple samples."""
)
P_import_seg.add_argument(
    "-c",
    "--chromosomes",
    help="""Mapping of chromosome indexes to names. Syntax: "from1:to1,from2:to2".
            Or use "human" for the preset: "23:X,24:Y,25:M".""",
)
P_import_seg.add_argument(
    "-p",
    "--prefix",
    help="""Prefix to add to chromosome names (e.g 'chr' to rename '8' in the SEG file
            to 'chr8' in the output).""",
)
P_import_seg.add_argument(
    "--from-log10",
    action="store_true",
    help="Convert base-10 logarithm values in the input to base-2 logs.",
)
P_import_seg.add_argument(
    "-d",
    "--output-dir",
    metavar="DIRECTORY",
    default=".",
    help="Output directory name.",
)
P_import_seg.set_defaults(func=_cmd_import_seg)


# import-theta ---------------------------------------------------------------

do_import_theta = public(importers.do_import_theta)


def _cmd_import_theta(args):
    """Convert THetA output to a BED-like, CNVkit-like tabular format.

    Equivalently, use the THetA results file to convert CNVkit .cns segments to
    integer copy number calls.
    """
    tumor_segs = read_cna(args.tumor_cns)
    for i, new_cns in enumerate(
        do_import_theta(tumor_segs, args.theta_results, args.ploidy)
    ):
        tabio.write(
            new_cns,
            os.path.join(args.output_dir, "%s-%d.cns" % (tumor_segs.sample_id, i + 1)),
        )


P_import_theta = AP_subparsers.add_parser(
    "import-theta", help=_cmd_import_theta.__doc__
)
P_import_theta.add_argument("tumor_cns")
P_import_theta.add_argument("theta_results")
P_import_theta.add_argument(
    "--ploidy",
    type=int,
    default=2,
    help="Ploidy of normal cells. [Default: %(default)d]",
)
P_import_theta.add_argument(
    "-d",
    "--output-dir",
    metavar="DIRECTORY",
    default=".",
    help="Output directory name.",
)
P_import_theta.set_defaults(func=_cmd_import_theta)


# import-rna ------------------------------------------------------------------

do_import_rna = public(import_rna.do_import_rna)


def _cmd_import_rna(args):
    """Convert a cohort of per-gene log2 ratios to CNVkit .cnr format."""
    all_data, cnrs = import_rna.do_import_rna(
        args.gene_counts,
        args.format,
        args.gene_resource,
        args.correlations,
        args.normal,
        args.do_gc,
        args.do_txlen,
        args.max_log2,
    )
    logging.info("Writing output files")
    if args.output:
        all_data.to_csv(args.output, sep="\t", index=True)
        logging.info("Wrote %s with %d rows", args.output, len(all_data))
    else:
        logging.info(all_data.describe(), file=sys.stderr)
    for cnr in cnrs:
        outfname = os.path.join(args.output_dir, cnr.sample_id + ".cnr")
        tabio.write(cnr, outfname, "tab")


P_import_rna = AP_subparsers.add_parser("import-rna", help=_cmd_import_rna.__doc__)
P_import_rna.add_argument(
    "gene_counts",
    nargs="+",
    metavar="FILES",
    help="""Tabular files with Ensembl gene ID and number of reads mapped to each gene,
            from RSEM or another transcript quantifier.""",
)
P_import_rna.add_argument(
    "-f",
    "--format",
    choices=("rsem", "counts"),
    default="counts",
    metavar="NAME",
    help="""Input format name: 'rsem' for RSEM gene-level read counts
            (*_rsem.genes.results), or 'counts' for generic 2-column gene IDs and their
            read counts (e.g. TCGA level 2 RNA expression).""",
)
P_import_rna.add_argument(
    "-g",
    "--gene-resource",
    metavar="FILE",
    required=True,
    help="Location of gene info table from Ensembl BioMart.",
)
P_import_rna.add_argument(
    "-c",
    "--correlations",
    metavar="FILE",
    help="""Correlation of each gene's copy number with expression. Output of
            cnv_expression_correlate.py.""",
)
P_import_rna.add_argument(
    "--max-log2",
    metavar="FLOAT",
    default=3.0,
    type=float,
    help="""Maximum log2 ratio in output. Observed values above this limit will be
            replaced with this value. [Default: %(default)s]""",
)
P_import_rna.add_argument(
    "-n",
    "--normal",
    nargs="+",
    default=[],
    help="""Normal samples (same format as `gene_counts`) to be used as a control to
            when normalizing and re-centering gene read depth ratios. All filenames
            following this option will be used.""",
)
P_import_rna.add_argument(
    "-d",
    "--output-dir",
    default=".",
    metavar="PATH",
    help="""Directory to write a CNVkit .cnr file for each input sample. [Default:
            %(default)s]""",
)
P_import_rna.add_argument(
    "-o", "--output", metavar="FILE", help="Output file name (summary table)."
)

P_import_rna_bias = P_import_rna.add_argument_group(
    "To disable specific automatic bias corrections"
)
P_import_rna_bias.add_argument(
    "--no-gc", dest="do_gc", action="store_false", help="Skip GC correction."
)
P_import_rna_bias.add_argument(
    "--no-txlen",
    dest="do_txlen",
    action="store_false",
    help="Skip transcript length correction.",
)

P_import_rna.set_defaults(func=_cmd_import_rna)


# export ----------------------------------------------------------------------

P_export = AP_subparsers.add_parser(
    "export", help="""Convert CNVkit output files to another format."""
)
P_export_subparsers = P_export.add_subparsers(
    help="Export formats (use with -h for more info)."
)


# BED special case: multiple samples's segments, like SEG
def _cmd_export_bed(args):
    """Convert segments to BED format.

    Input is a segmentation file (.cns) where, preferably, log2 ratios have
    already been adjusted to integer absolute values using the 'call' command.
    """
    bed_tables = []
    for segfname in args.segments:
        segments = read_cna(segfname)
        # ENH: args.sample_sex as a comma-separated list
        is_sample_female = verify_sample_sex(
            segments, args.sample_sex, args.male_reference, args.diploid_parx_genome
        )
        if args.sample_id:
            label = args.sample_id
        elif args.label_genes:
            label = None
        else:
            label = segments.sample_id
        tbl = export.export_bed(
            segments,
            args.ploidy,
            args.male_reference,
            args.diploid_parx_genome,
            is_sample_female,
            label,
            args.show,
        )
        bed_tables.append(tbl)
    table = pd.concat(bed_tables)
    write_dataframe(args.output, table, header=False)


P_export_bed = P_export_subparsers.add_parser("bed", help=_cmd_export_bed.__doc__)
P_export_bed.add_argument(
    "segments",
    nargs="+",
    help="""Segmented copy ratio data files (*.cns), the output of the 'segment' or
            'call' sub-commands.""",
)
P_export_bed.add_argument(
    "-i",
    "--sample-id",
    metavar="LABEL",
    help="""Identifier to write in the 4th column of the BED file.
            [Default: use the sample ID, taken from the file name]""",
)
P_export_bed.add_argument(
    "--label-genes",
    action="store_true",
    help="""Show gene names in the 4th column of the BED file. (This is a bad idea if >1
            input files are given.)""",
)
P_export_bed.add_argument(
    "--ploidy",
    type=int,
    default=2,
    help="Ploidy of the sample cells. [Default: %(default)d]",
)
P_export_bed.add_argument(
    "-x",
    "--sample-sex",
    "-g",
    "--gender",
    dest="sample_sex",
    choices=("m", "y", "male", "Male", "f", "x", "female", "Female"),
    help="""Specify the sample's chromosomal sex as male or female. (Otherwise guessed
            from X and Y coverage).""",
)
P_export_bed.add_argument(
    "--show",
    choices=("ploidy", "variant", "all"),
    default="ploidy",
    help="""Which segmented regions to show:
            'all' = all segment regions;
            'variant' = CNA regions with non-neutral copy number;
            'ploidy' = CNA regions with non-default ploidy.
            [Default: %(default)s]""",
)
P_export_bed.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Was a male reference used?  If so, expect half ploidy on chrX and chrY;
            otherwise, only chrY has half ploidy. In CNVkit, if a male reference was
            used, the "neutral" copy number (ploidy) of chrX is 1; chrY is haploid for
            either reference sex.""",
)
P_export_bed.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)
P_export_bed.set_defaults(func=_cmd_export_bed)


# SEG special case: segment coords don't match across samples
def _cmd_export_seg(args):
    """Convert segments to SEG format.

    Compatible with IGV and GenePattern.
    """
    table = export.export_seg(args.filenames, chrom_ids=args.enumerate_chroms)
    write_dataframe(args.output, table)


P_export_seg = P_export_subparsers.add_parser("seg", help=_cmd_export_seg.__doc__)
P_export_seg.add_argument(
    "filenames",
    nargs="+",
    help="""Segmented copy ratio data file(s) (*.cns), the output of the 'segment'
            sub-command.""",
)
P_export_seg.add_argument(
    "--enumerate-chroms",
    action="store_true",
    help="""Replace chromosome names with sequential integer IDs.""",
)
P_export_seg.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)
P_export_seg.set_defaults(func=_cmd_export_seg)


# VCF special case: only 1 sample, for now
def _cmd_export_vcf(args):
    """Convert segments to VCF format.

    Input is a segmentation file (.cns) where, preferably, log2 ratios have
    already been adjusted to integer absolute values using the 'call' command.
    """
    segarr = read_cna(args.segments)
    cnarr = read_cna(args.cnr) if args.cnr else None
    is_sample_female = verify_sample_sex(segarr, args.sample_sex, args.male_reference, args.diploid_parx_genome)
    header, body = export.export_vcf(
        segarr,
        args.ploidy,
        args.male_reference,
        args.diploid_parx_genome,
        is_sample_female,
        args.sample_id,
        cnarr,
    )
    write_text(args.output, header, body)


P_export_vcf = P_export_subparsers.add_parser("vcf", help=_cmd_export_vcf.__doc__)
P_export_vcf.add_argument(
    "segments",  # nargs='1',
    help="""Segmented copy ratio data file (*.cns), the output of the 'segment' or
            'call' sub-commands.""",
)
# ENH?: Incorporate left/right CI into .cns via 'segment' or 'segmetrics',
#   potentially calculated another way besides adjacent bin boundaries
P_export_vcf.add_argument(
    "--cnr",
    help="""Bin-level copy ratios (*.cnr). Used to indicate fuzzy boundaries for
            segments in the output VCF via the CIPOS and CIEND tags.""",
)
P_export_vcf.add_argument(
    "-i",
    "--sample-id",
    metavar="LABEL",
    help="""Sample name to write in the genotype field of the output VCF file.
            [Default: use the sample ID, taken from the file name]""",
)
P_export_vcf.add_argument(
    "--ploidy",
    type=int,
    default=2,
    help="Ploidy of the sample cells. [Default: %(default)d]",
)
P_export_vcf.add_argument(
    "-x",
    "--sample-sex",
    "-g",
    "--gender",
    dest="sample_sex",
    choices=("m", "y", "male", "Male", "f", "x", "female", "Female"),
    help="""Specify the sample's chromosomal sex as male or female.
            (Otherwise guessed from X and Y coverage).""",
)
P_export_vcf.add_argument(
    "-y",
    "--male-reference",
    "--haploid-x-reference",
    action="store_true",
    help="""Was a male reference used? If so, expect half ploidy on chrX and chrY;
            otherwise, only chrY has half ploidy. In CNVkit, if a male reference was
            used, the "neutral" copy number (ploidy) of chrX is 1; chrY is haploid for
            either reference sex.""",
)
P_export_vcf.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)
P_export_vcf.set_defaults(func=_cmd_export_vcf)
add_diploid_parx_genome(P_export_vcf)


# THetA special case: takes tumor .cns and normal .cnr or reference.cnn
def _cmd_export_theta(args):
    """Convert segments to THetA2 input file format (*.input)."""
    tumor_cn = read_cna(args.tumor_segment)
    normal_cn = read_cna(args.reference) if args.reference else None
    table = export.export_theta(tumor_cn, normal_cn)
    if not args.output:
        args.output = tumor_cn.sample_id + ".interval_count"
    table.to_csv(args.output, sep="\t", index=False)
    logging.info("Wrote %s", args.output)
    if args.vcf:
        variants = load_het_snps(
            args.vcf,
            args.sample_id,  # or tumor_cn.sample_id,
            args.normal_id,
            args.min_variant_depth,
            args.zygosity_freq,
        )
        if not len(variants):
            raise ValueError("VCF contains no usable SNV records")
        try:
            tumor_snps, normal_snps = export.export_theta_snps(variants)
        except ValueError as exc:
            raise ValueError(
                "VCF does not contain any tumor-normal paired samples"
            ) from exc
        for title, table in [("tumor", tumor_snps), ("normal", normal_snps)]:
            out_fname = f"{tumor_cn.sample_id}.{title}.snp_formatted.txt"
            table.to_csv(out_fname, sep="\t", index=False)
            logging.info("Wrote %s", out_fname)


P_export_theta = P_export_subparsers.add_parser("theta", help=_cmd_export_theta.__doc__)
P_export_theta.add_argument(
    "tumor_segment", help="""Tumor-sample segmentation file from CNVkit (.cns)."""
)
P_export_theta.add_argument(
    "-r",
    "--reference",
    help="""Reference copy number profile (.cnn), or normal-sample bin-level log2 copy
            ratios (.cnr). Use if the tumor_segment input file does not contain a
            "weight" column.""",
)
P_export_theta.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)

P_extheta_vcf = P_export_theta.add_argument_group(
    "To also output tables of SNP b-allele frequencies for THetA2"
)
P_extheta_vcf.add_argument(
    "-v",
    "--vcf",
    help="""VCF file containing SNVs observed in both the tumor and normal samples.
            Tumor sample ID should match the `tumor_segment` filename or be specified
            with -i/--sample-id.""",
)
P_extheta_vcf.add_argument(
    "-i",
    "--sample-id",
    help="""Specify the name of the tumor sample in the VCF (given with -v/--vcf).
            [Default: taken the tumor_segment file name]""",
)
P_extheta_vcf.add_argument(
    "-n", "--normal-id", help="Corresponding normal sample ID in the input VCF."
)
P_extheta_vcf.add_argument(
    "-m",
    "--min-variant-depth",
    type=int,
    default=20,
    help="""Minimum read depth for a SNP in the VCF to be counted.
            [Default: %(default)s]""",
)
P_extheta_vcf.add_argument(
    "-z",
    "--zygosity-freq",
    metavar="ALT_FREQ",
    nargs="?",
    type=float,
    const=0.25,
    help="""Ignore VCF's genotypes (GT field) and instead infer zygosity from allele
            frequencies. [Default if used without a number: %(const)s]""",
)

P_export_theta.set_defaults(func=_cmd_export_theta)


# Nexus "basic" special case: can only represent 1 sample
def _cmd_export_nb(args):
    """Convert bin-level log2 ratios to Nexus Copy Number "basic" format."""
    cnarr = read_cna(args.filename)
    table = export.export_nexus_basic(cnarr)
    write_dataframe(args.output, table)


P_export_nb = P_export_subparsers.add_parser("nexus-basic", help=_cmd_export_nb.__doc__)
P_export_nb.add_argument(
    "filename",
    help="""Log2 copy ratio data file (*.cnr), the output of the 'fix' sub-command.""",
)
P_export_nb.add_argument("-o", "--output", metavar="FILENAME", help="Output file name.")
P_export_nb.set_defaults(func=_cmd_export_nb)


# Nexus "Custom-OGT" special case: can only represent 1 sample
def _cmd_export_nbo(args):
    """Convert log2 ratios and b-allele freqs to Nexus "Custom-OGT" format."""
    cnarr = read_cna(args.filename)
    varr = load_het_snps(
        args.vcf,
        args.sample_id,
        args.normal_id,
        args.min_variant_depth,
        args.zygosity_freq,
    )
    table = export.export_nexus_ogt(cnarr, varr, args.min_weight)
    write_dataframe(args.output, table)


P_export_nbo = P_export_subparsers.add_parser("nexus-ogt", help=_cmd_export_nbo.__doc__)
P_export_nbo.add_argument(
    "filename",
    help="""Log2 copy ratio data file (*.cnr), the output of the 'fix' sub-command.""",
)
P_export_nbo.add_argument(
    "vcf",
    help="""VCF of SNVs for the same sample, to calculate b-allele frequencies.""",
)
P_export_nbo.add_argument(
    "-i",
    "--sample-id",
    help="""Specify the name of the sample in the VCF to use to extract b-allele
            frequencies.""",
)
P_export_nbo.add_argument(
    "-n", "--normal-id", help="Corresponding normal sample ID in the input VCF."
)
P_export_nbo.add_argument(
    "-m",
    "--min-variant-depth",
    type=int,
    default=20,
    help="""Minimum read depth for a SNV to be included in the b-allele frequency
            calculation. [Default: %(default)s]""",
)
P_export_nbo.add_argument(
    "-z",
    "--zygosity-freq",
    metavar="ALT_FREQ",
    nargs="?",
    type=float,
    const=0.25,
    help="""Ignore VCF's genotypes (GT field) and instead infer zygosity from allele
            frequencies. [Default if used without a number: %(const)s]""",
)
P_export_nbo.add_argument(
    "-w",
    "--min-weight",
    type=float,
    default=0.0,
    help="""Minimum weight (between 0 and 1) for a bin to be included in the output.
            [Default: %(default)s]""",
)
P_export_nbo.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)
P_export_nbo.set_defaults(func=_cmd_export_nbo)


# All else: export any number of .cnr or .cns files


def _cmd_export_cdt(args):
    """Convert log2 ratios to CDT format. Compatible with Java TreeView."""
    sample_ids = list(map(core.fbase, args.filenames))
    table = export.merge_samples(args.filenames)
    formatter = export.EXPORT_FORMATS["cdt"]
    outheader, outrows = formatter(sample_ids, table)
    write_tsv(args.output, outrows, colnames=outheader)


P_export_cdt = P_export_subparsers.add_parser("cdt", help=_cmd_export_cdt.__doc__)
P_export_cdt.add_argument(
    "filenames",
    nargs="+",
    help="""Log2 copy ratio data file(s) (*.cnr), the output of the 'fix'
            sub-command.""",
)
P_export_cdt.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)
P_export_cdt.set_defaults(func=_cmd_export_cdt)


def _cmd_export_jtv(args):
    """Convert log2 ratios to Java TreeView's native format."""
    sample_ids = list(map(core.fbase, args.filenames))
    table = export.merge_samples(args.filenames)
    formatter = export.EXPORT_FORMATS["jtv"]
    outheader, outrows = formatter(sample_ids, table)
    write_tsv(args.output, outrows, colnames=outheader)


P_export_jtv = P_export_subparsers.add_parser("jtv", help=_cmd_export_jtv.__doc__)
P_export_jtv.add_argument(
    "filenames",
    nargs="+",
    help="""Log2 copy ratio data file(s) (*.cnr), the output of the 'fix'
            sub-command.""",
)
P_export_jtv.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)
P_export_jtv.set_defaults(func=_cmd_export_jtv)


def _cmd_export_gistic(args):
    formatter = export.EXPORT_FORMATS["gistic"]
    outdf = formatter(args.filenames)
    write_dataframe(args.output, outdf)


P_export_gistic = P_export_subparsers.add_parser(
    "gistic", help=_cmd_export_gistic.__doc__
)
P_export_gistic.add_argument(
    "filenames",
    nargs="+",
    help="""Log2 copy ratio data file(s) (*.cnr), the output of the 'fix'
            sub-command.""",
)
P_export_gistic.add_argument(
    "-o", "--output", metavar="FILENAME", help="Output file name."
)
P_export_gistic.set_defaults(func=_cmd_export_gistic)


# version ---------------------------------------------------------------------


def print_version(_args):
    """Display this program's version."""
    print(__version__)


P_version = AP_subparsers.add_parser("version", help=print_version.__doc__)
P_version.set_defaults(func=print_version)


# _____________________________________________________________________________
# Shim for command-line execution


def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)
