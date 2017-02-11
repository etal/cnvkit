"""Command-line interface and corresponding API for CNVkit."""
# NB: argparse CLI definitions and API functions are interwoven:
#   "_cmd_*" handles I/O and arguments processing for the command
#   "do_*" runs the command's functionality as an API
from __future__ import absolute_import, division, print_function
from builtins import map, zip

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

# If running headless, use a suitable GUI-less plotting backend
if not os.environ.get('DISPLAY'):
    import matplotlib
    matplotlib.use("Agg", force=True)

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
pyplot.ioff()

from . import (access, antitarget, batch, call, core, coverage, descriptives,
               diagram, export, fix, heatmap, importers, metrics, parallel,
               reference, reports, scatter, segmentation, tabio, target)
from .cnary import CopyNumArray as _CNA
from .genome import GenomicArray as _GA
from ._version import __version__


__all__ = []
def public(fn):
    __all__.append(fn.__name__)
    return fn


AP = argparse.ArgumentParser(
        description="CNVkit, a command-line toolkit for copy number analysis.",
        epilog="Contact Eric Talevich <eric.talevich@ucsf.edu> for help.")
AP_subparsers = AP.add_subparsers(
        help="Sub-commands (use with -h for more info)")


# _____________________________________________________________________________
# Core pipeline

# batch -----------------------------------------------------------------------

def _cmd_batch(args):
    """Run the complete CNVkit pipeline on one or more BAM files."""
    # Validate/restrict options, beyond what argparse mutual exclusion can do
    if args.reference:
        bad_flags = [flag
                     for is_used, flag in (
                         (args.normal is not None,  '-n/--normal'),
                         (args.fasta,               '-f/--fasta'),
                         (args.targets,             '-t/--targets'),
                         (args.antitargets,         '-a/--antitargets'),
                         (args.access,              '-g/--access'),
                         (args.annotate,            '--annotate'),
                         (args.short_names,         '--short-names'),
                         (args.split,               '--split'),
                         (args.target_avg_size,     '--target-avg-size'),
                         (args.antitarget_avg_size, '--antitarget-avg-size'),
                         (args.antitarget_min_size, '--antitarget-min-size'),
                     ) if is_used]
        if bad_flags:
            sys.exit("If -r/--reference is given, options to construct a new " +
                     "reference (" + ", ".join(bad_flags) +
                     ") should not be used." +
                     "\n(See: cnvkit.py batch -h)")
    elif (args.normal is None or
          (args.method in ('hybrid', 'amplicon') and not args.targets)):
        sys.exit("Options -n/--normal and -t/--targets (at least) must be "
                 "given to build a new reference if -r/--reference is not used."
                 "\n(See: cnvkit.py batch -h)")

    if args.processes < 1:
        import multiprocessing
        args.processes = multiprocessing.cpu_count()

    # Ensure sample IDs are unique to avoid overwriting outputs
    seen_sids = {}
    for fname in (args.bam_files or []) + (args.normal or []):
        sid = core.fbase(fname)
        if sid in seen_sids:
            sys.exit("Duplicate sample ID %r (from %s and %s)"
                     % (sid, fname, seen_sids[sid]))
        seen_sids[sid] = fname

    if not args.reference:
        # Build a copy number reference; update (anti)targets upon request
        args.reference, args.targets, args.antitargets = batch.batch_make_reference(
            args.normal, args.targets, args.antitargets, args.male_reference,
            args.fasta, args.annotate, args.short_names, args.target_avg_size,
            args.access, args.antitarget_avg_size, args.antitarget_min_size,
            args.output_reference, args.output_dir, args.processes,
            args.count_reads, args.method)
    elif args.targets is None and args.antitargets is None:
        # Extract (anti)target BEDs from the given, existing CN reference
        ref_arr = tabio.read_cna(args.reference)
        targets, antitargets = reference.reference2regions(ref_arr)
        ref_pfx = os.path.join(args.output_dir, core.fbase(args.reference))
        args.targets = ref_pfx + '.target-tmp.bed'
        args.antitargets = ref_pfx + '.antitarget-tmp.bed'
        tabio.write(targets, args.targets, 'bed4')
        tabio.write(antitargets, args.antitargets, 'bed4')

    if args.bam_files:
        if args.processes == 1:
            procs_per_bam = 1
            logging.info("Running %d samples in serial", len(args.bam_files))
        else:
            procs_per_bam = max(1, args.processes // len(args.bam_files))
            logging.info("Running %d samples in %d processes "
                         "(that's %d processes per bam)",
                         len(args.bam_files), args.processes, procs_per_bam)

        with parallel.pick_pool(args.processes) as pool:
            for bam in args.bam_files:
                pool.submit(batch.batch_run_sample,
                            bam, args.targets, args.antitargets, args.reference,
                            args.output_dir, args.male_reference, args.scatter,
                            args.diagram, args.rlibpath, args.count_reads,
                            args.drop_low_coverage, args.method, procs_per_bam)


P_batch = AP_subparsers.add_parser('batch', help=_cmd_batch.__doc__)
P_batch.add_argument('bam_files', nargs='*',
        help="Mapped sequence reads (.bam)")
P_batch.add_argument('-m', '--method',
        choices=('hybrid', 'amplicon', 'wgs'), default='hybrid',
        help="""Sequencing protocol: hybridization capture ('hybrid'), targeted
                amplicon sequencing ('amplicon'), or whole genome sequencing
                ('wgs'). Determines whether and how to use antitarget bins.
                [Default: %(default)s]""")
P_batch.add_argument('-y', '--male-reference', action='store_true',
        help="""Use or assume a male reference (i.e. female samples will have +1
                log-CNR of chrX; otherwise male samples would have -1 chrX).""")
P_batch.add_argument('-c', '--count-reads', action='store_true',
        help="""Get read depths by counting read midpoints within each bin.
                (An alternative algorithm).""")
P_batch.add_argument("--drop-low-coverage", action='store_true',
        help="""Drop very-low-coverage bins before segmentation to avoid
                false-positive deletions in poor-quality tumor samples.""")
P_batch.add_argument('-p', '--processes',
        nargs='?', type=int, const=0, default=1,
        help="""Number of subprocesses used to running each of the BAM files in
                parallel. Without an argument, use the maximum number of
                available CPUs. [Default: process each BAM in serial]""")
P_batch.add_argument("--rlibpath",
        help="Path to an alternative site-library to use for R packages.")

# Reference-building options
P_batch_newref = P_batch.add_argument_group(
    "To construct a new copy number reference")
P_batch_newref.add_argument('-n', '--normal', nargs='*',
        help="""Normal samples (.bam) to construct the pooled reference.
                If this option is used but no files are given, a "flat"
                reference will be built.""")
P_batch_newref.add_argument('-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
P_batch_newref.add_argument('-t', '--targets', #required=True,
        help="Target intervals (.bed or .list)")
P_batch_newref.add_argument('-a', '--antitargets', #required=True,
        help="Antitarget intervals (.bed or .list)")
# For pre-processing targets
P_batch_newref.add_argument('--annotate',
        help="""Use gene models from this file to assign names to the target
                regions. Format: UCSC refFlat.txt or ensFlat.txt file
                (preferred), or BED, interval list, GFF, or similar.""")
P_batch_newref.add_argument('--short-names', action='store_true',
        help="Reduce multi-accession bait labels to be short and consistent.")
P_batch_newref.add_argument('--split', action='store_true', # DEPRECATED
        help=argparse.SUPPRESS)
P_batch_newref.add_argument('--target-avg-size', type=int,
        help="Average size of split target bins (results are approximate).")
# For antitargets:
P_batch_newref.add_argument('-g', '--access',
        help="""Regions of accessible sequence on chromosomes (.bed), as
                output by the 'access' command.""")
P_batch_newref.add_argument('--antitarget-avg-size', type=int,
        help="Average size of antitarget bins (results are approximate).")
P_batch_newref.add_argument('--antitarget-min-size', type=int,
        help="Minimum size of antitarget bins (smaller regions are dropped).")
P_batch_newref.add_argument('--output-reference',
        help="""Output filename/path for the new reference file being created.
                (If given, ignores the -o/--output-dir option and will write the
                file to the given path. Otherwise, \"reference.cnn\" will be
                created in the current directory or specified output directory.)
                """)

P_batch_oldref = P_batch.add_argument_group("To reuse an existing reference")
P_batch_oldref.add_argument('-r', '--reference', #required=True,
        help="Copy number reference file (.cnn).")

# Reporting options
P_batch_report = P_batch.add_argument_group("Output options")
P_batch_report.add_argument('-d', '--output-dir', default='.',
        help="Output directory.")
P_batch_report.add_argument('--scatter', action='store_true',
        help="Create a whole-genome copy ratio profile as a PDF scatter plot.")
P_batch_report.add_argument('--diagram', action='store_true',
        help="Create an ideogram of copy ratios on chromosomes as a PDF.")

P_batch.set_defaults(func=_cmd_batch)


# target ----------------------------------------------------------------------

do_target = public(target.do_target)


def _cmd_target(args):
    """Transform bait intervals into targets more suitable for CNVkit."""
    regions = tabio.read_auto(args.interval)
    regions = target.do_target(regions, args.annotate, args.short_names,
                               args.split, args.avg_size)
    tabio.write(regions, args.output, "bed4")


P_target = AP_subparsers.add_parser('target', help=_cmd_target.__doc__)
P_target.add_argument('interval',
        help="""BED or interval file listing the targeted regions.""")
P_target.add_argument('--annotate',
        help="""Use gene models from this file to assign names to the target
                regions. Format: UCSC refFlat.txt or ensFlat.txt file
                (preferred), or BED, interval list, GFF, or similar.""")
P_target.add_argument('--short-names', action='store_true',
        help="Reduce multi-accession bait labels to be short and consistent.")
P_target.add_argument('--split', action='store_true',
        help="Split large tiled intervals into smaller, consecutive targets.")
# Exons: [114--188==203==292--21750], mean=353 -> outlier=359, extreme=515
#   NV2:  [65--181==190==239--12630], mean=264 -> outlier=277, extreme=364
# Default avg_size chosen s.t. minimum bin size after split is ~= median
P_target.add_argument('-a', '--avg-size', type=int, default=200 / .75,
        help="""Average size of split target bins (results are approximate).
                [Default: %(default)s]""")
P_target.add_argument('-o', '--output', help="""Output file name.""")
P_target.set_defaults(func=_cmd_target)


# access ----------------------------------------------------------------------

do_access = public(access.do_access)


def _cmd_access(args):
    """List the locations of accessible sequence regions in a FASTA file."""
    access_arr = access.do_access(args.fa_fname, args.exclude,
                                  args.min_gap_size)
    tabio.write(access_arr, args.output, "bed3")


P_access = AP_subparsers.add_parser('access', help=_cmd_access.__doc__)
P_access.add_argument("fa_fname",
                help="Genome FASTA file name")
P_access.add_argument("-s", "--min-gap-size", type=int, default=5000,
                help="""Minimum gap size between accessible sequence
                regions.  Regions separated by less than this distance will
                be joined together. [Default: %(default)s]""")
P_access.add_argument("-x", "--exclude", action="append", default=[],
                help="""Additional regions to exclude, in BED format. Can be
                used multiple times.""")
P_access.add_argument("-o", "--output",
                type=argparse.FileType('w'), default=sys.stdout,
                help="Output file name")
P_access.set_defaults(func=_cmd_access)


# antitarget ------------------------------------------------------------------

do_antitarget = public(antitarget.do_antitarget)

def _cmd_antitarget(args):
    """Derive a background/antitarget BED file from a target BED file."""
    targets = tabio.read_auto(args.targets)
    access = tabio.read_auto(args.access)
    out_arr = antitarget.do_antitarget(targets, access, args.avg_size,
                                       args.min_size)
    if not args.output:
        base, ext = args.interval.rsplit('.', 1)
        args.output = base + '.antitarget.' + ext
    tabio.write(out_arr, args.output, "bed4")


P_anti = AP_subparsers.add_parser('antitarget', help=_cmd_antitarget.__doc__)
P_anti.add_argument('targets',
        help="""BED or interval file listing the targeted regions.""")
P_anti.add_argument('-g', '--access',
        help="""Regions of accessible sequence on chromosomes (.bed), as
                output by genome2access.py.""")
P_anti.add_argument('-a', '--avg-size', type=int, default=150000,
        help="""Average size of antitarget bins (results are approximate).
                [Default: %(default)s]""")
P_anti.add_argument('-m', '--min-size', type=int,
        help="""Minimum size of antitarget bins (smaller regions are dropped).
                [Default: 1/16 avg size, calculated]""")
P_anti.add_argument('-o', '--output', help="""Output file name.""")
P_anti.set_defaults(func=_cmd_antitarget)


# coverage --------------------------------------------------------------------

do_coverage = public(coverage.do_coverage)


def _cmd_coverage(args):
    """Calculate coverage in the given regions from BAM read depths."""
    pset = coverage.do_coverage(args.interval, args.bam_file, args.count,
                                args.min_mapq, args.processes)
    if not args.output:
        # Create an informative but unique name for the coverage output file
        bambase = core.fbase(args.bam_file)
        bedbase = core.fbase(args.interval)
        tgtbase = ('antitargetcoverage'
                   if 'anti' in bedbase.lower()
                   else 'targetcoverage')
        args.output = '%s.%s.cnn' % (bambase, tgtbase)
        if os.path.exists(args.output):
            args.output = '%s.%s.cnn' % (bambase, bedbase)
    core.ensure_path(args.output)
    tabio.write(pset, args.output)


P_coverage = AP_subparsers.add_parser('coverage', help=_cmd_coverage.__doc__)
P_coverage.add_argument('bam_file', help="Mapped sequence reads (.bam)")
P_coverage.add_argument('interval', help="Intervals (.bed or .list)")
P_coverage.add_argument('-c', '--count', action='store_true',
        help="""Get read depths by counting read midpoints within each bin.
                (An alternative algorithm).""")
P_coverage.add_argument('-q', '--min-mapq', type=int, default=0,
        help="""Minimum mapping quality score (phred scale 0-60) to count a read
                for coverage depth.  [Default: %(default)s]""")
P_coverage.add_argument('-o', '--output', help="""Output file name.""")
P_coverage.add_argument('-p', '--processes',
        nargs='?', type=int, const=0, default=1,
        help="""Number of subprocesses to calculate coverage in parallel.
                Without an argument, use the maximum number of available CPUs.
                [Default: use 1 process]""")
P_coverage.set_defaults(func=_cmd_coverage)


# reference -------------------------------------------------------------------

do_reference = public(reference.do_reference)
do_reference_flat = public(reference.do_reference_flat)


def _cmd_reference(args):
    """Compile a coverage reference from the given files (normal samples)."""
    usage_err_msg = ("Give .cnn samples OR targets and antitargets.")
    if args.targets and args.antitargets:
        # Flat refence
        assert not args.references, usage_err_msg
        ref_probes = reference.do_reference_flat(args.targets, args.antitargets,
                                                 args.fasta,
                                                 args.male_reference)
    elif args.references:
        # Pooled reference
        assert not args.targets and not args.antitargets, usage_err_msg
        filenames = []
        for path in args.references:
            if os.path.isdir(path):
                filenames.extend(os.path.join(path, f) for f in os.listdir(path)
                                 if f.endswith('targetcoverage.cnn'))
            else:
                filenames.append(path)
        targets = [f for f in filenames if 'antitarget' not in f]
        antitargets = [f for f in filenames if 'antitarget' in f]
        logging.info("Number of target and antitarget files: %d, %d",
                     len(targets), len(antitargets))
        female_samples = ((args.sample_sex.lower() not in ['y', 'm', 'male'])
                          if args.sample_sex else None)
        ref_probes = reference.do_reference(targets, antitargets, args.fasta,
                                            args.male_reference, female_samples,
                                            args.do_gc, args.do_edge,
                                            args.do_rmask)
    else:
        raise ValueError(usage_err_msg)

    ref_fname = args.output or "cnv_reference.cnn"
    core.ensure_path(ref_fname)
    tabio.write(ref_probes, ref_fname)


P_reference = AP_subparsers.add_parser('reference', help=_cmd_reference.__doc__)
P_reference.add_argument('references', nargs='*',
        help="""Normal-sample target or antitarget .cnn files, or the
                directory that contains them.""")
P_reference.add_argument('-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
P_reference.add_argument('-o', '--output', help="Output file name.")
P_reference.add_argument('-y', '--male-reference', action='store_true',
        help="""Create a male reference: shift female samples' chrX
                log-coverage by -1, so the reference chrX average is -1.
                Otherwise, shift male samples' chrX by +1, so the reference chrX
                average is 0.""")
P_reference.add_argument('-x', '--sample-sex', '-g', '--gender',
        dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the chromosomal sex of all given samples as male or
                female. (Default: guess each sample from coverage of X and Y
                chromosomes).""")

P_reference_flat = P_reference.add_argument_group(
    "To construct a generic, \"flat\" copy number reference with neutral "
    "expected coverage")
P_reference_flat.add_argument('-t', '--targets',
        help="Target intervals (.bed or .list)")
P_reference_flat.add_argument('-a', '--antitargets',
        help="Antitarget intervals (.bed or .list)")

P_reference_bias = P_reference.add_argument_group(
    "To disable specific automatic bias corrections")
P_reference_bias.add_argument('--no-gc', dest='do_gc', action='store_false',
        help="Skip GC correction.")
P_reference_bias.add_argument('--no-edge', dest='do_edge', action='store_false',
        help="Skip edge-effect correction.")
P_reference_bias.add_argument('--no-rmask', dest='do_rmask', action='store_false',
        help="Skip RepeatMasker correction.")
P_reference.set_defaults(func=_cmd_reference)


# fix -------------------------------------------------------------------------

do_fix = public(fix.do_fix)


def _cmd_fix(args):
    """Combine target and antitarget coverages and correct for biases.

    Adjust raw coverage data according to the given reference, correct potential
    biases and re-center.
    """
    # Verify that target and antitarget are from the same sample
    tgt_raw = tabio.read_cna(args.target)
    anti_raw = tabio.read_cna(args.antitarget)
    if len(anti_raw) and tgt_raw.sample_id != anti_raw.sample_id:
        raise ValueError("Sample IDs do not match:"
                         "'%s' (target) vs. '%s' (antitarget)"
                         % (tgt_raw.sample_id, anti_raw.sample_id))
    target_table = fix.do_fix(tgt_raw, anti_raw,
                              tabio.read_cna(args.reference),
                              args.do_gc, args.do_edge, args.do_rmask)
    tabio.write(target_table, args.output or tgt_raw.sample_id + '.cnr')


P_fix = AP_subparsers.add_parser('fix', help=_cmd_fix.__doc__)
P_fix.add_argument('target',
        help="Target coverage file (.targetcoverage.cnn).")
P_fix.add_argument('antitarget',
        help="Antitarget coverage file (.antitargetcoverage.cnn).")
P_fix.add_argument('reference',
        help="Reference coverage (.cnn).")
# P_fix.add_argument('--do-gc', action='store_true', default=True,
#         help="Do GC correction.")
# P_fix.add_argument('--do-edge', action='store_true',
#         help="Do edge-effect correction.")
# P_fix.add_argument('--do-size', action='store_true',
#         help="Do interval-size correction.")
P_fix.add_argument('--no-gc', dest='do_gc', action='store_false',
        help="Skip GC correction.")
P_fix.add_argument('--no-edge', dest='do_edge', action='store_false',
        help="Skip edge-effect correction.")
P_fix.add_argument('--no-rmask', dest='do_rmask', action='store_false',
        help="Skip RepeatMasker correction.")
P_fix.add_argument('-o', '--output',
        help="Output file name.")
P_fix.set_defaults(func=_cmd_fix)


# segment ---------------------------------------------------------------------

do_segmentation = public(segmentation.do_segmentation)


def _cmd_segment(args):
    """Infer copy number segments from the given coverage table."""
    cnarr = tabio.read_cna(args.filename)
    variants = load_het_snps(args.vcf, args.sample_id, args.normal_id,
                             args.min_variant_depth, args.zygosity_freq)
    results = segmentation.do_segmentation(cnarr, args.method, args.threshold,
                                           variants=variants,
                                           skip_low=args.drop_low_coverage,
                                           skip_outliers=args.drop_outliers,
                                           save_dataframe=bool(args.dataframe),
                                           rlibpath=args.rlibpath,
                                           processes=args.processes)
    if args.dataframe:
        segments, dframe = results
        with open(args.dataframe, 'w') as handle:
            handle.write(dframe)
        logging.info("Wrote %s", args.dataframe)
    else:
        segments = results
    tabio.write(segments, args.output or segments.sample_id + '.cns')


P_segment = AP_subparsers.add_parser('segment', help=_cmd_segment.__doc__)
P_segment.add_argument('filename',
        help="Bin-level log2 ratios (.cnr file), as produced by 'fix'.")
P_segment.add_argument('-o', '--output',
        help="Output table file name (CNR-like table of segments, .cns).")
P_segment.add_argument('-d', '--dataframe',
        help="""File name to save the raw R dataframe emitted by CBS or
                Fused Lasso. (Useful for debugging.)""")
P_segment.add_argument('-m', '--method',
        choices=('cbs', 'haar', 'flasso'), default='cbs',
        help="""Segmentation method (CBS, HaarSeg, or Fused Lasso).
                [Default: %(default)s]""")
P_segment.add_argument('-t', '--threshold', type=float,
        help="""Significance threshold (p-value or FDR, depending on method) to
                accept breakpoints during segmentation.""")
P_segment.add_argument("--drop-low-coverage", action='store_true',
        help="""Drop very-low-coverage bins before segmentation to avoid
                false-positive deletions in poor-quality tumor samples.""")
P_segment.add_argument("--drop-outliers",
        type=float, default=10, metavar="FACTOR",
        help="""Drop outlier bins more than this many multiples of the 95th
                quantile away from the average within a rolling window.
                Set to 0 for no outlier filtering.
                [Default: %(default)g]""")
P_segment.add_argument("--rlibpath",
        help="Path to an alternative site-library to use for R packages.")
P_segment.add_argument('-p', '--processes',
        nargs='?', type=int, const=0, default=1,
        help="""Number of subprocesses to segment in parallel.
                Give 0 or a negative value to use the maximum number
                of available CPUs. [Default: use 1 process]""")

P_segment_vcf = P_segment.add_argument_group(
    "To additionally segment SNP b-allele frequencies")
P_segment_vcf.add_argument('-v', '--vcf', metavar="FILENAME",
        help="""VCF file name containing variants for segmentation by allele
                frequencies.""")
P_segment_vcf.add_argument('-i', '--sample-id',
        help="""Specify the name of the sample in the VCF (-v/--vcf) to use for
                b-allele frequency extraction and as the default plot title.""")
P_segment_vcf.add_argument('-n', '--normal-id',
        help="""Corresponding normal sample ID in the input VCF (-v/--vcf).
                This sample is used to select only germline SNVs to plot
                b-allele frequencies.""")
P_segment_vcf.add_argument('--min-variant-depth', type=int, default=20,
        help="""Minimum read depth for a SNV to be displayed in the b-allele
                frequency plot. [Default: %(default)s]""")
P_segment_vcf.add_argument('-z', '--zygosity-freq',
        metavar='ALT_FREQ', nargs='?', type=float, const=0.25,
        help="""Ignore VCF's genotypes (GT field) and instead infer zygosity
                from allele frequencies.  [Default if used without a number:
                %(const)s]""")

P_segment.set_defaults(func=_cmd_segment)


# call ------------------------------------------------------------------------

do_call = public(call.do_call)


def _cmd_call(args):
    """Call copy number variants from segmented log2 ratios."""
    if args.purity and not 0.0 < args.purity <= 1.0:
        raise RuntimeError("Purity must be between 0 and 1.")

    cnarr = tabio.read_cna(args.filename)
    if args.center:
        cnarr.center_all(args.center)
    varr = load_het_snps(args.vcf, args.sample_id, args.normal_id,
                         args.min_variant_depth, args.zygosity_freq)
    is_sample_female = (verify_sample_sex(cnarr, args.sample_sex,
                                          args.male_reference)
                        if args.purity and args.purity < 1.0
                        else None)
    cnarr = call.do_call(cnarr, varr, args.method, args.ploidy, args.purity,
                         args.male_reference, is_sample_female, args.filters,
                         args.thresholds)
    tabio.write(cnarr, args.output or cnarr.sample_id + '.call.cns')


def csvstring(text):
    return tuple(map(float, text.split(",")))


P_call = AP_subparsers.add_parser('call', help=_cmd_call.__doc__)
P_call.add_argument('filename',
        help="Copy ratios (.cnr or .cns).")
P_call.add_argument("--center",
        choices=('mean', 'median', 'mode', 'biweight'),
        help="""Re-center the log2 ratio values using this estimate of the
                center or average value.""")
P_call.add_argument('--filter', action='append', default=[], dest='filters',
        choices=('ampdel', 'cn', 'ci', 'sem', # 'bic'
                ),
        help="""Merge segments flagged by the specified filter(s) with the
                adjacent segment(s).""")
P_call.add_argument('-m', '--method',
        choices=('threshold', 'clonal', 'none'), default='threshold',
        help="""Calling method. [Default: %(default)s]""")
P_call.add_argument('-t', '--thresholds',
        type=csvstring, default="-1.1,-0.25,0.2,0.7",
        help="""Hard thresholds for calling each integer copy number, separated
                by commas. Use the '=' sign on the command line, e.g.: -t=-1,0,1
                [Default: %(default)s]""")
P_call.add_argument("--ploidy", type=int, default=2,
        help="Ploidy of the sample cells. [Default: %(default)d]")
P_call.add_argument("--purity", type=float,
        help="Estimated tumor cell fraction, a.k.a. purity or cellularity.")
P_call.add_argument('-x', '--sample-sex', '-g', '--gender', dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the sample's chromosomal sex as male or female.
                (Otherwise guessed from X and Y coverage).""")
P_call.add_argument('-y', '--male-reference', action='store_true',
        help="""Was a male reference used?  If so, expect half ploidy on
                chrX and chrY; otherwise, only chrY has half ploidy.  In CNVkit,
                if a male reference was used, the "neutral" copy number (ploidy)
                of chrX is 1; chrY is haploid for either reference sex.""")
P_call.add_argument('-o', '--output',
        help="Output table file name (CNR-like table of segments, .cns).")

P_call_vcf = P_call.add_argument_group(
    "To additionally process SNP b-allele frequencies for allelic copy number")
P_call_vcf.add_argument('-v', '--vcf', metavar="FILENAME",
        help="""VCF file name containing variants for calculation of b-allele
                frequencies.""")
P_call_vcf.add_argument('-i', '--sample-id',
        help="""Name of the sample in the VCF (-v/--vcf) to use for b-allele
                frequency extraction.""")
P_call_vcf.add_argument('-n', '--normal-id',
        help="""Corresponding normal sample ID in the input VCF (-v/--vcf).
                This sample is used to select only germline SNVs to calculate
                b-allele frequencies.""")
P_call_vcf.add_argument('--min-variant-depth', type=int, default=20,
        help="""Minimum read depth for a SNV to be used in the b-allele
                frequency calculation. [Default: %(default)s]""")
P_call_vcf.add_argument('-z', '--zygosity-freq',
        metavar='ALT_FREQ', nargs='?', type=float, const=0.25,
        help="""Ignore VCF's genotypes (GT field) and instead infer zygosity
                from allele frequencies.  [Default if used without a number:
                %(const)s]""")

P_call.set_defaults(func=_cmd_call)


# _____________________________________________________________________________
# Plots and graphics

# diagram ---------------------------------------------------------------------

def _cmd_diagram(args):
    """Draw copy number (log2 coverages, CBS calls) on chromosomes as a diagram.

    If both the raw probes and segments are given, show them side-by-side on
    each chromosome (segments on the left side, probes on the right side).
    """
    if not args.filename and not args.segment:
        raise ValueError("Must specify a filename as an argument or with "
                         "the '-s' option, or both. You did neither.")

    cnarr = tabio.read_cna(args.filename) if args.filename else None
    segarr = tabio.read_cna(args.segment) if args.segment else None
    if args.adjust_xy:
        is_sample_female = verify_sample_sex(cnarr or segarr, args.sample_sex,
                                             args.male_reference)
        if cnarr:
            cnarr = cnarr.shift_xx(args.male_reference, is_sample_female)
        if segarr:
            segarr = segarr.shift_xx(args.male_reference, is_sample_female)
    outfname = diagram.create_diagram(cnarr, segarr, args.threshold,
                                      args.min_probes, args.output)
    logging.info("Wrote %s", outfname)


P_diagram = AP_subparsers.add_parser('diagram', help=_cmd_diagram.__doc__)
P_diagram.add_argument('filename', nargs='?',
        help="""Processed coverage data file (*.cnr), the output of the
                'fix' sub-command.""")
P_diagram.add_argument('-s', '--segment',
        help="Segmentation calls (.cns), the output of the 'segment' command.")
P_diagram.add_argument('-t', '--threshold', type=float, default=0.5,
        help="""Copy number change threshold to label genes.
                [Default: %(default)s]""")
P_diagram.add_argument('-m', '--min-probes', type=int, default=3,
        help="""Minimum number of covered probes to label a gene.
                [Default: %(default)d]""")
P_diagram.add_argument('-y', '--male-reference', action='store_true',
        help="""Assume inputs are already corrected against a male
                reference (i.e. female samples will have +1 log-CNR of
                chrX; otherwise male samples would have -1 chrX).""")
P_diagram.add_argument('-x', '--sample-sex', '-g', '--gender',
        dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the sample's chromosomal sex as male or female.
                (Otherwise guessed from X and Y coverage).""")
P_diagram.add_argument('--no-shift-xy', dest='adjust_xy', action='store_false',
        help="Don't adjust the X and Y chromosomes according to sample sex.")
P_diagram.add_argument('-o', '--output',
        help="Output PDF file name.")
P_diagram.set_defaults(func=_cmd_diagram)


# scatter ---------------------------------------------------------------------

do_scatter = public(scatter.do_scatter)


def _cmd_scatter(args):
    """Plot probe log2 coverages and segmentation calls together."""
    cnarr = tabio.read_cna(args.filename, sample_id=args.sample_id
                          ) if args.filename else None
    segarr = tabio.read_cna(args.segment, sample_id=args.sample_id
                           ) if args.segment else None
    varr = load_het_snps(args.vcf, args.sample_id, args.normal_id,
                         args.min_variant_depth, args.zygosity_freq)
    if args.range_list:
        with PdfPages(args.output) as pdf_out:
            for region in tabio.read_auto(args.range_list).coords():
                try:
                    scatter.do_scatter(cnarr, segarr, varr, region, None,
                                       args.background_marker, args.trend,
                                       args.width, args.y_min, args.y_max,
                                       ("%s %s" % (args.title,
                                                   region.chromosome)
                                        if args.title else None),
                                       args.segment_color)
                except ValueError as exc:
                    # Probably no bins in the selected region
                    logging.warn("Not plotting region %r: %s",
                                 _GA.row2label(region), exc)
                pdf_out.savefig()
                pyplot.close()
    else:
        scatter.do_scatter(cnarr, segarr, varr, args.chromosome, args.gene,
                           args.background_marker, args.trend, args.width,
                           args.y_min, args.y_max, args.title,
                           args.segment_color)
        if args.output:
            oformat = os.path.splitext(args.output)[-1].replace(".", "")
            pyplot.savefig(args.output, format=oformat, bbox_inches="tight")
            logging.info("Wrote %s", args.output)
        else:
            pyplot.show()


P_scatter = AP_subparsers.add_parser('scatter', help=_cmd_scatter.__doc__)
P_scatter.add_argument('filename', nargs="?",
        help="""Processed bin-level copy ratios (*.cnr), the output
                of the 'fix' sub-command.""")
P_scatter.add_argument('-s', '--segment', metavar="FILENAME",
        help="Segmentation calls (.cns), the output of the 'segment' command.")
P_scatter.add_argument('-c', '--chromosome', metavar="RANGE",
        help="""Chromosome or chromosomal range, e.g. 'chr1' or
                'chr1:2333000-2444000', to display. If a range is given,
                all targeted genes in this range will be shown, unless
                -g/--gene is also given.""")
P_scatter.add_argument('-g', '--gene',
        help="Name of gene or genes (comma-separated) to display.")
P_scatter.add_argument('-l', '--range-list',
        help="""File listing the chromosomal ranges to display, as BED, interval
                list or 'chr:start-end' text. Creates focal plots similar to
                -c/--chromosome for each listed region, combined into a
                multi-page PDF.  The output filename must also be
                specified (-o/--output).""")
P_scatter.add_argument('-w', '--width', type=float, default=1e6,
        help="""Width of margin to show around the selected gene or chromosomal
                region (-g/--gene or -c/--chromosome). [Default: %(default)d]
                """)
P_scatter.add_argument('-o', '--output', metavar="FILENAME",
        help="Output PDF file name.")

P_scatter_aes = P_scatter.add_argument_group("Plot aesthetics")
P_scatter_aes.add_argument('-b', '--background-marker', metavar='CHARACTER',
        default=None,
        help="""Plot antitargets using this symbol when plotting in a selected
                chromosomal region (-g/--gene or -c/--chromosome).
                [Default: same as targets]""")
P_scatter_aes.add_argument('--segment-color', default=scatter.SEG_COLOR,
        help="""Plot segment lines in this color. Value can be any string
                accepted by matplotlib, e.g. 'red' or '#CC0000'.""")
P_scatter_aes.add_argument('--title',
        help="Plot title. [Default: sample ID, from filename or -i]")
P_scatter_aes.add_argument('-t', '--trend', action='store_true',
        help="Draw a smoothed local trendline on the scatter plot.")
P_scatter_aes.add_argument('--y-max', type=float, help="y-axis upper limit.")
P_scatter_aes.add_argument('--y-min', type=float, help="y-axis lower limit.")

P_scatter_vcf = P_scatter.add_argument_group(
    "To plot SNP b-allele frequencies")
P_scatter_vcf.add_argument('-v', '--vcf', metavar="FILENAME",
        help="""VCF file name containing variants to plot for SNV b-allele
                frequencies.""")
P_scatter_vcf.add_argument('-i', '--sample-id',
        help="""Name of the sample in the VCF to use for b-allele frequency
                extraction and as the default plot title.""")
P_scatter_vcf.add_argument('-n', '--normal-id',
        help="""Corresponding normal sample ID in the input VCF. This sample is
                used to select only germline SNVs to plot.""")
P_scatter_vcf.add_argument('-m', '--min-variant-depth', type=int, default=20,
        help="""Minimum read depth for a SNV to be used in the b-allele
                frequency calculation. [Default: %(default)s]""")
P_scatter_vcf.add_argument('-z', '--zygosity-freq',
        metavar='ALT_FREQ', nargs='?', type=float, const=0.25,
        help="""Ignore VCF's genotypes (GT field) and instead infer zygosity
                from allele frequencies.  [Default if used without a number:
                %(const)s]""")

P_scatter.set_defaults(func=_cmd_scatter)


# heatmap ---------------------------------------------------------------------

do_heatmap = public(heatmap.do_heatmap)


def _cmd_heatmap(args):
    """Plot copy number for multiple samples as a heatmap."""
    cnarrs = []
    for fname in args.filenames:
        cnarr = tabio.read_cna(fname)
        if args.adjust_xy:
            is_sample_female = verify_sample_sex(cnarr, args.sample_sex,
                                                 args.male_reference)
            cnarr = cnarr.shift_xx(args.male_reference, is_sample_female)
        cnarrs.append(cnarr)
    heatmap.do_heatmap(cnarrs, args.chromosome, args.desaturate)
    if args.output:
        oformat = os.path.splitext(args.output)[-1].replace(".", "")
        pyplot.savefig(args.output, format=oformat, bbox_inches="tight")
        logging.info("Wrote %s", args.output)
    else:
        pyplot.show()


P_heatmap = AP_subparsers.add_parser('heatmap', help=_cmd_heatmap.__doc__)
P_heatmap.add_argument('filenames', nargs='+',
        help="Sample coverages as raw probes (.cnr) or segments (.cns).")
P_heatmap.add_argument('-c', '--chromosome',
        help="""Chromosome (e.g. 'chr1') or chromosomal range (e.g.
                'chr1:2333000-2444000') to display. If a range is given,
                all targeted genes in this range will be shown, unless
                '--gene'/'-g' is already given.""")
# P_heatmap.add_argument('-g', '--gene',
#         help="Name of gene to display.")
P_heatmap.add_argument('-d', '--desaturate', action='store_true',
        help="Tweak color saturation to focus on significant changes.")
P_heatmap.add_argument('-y', '--male-reference', action='store_true',
        help="""Assume inputs are already corrected against a male
                reference (i.e. female samples will have +1 log-CNR of
                chrX; otherwise male samples would have -1 chrX).""")
P_heatmap.add_argument('-x', '--sample-sex', '-g', '--gender',
        dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the chromosomal sex of all given samples as male or
                female. (Default: guess each sample from coverage of X and Y
                chromosomes).""")
P_heatmap.add_argument('--no-shift-xy', dest='adjust_xy', action='store_false',
        help="Don't adjust the X and Y chromosomes according to sample sex.")
P_heatmap.add_argument('-o', '--output',
        help="Output PDF file name.")
P_heatmap.set_defaults(func=_cmd_heatmap)


# _____________________________________________________________________________
# Tabular outputs

# breaks ----------------------------------------------------------------------

def _cmd_breaks(args):
    """List the targeted genes in which a copy number breakpoint occurs."""
    cnarr = tabio.read_cna(args.filename)
    segarr = tabio.read_cna(args.segment)
    bpoints = do_breaks(cnarr, segarr, args.min_probes)
    logging.info("Found %d gene breakpoints", len(bpoints))
    core.write_dataframe(args.output, bpoints)


@public
def do_breaks(probes, segments, min_probes=1):
    """List the targeted genes in which a copy number breakpoint occurs."""
    intervals = reports.get_gene_intervals(probes)
    bpoints = reports.get_breakpoints(intervals, segments, min_probes)
    return pd.DataFrame.from_records(bpoints,
                                     columns=['gene', 'chromosome',
                                              'location', 'change',
                                              'probes_left', 'probes_right'])


P_breaks = AP_subparsers.add_parser('breaks', help=_cmd_breaks.__doc__)
P_breaks.add_argument('filename',
        help="""Processed sample coverage data file (*.cnr), the output
                of the 'fix' sub-command.""")
P_breaks.add_argument('segment',
        help="Segmentation calls (.cns), the output of the 'segment' command).")
P_breaks.add_argument('-m', '--min-probes', type=int, default=1,
        help="""Minimum number of within-gene probes on both sides of a
                breakpoint to report it. [Default: %(default)d]""")
P_breaks.add_argument('-o', '--output',
        help="Output table file name.")
P_breaks.set_defaults(func=_cmd_breaks)


# gainloss --------------------------------------------------------------------

def _cmd_gainloss(args):
    """Identify targeted genes with copy number gain or loss."""
    cnarr = tabio.read_cna(args.filename)
    segarr = tabio.read_cna(args.segment) if args.segment else None
    is_sample_female = verify_sample_sex(cnarr, args.sample_sex,
                                         args.male_reference)
    gainloss = do_gainloss(cnarr, segarr, args.threshold,
                           args.min_probes, args.drop_low_coverage,
                           args.male_reference, is_sample_female)
    logging.info("Found %d gene-level gains and losses", len(gainloss))
    core.write_dataframe(args.output, gainloss)


@public
def do_gainloss(cnarr, segments=None, threshold=0.2, min_probes=3,
                skip_low=False, male_reference=False, is_sample_female=None):
    """Identify targeted genes with copy number gain or loss."""
    if is_sample_female is None:
        is_sample_female = cnarr.guess_xx(male_reference=male_reference)
    cnarr = cnarr.shift_xx(male_reference, is_sample_female)
    if segments:
        segments = segments.shift_xx(male_reference, is_sample_female)
        gainloss = reports.gainloss_by_segment(cnarr, segments, threshold,
                                               skip_low)
    else:
        gainloss = reports.gainloss_by_gene(cnarr, threshold, skip_low)
    gainloss = list(gainloss)
    columns = (gainloss[0].index if len(gainloss) else cnarr._required_columns)
    columns = ["gene"] + [col for col in columns if col != "gene"]
    gainloss = pd.DataFrame.from_records(gainloss).reindex(columns=columns)
    if min_probes and len(gainloss):
        gainloss = gainloss[gainloss.probes >= min_probes]
    return gainloss


P_gainloss = AP_subparsers.add_parser('gainloss', help=_cmd_gainloss.__doc__)
P_gainloss.add_argument('filename',
        help="""Processed sample coverage data file (*.cnr), the output
                of the 'fix' sub-command.""")
P_gainloss.add_argument('-s', '--segment',
        help="Segmentation calls (.cns), the output of the 'segment' command).")
P_gainloss.add_argument('-t', '--threshold', type=float, default=0.2,
        help="""Copy number change threshold to report a gene gain/loss.
                [Default: %(default)s]""")
P_gainloss.add_argument('-m', '--min-probes', type=int, default=3,
        help="""Minimum number of covered probes to report a gain/loss.
                [Default: %(default)d]""")
P_gainloss.add_argument("--drop-low-coverage", action='store_true',
        help="""Drop very-low-coverage bins before segmentation to avoid
                false-positive deletions in poor-quality tumor samples.""")
P_gainloss.add_argument('-y', '--male-reference', action='store_true',
        help="""Assume inputs are already corrected against a male
                reference (i.e. female samples will have +1 log-coverage of
                chrX; otherwise male samples would have -1 chrX).""")
P_gainloss.add_argument('-x', '--sample-sex', '-g', '--gender',
        dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the sample's chromosomal sex as male or female.
                (Otherwise guessed from X and Y coverage).""")
P_gainloss.add_argument('-o', '--output',
        help="Output table file name.")
P_gainloss.set_defaults(func=_cmd_gainloss)


# sex/gender ------------------------------------------------------------------

def _cmd_sex(args):
    """Guess samples' sex from the relative coverage of chromosomes X and Y."""
    cnarrs = (tabio.read_cna(fname) for fname in args.filenames)
    table = do_sex(cnarrs, args.male_reference)
    core.write_dataframe(args.output, table, header=True)


@public
def do_sex(cnarrs, is_male_reference):
    """Guess samples' sex from the relative coverage of chromosomes X and Y."""
    def strsign(num):
        if num > 0:
            return "+%.3g" % num
        return "%.3g" % num

    def guess_and_format(cna):
        is_xy, stats = cna.compare_sex_chromosomes(is_male_reference)
        return (cna.meta["filename"] or cna.sample_id,
                ("Female", "Male")[is_xy],
                strsign(stats['chrx_ratio']),
                strsign(stats['chry_ratio']))

    rows = (guess_and_format(cna) for cna in cnarrs)
    columns = ["sample", "sex", "X_logratio", "Y_logratio"]
    return pd.DataFrame.from_records(rows, columns=columns)


P_sex = AP_subparsers.add_parser('sex', help=_cmd_sex.__doc__)
P_sex.add_argument('filenames', nargs='+',
        help="Copy number or copy ratio files (*.cnn, *.cnr).")
P_sex.add_argument('-y', '--male-reference', action='store_true',
        help="""Assume inputs are already normalized to a male reference
                (i.e. female samples will have +1 log-coverage of chrX;
                otherwise male samples would have -1 chrX).""")
P_sex.add_argument('-o', '--output',
        help="Output table file name.")
P_sex.set_defaults(func=_cmd_sex)

# Shims
AP_subparsers._name_parser_map['gender'] = P_sex
do_gender = public(do_sex)


# metrics ---------------------------------------------------------------------

do_metrics = public(metrics.do_metrics)


def _cmd_metrics(args):
    """Compute coverage deviations and other metrics for self-evaluation."""
    if (len(args.cnarrays) > 1 and
        args.segments and len(args.segments) > 1 and
        len(args.cnarrays) != len(args.segments)):
        raise ValueError("Number of coverage/segment filenames given must be "
                         "equal, if more than 1 segment file is given.")

    cnarrs = map(tabio.read_cna, args.cnarrays)
    if args.segments:
        args.segments = map(tabio.read_cna, args.segments)
    table = metrics.do_metrics(cnarrs, args.segments, args.drop_low_coverage)
    core.write_dataframe(args.output, table)


P_metrics = AP_subparsers.add_parser('metrics', help=_cmd_metrics.__doc__)
P_metrics.add_argument('cnarrays', nargs='+',
        help="""One or more bin-level coverage data files (*.cnn, *.cnr).""")
P_metrics.add_argument('-s', '--segments', nargs='+',
        help="""One or more segmentation data files (*.cns, output of the
                'segment' command).  If more than one file is given, the number
                must match the coverage data files, in which case the input
                files will be paired together in the given order. Otherwise, the
                same segments will be used for all coverage files.""")
P_metrics.add_argument("--drop-low-coverage", action='store_true',
        help="""Drop very-low-coverage bins before calculations to reduce
                negative "fat tail" of bin log2 values in poor-quality
                tumor samples.""")
P_metrics.add_argument('-o', '--output',
        help="Output table file name.")
P_metrics.set_defaults(func=_cmd_metrics)


# segmetrics ------------------------------------------------------------------

def _cmd_segmetrics(args):
    """Compute segment-level metrics from bin-level log2 ratios."""
    if not 0.0 < args.alpha <= 1.0:
        raise RuntimeError("alpha must be between 0 and 1.")

    from scipy.stats import sem
    # silence sem's "Degrees of freedom <= 0 for slice"; NaN is OK
    import warnings
    warnings.simplefilter("ignore", RuntimeWarning)

    stats = {
        'mean': np.mean,
        'median': np.median,
        'mode': descriptives.modal_location,
        'stdev': np.std,
        'mad':  descriptives.median_absolute_deviation,
        'mse':  descriptives.mean_squared_error,
        'iqr':  descriptives.interquartile_range,
        'bivar': descriptives.biweight_midvariance,
        'sem': sem,
        'ci': lambda x: metrics.confidence_interval_bootstrap(x, args.alpha,
                                                              args.bootstrap),
        'pi': lambda x: metrics.prediction_interval(x, args.alpha),
    }
    if not any(getattr(args, name) for name in stats):
        logging.info("No stats specified")
        return

    # Calculate all metrics
    cnarr = tabio.read_cna(args.cnarray)
    if args.drop_low_coverage:
        cnarr = cnarr.drop_low_coverage()
    segarr = tabio.read_cna(args.segments)
    segments, segbins = zip(*cnarr.by_ranges(segarr))
    # Measures of location
    for statname in ("mean", "median", "mode"):
        if getattr(args, statname):
            func = stats[statname]
            segarr[statname] = np.asfarray([func(sb.log2) for sb in segbins])
    # Measures of spread
    deviations = [sb.log2 - seg.log2 for seg, sb in zip(segments, segbins)]
    for statname in ("stdev", "sem", "mad", "mse", "iqr", "bivar"):
        if getattr(args, statname):
            func = stats[statname]
            segarr[statname] = np.asfarray([func(d) for d in deviations])

    # Interval calculations
    if args.ci:
        segarr["ci_lo"], segarr["ci_hi"] = _segmetric_interval(segarr, cnarr,
                                                               stats['ci'])
    if args.pi:
        segarr["pi_lo"], segarr["pi_hi"] = _segmetric_interval(segarr, cnarr,
                                                               stats['pi'])

    tabio.write(segarr, args.output or segarr.sample_id + ".segmetrics.cns")


def _segmetric_interval(segarr, cnarr, func):
    """Compute a stat that yields intervals (low & high values)"""
    out_vals_lo =  np.repeat(np.nan, len(segarr))
    out_vals_hi = np.repeat(np.nan, len(segarr))
    for i, (_segment, bins) in enumerate(cnarr.by_ranges(segarr)):
        if len(bins):
            out_vals_lo[i], out_vals_hi[i] = func(bins)
    return out_vals_lo, out_vals_hi


P_segmetrics = AP_subparsers.add_parser('segmetrics', help=_cmd_segmetrics.__doc__)
P_segmetrics.add_argument('cnarray',
        help="""Bin-level copy ratio data file (*.cnn, *.cnr).""")
P_segmetrics.add_argument('-s', '--segments', required=True,
        help="Segmentation data file (*.cns, output of the 'segment' command).")
P_segmetrics.add_argument("--drop-low-coverage", action='store_true',
        help="""Drop very-low-coverage bins before calculations to avoid
                negative bias in poor-quality tumor samples.""")
P_segmetrics.add_argument('-o', '--output',
        help="Output table file name.")

P_segmetrics_stats = P_segmetrics.add_argument_group(
    "Statistics available")
P_segmetrics_stats.add_argument('--mean', action='store_true',
        help="Mean log2 value (unweighted).")
P_segmetrics_stats.add_argument('--median', action='store_true',
        help="Median.")
P_segmetrics_stats.add_argument('--mode', action='store_true',
        help="Mode (i.e. peak density of log2 values).")
P_segmetrics_stats.add_argument('--stdev', action='store_true',
        help="Standard deviation.")
P_segmetrics_stats.add_argument('--sem', action='store_true',
        help="Standard error of the mean.")
P_segmetrics_stats.add_argument('--mad', action='store_true',
        help="Median absolute deviation (standardized).")
P_segmetrics_stats.add_argument('--mse', action='store_true',
        help="Mean squared error.")
P_segmetrics_stats.add_argument('--iqr', action='store_true',
        help="Inter-quartile range.")
P_segmetrics_stats.add_argument('--bivar', action='store_true',
        help="Tukey's biweight midvariance.")
P_segmetrics_stats.add_argument('--ci', action='store_true',
        help="Confidence interval (by bootstrap).")
P_segmetrics_stats.add_argument('--pi', action='store_true',
        help="Prediction interval.")
P_segmetrics_stats.add_argument('-a', '--alpha', type=float, default=.05,
        help="""Level to estimate confidence and prediction intervals;
                use with --ci and --pi. [Default: %(default)s]""")
P_segmetrics_stats.add_argument('-b', '--bootstrap', type=int, default=100,
        help="""Number of bootstrap iterations to estimate confidence interval;
                use with --ci. [Default: %(default)d]""")
P_segmetrics.set_defaults(func=_cmd_segmetrics)


# _____________________________________________________________________________
# Other I/O and compatibility

# import-picard ---------------------------------------------------------------

def _cmd_import_picard(args):
    """Convert Picard CalculateHsMetrics tabular output to CNVkit .cnn files.

    The input file is generated by the PER_TARGET_COVERAGE option in the
    CalculateHsMetrics script in Picard tools.
    """
    for fname in importers.find_picard_files(args.targets):
        cnarr = tabio.read(fname, "picardhs")
        outfname = os.path.basename(fname)[:-4] + '.cnn'
        if args.output_dir:
            if not os.path.isdir(args.output_dir):
                os.mkdir(args.output_dir)
                logging.info("Created directory %s", args.output_dir)
            outfname = os.path.join(args.output_dir, outfname)
        tabio.write(cnarr, outfname)


P_import_picard = AP_subparsers.add_parser('import-picard',
        help=_cmd_import_picard.__doc__)
P_import_picard.add_argument('targets', nargs='*', default=['.'],
        help="""Sample coverage .csv files (target and antitarget), or the
                directory that contains them.""")
P_import_picard.add_argument('-d', '--output-dir', default='.',
        help="Output directory name.")
P_import_picard.set_defaults(func=_cmd_import_picard)


# import-seg ------------------------------------------------------------------

def _cmd_import_seg(args):
    """Convert a SEG file to CNVkit .cns files."""
    if args.chromosomes:
        if args.chromosomes == 'human':
            chrom_names = {'23': 'X', '24': 'Y', '25': 'M'}
        else:
            chrom_names = dict(kv.split(':')
                               for kv in args.chromosomes.split(','))
    else:
        chrom_names = args.chromosomes

    for sid, segtable in tabio.seg.parse_seg(args.segfile, chrom_names, args.prefix,
                                             args.from_log10):
        segarr = _CNA(segtable, {"sample_id": sid})
        tabio.write(segarr, os.path.join(args.output_dir, sid + '.cns'))


P_import_seg = AP_subparsers.add_parser('import-seg',
        help=_cmd_import_seg.__doc__)
P_import_seg.add_argument('segfile',
        help="""Input file in SEG format. May contain multiple samples.""")
P_import_seg.add_argument('-c', '--chromosomes',
        help="""Mapping of chromosome indexes to names. Syntax:
                "from1:to1,from2:to2". Or use "human" for the preset:
                "23:X,24:Y,25:M".""")
P_import_seg.add_argument('-p', '--prefix',
        help="""Prefix to add to chromosome names (e.g 'chr' to rename '8' in
                the SEG file to 'chr8' in the output).""")
P_import_seg.add_argument('--from-log10', action='store_true',
        help="Convert base-10 logarithm values in the input to base-2 logs.")
P_import_seg.add_argument('-d', '--output-dir', default='.',
        help="Output directory name.")
P_import_seg.set_defaults(func=_cmd_import_seg)


# import-theta ---------------------------------------------------------------

do_import_theta = public(importers.do_import_theta)


def _cmd_import_theta(args):
    """Convert THetA output to a BED-like, CNVkit-like tabular format.

    Equivalently, use the THetA results file to convert CNVkit .cns segments to
    integer copy number calls.
    """
    tumor_segs = tabio.read_cna(args.tumor_cns)
    for i, new_cns in enumerate(do_import_theta(tumor_segs, args.theta_results,
                                                args.ploidy)):
        tabio.write(new_cns,
                    os.path.join(args.output_dir,
                                 "%s-%d.cns" % (tumor_segs.sample_id, i + 1)))


P_import_theta = AP_subparsers.add_parser('import-theta',
        help=_cmd_import_theta.__doc__)
P_import_theta.add_argument("tumor_cns")
P_import_theta.add_argument("theta_results")
P_import_theta.add_argument("--ploidy", type=int, default=2,
        help="Ploidy of normal cells. [Default: %(default)d]")
P_import_theta.add_argument('-d', '--output-dir', default='.',
        help="Output directory name.")
P_import_theta.set_defaults(func=_cmd_import_theta)


# export ----------------------------------------------------------------------

P_export = AP_subparsers.add_parser('export',
        help="""Convert CNVkit output files to another format.""")
P_export_subparsers = P_export.add_subparsers(
        help="Export formats (use with -h for more info).")


# BED special case: multiple samples's segments, like SEG
def _cmd_export_bed(args):
    """Convert segments to BED format.

    Input is a segmentation file (.cns) where, preferably, log2 ratios have
    already been adjusted to integer absolute values using the 'call' command.
    """
    bed_tables = []
    for segfname in args.segments:
        segments = tabio.read_cna(segfname)
        # ENH: args.sample_sex as a comma-separated list
        is_sample_female = verify_sample_sex(segments, args.sample_sex,
                                             args.male_reference)
        tbl = export.export_bed(segments, args.ploidy,
                                args.male_reference, is_sample_female,
                                args.sample_id or segments.sample_id,
                                args.show)
        bed_tables.append(tbl)
    table = pd.concat(bed_tables)
    core.write_dataframe(args.output, table, header=False)

P_export_bed = P_export_subparsers.add_parser('bed',
        help=_cmd_export_bed.__doc__)
P_export_bed.add_argument('segments', nargs='+',
        help="""Segmented copy ratio data files (*.cns), the output of the
                'segment' or 'call' sub-commands.""")
P_export_bed.add_argument("-i", "--sample-id", metavar="LABEL",
        help="""Identifier to write in the 4th column of the BED file.
                [Default: use the sample ID, taken from the file name]""")
P_export_bed.add_argument("--ploidy", type=int, default=2,
        help="Ploidy of the sample cells. [Default: %(default)d]")
P_export_bed.add_argument('-x', '--sample-sex', '-g', '--gender',
        dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the sample's chromosomal sex as male or female.
                (Otherwise guessed from X and Y coverage).""")
P_export_bed.add_argument("--show",
        choices=('ploidy', 'variant', 'all'), default="ploidy",
        help="""Which segmented regions to show:
                'all' = all segment regions;
                'variant' = CNA regions with non-neutral copy number;
                'ploidy' = CNA regions with non-default ploidy.
                [Default: %(default)s]""")
P_export_bed.add_argument("-y", "--male-reference", action="store_true",
        help="""Was a male reference used?  If so, expect half ploidy on
                chrX and chrY; otherwise, only chrY has half ploidy.  In CNVkit,
                if a male reference was used, the "neutral" copy number (ploidy)
                of chrX is 1; chrY is haploid for either reference sex.""")
P_export_bed.add_argument('-o', '--output', help="Output file name.")
P_export_bed.set_defaults(func=_cmd_export_bed)


# SEG special case: segment coords don't match across samples
def _cmd_export_seg(args):
    """Convert segments to SEG format.

    Compatible with IGV and GenePattern.
    """
    table = export.export_seg(args.filenames)
    core.write_dataframe(args.output, table)

P_export_seg = P_export_subparsers.add_parser('seg',
        help=_cmd_export_seg.__doc__)
P_export_seg.add_argument('filenames', nargs='+',
        help="""Segmented copy ratio data file(s) (*.cns), the output of the
                'segment' sub-command.""")
P_export_seg.add_argument('-o', '--output', help="Output file name.")
P_export_seg.set_defaults(func=_cmd_export_seg)


# VCF special case: only 1 sample, for now
def _cmd_export_vcf(args):
    """Convert segments to VCF format.

    Input is a segmentation file (.cns) where, preferably, log2 ratios have
    already been adjusted to integer absolute values using the 'call' command.
    """
    segments = tabio.read_cna(args.segments)
    is_sample_female = verify_sample_sex(segments, args.sample_sex,
                                         args.male_reference)
    header, body = export.export_vcf(segments, args.ploidy, args.male_reference,
                                     is_sample_female, args.sample_id)
    core.write_text(args.output, header, body)

P_export_vcf = P_export_subparsers.add_parser('vcf',
        help=_cmd_export_vcf.__doc__)
P_export_vcf.add_argument('segments', #nargs='1',
        help="""Segmented copy ratio data file (*.cns), the output of the
                'segment' or 'call' sub-commands.""")
P_export_vcf.add_argument("-i", "--sample-id", metavar="LABEL",
        help="""Sample name to write in the genotype field of the output VCF file.
                [Default: use the sample ID, taken from the file name]""")
P_export_vcf.add_argument("--ploidy", type=int, default=2,
        help="Ploidy of the sample cells. [Default: %(default)d]")
P_export_vcf.add_argument('-x', '--sample-sex', '-g', '--gender',
        dest='sample_sex',
        choices=('m', 'y', 'male', 'Male', 'f', 'x', 'female', 'Female'),
        help="""Specify the sample's chromosomal sex as male or female.
                (Otherwise guessed from X and Y coverage).""")
P_export_vcf.add_argument("-y", "--male-reference", action="store_true",
        help="""Was a male reference used?  If so, expect half ploidy on
                chrX and chrY; otherwise, only chrY has half ploidy.  In CNVkit,
                if a male reference was used, the "neutral" copy number (ploidy)
                of chrX is 1; chrY is haploid for either reference sex.""")
P_export_vcf.add_argument('-o', '--output', help="Output file name.")
P_export_vcf.set_defaults(func=_cmd_export_vcf)


# THetA special case: takes tumor .cns and normal .cnr or reference.cnn
def _cmd_export_theta(args):
    """Convert segments to THetA2 input file format (*.input)."""
    tumor_cn = tabio.read_cna(args.tumor_segment)
    normal_cn = (tabio.read_cna(args.reference)
                 if args.reference else None)
    table = export.export_theta(tumor_cn, normal_cn)
    if not args.output:
        args.output = tumor_cn.sample_id + ".interval_count"
    table.to_csv(args.output, sep='\t', index=False)
    logging.info("Wrote %s", args.output)
    if args.vcf:
        variants = load_het_snps(args.vcf,
                                 args.sample_id,  # or tumor_cn.sample_id,
                                 args.normal_id, args.min_variant_depth,
                                 args.zygosity_freq)
        if not len(variants):
            raise ValueError("VCF contains no usable SNV records")
        try:
            tumor_snps, normal_snps = export.export_theta_snps(variants)
        except ValueError:
            raise ValueError("VCF does not contain any tumor/normal paired "
                             "samples")
        for title, table in [("tumor", tumor_snps), ("normal", normal_snps)]:
            out_fname = "{}.{}.snp_formatted.txt".format(tumor_cn.sample_id, title)
            table.to_csv(out_fname, sep='\t', index=False)
            logging.info("Wrote %s", out_fname)

P_export_theta = P_export_subparsers.add_parser('theta',
        help=_cmd_export_theta.__doc__)
P_export_theta.add_argument('tumor_segment',
        help="""Tumor-sample segmentation file from CNVkit (.cns).""")
P_export_theta.add_argument('-r', '--reference',
        help="""Reference copy number profile (.cnn), or normal-sample bin-level
                log2 copy ratios (.cnr). Use if the tumor_segment input file
                does not contain a "weight" column.""")
P_export_theta.add_argument('-o', '--output', help="Output file name.")

P_extheta_vcf = P_export_theta.add_argument_group(
    "To also output tables of SNP b-allele frequencies for THetA2")
P_extheta_vcf.add_argument('-v', '--vcf',
        help="""VCF file containing SNVs observed in both the tumor and normal
                samples. Tumor sample ID should match the `tumor_segment`
                filename or be specified with -i/--sample-id.""")
P_extheta_vcf.add_argument('-i', '--sample-id',
        help="""Specify the name of the tumor sample in the VCF (given with
                -v/--vcf). [Default: taken the tumor_segment file name]""")
P_extheta_vcf.add_argument('-n', '--normal-id',
        help="Corresponding normal sample ID in the input VCF.")
P_extheta_vcf.add_argument('-m', '--min-variant-depth', type=int, default=20,
        help="""Minimum read depth for a SNP in the VCF to be counted.
                [Default: %(default)s]""")
P_extheta_vcf.add_argument('-z', '--zygosity-freq',
        metavar='ALT_FREQ', nargs='?', type=float, const=0.25,
        help="""Ignore VCF's genotypes (GT field) and instead infer zygosity
                from allele frequencies.  [Default if used without a number:
                %(const)s]""")

P_export_theta.set_defaults(func=_cmd_export_theta)


# Nexus "basic" special case: can only represent 1 sample
def _cmd_export_nb(args):
    """Convert bin-level log2 ratios to Nexus Copy Number "basic" format."""
    cnarr = tabio.read_cna(args.filename)
    table = export.export_nexus_basic(cnarr)
    core.write_dataframe(args.output, table)

P_export_nb = P_export_subparsers.add_parser('nexus-basic',
        help=_cmd_export_nb.__doc__)
P_export_nb.add_argument('filename',
        help="""Log2 copy ratio data file (*.cnr), the output of the 'fix'
                sub-command.""")
P_export_nb.add_argument('-o', '--output', help="Output file name.")
P_export_nb.set_defaults(func=_cmd_export_nb)


# Nexus "Custom-OGT" special case: can only represent 1 sample
def _cmd_export_nbo(args):
    """Convert log2 ratios and b-allele freqs to Nexus "Custom-OGT" format."""
    cnarr = tabio.read_cna(args.filename)
    varr = load_het_snps(args.vcf, args.sample_id, args.normal_id,
                         args.min_variant_depth, args.zygosity_freq)
    table = export.export_nexus_ogt(cnarr, varr, args.min_weight)
    core.write_dataframe(args.output, table)

P_export_nbo = P_export_subparsers.add_parser('nexus-ogt',
        help=_cmd_export_nbo.__doc__)
P_export_nbo.add_argument('filename',
        help="""Log2 copy ratio data file (*.cnr), the output of the 'fix'
                sub-command.""")
P_export_nbo.add_argument('vcf',
        help="""VCF of SNVs for the same sample, to calculate b-allele
                frequencies.""")
P_export_nbo.add_argument('-i', '--sample-id',
        help="""Specify the name of the sample in the VCF to use to extract
                b-allele frequencies.""")
P_export_nbo.add_argument('-n', '--normal-id',
        help="Corresponding normal sample ID in the input VCF.")
P_export_nbo.add_argument('-m', '--min-variant-depth', type=int, default=20,
        help="""Minimum read depth for a SNV to be included in the b-allele
                frequency calculation. [Default: %(default)s]""")
P_export_nbo.add_argument('-z', '--zygosity-freq',
        metavar='ALT_FREQ', nargs='?', type=float, const=0.25,
        help="""Ignore VCF's genotypes (GT field) and instead infer zygosity
                from allele frequencies.  [Default if used without a number:
                %(const)s]""")
P_export_nbo.add_argument('-w', '--min-weight', type=float, default=0.0,
        help="""Minimum weight (between 0 and 1) for a bin to be included in
                the output. [Default: %(default)s]""")
P_export_nbo.add_argument('-o', '--output', help="Output file name.")
P_export_nbo.set_defaults(func=_cmd_export_nbo)


# All else: export any number of .cnr or .cns files

for fmt_key, fmt_descr in (
    ('cdt', "Convert log2 ratios to CDT format. Compatible with Java TreeView."),
    ('jtv', "Convert log2 ratios to Java TreeView's native format."),
    # Not implemented yet:
    # 'gct' (GenePattern).
):
    def _cmd_export_simple(args):
        sample_ids = list(map(core.fbase, args.filenames))
        table = export.merge_samples(args.filenames)
        formatter = export.EXPORT_FORMATS[fmt_key]
        outheader, outrows = formatter(sample_ids, table)
        core.write_tsv(args.output, outrows, colnames=outheader)

    P_export_simple = P_export_subparsers.add_parser(fmt_key, help=fmt_descr)
    P_export_simple.add_argument('filenames', nargs='+',
            help="""Log2 copy ratio data file(s) (*.cnr), the output of the
                    'fix' sub-command.""")
    P_export_simple.add_argument('-o', '--output', help="Output file name.")
    P_export_simple.set_defaults(func=_cmd_export_simple)


# version ---------------------------------------------------------------------

def print_version(_args):
    """Display this program's version."""
    print(__version__)


P_version = AP_subparsers.add_parser('version', help=print_version.__doc__)
P_version.set_defaults(func=print_version)


# _____________________________________________________________________________
# Functions reused within this module

def load_het_snps(vcf_fname, sample_id, normal_id, min_variant_depth,
                  zygosity_freq, tumor_boost=False):
    if vcf_fname is None:
        return None
    varr = tabio.read(vcf_fname, 'vcf',
                      sample_id=sample_id,
                      normal_id=normal_id,
                      min_depth=min_variant_depth,
                      skip_somatic=True)
    if (zygosity_freq is None and 'n_zygosity' in varr and
        not varr['n_zygosity'].any()):
        # Mutect2 sets all normal genotypes to 0/0 -- work around it
        logging.warn("VCF normal sample's genotypes are all 0/0 or missing; "
                     "inferring genotypes from allele frequency instead")
        zygosity_freq = 0.25
    if zygosity_freq is not None:
        varr = varr.zygosity_from_freq(zygosity_freq, 1 - zygosity_freq)
    if 'n_zygosity' in varr:
        # Infer & drop (more) somatic loci based on genotype
        somatic_idx = (varr['zygosity'] != 0.0) & (varr['n_zygosity'] == 0.0)
        if somatic_idx.any() and not somatic_idx.all():
            logging.info("Skipping %d additional somatic record based on "
                         "T/N genotypes", somatic_idx.sum())
        varr = varr[~somatic_idx]
    orig_len = len(varr)
    varr = varr.heterozygous()
    logging.info("Kept %d heterozygous of %d VCF records",
                 len(varr), orig_len)
    # TODO use/explore tumor_boost option
    if tumor_boost:
        varr['alt_freq'] = varr.tumor_boost()
    return varr


def verify_sample_sex(cnarr, sex_arg, is_male_reference):
    is_sample_female = cnarr.guess_xx(is_male_reference, verbose=False)
    if sex_arg:
        is_sample_female_given = (sex_arg.lower() not in ['y', 'm', 'male'])
        if is_sample_female != is_sample_female_given:
            logging.warn("Sample sex specified as %s "
                         "but chromosomal X/Y ploidy looks like %s",
                         "female" if is_sample_female_given else "male",
                         "female" if is_sample_female else "male")
            is_sample_female = is_sample_female_given
    logging.info("Treating sample %s as %s",
                 cnarr.sample_id or '',
                 "female" if is_sample_female else "male")
    return is_sample_female


# _____________________________________________________________________________
# Shim for command-line execution

def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)
