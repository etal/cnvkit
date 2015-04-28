"""Command-line interface and corresponding API for CNVkit."""
# NB: argparse CLI definitions and API functions are interwoven:
#   "_cmd_*" handles I/O and arguments processing for the command
#   "do_*" runs the command's functionality as an API
from __future__ import absolute_import, division
import argparse
import collections
import math
import multiprocessing
import os
import sys

import numpy
import matplotlib as mpl
mpl.use("Agg", force=True)
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
# pyplot.ioff()


from Bio._py3k import map, range, zip
iteritems = (dict.iteritems if sys.version_info[0] < 3 else dict.items)

from . import (core, target, antitarget, coverage, fix, metrics, reference,
               reports, export, importers, segmentation,
               params, ngfrills, plots)
from .ngfrills import echo
from .cnarray import CopyNumArray as CNA


AP = argparse.ArgumentParser(description=__doc__,
        epilog="Contact Eric Talevich <eric.talevich@ucsf.edu> for help.")
AP_subparsers = AP.add_subparsers(
        help="Sub-commands (use with -h for more info)")


# _____________________________________________________________________________
# Core pipeline

class SerialPool(object):
    """Mimic the multiprocessing.Pool interface, but run in serial."""

    def __init__(self):
        pass

    def apply_async(self, func, args):
        """Just call the function."""
        func(*args)

    # No-ops to mimic Pool
    def close(self): pass
    def join(self): pass


def pick_pool(nprocs):
    if nprocs == 1:
        return SerialPool()
    return multiprocessing.Pool(nprocs)


# batch -----------------------------------------------------------------------

def _cmd_batch(args):
    """Run the complete CNVkit pipeline on one or more BAM files."""
    # Mutual exclusion of -r/-n
    if not ((args.normal is not None) ^ bool(args.reference)):
        raise ValueError("One of the arguments -n/--normal or -r/--reference "
                         "is required.")

    if args.processes < 1:
        args.processes = multiprocessing.cpu_count()

    if not args.reference:
        # Need target BED file to build a new reference
        if not args.targets:
            raise ValueError("Argument -t/--target is required.")
        # Build a copy number reference; update (anti)targets upon request
        args.reference, args.targets, args.antitargets = batch_make_reference(
            args.normal, args.targets, args.antitargets, args.male_reference,
            args.fasta, args.annotate, args.short_names, args.split,
            args.target_avg_size, args.access, args.antitarget_avg_size,
            args.antitarget_min_size, args.output_reference, args.output_dir,
            args.processes, args.count_reads)
    elif args.targets is None and args.antitargets is None:
        # Extract (anti)target BEDs from the given, existing CN reference
        ref_arr = CNA.read(args.reference)
        target_coords, antitarget_coords = reference.reference2regions(ref_arr)
        ref_pfx = os.path.join(args.output_dir, core.fbase(args.reference))
        args.targets = ref_pfx + '.target-tmp.bed'
        args.antitargets = ref_pfx + '.antitarget-tmp.bed'
        core.write_tsv(args.targets, target_coords)
        core.write_tsv(args.antitargets, antitarget_coords)

    if args.bam_files:
        echo("Running", len(args.bam_files), "samples in",
             (("%d processes" % args.processes)
              if args.processes > 1 else "serial"))

        pool = pick_pool(args.processes)
        for bam in args.bam_files:
            pool.apply_async(batch_run_sample,
                             (bam, args.targets, args.antitargets, args.reference,
                              args.output_dir, args.male_reference, args.scatter,
                              args.diagram, args.rlibpath, args.count_reads))
        pool.close()
        pool.join()


def batch_make_reference(normal_bams, target_bed, antitarget_bed, male_reference,
                         fasta, annotate, short_names, split, target_avg_size,
                         access, antitarget_avg_size, antitarget_min_size,
                         output_reference, output_dir, processes, by_count):
    """Build the CN reference from normal samples, targets and antitargets."""
    # To make temporary filenames for processed targets or antitargets
    tgt_name_base, tgt_name_ext = os.path.splitext(os.path.basename(target_bed))
    if output_dir:
        tgt_name_base = os.path.join(output_dir, tgt_name_base)

    if annotate or short_names or split:
        # Pre-process baits/targets
        new_target_fname = tgt_name_base + '.target' + tgt_name_ext
        do_targets(target_bed, new_target_fname,
                    annotate, short_names, split,
                    **({'avg_size': target_avg_size}
                        if split and target_avg_size
                        else {}))
        target_bed = new_target_fname

    if not antitarget_bed:
        # Build antitarget BED from the given targets
        anti_kwargs = {}
        if access:
            anti_kwargs['access_bed'] = access
        if antitarget_avg_size:
            anti_kwargs['avg_bin_size'] = antitarget_avg_size
        if antitarget_min_size:
            anti_kwargs['min_bin_size'] = antitarget_min_size
        anti_rows = do_antitarget(target_bed, **anti_kwargs)
        # Devise a temporary antitarget filename
        antitarget_bed = tgt_name_base + '.antitarget' + tgt_name_ext
        with ngfrills.safe_write(antitarget_bed, False) as anti_file:
            i = 0
            for i, row in enumerate(anti_rows):
                anti_file.write("\t".join(row) + '\n')
            echo("Wrote", antitarget_bed,
                    "with", i + 1, "background intervals")

    if len(normal_bams) == 0:
        echo("Building a flat reference...")
        ref_arr = do_reference_flat(target_bed, antitarget_bed, fasta,
                                    male_reference)
    else:
        echo("Building a copy number reference from normal samples...")
        target_fnames = []
        antitarget_fnames = []
        # Run coverage on all normals
        pool = pick_pool(processes)
        for nbam in normal_bams:
            sample_id = core.fbase(nbam)
            sample_pfx = os.path.join(output_dir, sample_id)
            tgt_fname = sample_pfx + '.targetcoverage.cnn'
            pool.apply_async(batch_write_coverage,
                             (target_bed, nbam, tgt_fname, by_count))
            target_fnames.append(tgt_fname)
            anti_fname = sample_pfx + '.antitargetcoverage.cnn'
            pool.apply_async(batch_write_coverage,
                             (antitarget_bed, nbam, anti_fname, by_count))
            antitarget_fnames.append(anti_fname)
        pool.close()
        pool.join()
        # Build reference from *.cnn
        ref_arr = do_reference(target_fnames, antitarget_fnames, fasta,
                               male_reference)
    if not output_reference:
        output_reference = os.path.join(output_dir, "reference.cnn")
    ngfrills.ensure_path(output_reference)
    ref_arr.write(output_reference)

    return output_reference, target_bed, antitarget_bed


def batch_write_coverage(bed_fname, bam_fname, out_fname, by_count):
    """Run coverage on one sample, write to file."""
    cnarr = do_coverage(bed_fname, bam_fname, by_count)
    cnarr.write(out_fname)


def batch_run_sample(bam_fname, target_bed, antitarget_bed, ref_fname,
                     output_dir, male_reference=False, scatter=False,
                     diagram=False, rlibpath=None, by_count=False):
    """Run the pipeline on one BAM file."""
    # ENH - return probes, segments (cnarr, segarr)
    echo("Running the CNVkit pipeline on", bam_fname, "...")
    sample_id = core.fbase(bam_fname)
    sample_pfx = os.path.join(output_dir, sample_id)

    raw_tgt = do_coverage(target_bed, bam_fname, by_count)
    raw_tgt.write(sample_pfx + '.targetcoverage.cnn')

    raw_anti = do_coverage(antitarget_bed, bam_fname, by_count)
    raw_anti.write(sample_pfx + '.antitargetcoverage.cnn')

    cnarr = do_fix(raw_tgt, raw_anti, CNA.read(ref_fname))
    cnarr.write(sample_pfx + '.cnr')

    echo("Segmenting", sample_pfx + '.cnr ...')
    segments = segmentation.do_segmentation(sample_pfx + '.cnr', False, 'cbs',
                                            rlibpath)
    segments.write(sample_pfx + '.cns')

    if scatter:
        do_scatter(cnarr, segments)
        pyplot.savefig(sample_pfx + '-scatter.pdf', format='pdf',
                       bbox_inches="tight")
        echo("Wrote", sample_pfx + '-scatter.pdf')

    if diagram:
        from cnvlib import diagram
        outfname = sample_pfx + '-diagram.pdf'
        diagram.create_diagram(cnarr, segments, 0.6, outfname, male_reference)
        echo("Wrote", outfname)


P_batch = AP_subparsers.add_parser('batch', help=_cmd_batch.__doc__)
P_batch.add_argument('bam_files', nargs='*',
        help="Mapped sequence reads (.bam)")
P_batch.add_argument('-y', '--male-reference', action='store_true',
        help="""Use or assume a male reference (i.e. female samples will have +1
                log-CNR of chrX; otherwise male samples would have -1 chrX).""")
P_batch.add_argument('-c', '--count-reads', action='store_true',
        help="""Get read depths by counting read midpoints within each bin.
                (An alternative algorithm).""")
P_batch.add_argument('-p', '--processes', type=int, default=1,
        help="""Number of subprocesses used to running each of the BAM files in
                parallel. Give 0 or a negative value to use the maximum number
                of available CPUs. [Default: process each BAM in serial]""")
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
        help="""UCSC refFlat.txt or ensFlat.txt file for the reference genome.
                Pull gene names from this file and assign them to the target
                regions.""")
P_batch_newref.add_argument('--short-names', action='store_true',
        help="Reduce multi-accession bait labels to be short and consistent.")
P_batch_newref.add_argument('--split', action='store_true',
        help="Split large tiled intervals into smaller, consecutive targets.")
P_batch_newref.add_argument('--target-avg-size', type=int,
        help="Average size of split target bins (results are approximate).")
# For antitargets:
P_batch_newref.add_argument('-g', '--access',
        help="""Regions of accessible sequence on chromosomes (.bed), as
                output by genome2access.py.""")
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
        help="Create a diagram of copy ratios on chromosomes as a PDF.")

P_batch.set_defaults(func=_cmd_batch)


# target ----------------------------------------------------------------------

def _cmd_target(args):
    """Transform bait intervals into targets more suitable for CNVkit."""
    do_targets(args.interval, args.output,
               args.annotate, args.short_names, args.split, args.avg_size)


def do_targets(bed_fname, out_fname, annotate=None, do_short_names=False,
               do_split=False, avg_size=200/.75):
    """Transform bait intervals into targets more suitable for CNVkit."""
    bed_rows = ngfrills.parse_regions(bed_fname, False,
                                      keep_strand=bool(annotate))
    if annotate:
        ngfrills.echo("Applying annotations as target names")
        bed_rows = target.add_refflat_names(bed_rows, annotate)
    if do_short_names:
        ngfrills.echo("Shortening interval labels")
        bed_rows = target.shorten_labels(bed_rows)
    if do_split:
        ngfrills.echo("Splitting large targets")
        bed_rows = target.split_targets(bed_rows, avg_size)
    # Output with logging
    with ngfrills.safe_write(out_fname or sys.stdout, False) as outfile:
        i = 0
        for i, row in enumerate(sorted(bed_rows,
                                       key=lambda r: (core.sorter_chrom(r[0]),
                                                      r[1]))):
            outfile.write("\t".join(map(str, row)) + '\n')
    if out_fname:
        ngfrills.echo("Wrote", out_fname, "with", i + 1, "target intervals")


P_target = AP_subparsers.add_parser('target', help=_cmd_target.__doc__)
P_target.add_argument('interval',
        help="""BED or interval file listing the targeted regions.""")
P_target.add_argument('--annotate',
        help="""UCSC refFlat.txt or ensFlat.txt file for the reference genome.
                Pull gene names from this file and assign them to the target
                regions.""")
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


# antitarget ------------------------------------------------------------------

def _cmd_antitarget(args):
    """Derive a background/antitarget BED file from a target BED file."""
    out_rows = do_antitarget(args.interval, args.access,
                             args.avg_size, args.min_size)
    if not args.output:
        base, ext = args.interval.rsplit('.', 1)
        args.output = base + '.antitarget.' + ext
    with ngfrills.safe_write(args.output, False) as outfile:
        i = -1
        for i, row in enumerate(out_rows):
            outfile.write("\t".join(row) + '\n')
        echo("Wrote", args.output, "with", i + 1, "background intervals")


def do_antitarget(target_bed, access_bed=None, avg_bin_size=150000,
                  min_bin_size=None):
    """Derive a background/antitarget BED file from a target BED file."""
    if not min_bin_size:
        min_bin_size = 2 * int(avg_bin_size * (2 ** params.MIN_BIN_COVERAGE))
    background_regions = antitarget.get_background(target_bed, access_bed,
                                                   avg_bin_size, min_bin_size)

    # Sniff the number of columns in the BED/interval file
    # Output will try to match this format
    ncols = ngfrills.sniff_num_columns(target_bed)
    if ncols == 3:
        # Just the coordinates
        suffix = []
    elif ncols == 5:
        # Standard interval: coordinates, strand, name
        suffix = ['+', 'Background']
    else:
        # Full or extended BED, or whatever
        suffix = ['Background']

    for row in background_regions:
        yield list(map(str, row)) + suffix


P_anti = AP_subparsers.add_parser('antitarget', help=_cmd_antitarget.__doc__)
P_anti.add_argument('interval',
        help="""BED or interval file listing the targeted regions.""")
P_anti.add_argument('-g', '--access',
        help="""Regions of accessible sequence on chromosomes (.bed), as
                output by genome2access.py.""")
P_anti.add_argument('-a', '--avg-size', type=int, default=150000,
        help="""Average size of antitarget bins (results are approximate).
                [Default: %(default)s]""")
P_anti.add_argument('-m', '--min-size', type=int,
        help="Minimum size of antitarget bins (smaller regions are dropped).")
P_anti.add_argument('-o', '--output', help="""Output file name.""")
P_anti.set_defaults(func=_cmd_antitarget)


# coverage --------------------------------------------------------------------

def _cmd_coverage(args):
    """Calculate coverage in the given regions from BAM read depths."""
    pset = do_coverage(args.interval, args.bam_file, args.count)
    if not args.output:
        # Create an informative but unique name for the coverage output file
        bambase = core.fbase(args.bam_file)
        bedbase = core.rbase(args.interval)
        tgtbase = ('antitargetcoverage'
                   if 'anti' in bedbase.lower()
                   else 'targetcoverage')
        args.output = '%s.%s.cnn' % (bambase, tgtbase)
        if os.path.exists(args.output):
            args.output = '%s.%s.cnn' % (bambase, bedbase)
    ngfrills.ensure_path(args.output)
    pset.write(args.output)


def do_coverage(bed_fname, bam_fname, by_count=False):
    """Calculate coverage in the given regions from BAM read depths."""
    if not ngfrills.ensure_bam_sorted(bam_fname):
        raise RuntimeError("BAM file %s must be sorted by coordinates"
                            % bam_fname)
    ngfrills.ensure_bam_index(bam_fname)
    # ENH: count importers.TOO_MANY_NO_COVERAGE & warn
    cnarr = coverage.interval_coverages(bed_fname, bam_fname, by_count)
    if len(cnarr):
        cnarr.center_all()
    return cnarr


P_coverage = AP_subparsers.add_parser('coverage', help=_cmd_coverage.__doc__)
P_coverage.add_argument('bam_file', help="Mapped sequence reads (.bam)")
P_coverage.add_argument('interval', help="Intervals (.bed or .list)")
P_coverage.add_argument('-c', '--count', action='store_true',
        help="""Get read depths by counting read midpoints within each bin.
                (An alternative algorithm).""")
P_coverage.add_argument('-o', '--output', help="""Output file name.""")
P_coverage.set_defaults(func=_cmd_coverage)


# reference -------------------------------------------------------------------

def _cmd_reference(args):
    """Compile a coverage reference from the given files (normal samples)."""
    usage_err_msg = ("Give .cnn samples OR targets and antitargets.")
    if args.targets and args.antitargets:
        # Flat refence
        assert not args.references, usage_err_msg
        ref_probes = do_reference_flat(args.targets, args.antitargets,
                                       args.fasta, args.male_reference)
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
        echo("Number of target and antitarget files: %d, %d"
            % (len(targets), len(antitargets)))
        ref_probes = do_reference(targets, antitargets, args.fasta,
                                  args.male_reference)
    else:
        raise ValueError(usage_err_msg)

    ref_fname = args.output or "cnv_reference.cnn"
    ngfrills.ensure_path(ref_fname)
    ref_probes.write(ref_fname)


def do_reference(target_fnames, antitarget_fnames, fa_fname=None,
                 male_reference=False):
    """Compile a coverage reference from the given files (normal samples)."""
    core.assert_equal("Unequal number of target and antitarget files given",
                      targets=len(target_fnames),
                      antitargets=len(antitarget_fnames))
    # Calculate & save probe centers
    ref_probes = reference.combine_probes(target_fnames, bool(fa_fname),
                                          male_reference)
    ref_probes.extend(reference.combine_probes(antitarget_fnames,
                                               bool(fa_fname), male_reference))
    ref_probes.sort()
    ref_probes.center_all()
    # Calculate GC and RepeatMasker content for each probe's genomic region
    if fa_fname:
        gc, rmask = reference.get_fasta_stats(ref_probes, fa_fname)
        ref_probes['gc'] = gc
        ref_probes['rmask'] = rmask
    else:
        echo("No FASTA reference genome provided; skipping GC, RM calculations")
    reference.warn_bad_probes(ref_probes)
    return ref_probes


def do_reference_flat(target_list, antitarget_list, fa_fname=None,
                      male_reference=False):
    """Compile a neutral-coverage reference from the given intervals.

    Combines the intervals, shifts chrX values if requested, and calculates GC
    and RepeatMasker content from the genome FASTA sequence.
    """
    ref_probes = reference.bed2probes(target_list)
    ref_probes.extend(reference.bed2probes(antitarget_list))
    ref_probes.sort()
    # Set sex chromosomes by "reference" gender
    chr_x = core.guess_chr_x(ref_probes)
    chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
    if male_reference:
        ref_probes['coverage'][(ref_probes.chromosome == chr_x) |
                               (ref_probes.chromosome == chr_y)] = -1.0
    else:
        ref_probes['coverage'][ref_probes.chromosome == chr_y] = -1.0
    # Calculate GC and RepeatMasker content for each probe's genomic region
    if fa_fname:
        gc, rmask = reference.get_fasta_stats(ref_probes, fa_fname)
        ref_probes['gc'] = gc
        ref_probes['rmask'] = rmask
        reference.warn_bad_probes(ref_probes)
    else:
        echo("No FASTA reference genome provided; skipping GC, RM calculations")
    return ref_probes


P_reference = AP_subparsers.add_parser('reference', help=_cmd_reference.__doc__)
P_reference.add_argument('references', nargs='*',
        help="""Normal-sample target or antitarget .cnn files, or the
                directory that contains them.""")
P_reference.add_argument('-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
P_reference.add_argument('-t', '--targets',
        help="Target intervals (.bed or .list)")
P_reference.add_argument('-a', '--antitargets',
        help="Antitarget intervals (.bed or .list)")
P_reference.add_argument('-y', '--male-reference', action='store_true',
        help="""Create a male reference: shift female samples' chrX
                log-coverage by -1, so the reference chrX average is -1.
                Otherwise, shift male samples' chrX by +1, so the reference chrX
                average is 0.""")
P_reference.add_argument('-o', '--output', help="Output file name.")
P_reference.set_defaults(func=_cmd_reference)


# fix -------------------------------------------------------------------------

def _cmd_fix(args):
    """Combine target and antitarget coverages and correct for biases.

    Adjust raw coverage data according to the given reference, correct potential
    biases and re-center.
    """
    # Verify that target and antitarget are from the same sample
    tgt_raw = CNA.read(args.target)
    anti_raw = CNA.read(args.antitarget)
    if tgt_raw.sample_id != anti_raw.sample_id:
        raise ValueError("Sample IDs do not match:"
                         "'%s' (target) vs. '%s' (antitarget)"
                         % (tgt_raw.sample_id, anti_raw.sample_id))
    target_table = do_fix(tgt_raw, anti_raw, CNA.read(args.reference),
                          args.do_gc, args.do_edge, args.do_rmask)
    target_table.write(args.output or tgt_raw.sample_id + '.cnr')


def do_fix(target_raw, antitarget_raw, reference,
           do_gc=True, do_edge=True, do_rmask=True):
    """Combine target and antitarget coverages and correct for biases."""
    # Load, recenter and GC-correct target & antitarget probes separately
    echo("Processing target:", target_raw.sample_id)
    cnarr = fix.load_adjust_coverages(target_raw, reference,
                                      do_gc, do_edge, False)
    echo("Processing antitarget:", antitarget_raw.sample_id)
    anti_cnarr = fix.load_adjust_coverages(antitarget_raw, reference,
                                           do_gc, False, do_rmask)
    # Merge target and antitarget & sort probes by chromosomal location
    cnarr.extend(anti_cnarr)
    cnarr.sort()
    cnarr.center_all()
    # Determine weights for each bin (used in segmentation)
    return fix.apply_weights(cnarr, reference)


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

def _cmd_segment(args):
    """Infer copy number segments from the given coverage table."""
    if args.dataframe:
        segments, dframe = segmentation.do_segmentation(args.filename, True,
                                                         args.method, args.rlibpath)
        with open(args.dataframe, 'w') as handle:
            handle.write(dframe)
        echo("Wrote", args.dataframe)
    else:
        segments = segmentation.do_segmentation(args.filename, False,
                                                args.method, args.rlibpath)
    segments.write(args.output or segments.sample_id + '.cns')


P_segment = AP_subparsers.add_parser('segment', help=_cmd_segment.__doc__)
P_segment.add_argument('filename',
        help="Coverage file (.cnr), as produced by 'fix'.")
P_segment.add_argument('-o', '--output',
        help="Output table file name (CNR-like table of segments, .cns).")
P_segment.add_argument('-d', '--dataframe',
        help="""File name to save the raw R dataframe emitted by CBS or
                Fused Lasso. (Useful for debugging.)""")
P_segment.add_argument('-m', '--method',
        choices=('cbs', 'haar', 'flasso'), default='cbs',
        help="""Segmentation method (CBS, HaarSeg, or Fused Lasso).
                [Default: %(default)s]""")
P_segment.add_argument("--rlibpath",
                       help="Path to an alternative site-library to use for R packages.")
P_segment.set_defaults(func=_cmd_segment)


# cbs -------------------------------------------------------------------------

def _cmd_cbs(args):
    """[DEPRECATED - use 'segment' instead]

    Run circular binary segmentation (CBS) on the given coverage table.
    """
    if args.dataframe:
        segments, seg_out = segmentation.do_segmentation(args.filename, True,
                                                         'cbs', args.rlibpath)
        with ngfrills.safe_write(args.dataframe) as handle:
            handle.write(seg_out)
    else:
        segments = segmentation.do_segmentation(args.filename, False, 'cbs')
    segments.write(args.output or segments.sample_id + '.cns')


P_cbs = AP_subparsers.add_parser('cbs', help=_cmd_cbs.__doc__)
P_cbs.add_argument('filename',
        help="Coverage file (.cnr), as produced by 'fix'.")
P_cbs.add_argument('-o', '--output',
        help="Output table file name (CNR-like table of segments, .cns).")
P_cbs.add_argument('-d', '--dataframe',
        help="Output filename for the unaltered dataframe emitted by CBS.")
P_cbs.add_argument("--rlibpath",
                   help="Path to an alternative site-library to use for R packages.")
P_cbs.set_defaults(func=_cmd_cbs)


# _____________________________________________________________________________
# Plots and graphics

# diagram ---------------------------------------------------------------------

def _cmd_diagram(args):
    """Draw copy number (log2 coverages, CBS calls) on chromosomes as a diagram.

    If both the raw probes and segments are given, show them side-by-side on
    each chromosome (segments on the left side, probes on the right side).
    """
    from cnvlib import diagram
    cnarr = CNA.read(args.filename) if args.filename else None
    segarr = CNA.read(args.segment) if args.segment else None
    outfname = diagram.create_diagram(cnarr, segarr, args.threshold,
                                      args.output, args.male_reference)
    echo("Wrote", outfname)


P_diagram = AP_subparsers.add_parser('diagram', help=_cmd_diagram.__doc__)
P_diagram.add_argument('filename', nargs='?',
        help="""Processed coverage data file (*.cnr), the output of the
                'fix' sub-command.""")
P_diagram.add_argument('-s', '--segment',
        help="Segmentation calls (.cns), the output of the 'segment' command.")
P_diagram.add_argument('-t', '--threshold', type=float, default=0.6,
        help="""Copy number change threshold to label genes.
                [Default: %(default)s]""")
P_diagram.add_argument('-y', '--male-reference', action='store_true',
        help="""Assume inputs are already corrected against a male
                reference (i.e. female samples will have +1 log-CNR of
                chrX; otherwise male samples would have -1 chrX).""")
P_diagram.add_argument('-o', '--output',
        help="Output PDF file name.")
P_diagram.set_defaults(func=_cmd_diagram)


# scatter ---------------------------------------------------------------------

def _cmd_scatter(args):
    """Plot probe log2 coverages and segmentation calls together."""
    pset_cvg = CNA.read(args.filename, args.sample_id)
    pset_seg = CNA.read(args.segment) if args.segment else None
    if args.range_list:
        with PdfPages(args.output) as pdf_out:
            for chrom, start, end in ngfrills.parse_regions(args.range_list,
                                                            True):
                region = "{}:{}-{}".format(chrom, start, end)
                do_scatter(pset_cvg, pset_seg, args.vcf, args.chromosome, args.gene,
                        region, args.background_marker, args.trend, args.width,
                        args.sample_id, args.min_variant_depth)
                pyplot.title(region)
                pdf_out.savefig()
                pyplot.close()
    else:
        do_scatter(pset_cvg, pset_seg, args.vcf,
                args.chromosome, args.gene, args.range,
                args.background_marker, args.trend, args.width,
                args.sample_id, args.min_variant_depth)
        if args.output:
            pyplot.savefig(args.output, format='pdf', bbox_inches="tight")
            echo("Wrote", args.output)
        else:
            pyplot.show()

def do_scatter(pset_cvg, pset_seg=None, vcf_fname=None,
               show_chromosome=None, show_gene=None, show_range=None,
               background_marker=None, do_trend=False, window_width=1e6,
               sample_id=None, min_variant_depth=20):
    """Plot probe log2 coverages and CBS calls together."""
    if not show_gene and not show_range and not show_chromosome:
        # Plot all chromosomes, concatenated on one plot
        PAD = 1e7
        if vcf_fname:
            # Lay out top 3/5 for the CN scatter, bottom 2/5 for LOH plot
            axgrid = pyplot.GridSpec(5, 1, hspace=.85)
            axis = pyplot.subplot(axgrid[:3])
            axis2 = pyplot.subplot(axgrid[3:], sharex=axis)
            chrom_snvs = ngfrills.load_vcf(vcf_fname, min_variant_depth,
                                           sample_id=sample_id)
            # Place chromosome labels between the CNR and LOH plots
            axis2.tick_params(labelbottom=False)
            chrom_sizes = plots.chromosome_sizes(pset_cvg)
            plots.plot_loh(axis2, chrom_snvs, chrom_sizes, do_trend, PAD)
        else:
            _fig, axis = pyplot.subplots()
        axis.set_title(pset_cvg.sample_id)
        plots.plot_genome(axis, pset_cvg, pset_seg, PAD, do_trend)
    else:
        # Plot a specified region on one chromosome
        window_coords = None
        chrom = None
        genes = []
        if show_gene:
            gene_names = show_gene.split(',')
            # Scan for probes matching the specified gene
            gene_coords = plots.gene_coords_by_name(pset_cvg, gene_names)
            if not len(gene_coords) == 1:
                raise ValueError("Genes %s are split across chromosomes %s"
                                 % (show_gene, gene_coords.keys()))
            chrom, genes = gene_coords.popitem()
            genes.sort()
            # Set the display window to the selected genes +/- a margin
            window_coords = (genes[0][0] - window_width,
                             genes[-1][1] + window_width)

        if show_range:
            if chrom:
                if not show_range.startswith(chrom):
                    raise ValueError("Selected genes are on chromosome " +
                                     chrom + " but specified range is " +
                                     show_range)
            chrom, start, end = plots.parse_range(show_range)
            if window_coords:
                if start > window_coords[0] or end < window_coords[1]:
                    raise ValueError("Selected gene range " + chrom +
                                     (":%d-%d" % window_coords) +
                                     " is outside specified range " +
                                     show_range)
            if not genes:
                genes = plots.gene_coords_by_range(pset_cvg, chrom,
                                                   start, end)[chrom]
            if not genes and window_width > (end - start) / 10.0:
                # No genes in the selected region, so highlight the region
                # itself (unless the selection is ~entire displayed window)
                genes = [(start, end, "Selection")]
            window_coords = (start - window_width, end + window_width)

        if show_chromosome:
            if chrom:
                # Confirm that the selected chromosomes match
                core.assert_equal("Chromosome also selected by region or gene "
                                  "does not match",
                                  **{"chromosome": show_chromosome,
                                     "region or gene": chrom})
            else:
                chrom = show_chromosome
            if window_coords and not (show_gene and show_range):
                # Reset window to show the whole chromosome
                # (unless range and gene were both specified)
                window_coords = None

        # Prune plotted elements to the selected region
        if window_coords:
            # Show the selected region
            sel_probes = pset_cvg.in_range(chrom, *window_coords)
            sel_seg = (pset_seg.in_range(chrom, *window_coords, trim=True)
                       if pset_seg else None)
        else:
            # Show the whole chromosome
            sel_probes = pset_cvg.in_range(chrom)
            sel_seg = (pset_seg.in_range(chrom)
                       if pset_seg else None)

        echo("Showing", len(sel_probes), "probes and", len(genes),
             "selected genes in range",
             (chrom + ":%d-%d" % window_coords if window_coords else chrom))

        # Similarly for SNV allele freqs, if given
        if vcf_fname:
            # Lay out top 3/5 for the CN scatter, bottom 2/5 for LOH plot
            axgrid = pyplot.GridSpec(5, 1, hspace=.5)
            axis = pyplot.subplot(axgrid[:3])
            axis2 = pyplot.subplot(axgrid[3:], sharex=axis)
            chrom_snvs = ngfrills.load_vcf(vcf_fname, min_variant_depth,
                                           sample_id=sample_id)
            # Plot LOH for only the selected region
            snv_x_y = [(pos * plots.MB, abs(altfreq - .5) + 0.5)
                       for pos, zyg, altfreq in chrom_snvs[chrom]]
            if window_coords:
                snv_x_y = [(x, y) for x, y in snv_x_y
                           if (window_coords[0] * plots.MB <= x <=
                               window_coords[1] * plots.MB)]
            snv_x, snv_y = (zip(*snv_x_y) if snv_x_y else ([], []))
            axis2.set_ylim(0.5, 1.0)
            axis2.set_ylabel("VAF")
            axis2.scatter(snv_x, snv_y, color='#808080', alpha=0.3)
            axis2.set_xlabel("Position (Mb)")
            axis2.get_yaxis().tick_left()
            axis2.get_xaxis().tick_top()
            axis2.tick_params(which='both', direction='out',
                              labelbottom=False, labeltop=False)
        else:
            _fig, axis = pyplot.subplots()
            axis.set_xlabel("Position (Mb)")

        # Plot CNVs
        plots.plot_chromosome(axis, sel_probes, sel_seg, chrom,
                              pset_cvg.sample_id, genes,
                              background_marker=background_marker,
                              do_trend=do_trend)

P_scatter = AP_subparsers.add_parser('scatter', help=_cmd_scatter.__doc__)
P_scatter.add_argument('filename',
        help="""Processed coverage sample data file (*.cnr), the output
                of the 'fix' sub-command.""")
P_scatter.add_argument('-s', '--segment',
        help="Segmentation calls (.cns), the output of the 'segment' command.")
P_scatter.add_argument('-c', '--chromosome',
        help="Display only the specified chromosome.")
P_scatter.add_argument('-g', '--gene',
        help="Name of gene or genes (comma-separated) to display.")
P_scatter.add_argument('-r', '--range',
        help="""Chromosomal range to display (e.g. chr1:2333000-2444000).
                All targeted genes in this range will be shown, unless
                '--gene'/'-g' is already given.""")
P_scatter.add_argument('-l', '--range-list',
        help="""File listing the chromosomal ranges to display, as BED, interval
                list or "chr:start-end" text. Creates focal plots similar to
                --range for each listed region, combined into a multi-page PDF.
                The output filename must also be specified (-o/--output).""")
P_scatter.add_argument("-i", "--sample-id",
        help="""Specify the name of the sample to show in plot title and use for
                LOH analysis.""")
P_scatter.add_argument('-b', '--background-marker', default=None,
        help="""Plot antitargets with this symbol, in zoomed/selected regions.
                [Default: same as targets]""")
P_scatter.add_argument('-t', '--trend', action='store_true',
        help="Draw a smoothed local trendline on the scatter plot.")
P_scatter.add_argument('-v', '--vcf',
        help="""VCF file name containing variants to plot for LOH.""")
P_scatter.add_argument('-m', '--min-variant-depth', type=int, default=20,
        help="""Minimum read depth for a variant to be displayed in the LOH
                plot. [Default: %(default)s]""")
P_scatter.add_argument('-w', '--width', type=float, default=1e6,
        help="""Width of margin to show around the selected gene or region
                on the chromosome (use with --gene or --region).
                [Default: %(default)d]""")
P_scatter.add_argument('-o', '--output',
        help="Output table file name.")
P_scatter.set_defaults(func=_cmd_scatter)


# loh -------------------------------------------------------------------------

def _cmd_loh(args):
    """Plot allelic frequencies at each variant position in a VCF file.

    Divergence from 0.5 indicates loss of heterozygosity in a tumor sample.
    """
    create_loh(args.variants, args.min_depth, args.trend, args.sample_id)
    if args.output:
        pyplot.savefig(args.output, format='pdf', bbox_inches="tight")
    else:
        pyplot.show()


def create_loh(variants, min_depth=20, do_trend=False, sample_id=None):
    """Plot allelic frequencies at each variant position in a VCF file."""
    # TODO - accept "normal" sample ID (default: look for "Normal")
    # if so: find het vars in normal (#1), look those up in tumor (#2)
    #   -> plot only those ones

    # TODO - if given segments (args.segment),
    #   do test for significance of shift (abs diff of alt from .5)
    #       between segment and all other points
    #   if significant:
    #       colorize those points
    #       draw a CBS-like bar at mean level in that segment

    _fig, axis = pyplot.subplots()
    axis.set_title("Variant allele frequencies: %s" % variants[0])
    chrom_snvs = ngfrills.load_vcf(variants[0], min_depth, sample_id=sample_id)
    chrom_sizes = collections.OrderedDict()
    for chrom in sorted(chrom_snvs, key=core.sorter_chrom):
        chrom_sizes[chrom] = max(v[0] for v in chrom_snvs[chrom])

    PAD = 2e7
    plots.plot_loh(axis, chrom_snvs, chrom_sizes, do_trend, PAD)


P_loh = AP_subparsers.add_parser('loh', help=_cmd_loh.__doc__)
P_loh.add_argument('variants', nargs='+',
        help="Sample variants in VCF format.")
P_loh.add_argument('-s', '--segment',
        help="Segmentation calls (.cns), the output of the 'segment' command.")
P_loh.add_argument('-m', '--min-depth', type=int, default=20,
        help="""Minimum read depth for a variant to be displayed.
                [Default: %(default)s]""")
P_loh.add_argument("-i", "--sample-id",
        help="Sample name to use for LOH calculations from the input VCF.")
P_loh.add_argument('-t', '--trend', action='store_true',
        help="Draw a smoothed local trendline on the scatter plot.")
P_loh.add_argument('-o', '--output',
        help="Output PDF file name.")
P_loh.set_defaults(func=_cmd_loh)


# heatmap ---------------------------------------------------------------------

def _cmd_heatmap(args):
    """Plot copy number for multiple samples as a heatmap."""
    create_heatmap(args.filenames, args.chromosome, args.desaturate)
    if args.output:
        pyplot.savefig(args.output, format='pdf', bbox_inches="tight")
        echo("Wrote", args.output)
    else:
        pyplot.show()


def create_heatmap(filenames, show_chromosome=None, do_desaturate=False):
    """Plot copy number for multiple samples as a heatmap."""
    # ENH - see the zip magic in _cmd_format
    # Also, for more efficient plotting:
    # http://matplotlib.org/examples/api/span_regions.html

    _fig, axis = pyplot.subplots()

    # List sample names on the y-axis
    axis.set_yticks([i + 0.5 for i in range(len(filenames))])
    axis.set_yticklabels(list(map(core.fbase, filenames)))
    axis.invert_yaxis()
    axis.set_ylabel("Samples")
    axis.set_axis_bgcolor('#DDDDDD')

    # Group each file's probes/segments by chromosome
    sample_data = [collections.defaultdict(list) for _f in filenames]
    for i, fname in enumerate(filenames):
        pset = CNA.read(fname)
        for chrom, subpset in pset.by_chromosome():
            sample_data[i][chrom] = list(zip(subpset['start'], subpset['end'],
                                             subpset['coverage']))

    # Calculate the size (max endpoint value) of each chromosome
    chrom_sizes = {}
    for row in sample_data:
        for chrom, data in iteritems(row):
            max_posn = max(coord[1] for coord in data)
            chrom_sizes[chrom] = max(max_posn, chrom_sizes.get(chrom, 0))
    chrom_sizes = collections.OrderedDict(sorted(iteritems(chrom_sizes),
                                                 key=core.sorter_chrom_at(0)))

    # Closes over do_desaturate, axis
    def plot_rect(y_idx, x_start, x_end, cvg):
        """Draw a rectangle in the given coordinates and color."""
        x_coords = (x_start, x_start, x_end + 1, x_end + 1)
        y_coords = (y_idx, y_idx + 1, y_idx + 1, y_idx)
        rgbcolor = plots.cvg2rgb(cvg, do_desaturate)
        axis.fill(x_coords, y_coords, color=rgbcolor)

    if show_chromosome:
        # Lay out only the selected chromosome
        chrom_offsets = {show_chromosome: 0.0}
        # Set x-axis the chromosomal positions (in Mb), title as the chromosome
        axis.set_xlim(0, chrom_sizes[show_chromosome] * plots.MB)
        axis.set_title(show_chromosome)
        axis.set_xlabel("Position (Mb)")
        axis.tick_params(which='both', direction='out')
        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()
        # Plot the individual probe/segment coverages
        for i, row in enumerate(sample_data):
            for start, end, cvg in row[show_chromosome]:
                plot_rect(i, start * plots.MB, end * plots.MB, cvg)

    else:
        # Lay out chromosome dividers and x-axis labels
        # (Just enough padding to avoid overlap with the divider line)
        chrom_offsets = plots.plot_x_dividers(axis, chrom_sizes, 1)
        # Plot the individual probe/segment coverages
        for i, row in enumerate(sample_data):
            for chrom, curr_offset in iteritems(chrom_offsets):
                for start, end, cvg in row[chrom]:
                    plot_rect(i, start + curr_offset, end + curr_offset, cvg)


P_heatmap = AP_subparsers.add_parser('heatmap', help=_cmd_heatmap.__doc__)
P_heatmap.add_argument('filenames', nargs='+',
        help="Sample coverages as raw probes (.cnr) or segments (.cns).")
P_heatmap.add_argument('-c', '--chromosome',
        help="Name of chromosome to display. [Default: show them all]")
# P_heatmap.add_argument('-g', '--gene',
#         help="Name of gene to display.")
# P_heatmap.add_argument('-r', '--range',
#         help="Chromosomal range to display (e.g. chr1:123000-456000).")
P_heatmap.add_argument('-d', '--desaturate', action='store_true',
        help="Tweak color saturation to focus on significant changes.")
P_heatmap.add_argument('-o', '--output',
        help="Output PDF file name.")
P_heatmap.set_defaults(func=_cmd_heatmap)


# _____________________________________________________________________________
# Tabular outputs

# breaks ----------------------------------------------------------------------

def _cmd_breaks(args):
    """List the targeted genes in which a copy number breakpoint occurs."""
    pset_cvg = CNA.read(args.filename)
    pset_seg = CNA.read(args.segment)
    bpoints = do_breaks(pset_cvg, pset_seg, args.min_probes)
    echo("Found", len(bpoints), "gene breakpoints")
    core.write_tsv(args.output, bpoints,
                   colnames=['Gene', 'Chrom.', 'Location', 'Change',
                             'ProbesLeft', 'ProbesRight'])


def do_breaks(probes, segments, min_probes=1):
    """List the targeted genes in which a copy number breakpoint occurs."""
    intervals = reports.get_gene_intervals(probes)
    bpoints = reports.get_breakpoints(intervals, segments, min_probes)
    return bpoints


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
    pset = CNA.read(args.filename)
    segs = CNA.read(args.segment) if args.segment else None
    gainloss = do_gainloss(pset, segs, args.male_reference, args.threshold,
                           args.min_probes)
    echo("Found", len(gainloss), "gene-level gains and losses")
    core.write_tsv(args.output, gainloss,
                   colnames=['Gene', 'Chrom.', 'Start', 'End', 'Coverage',
                             'Probes'])


def do_gainloss(probes, segments=None, male_reference=False,
                threshold=0.6, min_probes=1):
    """Identify targeted genes with copy number gain or loss."""
    probes = core.shift_xx(probes, male_reference)
    gainloss = []
    if segments:
        segments = core.shift_xx(segments, male_reference)
        for segment, subprobes in probes.by_segment(segments):
            if abs(segment['coverage']) >= threshold:
                for (gene, chrom, start, end, coverages
                    ) in reports.group_by_genes(subprobes):
                    gainloss.append((gene, chrom, start, end,
                                     segment['coverage'], len(coverages)))
    else:
        for gene, chrom, start, end, coverages in reports.group_by_genes(probes):
            gene_coverage = numpy.median(coverages)
            if abs(gene_coverage) >= threshold:
                gainloss.append((gene, chrom, start, end,
                                gene_coverage, len(coverages)))
    return [row for row in gainloss if row[5] >= min_probes]


P_gainloss = AP_subparsers.add_parser('gainloss', help=_cmd_gainloss.__doc__)
P_gainloss.add_argument('filename',
        help="""Processed sample coverage data file (*.cnr), the output
                of the 'fix' sub-command.""")
P_gainloss.add_argument('-s', '--segment',
        help="Segmentation calls (.cns), the output of the 'segment' command).")
P_gainloss.add_argument('-t', '--threshold', type=float, default=0.6,
        help="""Copy number change threshold to report a gene gain/loss.
                [Default: %(default)s]""")
P_gainloss.add_argument('-m', '--min-probes', type=int, default=1,
        help="""Minimum number of covered probes to report a gain/loss.
                [Default: %(default)d]""")
P_gainloss.add_argument('-y', '--male-reference', action='store_true',
        help="""Assume inputs are already corrected against a male
                reference (i.e. female samples will have +1 log-coverage of
                chrX; otherwise male samples would have -1 chrX).""")
P_gainloss.add_argument('-o', '--output',
        help="Output table file name.")
P_gainloss.set_defaults(func=_cmd_gainloss)


# gender ----------------------------------------------------------------------

def _cmd_gender(args):
    """Guess samples' gender from the relative coverage of chromosome X."""
    outrows = []
    for fname in args.targets:
        rel_chrx_cvg = core.get_relative_chrx_cvg(CNA.read(fname))
        if args.male_reference:
            is_xx = (rel_chrx_cvg >= 0.7)
        else:
            is_xx = (rel_chrx_cvg >= -0.4)
        outrows.append((fname,
                        ("Female" if is_xx else "Male"),
                        "%s%.3g" % ('+' if rel_chrx_cvg > 0 else '',
                                    rel_chrx_cvg)))
    core.write_tsv(args.output, outrows)


P_gender = AP_subparsers.add_parser('gender', help=_cmd_gender.__doc__)
P_gender.add_argument('targets', nargs='+',
        help="Copy number or copy ratio files (*.cnn, *.cnr).")
P_gender.add_argument('-y', '--male-reference', action='store_true',
        help="""Assume inputs are already corrected against a male
                reference (i.e. female samples will have +1 log-coverage of
                chrX; otherwise male samples would have -1 chrX).""")
P_gender.add_argument('-o', '--output',
        help="Output table file name.")
P_gender.set_defaults(func=_cmd_gender)


# metrics ---------------------------------------------------------------------

def _cmd_metrics(args):
    """Compute coverage deviations and other metrics for self-evaluation.
    """
    if (len(args.coverages) > 1 and len(args.segments) > 1 and
        len(args.coverages) != len(args.segments)):
        raise ValueError("Number of coverage/segment filenames given must be "
                         "equal, if more than 1 segment file is given.")

    # Repeat a single segment file to match the number of coverage files
    if len(args.coverages) > 1 and len(args.segments) == 1:
        args.segments = [args.segments[0] for _i in range(len(args.coverages))]

    # Calculate all metrics
    outrows = []
    for probes_fname, segs_fname in zip(args.coverages, args.segments):
        probes = CNA.read(probes_fname)
        segments = CNA.read(segs_fname)
        values = metrics.ests_of_scale(
            metrics.probe_deviations_from_segments(probes, segments))
        outrows.append([core.rbase(probes_fname), len(segments)] +
                       ["%.7f" % val for val in values])

    if len(outrows) == 1:
        # Plain-text output for one sample
        sample_id, nseg, stdev, mad, iqr, biweight = outrows[0]
        with ngfrills.safe_write(args.output or sys.stdout) as handle:
            handle.write("Sample: %s\n" % sample_id)
            handle.write("Number of called segments: %d\n" % nseg)
            handle.write("Deviation of probe coverages from segment calls:\n")
            handle.write("  Standard deviation = %s\n" % stdev)
            handle.write("  Median absolute deviation = %s\n" % mad)
            handle.write("  Interquartile range = %s\n" % iqr)
            handle.write("  Biweight midvariance = %s\n" % biweight)
    else:
        # Tabular output for multiple samples
        core.write_tsv(args.output, outrows,
                       colnames=("sample", "segments", "stdev", "mad", "iqr",
                                 "bivar"))


P_metrics = AP_subparsers.add_parser('metrics', help=_cmd_metrics.__doc__)
P_metrics.add_argument('coverages', nargs='+',
        help="""One or more coverage data files (*.cnn, *.cnr).""")
P_metrics.add_argument('-s', '--segments', nargs='+',
        help="""One or more segmentation data files (*.cns, output of the
                'segment' command).  If more than one file is given, the number
                must match the coverage data files, in which case the input
                files will be paired together in the given order. Otherwise, the
                same segments will be used for all coverage files.""")
P_metrics.add_argument('-o', '--output',
        help="Output table file name.")
P_metrics.set_defaults(func=_cmd_metrics)


# _____________________________________________________________________________
# Other I/O and compatibility

# import-picard ---------------------------------------------------------------

def _cmd_import_picard(args):
    """Convert Picard CalculateHsMetrics tabular output to CNVkit .cnn files."""
    for fname in importers.find_picard_files(args.targets):
        cnarr = importers.load_targetcoverage_csv(fname)
        outfname = os.path.basename(fname)[:-4] + '.cnn'
        if args.output_dir:
            if not os.path.isdir(args.output_dir):
                os.mkdir(args.output_dir)
                echo("Created directory", args.output_dir)
            outfname = os.path.join(args.output_dir, outfname)
        cnarr.write(outfname)


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

    for segset in importers.import_seg(args.segfile, chrom_names, args.prefix,
                                       args.from_log10):
        segset.write(os.path.join(args.output_dir, segset.sample_id + '.cns'))


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

def _cmd_import_theta(args):
    """Convert THetA output to a BED-like, CNVkit-like tabular format.

    Equivalently, use the THetA results file to convert CNVkit .cns segments to
    integer copy number calls.
    """
    tumor_segs = CNA.read(args.tumor_cns)
    theta = importers.parse_theta_results(args.theta_results)
    for i, copies in enumerate(theta['C']):
        # Replace segment values with these integers
        # Drop any segments where the C value is None
        new_segs = []
        for seg, ncop in zip(tumor_segs.copy(), copies):
            if ncop is None:
                continue
            seg["coverage"] = math.log((ncop or 0.5) / args.ploidy, 2)
            new_segs.append(seg)
        new_cns = tumor_segs.to_rows(new_segs)
        new_cns.write(os.path.join(args.output_dir,
                                   "%s-%d.cns" % (tumor_segs.sample_id, i + 1)))


P_import_theta = AP_subparsers.add_parser('import-theta',
        help=_cmd_import_theta.__doc__)
P_import_theta.add_argument("tumor_cns")
P_import_theta.add_argument("theta_results")
P_import_theta.add_argument("--ploidy", type=int, default=2,
        help="Ploidy of the sample cells. [Default: %(default)d]")
P_import_theta.add_argument('-d', '--output-dir', default='.',
        help="Output directory name.")
P_import_theta.set_defaults(func=_cmd_import_theta)


# export ----------------------------------------------------------------------


P_export = AP_subparsers.add_parser('export',
        help="""Convert CNVkit output files to another format.""")
P_export_subparsers = P_export.add_subparsers(
        help="Export formats (use with -h for more info).")


# SEG special case: segment coords don't match across samples
def _cmd_export_seg(args):
    """Convert segments to SEG format.

    Compatible with IGV and GenePattern.
    """
    outheader, outrows = export.export_seg(args.filenames)
    core.write_tsv(args.output, outrows, colnames=outheader)

P_export_seg = P_export_subparsers.add_parser('seg', help=_cmd_export_seg.__doc__)
P_export_seg.add_argument('filenames', nargs='+',
        help="""Segmented copy ratio data file(s) (*.cns), the output of the
                'segment' sub-command.""")
P_export_seg.add_argument('-o', '--output', help="Output file name.")
P_export_seg.set_defaults(func=_cmd_export_seg)


# Nexus "basic" special case: can only represent 1 sample
def _cmd_export_nb(args):
    """Convert log2 copy ratios to Nexus Copy Number "basic" format."""
    outheader, outrows = export.export_nexus_basic(args.filename)
    core.write_tsv(args.output, outrows, colnames=outheader)

P_export_nb = P_export_subparsers.add_parser('nexus-basic',
        help=_cmd_export_nb.__doc__)
P_export_nb.add_argument('filename',
        help="""Log2 copy ratio data file (*.cnr), the output of the 'fix'
                sub-command.""")
P_export_nb.add_argument('-o', '--output', help="Output file name.")
P_export_nb.set_defaults(func=_cmd_export_nb)


# FreeBayes special case: multiple samples's segments, like SEG
def _cmd_export_fb(args):
    """Convert segments to FreeBayes --cnv-map format (BED-like).

    Generates an input file for use with FreeBayes's --cnv-map option.

    Input is a segmentation file (.cns). This may be imported from THetA to
    account for normal-cell contamination and subclonal tumor cell population
    structure; otherwise, the --purity argument provides a simpler adjustment.
    """
    outheader, outrows = export.export_freebayes(args.segments, args)
    core.write_tsv(args.output, outrows, colnames=outheader)

P_export_fb = P_export_subparsers.add_parser('freebayes',
        help=_cmd_export_fb.__doc__)
P_export_fb.add_argument('segments', nargs='+',
        help="""Segmented copy ratio data files (*.cns), the output of the
                'segment' sub-command.""")
P_export_fb.add_argument("-i", "--sample-id",
        help="Sample name, as FreeBayes should see it.")
# Arguments that could be shared across 'export'
P_export_fb.add_argument("--ploidy", type=int, default=2,
        help="Ploidy of the sample cells. [Default: %(default)d]")
P_export_fb.add_argument("--purity", type=float,
        help="Estimated tumor cell purity or cellularity.")
P_export_fb.add_argument("-g", "--gender",
        choices=('m', 'male', 'Male', 'f', 'female', 'Female'),
        help="""Specify the sample's gender as male or female. (Otherwise
                guessed from chrX copy number).""")
P_export_fb.add_argument("-y", "--male-reference", action="store_true",
        help="""Was a male reference used?  If so, expect half ploidy on
                chrX and chrY; otherwise, only chrY has half ploidy.  In CNVkit,
                if a male reference was used, the "neutral" copy number (ploidy)
                of chrX is 1; chrY is haploid for either gender reference.""")
# /
P_export_fb.add_argument('-o', '--output', help="Output file name.")
P_export_fb.set_defaults(func=_cmd_export_fb)


# THetA special case: takes tumor .cns and normal .cnr or reference.cnn
def _cmd_export_theta(args):
    """Convert segments to THetA2 input file format (*.input)."""
    outheader, outrows = export.export_theta(args.tumor_segment,
                                             args.normal_reference)
    # if not args.output:
    #     args.output = tumor_segs.sample_id + ".input"
    core.write_tsv(args.output, outrows, colnames=outheader)

P_export_theta = P_export_subparsers.add_parser('theta',
        help=_cmd_export_theta.__doc__)
P_export_theta.add_argument('tumor_segment',
        help="""Tumor-sample segmentation file from CNVkit (.cns).""")
P_export_theta.add_argument('normal_reference',
        help="""Reference copy number profile (.cnn), or normal-sample bin-level
                log2 copy ratios (.cnr).""")
P_export_theta.add_argument('-o', '--output', help="Output file name.")
P_export_theta.set_defaults(func=_cmd_export_theta)


# All else: export any number of .cnr or .cns files

for fmt_key, fmt_descr in (
    ('cdt', "Convert log2 ratios to CDT format. Compatible with Java TreeView."),
    ('jtv', "Convert log2 ratios to Java TreeView's native format."),
    # Not implemented yet:
    # 'multi' (Nexus Copy Number "multi1")
    # 'gct' (GenePattern).
):
    def _cmd_export_simple(args):
        sample_ids = list(map(core.fbase, args.filenames))
        rows = export.merge_samples(args.filenames)
        formatter = export.EXPORT_FORMATS[fmt_key]
        outheader, outrows = formatter(sample_ids, rows)
        core.write_tsv(args.output, outrows, colnames=outheader)

    P_export_simple = P_export_subparsers.add_parser(fmt_key, help=fmt_descr)
    P_export_simple.add_argument('filenames', nargs='+',
            help="""Log2 copy ratio data file(s) (*.cnr), the output of the
                    'fix' sub-command.""")
    P_export_simple.add_argument('-o', '--output', help="Output file name.")
    P_export_simple.set_defaults(func=_cmd_export_simple)



# _____________________________________________________________________________
# Shim for command-line execution

def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)
