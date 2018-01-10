"""The 'batch' command."""
from __future__ import absolute_import, division, print_function

import logging
import os

from matplotlib import pyplot
from skgenome import tabio, GenomicArray as GA

from . import (access, antitarget, autobin, core, coverage, diagram, fix,
               parallel, reference, scatter, segmentation, target)
from .cmdutil import read_cna


def batch_make_reference(normal_bams, target_bed, antitarget_bed,
                         male_reference, fasta, annotate, short_names,
                         target_avg_size, access_bed, antitarget_avg_size,
                         antitarget_min_size, output_reference, output_dir,
                         processes, by_count, method):
    """Build the CN reference from normal samples, targets and antitargets."""
    if method in ("wgs", "amplicon"):
        if antitarget_bed:
            raise ValueError("%r protocol: antitargets should not be "
                             "given/specified." % method)
        if access_bed and target_bed and access_bed != target_bed:
            raise ValueError("%r protocol: targets and access should not be "
                             "different." % method)

    bait_arr = None
    if method == "wgs":
        if not annotate:
            # TODO check if target_bed has gene names
            logging.warning("WGS protocol: recommend '--annotate' option "
                            "(e.g. refFlat.txt) to help locate genes "
                            "in output files.")
        access_arr = None
        if not target_bed:
            # TODO - drop weird contigs before writing, see antitargets.py
            if access_bed:
                target_bed = access_bed
            elif fasta:
                # Run 'access' on the fly
                access_arr = access.do_access(fasta)
                # Take filename base from FASTA, lacking any other clue
                target_bed = os.path.splitext(os.path.basename(fasta)
                                             )[0] + ".bed"
                tabio.write(access_arr, target_bed, "bed3")
            else:
                raise ValueError("WGS protocol: need to provide --targets, "
                                 "--access, or --fasta options.")

        # Tweak default parameters
        if not target_avg_size:
            if normal_bams:
                # Calculate bin size from .bai & access
                if fasta and not access_arr:
                    # Calculate wgs depth from all
                    # sequencing-accessible area (it doesn't take that long
                    # compared to WGS coverage); user-provided access might be
                    # something else that excludes a significant number of
                    # mapped reads.
                    access_arr = access.do_access(fasta)
                if access_arr:
                    autobin_args = ['wgs', None, access_arr]
                else:
                    # Don't assume the given targets/access covers the whole
                    # genome; use autobin sampling to estimate bin size, as we
                    # do for amplicon
                    bait_arr = tabio.read_auto(target_bed)
                    autobin_args = ['amplicon', bait_arr]
                # Choose median-size normal bam or tumor bam
                bam_fname = autobin.midsize_file(normal_bams)
                (wgs_depth, target_avg_size), _ = autobin.do_autobin(
                    bam_fname, *autobin_args, bp_per_bin=50000.)
                logging.info("WGS average depth %.2f --> using bin size %d",
                             wgs_depth, target_avg_size)
            else:
                # This bin size is OK down to 10x
                target_avg_size = 5000

    # To make temporary filenames for processed targets or antitargets
    tgt_name_base, _tgt_ext = os.path.splitext(os.path.basename(target_bed))
    if output_dir:
        tgt_name_base = os.path.join(output_dir, tgt_name_base)

    # Pre-process baits/targets
    new_target_fname = tgt_name_base + '.target.bed'
    if bait_arr is None:
        bait_arr = tabio.read_auto(target_bed)
    target_arr = target.do_target(bait_arr, annotate, short_names, True,
                                  **({'avg_size': target_avg_size}
                                     if target_avg_size
                                     else {}))
    tabio.write(target_arr, new_target_fname, 'bed4')
    target_bed = new_target_fname

    if not antitarget_bed:
        # Devise a temporary antitarget filename
        antitarget_bed = tgt_name_base + '.antitarget.bed'
        if method == "hybrid":
            # Build antitarget BED from the given targets
            anti_kwargs = {}
            if access_bed:
                anti_kwargs['access'] = tabio.read_auto(access_bed)
            if antitarget_avg_size:
                anti_kwargs['avg_bin_size'] = antitarget_avg_size
            if antitarget_min_size:
                anti_kwargs['min_bin_size'] = antitarget_min_size
            anti_arr = antitarget.do_antitarget(target_arr, **anti_kwargs)
        else:
            # No antitargets for wgs, amplicon
            anti_arr = GA([])
        tabio.write(anti_arr, antitarget_bed, "bed4")

    if len(normal_bams) == 0:
        logging.info("Building a flat reference...")
        ref_arr = reference.do_reference_flat(target_bed, antitarget_bed, fasta,
                                              male_reference)
    else:
        logging.info("Building a copy number reference from normal samples...")
        # Run coverage on all normals
        with parallel.pick_pool(processes) as pool:
            tgt_futures = []
            anti_futures = []
            procs_per_cnn = max(1, processes // (2 * len(normal_bams)))
            for nbam in normal_bams:
                sample_id = core.fbase(nbam)
                sample_pfx = os.path.join(output_dir, sample_id)
                tgt_futures.append(
                    pool.submit(batch_write_coverage,
                                target_bed, nbam,
                                sample_pfx + '.targetcoverage.cnn',
                                by_count, procs_per_cnn))
                anti_futures.append(
                    pool.submit(batch_write_coverage,
                                antitarget_bed, nbam,
                                sample_pfx + '.antitargetcoverage.cnn',
                                by_count, procs_per_cnn))

        target_fnames = [tf.result() for tf in tgt_futures]
        antitarget_fnames = [af.result() for af in anti_futures]
        # Build reference from *.cnn
        ref_arr = reference.do_reference(target_fnames, antitarget_fnames,
                                         fasta, male_reference, None,
                                         do_gc=True,
                                         do_edge=(method == "hybrid"),
                                         do_rmask=True)
    if not output_reference:
        output_reference = os.path.join(output_dir, "reference.cnn")
    core.ensure_path(output_reference)
    tabio.write(ref_arr, output_reference)
    return output_reference, target_bed, antitarget_bed


def batch_write_coverage(bed_fname, bam_fname, out_fname, by_count, processes):
    """Run coverage on one sample, write to file."""
    cnarr = coverage.do_coverage(bed_fname, bam_fname, by_count, 0, processes)
    tabio.write(cnarr, out_fname)
    return out_fname


def batch_run_sample(bam_fname, target_bed, antitarget_bed, ref_fname,
                     output_dir, male_reference, plot_scatter, plot_diagram,
                     rlibpath, by_count, skip_low, method, processes):
    """Run the pipeline on one BAM file."""
    # ENH - return probes, segments (cnarr, segarr)
    logging.info("Running the CNVkit pipeline on %s ...", bam_fname)
    sample_id = core.fbase(bam_fname)
    sample_pfx = os.path.join(output_dir, sample_id)

    raw_tgt = coverage.do_coverage(target_bed, bam_fname, by_count, 0,
                                   processes)
    tabio.write(raw_tgt, sample_pfx + '.targetcoverage.cnn')

    raw_anti = coverage.do_coverage(antitarget_bed, bam_fname, by_count, 0,
                                    processes)
    tabio.write(raw_anti, sample_pfx + '.antitargetcoverage.cnn')

    cnarr = fix.do_fix(raw_tgt, raw_anti, read_cna(ref_fname),
                       do_gc=True, do_edge=(method == "hybrid"), do_rmask=True)
    tabio.write(cnarr, sample_pfx + '.cnr')

    logging.info("Segmenting %s.cnr ...", sample_pfx)
    segments = segmentation.do_segmentation(cnarr, 'cbs',
                                            rlibpath=rlibpath,
                                            skip_low=skip_low,
                                            processes=processes,
                                            **({'threshold': 1e-6}
                                               if method == 'wgs'
                                               else {}))
    tabio.write(segments, sample_pfx + '.cns')

    if plot_scatter:
        scatter.do_scatter(cnarr, segments)
        pyplot.savefig(sample_pfx + '-scatter.pdf', format='pdf',
                       bbox_inches="tight")
        logging.info("Wrote %s-scatter.pdf", sample_pfx)

    if plot_diagram:
        is_xx = cnarr.guess_xx(male_reference)
        outfname = sample_pfx + '-diagram.pdf'
        diagram.create_diagram(cnarr.shift_xx(male_reference, is_xx),
                               segments.shift_xx(male_reference, is_xx),
                               0.5, 3, outfname)
        logging.info("Wrote %s", outfname)
