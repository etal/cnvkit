"""Supporting functions for the 'reference' command."""
from __future__ import absolute_import, division, print_function
from builtins import map, zip

import collections
import logging

import numpy as np
import pyfaidx
from skgenome import tabio, GenomicArray as GA

from . import core, fix, descriptives, params
from .cmdutil import read_cna
from .cnary import CopyNumArray as CNA


def do_reference(target_fnames, antitarget_fnames=None, fa_fname=None,
                 male_reference=False, female_samples=None,
                 do_gc=True, do_edge=True, do_rmask=True):
    """Compile a coverage reference from the given files (normal samples)."""
    if antitarget_fnames:
        core.assert_equal("Unequal number of target and antitarget files given",
                          targets=len(target_fnames),
                          antitargets=len(antitarget_fnames))
    if not fa_fname:
        logging.info("No FASTA reference genome provided; "
                     "skipping GC, RM calculations")

    if female_samples is None:
        # NB: Antitargets are usually preferred for inferring sex, but might be
        # empty files, in which case no inference can be done. Since targets are
        # guaranteed to exist, infer from those first, then replace those
        # values where antitargets are suitable.
        sexes = infer_sexes(target_fnames, male_reference)
        if antitarget_fnames:
            a_sexes = infer_sexes(antitarget_fnames, male_reference)
            for sid, a_is_xx in a_sexes.items():
                t_is_xx = sexes.get(sid)
                if t_is_xx is None:
                    sexes[sid] = a_is_xx
                elif t_is_xx != a_is_xx and a_is_xx is not None:
                    logging.warning("Sample %s chromosomal X/Y ploidy looks "
                                    "like %s in targets but %s in antitargets; "
                                    "preferring antitargets",
                                    sid,
                                    "female" if t_is_xx else "male",
                                    "female" if a_is_xx else "male")
                    sexes[sid] = a_is_xx
    else:
        sexes = collections.defaultdict(lambda: female_samples)

    # Calculate & save probe centers
    ref_probes = combine_probes(target_fnames, fa_fname,
                                male_reference, sexes, True,
                                do_gc, do_edge, False)
    if antitarget_fnames:
        ref_probes.add(combine_probes(antitarget_fnames, fa_fname,
                                      male_reference, sexes, False,
                                      do_gc, False, do_rmask))
    ref_probes.center_all(skip_low=True)
    ref_probes.sort_columns()
    warn_bad_bins(ref_probes)
    return ref_probes


def do_reference_flat(targets, antitargets=None, fa_fname=None,
                      male_reference=False):
    """Compile a neutral-coverage reference from the given intervals.

    Combines the intervals, shifts chrX values if requested, and calculates GC
    and RepeatMasker content from the genome FASTA sequence.
    """
    ref_probes = bed2probes(targets)
    if antitargets:
        ref_probes.add(bed2probes(antitargets))
    # Set sex chromosomes by "reference" sex
    ref_probes['log2'] = ref_probes.expect_flat_log2(male_reference)
    ref_probes['depth'] = np.exp2(ref_probes['log2'])  # Shim
    # Calculate GC and RepeatMasker content for each probe's genomic region
    if fa_fname:
        gc, rmask = get_fasta_stats(ref_probes, fa_fname)
        ref_probes['gc'] = gc
        ref_probes['rmask'] = rmask
        # warn_bad_bins(ref_probes)
    else:
        logging.info("No FASTA reference genome provided; "
                     "skipping GC, RM calculations")
    ref_probes.sort_columns()
    return ref_probes


def bed2probes(bed_fname):
    """Create a neutral-coverage CopyNumArray from a file of regions."""
    regions = tabio.read_auto(bed_fname)
    table = regions.data.loc[:, ("chromosome", "start", "end")]
    table["gene"] = (regions.data["gene"] if "gene" in regions.data else '-')
    table["log2"] = 0.0
    table["spread"] = 0.0
    return CNA(table, {"sample_id": core.fbase(bed_fname)})


def infer_sexes(cnn_fnames, is_male_reference):
    """Map sample IDs to inferred chromosomal sex, where possible.

    For samples where the source file is empty or does not include either sex
    chromosome, that sample ID will not be in the returned dictionary.
    """
    sexes = {}
    for fname in cnn_fnames:
        cnarr = read_cna(fname)
        if cnarr:
            is_xx = cnarr.guess_xx(is_male_reference)
            if is_xx is not None:
                sexes[cnarr.sample_id] = is_xx
    return sexes


def combine_probes(filenames, fa_fname, is_male_reference, sexes, skip_low,
                   fix_gc, fix_edge, fix_rmask):
    """Calculate the median coverage of each bin across multiple samples.

    Parameters
    ----------
    filenames : list
        List of string filenames, corresponding to targetcoverage.cnn and
        antitargetcoverage.cnn files, as generated by 'coverage' or
        'import-picard'.
    fa_fname : str
        Reference genome sequence in FASTA format, used to extract GC and
        RepeatMasker content of each genomic bin.
    is_male_reference : bool
    skip_low : bool
    fix_gc : bool
    fix_edge : bool
    fix_rmask : bool

    Returns
    -------
    CopyNumArray
        One object summarizing the coverages of the input samples, including
        each bin's "average" coverage, "spread" of coverages, and GC content.
    """
    columns = {}

    # Load coverage from target/antitarget files
    logging.info("Loading %s", filenames[0])
    cnarr1 = read_cna(filenames[0])
    if not len(cnarr1):
        # Just create an empty array with the right columns
        col_names = ['chromosome', 'start', 'end', 'gene', 'log2', 'depth']
        if 'gc' in cnarr1 or fa_fname:
            col_names.append('gc')
        if fa_fname:
            col_names.append('rmask')
        col_names.append('spread')
        return CNA.from_rows([], col_names, {'sample_id': "reference"})

    # Calculate GC and RepeatMasker content for each probe's genomic region
    if fa_fname and (fix_rmask or fix_gc):
        gc, rmask = get_fasta_stats(cnarr1, fa_fname)
        if fix_gc:
            columns['gc'] = gc
        if fix_rmask:
            columns['rmask'] = rmask
    elif 'gc' in cnarr1 and fix_gc:
        # Reuse .cnn GC values if they're already stored (via import-picard)
        gc = cnarr1['gc']
        columns['gc'] = gc

    # Make the sex-chromosome coverages of male and female samples compatible
    is_chr_x = (cnarr1.chromosome == cnarr1._chr_x_label)
    is_chr_y = (cnarr1.chromosome == cnarr1._chr_y_label)
    flat_coverage = cnarr1.expect_flat_log2(is_male_reference)
    def shift_sex_chroms(cnarr):
        """Shift sample X and Y chromosomes to match the reference sex.

        Reference values:
            XY: chrX -1, chrY -1
            XX: chrX 0, chrY -1

        Plan:
          chrX:
            xx sample, xx ref: 0    (from 0)
            xx sample, xy ref: -= 1 (from -1)
            xy sample, xx ref: += 1 (from 0)    +1
            xy sample, xy ref: 0    (from -1)   +1
          chrY:
            xx sample, xx ref: = -1 (from -1)
            xx sample, xy ref: = -1 (from -1)
            xy sample, xx ref: 0    (from -1)   +1
            xy sample, xy ref: 0    (from -1)   +1

        """
        is_xx = sexes.get(cnarr.sample_id)
        cnarr['log2'] += flat_coverage
        if is_xx:
            # chrX has same ploidy as autosomes; chrY is just unusable noise
            cnarr[is_chr_y, 'log2'] = -1.0  # np.nan is worse
        else:
            # 1/2 #copies of each sex chromosome
            cnarr[is_chr_x | is_chr_y, 'log2'] += 1.0

    edge_bias = fix.get_edge_bias(cnarr1, params.INSERT_SIZE)
    def bias_correct_coverage(cnarr):
        """Perform bias corrections on the sample."""
        cnarr.center_all(skip_low=skip_low)
        shift_sex_chroms(cnarr)
        # Skip bias corrections if most bins have no coverage (e.g. user error)
        if (cnarr['log2'] > params.NULL_LOG2_COVERAGE - params.MIN_REF_COVERAGE
           ).sum() <= len(cnarr) // 2:
            logging.warning("WARNING: most bins have no or very low coverage; "
                            "check that the right BED file was used")
        else:
            if 'gc' in columns and fix_gc:
                logging.info("Correcting for GC bias...")
                cnarr = fix.center_by_window(cnarr, .1, columns['gc'])
            if 'rmask' in columns and fix_rmask:
                logging.info("Correcting for RepeatMasker bias...")
                cnarr = fix.center_by_window(cnarr, .1, columns['rmask'])
            if fix_edge:
                logging.info("Correcting for density bias...")
                cnarr = fix.center_by_window(cnarr, .1, edge_bias)
        return cnarr['log2']

    # Pseudocount of 1 "flat" sample
    all_depths = [cnarr1['depth'] if 'depth' in cnarr1
                  else np.exp2(cnarr1['log2'])]
    all_coverages = [flat_coverage, bias_correct_coverage(cnarr1)]
    for fname in filenames[1:]:
        logging.info("Loading target %s", fname)
        cnarrx = read_cna(fname)
        # Bin information should match across all files
        if not (len(cnarr1) == len(cnarrx)
                and (cnarr1.chromosome == cnarrx.chromosome).all()
                and (cnarr1.start == cnarrx.start).all()
                and (cnarr1.end == cnarrx.end).all()
                and (cnarr1['gene'] == cnarrx['gene']).all()):
            raise RuntimeError("%s bins do not match those in %s"
                               % (fname, filenames[0]))
        all_depths.append(cnarrx['depth'] if 'depth' in cnarrx
                          else np.exp2(cnarrx['log2']))
        all_coverages.append(bias_correct_coverage(cnarrx))
    all_coverages = np.vstack(all_coverages)

    logging.info("Calculating average bin coverages")
    cvg_centers = np.apply_along_axis(descriptives.biweight_location, 0,
                                      all_coverages)
    depth_centers = np.apply_along_axis(descriptives.biweight_location, 0,
                                        np.vstack(all_depths))
    logging.info("Calculating bin spreads")
    spreads = np.array([descriptives.biweight_midvariance(a, initial=i)
                        for a, i in zip(all_coverages.T, cvg_centers)])
    columns.update({
        'chromosome': cnarr1.chromosome,
        'start': cnarr1.start,
        'end': cnarr1.end,
        'gene': cnarr1['gene'],
        'log2': cvg_centers,
        'depth': depth_centers,
        'spread': spreads,
    })
    return CNA.from_columns(columns, {'sample_id': "reference"})


def warn_bad_bins(cnarr, max_name_width=50):
    """Warn about target bins where coverage is poor.

    Prints a formatted table to stderr.
    """
    bad_bins = cnarr[fix.mask_bad_bins(cnarr)]
    fg_index = ~bad_bins['gene'].isin(params.ANTITARGET_ALIASES)
    fg_bad_bins = bad_bins[fg_index]
    if len(fg_bad_bins) > 0:
        bad_pct = (100 * len(fg_bad_bins)
                   / sum(~cnarr['gene'].isin(params.ANTITARGET_ALIASES)))
        logging.info("Targets: %d (%s) bins failed filters:",
                     len(fg_bad_bins), "%.4f" % bad_pct + '%')
        gene_cols = min(max_name_width, max(map(len, fg_bad_bins['gene'])))
        labels = fg_bad_bins.labels()
        chrom_cols = max(labels.apply(len))
        last_gene = None
        for label, probe in zip(labels, fg_bad_bins):
            if probe.gene == last_gene:
                gene = '  "'
            else:
                gene = probe.gene
                last_gene = gene
            if len(gene) > max_name_width:
                gene = gene[:max_name_width-3] + '...'
            if 'rmask' in cnarr:
                logging.info("  %s  %s  log2=%.3f  spread=%.3f  rmask=%.3f",
                             gene.ljust(gene_cols), label.ljust(chrom_cols),
                             probe.log2, probe.spread, probe.rmask)
            else:
                logging.info("  %s  %s  log2=%.3f  spread=%.3f",
                             gene.ljust(gene_cols), label.ljust(chrom_cols),
                             probe.log2, probe.spread)

    # Count the number of BG bins dropped, too (names are all "Antitarget")
    bg_bad_bins = bad_bins[~fg_index]
    if len(bg_bad_bins) > 0:
        bad_pct = (100 * len(bg_bad_bins)
                   / sum(cnarr['gene'].isin(params.ANTITARGET_ALIASES)))
        logging.info("Antitargets: %d (%s) bins failed filters",
                     len(bg_bad_bins), "%.4f" % bad_pct + '%')


def get_fasta_stats(cnarr, fa_fname):
    """Calculate GC and RepeatMasker content of each bin in the FASTA genome."""
    logging.info("Calculating GC and RepeatMasker content in %s ...", fa_fname)
    gc_rm_vals = [calculate_gc_lo(subseq)
                  for subseq in fasta_extract_regions(fa_fname, cnarr)]
    gc_vals, rm_vals = zip(*gc_rm_vals)
    return np.asfarray(gc_vals), np.asfarray(rm_vals)


def calculate_gc_lo(subseq):
    """Calculate the GC and lowercase (RepeatMasked) content of a string."""
    cnt_at_lo = subseq.count('a') + subseq.count('t')
    cnt_at_up = subseq.count('A') + subseq.count('T')
    cnt_gc_lo = subseq.count('g') + subseq.count('c')
    cnt_gc_up = subseq.count('G') + subseq.count('C')
    tot = float(cnt_gc_up + cnt_gc_lo + cnt_at_up + cnt_at_lo)
    if not tot:
        return 0.0, 0.0
    frac_gc = (cnt_gc_lo + cnt_gc_up) / tot
    frac_lo = (cnt_at_lo + cnt_gc_lo) / tot
    return frac_gc, frac_lo


def fasta_extract_regions(fa_fname, intervals):
    """Extract an iterable of regions from an indexed FASTA file.

    Input: FASTA file name; iterable of (seq_id, start, end) (1-based)
    Output: iterable of string sequences.
    """
    with pyfaidx.Fasta(fa_fname, as_raw=True) as fa_file:
        for chrom, subarr in intervals.by_chromosome():
            logging.info("Extracting sequences from chromosome %s", chrom)
            for _chrom, start, end in subarr.coords():
                yield fa_file[_chrom][int(start):int(end)]


def reference2regions(refarr):
    """Split reference into target and antitarget regions."""
    is_bg = (refarr['gene'].isin(params.ANTITARGET_ALIASES))
    regions = GA(refarr.data.loc[:, ('chromosome', 'start', 'end', 'gene')],
                 {'sample_id': 'reference'})
    targets = regions[~is_bg]
    antitargets = regions[is_bg]
    return targets, antitargets
