"""The 'batch' command."""

from __future__ import annotations
import logging
import os

from matplotlib import pyplot
from skgenome import tabio, GenomicArray as GA

from . import (
    access,
    antitarget,
    autobin,
    bintest,
    call,
    core,
    coverage,
    diagram,
    fix,
    parallel,
    reference,
    scatter,
    segmentation,
    segmetrics,
    target,
)
from .cmdutil import read_cna


def batch_make_reference(
    normal_fnames: list[str],
    target_bed: str | None,
    antitarget_bed: str | None,
    is_haploid_x_reference: bool,
    diploid_parx_genome: str | None,
    fasta: str | None,
    annotate: str | None,
    short_names: bool,
    target_avg_size: int,
    access_bed: str | None,
    antitarget_avg_size: int | None,
    antitarget_min_size: int | None,
    output_reference: str | None,
    output_dir: str,
    processes: int,
    by_count: bool,
    min_mapq: int,
    method: str,
    do_cluster: bool,
) -> tuple[str, str, str]:
    """Build a complete copy number reference from normal samples.

    This is the high-level workflow function used by the `batch` command to
    create a reference. It coordinates target/antitarget preparation, coverage
    calculation across normal samples, and reference construction.

    The function handles three sequencing methods: hybrid capture (default),
    whole genome sequencing (WGS), and targeted amplicon sequencing.

    Parameters
    ----------
    normal_fnames : list of str
        Paths to BAM or .bed.gz files from normal/control samples. If empty, creates a
        "flat" reference with uniform expected coverage.
    target_bed : str, optional
        Path to BED file defining target/baited regions. Required for hybrid
        capture and amplicon. Optional for WGS (can be auto-generated).
    antitarget_bed : None
        Path to pre-computed antitarget BED file. If None (typical), antitargets
        are automatically generated for hybrid capture.
    is_haploid_x_reference : bool
        True if reference samples are male (haploid X chromosome).
    diploid_parx_genome : str, optional
        Reference genome name (e.g., 'hg19', 'hg38') for pseudo-autosomal region
        handling on X/Y chromosomes.
    fasta : str
        Path to reference genome FASTA file. Required for GC/RepeatMasker
        calculations and for auto-generating WGS targets.
    annotate : str, optional
        Path to gene annotation file (e.g., refFlat.txt) to add gene names
        to targets. Recommended for WGS.
    short_names : bool
        Shorten target names to gene symbols when possible.
    target_avg_size : int
        Target average bin size in base pairs. If 0 or None, automatically
        determined. For WGS, defaults to 5000 bp bins.
    access_bed : None
        Path to sequencing-accessible regions BED file (e.g., from `access`
        command). Used for WGS to exclude low-mappability regions.
    antitarget_avg_size : int, optional
        Average size for antitarget bins in bp (hybrid capture only).
    antitarget_min_size : int, optional
        Minimum size for antitarget bins in bp (hybrid capture only).
    output_reference : None
        Path for output reference file. If None, writes to
        "{output_dir}/reference.cnn".
    output_dir : str
        Directory for output files (intermediate coverage files, reference).
    processes : int
        Number of parallel processes for coverage calculation.
        0 = use all available CPUs.
    by_count : bool
        Calculate coverage by read count instead of read depth.
    method : str
        Sequencing protocol: 'hybrid' (default), 'wgs', or 'amplicon'.
        Determines target/antitarget handling and default parameters.
    do_cluster : bool
        Apply hierarchical clustering to identify and separate reference
        sample subgroups (useful for heterogeneous normal cohorts).

    Returns
    -------
    tuple of (str, str, str)
        Paths to: (reference .cnn file, targets .bed file, antitargets .bed file)

    Notes
    -----
    The reference creation workflow:

    1. **Target preparation** (protocol-dependent):

       - **Hybrid capture**: Uses provided target_bed, generates antitargets
       - **WGS**: Auto-generates genome-wide bins, no antitargets
       - **Amplicon**: Uses provided target_bed, no antitargets

    2. **Coverage calculation** (parallel across normal samples):

       - Runs `coverage` command on each normal BAM
       - Creates .targetcoverage.cnn and .antitargetcoverage.cnn files
       - Parallelized across samples and within each coverage calculation

    3. **Reference construction**:

       - Pools coverage across normal samples
       - Calculates GC content, RepeatMasker, edge effects from FASTA
       - Optionally clusters samples and creates sub-references
       - Outputs reference.cnn file

    **Flat reference**: If no normal BAMs provided, creates a uniform reference
    assuming equal coverage in all bins. Less accurate but useful when no
    normals are available.

    **Protocol differences**:

    - **hybrid**: Uses both targets and antitargets, edge correction enabled
    - **wgs**: Large genome-wide bins, no antitargets, no edge correction
    - **amplicon**: Uses only targets, no antitargets, no edge correction

    Raises
    ------
    ValueError
        If method='wgs' or 'amplicon' but antitargets are specified.
        If method='wgs' and neither targets, access, nor fasta provided.
        If method='wgs' and targets != access (when both given).

    See Also
    --------
    do_reference : Creates reference from coverage files
    do_reference_flat : Creates flat reference without normal samples
    batch_write_coverage : Calculates coverage for a single BAM
    autobin : Automatically determines optimal bin size

    Examples
    --------
    Hybrid capture with normal samples:
    >>> ref, tgt, anti = batch_make_reference(
    ...     normal_fnames=['N1.bam', 'N2.bam', 'N3.bam'],
    ...     target_bed='baits.bed',
    ...     antitarget_bed=None,  # Auto-generate
    ...     is_haploid_x_reference=False,
    ...     diploid_parx_genome='hg38',
    ...     fasta='hg38.fa',
    ...     annotate='refFlat.txt',
    ...     short_names=True,
    ...     target_avg_size=0,  # Auto-determine
    ...     access_bed=None,
    ...     antitarget_avg_size=None,
    ...     antitarget_min_size=None,
    ...     output_reference=None,
    ...     output_dir='output/',
    ...     processes=8,
    ...     by_count=False,
    ...     method='hybrid',
    ...     do_cluster=False
    ... )

    WGS without normal samples (flat reference):
    >>> ref, tgt, anti = batch_make_reference(
    ...     normal_fnames=[],
    ...     target_bed=None,  # Auto-generate from FASTA
    ...     antitarget_bed=None,
    ...     is_haploid_x_reference=True,
    ...     diploid_parx_genome='hg19',
    ...     fasta='hg19.fa',
    ...     annotate='refFlat.txt',
    ...     short_names=False,
    ...     target_avg_size=5000,
    ...     access_bed='access-5k.hg19.bed',
    ...     antitarget_avg_size=None,
    ...     antitarget_min_size=None,
    ...     output_reference='flat_reference.cnn',
    ...     output_dir='wgs_ref/',
    ...     processes=1,
    ...     by_count=False,
    ...     method='wgs',
    ...     do_cluster=False
    ... )
    """
    if method in ("wgs", "amplicon"):
        if antitarget_bed:
            raise ValueError(
                r"{method!r} protocol: antitargets should not be given/specified."
            )
        if access_bed and target_bed and access_bed != target_bed:
            raise ValueError(
                f"{method!r} protocol: targets and access should not be different."
            )

    bait_arr = None
    if method == "wgs":
        if not annotate:
            # TODO check if target_bed has gene names
            logging.warning(
                "WGS protocol: recommend '--annotate' option "
                "(e.g. refFlat.txt) to help locate genes "
                "in output files."
            )
        access_arr = None
        if not target_bed:
            # TODO - drop weird contigs before writing, see antitargets.py
            if access_bed:
                target_bed = access_bed
            elif fasta:
                # Run 'access' on the fly
                access_arr = access.do_access(fasta)
                # Take filename base from FASTA, lacking any other clue
                target_bed = os.path.splitext(os.path.basename(fasta))[0] + ".bed"
                if output_dir:
                    target_bed = os.path.join(output_dir, target_bed)

                tabio.write(access_arr, target_bed, "bed3")
            else:
                raise ValueError(
                    "WGS protocol: need to provide --targets, --access, or "
                    "--fasta options."
                )

        # Tweak default parameters
        if not target_avg_size:
            if normal_fnames:
                # Calculate bin size from .bai & access
                if fasta and not access_arr:
                    # Calculate wgs depth from all
                    # sequencing-accessible area (it doesn't take that long
                    # compared to WGS coverage); user-provided access might be
                    # something else that excludes a significant number of
                    # mapped reads.
                    access_arr = access.do_access(fasta)
                if access_arr:
                    autobin_args = ["wgs", None, access_arr]
                else:
                    # Don't assume the given targets/access covers the whole
                    # genome; use autobin sampling to estimate bin size, as we
                    # do for amplicon
                    bait_arr = tabio.read_auto(target_bed)
                    autobin_args = ["amplicon", bait_arr]
                # Choose median-size normal or tumor sample file
                bam_fname = autobin.midsize_file(normal_fnames)
                (wgs_depth, target_avg_size), _ = autobin.do_autobin(  # type: ignore[assignment,misc]
                    bam_fname,
                    *autobin_args,  # type: ignore[arg-type]
                    bp_per_bin=50000.0,
                    fasta=fasta,
                )
                logging.info(
                    "WGS average depth %.2f --> using bin size %d",
                    wgs_depth,
                    target_avg_size,
                )
            else:
                # This bin size is OK down to 10x
                target_avg_size = 5000

    # To make temporary filenames for processed targets or antitargets
    assert target_bed is not None
    tgt_name_base, _tgt_ext = os.path.splitext(os.path.basename(target_bed))
    if output_dir:
        tgt_name_base = os.path.join(output_dir, tgt_name_base)

    # Pre-process baits/targets
    new_target_fname = tgt_name_base + ".target.bed"
    if bait_arr is None:
        bait_arr = tabio.read_auto(target_bed)
    target_arr = target.do_target(
        bait_arr,
        annotate,
        short_names,
        True,
        **({"avg_size": target_avg_size} if target_avg_size else {}),
    )
    tabio.write(target_arr, new_target_fname, "bed4")
    target_bed = new_target_fname

    if not antitarget_bed:
        # Devise a temporary antitarget filename
        antitarget_bed = tgt_name_base + ".antitarget.bed"
        if method == "hybrid":
            # Build antitarget BED from the given targets
            anti_kwargs: dict = {}
            if access_bed:
                anti_kwargs["access"] = tabio.read_auto(access_bed)
            if antitarget_avg_size:
                anti_kwargs["avg_bin_size"] = antitarget_avg_size
            if antitarget_min_size:
                anti_kwargs["min_bin_size"] = antitarget_min_size
            anti_arr = antitarget.do_antitarget(target_arr, **anti_kwargs)  # type: ignore[arg-type]
        else:
            # No antitargets for wgs, amplicon
            anti_arr = GA([])
        tabio.write(anti_arr, antitarget_bed, "bed4")

    assert antitarget_bed is not None
    if len(normal_fnames) == 0:
        logging.info("Building a flat reference...")
        ref_arr = reference.do_reference_flat(
            target_bed,
            antitarget_bed,
            fasta,
            is_haploid_x_reference,
            diploid_parx_genome,
        )
    else:
        logging.info("Building a copy number reference from normal samples...")
        # Run coverage on all normals
        with parallel.pick_pool(processes) as pool:
            tgt_futures = []
            anti_futures = []
            procs_per_cnn = max(1, processes // (2 * len(normal_fnames)))
            for nbam in normal_fnames:
                sample_id = core.fbase(nbam)
                sample_pfx = os.path.join(output_dir, sample_id)
                tgt_futures.append(
                    pool.submit(
                        batch_write_coverage,
                        target_bed,
                        nbam,
                        sample_pfx + ".targetcoverage.cnn",
                        by_count,
                        min_mapq,
                        procs_per_cnn,
                        fasta,
                    )
                )
                anti_futures.append(
                    pool.submit(
                        batch_write_coverage,
                        antitarget_bed,
                        nbam,
                        sample_pfx + ".antitargetcoverage.cnn",
                        by_count,
                        min_mapq,
                        procs_per_cnn,
                        fasta,
                    )
                )

        target_fnames = [tf.result() for tf in tgt_futures]
        antitarget_fnames = [af.result() for af in anti_futures]
        # Build reference from *.cnn
        ref_arr = reference.do_reference(
            target_fnames,
            antitarget_fnames,
            fasta,
            is_haploid_x_reference,
            diploid_parx_genome,
            None,
            do_gc=True,
            do_edge=(method == "hybrid"),
            do_rmask=True,
            do_cluster=do_cluster,
        )
    if not output_reference:
        output_reference = os.path.join(output_dir, "reference.cnn")
    core.ensure_path(output_reference)
    tabio.write(ref_arr, output_reference)
    return output_reference, target_bed, antitarget_bed


def batch_write_coverage(
    bed_fname: str,
    sample_fname: str,
    out_fname: str,
    by_count: bool,
    min_mapq: int,
    processes: int,
    fasta: str | None,
) -> str:
    """Run coverage on one sample (BAM or bedGraph), write to file."""
    cnarr = coverage.do_coverage(
        bed_fname, sample_fname, by_count, min_mapq, processes, fasta
    )
    tabio.write(cnarr, out_fname)
    return out_fname


def batch_run_sample(
    sample_fname: str,
    target_bed: str,
    antitarget_bed: str,
    ref_fname: str,
    output_dir: str,
    is_haploid_x_reference: bool,
    diploid_parx_genome: str | None,
    plot_scatter: bool,
    plot_diagram: bool,
    rscript_path: str,
    by_count: bool,
    min_mapq: int,
    skip_low: bool,
    seq_method: str,
    segment_method: str,
    processes: int,
    do_cluster: bool,
    fasta: str | None = None,
) -> None:
    """Run the pipeline on one sample (BAM or bedGraph file)."""
    # ENH - return probes, segments (cnarr, segarr)
    logging.info("Running the CNVkit pipeline on %s ...", sample_fname)
    sample_id = core.fbase(sample_fname)
    sample_pfx = os.path.join(output_dir, sample_id)

    raw_tgt = coverage.do_coverage(
        target_bed, sample_fname, by_count, min_mapq, processes, fasta
    )
    tabio.write(raw_tgt, sample_pfx + ".targetcoverage.cnn")

    raw_anti = coverage.do_coverage(
        antitarget_bed, sample_fname, by_count, min_mapq, processes, fasta
    )
    tabio.write(raw_anti, sample_pfx + ".antitargetcoverage.cnn")

    cnarr = fix.do_fix(
        raw_tgt,
        raw_anti,
        read_cna(ref_fname),
        diploid_parx_genome,
        do_gc=True,
        do_edge=(seq_method == "hybrid"),
        do_rmask=True,
        do_cluster=do_cluster,
    )
    tabio.write(cnarr, sample_pfx + ".cnr")

    logging.info("Segmenting %s.cnr ...", sample_pfx)
    segments = segmentation.do_segmentation(
        cnarr,
        segment_method,
        diploid_parx_genome,
        rscript_path=rscript_path,
        skip_low=skip_low,
        processes=processes,
        **({"threshold": 1e-6} if seq_method == "wgs" else {}),  # type: ignore[arg-type]
    )

    logging.info("Post-processing %s.cns ...", sample_pfx)
    # TODO/ENH take centering shift & apply to .cnr for use in segmetrics
    seg_metrics = segmetrics.do_segmetrics(
        cnarr,
        segments,  # type: ignore[arg-type]
        interval_stats=["ci"],
        alpha=0.5,
        smoothed=True,
        skip_low=skip_low,
    )
    tabio.write(seg_metrics, sample_pfx + ".cns")

    # Remove likely false-positive breakpoints
    seg_call = call.do_call(seg_metrics, method="none", filters=["ci"])
    # Calculate another segment-level test p-value
    seg_alltest = segmetrics.do_segmetrics(
        cnarr, seg_call, location_stats=["p_ttest"], skip_low=skip_low
    )
    # Finally, assign absolute copy number values to each segment
    seg_alltest.center_all("median", diploid_parx_genome=diploid_parx_genome)
    seg_final = call.do_call(seg_alltest, method="threshold")
    tabio.write(seg_final, sample_pfx + ".call.cns")

    # Test for single-bin CNVs separately
    seg_bintest = bintest.do_bintest(cnarr, seg_call, target_only=True)
    tabio.write(seg_bintest, sample_pfx + ".bintest.cns")

    if plot_scatter:
        scatter.do_scatter(cnarr, seg_final)
        pyplot.savefig(sample_pfx + "-scatter.png", format="png", bbox_inches="tight")
        logging.info("Wrote %s-scatter.png", sample_pfx)

    if plot_diagram:
        is_xx = cnarr.guess_xx(is_haploid_x_reference, diploid_parx_genome)
        outfname = sample_pfx + "-diagram.pdf"
        diagram.create_diagram(
            cnarr.shift_xx(is_haploid_x_reference, is_xx, diploid_parx_genome),
            seg_final.shift_xx(is_haploid_x_reference, is_xx, diploid_parx_genome),
            0.5,
            3,
            outfname,
        )
        logging.info("Wrote %s", outfname)
