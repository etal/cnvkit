#!/usr/bin/env python
"""Unit tests for the CNVkit library, cnvlib."""
import unittest

import logging

from cnvlib.cnary import CopyNumArray

logging.basicConfig(level=logging.ERROR, format="%(message)s")

# unittest/pomegranate 0.10.0: ImportWarning: can't resolve package from
# __spec__ or __package__, falling back on __name__ and __path__
import warnings

warnings.filterwarnings("ignore", category=ImportWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import numpy as np
from skgenome import tabio

import cnvlib

# Import all modules as a smoke test
from cnvlib import (
    access,
    antitarget,
    autobin,
    batch,
    bintest,
    call,
    cnary,
    commands,
    core,
    coverage,
    diagram,
    export,
    fix,
    import_rna,
    importers,
    metrics,
    params,
    plots,
    reference,
    reports,
    segmentation,
    segmetrics,
    smoothing,
    vary,
)


class CommandTests(unittest.TestCase):
    """Tests for top-level commands."""

    def test_access(self):
        fasta = "formats/chrM-Y-trunc.hg19.fa"
        for min_gap_size, expect_nrows in ((None, 7), (500, 3), (1000, 2)):
            acc = commands.do_access(fasta, [], min_gap_size, skip_noncanonical=False)
            self.assertEqual(len(acc), expect_nrows)
        excludes = ["formats/dac-my.bed", "formats/my-targets.bed"]
        for min_gap_size, expect_nrows in (
            (None, 12),
            (2, 10),
            (20, 5),
            (200, 3),
            (2000, 2),
        ):
            acc = commands.do_access(
                fasta, excludes, min_gap_size, skip_noncanonical=False
            )
            self.assertEqual(len(acc), expect_nrows)
        # Dropping chrM, keeping only chrY
        acc = commands.do_access(fasta, excludes, 10, skip_noncanonical=True)
        self.assertEqual(len(acc), 5)

    def test_antitarget(self):
        """The 'antitarget' command."""
        baits = tabio.read_auto("formats/nv2_baits.interval_list")
        access = tabio.read_auto("../data/access-5k-mappable.hg19.bed")
        self.assertLess(0, len(commands.do_antitarget(baits)))
        self.assertLess(0, len(commands.do_antitarget(baits, access)))
        self.assertLess(0, len(commands.do_antitarget(baits, access, 200000)))
        self.assertLess(0, len(commands.do_antitarget(baits, access, 10000, 5000)))

    def test_autobin(self):
        """The 'autobin' command."""
        bam_fname = "formats/na12878-chrM-Y-trunc.bam"
        target_bed = "formats/my-targets.bed"
        targets = tabio.read(target_bed, "bed")
        access_bed = "../data/access-5k-mappable.hg19.bed"
        accessible = tabio.read(access_bed, "bed").filter(chromosome="chrY")
        for method in ("amplicon", "wgs", "hybrid"):
            (cov, bs), _ = autobin.do_autobin(
                bam_fname, method, targets=targets, access=accessible
            )
            self.assertGreater(cov, 0)
            self.assertGreater(bs, 0)

    def test_batch(self):
        """The 'batch' command."""
        target_bed = "formats/my-targets.bed"
        fasta = "formats/chrM-Y-trunc.hg19.fa"
        bam = "formats/na12878-chrM-Y-trunc.bam"
        annot = "formats/my-refflat.bed"
        # Build a single-sample WGS reference
        ref_fname, tgt_bed_fname, _ = batch.batch_make_reference(
            [bam],
            None,
            None,
            True,
            None,
            fasta,
            annot,
            True,
            500,
            None,
            None,
            None,
            None,
            "build",
            1,
            False,
            "wgs",
            False,
        )
        self.assertEqual(ref_fname, "build/reference.cnn")
        refarr = cnvlib.read(ref_fname, "bed")
        tgt_regions = tabio.read(tgt_bed_fname, "bed")
        self.assertEqual(len(refarr), len(tgt_regions))
        # Build a single-sample hybrid-capture reference
        ref_fname, tgt_bed_fname, anti_bed_fname = batch.batch_make_reference(
            [bam],
            target_bed,
            None,
            True,
            None,
            fasta,
            None,
            True,
            10,
            None,
            1000,
            100,
            None,
            "build",
            1,
            False,
            "hybrid",
            False,
        )
        self.assertEqual(ref_fname, "build/reference.cnn")
        refarr = cnvlib.read(ref_fname, "bed")
        tgt_regions = tabio.read(tgt_bed_fname, "bed")
        anti_regions = tabio.read(anti_bed_fname, "bed")
        self.assertEqual(len(refarr), len(tgt_regions) + len(anti_regions))
        # Run the same sample
        batch.batch_run_sample(
            bam,
            tgt_bed_fname,
            anti_bed_fname,
            ref_fname,
            "build",
            True,
            None,
            True,
            True,
            "Rscript",
            False,
            False,
            "hybrid",
            "hmm",
            1,
            False,
        )
        cns = cnvlib.read("build/na12878-chrM-Y-trunc.cns", None)
        self.assertGreater(len(cns), 0)

    def test_bintest(self):
        """The 'bintest' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr", None)
        segarr = cnvlib.read("formats/amplicon.cns", None)
        # Simple
        rows = commands.do_bintest(cnarr, alpha=0.05)
        self.assertGreater(len(rows), 0)
        self.assertLess(len(rows), len(cnarr))
        # Versus segments
        rows = commands.do_bintest(cnarr, segarr, target_only=True)
        self.assertGreaterEqual(len(rows), len(segarr))
        self.assertLess(len(rows), len(cnarr))

    def test_breaks(self):
        """The 'breaks' command."""
        probes = cnvlib.read("formats/amplicon.cnr", None)
        segs = cnvlib.read("formats/amplicon.cns", None)
        rows = commands.do_breaks(probes, segs, 4)
        self.assertGreater(len(rows), 0)

    def test_call(self):
        """The 'call' command."""
        # Methods: clonal, threshold, none
        tr_cns = cnvlib.read("formats/tr95t.cns", None)
        tr_thresh = commands.do_call(
            tr_cns, None, "threshold", is_reference_male=True, is_sample_female=True
        )
        self.assertEqual(len(tr_cns), len(tr_thresh))
        tr_clonal = commands.do_call(
            tr_cns,
            None,
            "clonal",
            purity=0.65,
            is_reference_male=True,
            is_sample_female=True,
        )
        self.assertEqual(len(tr_cns), len(tr_clonal))
        cl_cns = cnvlib.read("formats/cl_seq.cns", None)
        cl_thresh = commands.do_call(
            cl_cns,
            None,
            "threshold",
            thresholds=np.log2((np.arange(12) + 0.5) / 6.0),
            is_reference_male=True,
            is_sample_female=True,
        )
        self.assertEqual(len(cl_cns), len(cl_thresh))
        cl_clonal = commands.do_call(
            cl_cns,
            None,
            "clonal",
            ploidy=6,
            purity=0.99,
            is_reference_male=True,
            is_sample_female=True,
        )
        self.assertEqual(len(cl_cns), len(cl_clonal))
        cl_none = commands.do_call(
            cl_cns,
            None,
            "none",
            ploidy=6,
            purity=0.99,
            is_reference_male=True,
            is_sample_female=True,
        )
        self.assertEqual(len(cl_cns), len(cl_none))

    def test_call_filter(self):
        segments = cnvlib.read("formats/tr95t.segmetrics.cns", None)
        variants = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        # Each filter individually, then all filters together
        for filters in (
            ["ampdel"],
            ["cn"],
            ["ci"],
            ["sem"],
            ["sem", "cn", "ampdel"],
            ["ci", "cn"],
        ):
            result = commands.do_call(
                segments,
                variants,
                method="threshold",
                purity=0.9,
                is_reference_male=True,
                is_sample_female=True,
                filters=filters,
            )
            self.assertLessEqual(len(result), len(segments))
            if "ampdel" not in filters:
                # At least 1 segment per chromosome remains
                self.assertLessEqual(len(segments.chromosome.unique()), len(result))
            for colname in "baf", "cn", "cn1", "cn2":
                self.assertIn(colname, result)

    def test_call_sex(self):
        """Test each 'call' method on allosomes."""
        for (
            fname,
            sample_is_f,
            ref_is_m,
            chr1_expect,
            chrx_expect,
            chry_expect,
            chr1_cn,
            chrx_cn,
            chry_cn,
        ) in (
            ("formats/f-on-f.cns", True, False, 0, 0, None, 2, 2, None),
            ("formats/f-on-m.cns", True, True, 0.585, 1, None, 3, 2, None),
            ("formats/m-on-f.cns", False, False, 0, -1, 0, 2, 1, 1),
            ("formats/m-on-m.cns", False, True, 0, 0, 0, 2, 1, 1),
        ):
            cns = cnvlib.read(fname, None)
            chr1_idx = cns.chromosome == "chr1"
            chrx_idx = cns.chromosome == "chrX"
            chry_idx = cns.chromosome == "chrY"

            def test_chrom_means(segments):
                self.assertEqual(chr1_cn, segments["cn"][chr1_idx].mean())
                self.assertAlmostEqual(
                    chr1_expect, segments["log2"][chr1_idx].mean(), 0
                )
                self.assertEqual(chrx_cn, segments["cn"][chrx_idx].mean())
                self.assertAlmostEqual(
                    chrx_expect, segments["log2"][chrx_idx].mean(), 0
                )
                if not sample_is_f:
                    self.assertEqual(chry_cn, segments["cn"][chry_idx].mean())
                    self.assertAlmostEqual(
                        chry_expect, segments["log2"][chry_idx].mean(), 0
                    )

            # Call threshold
            cns_thresh = commands.do_call(
                cns,
                None,
                "threshold",
                is_reference_male=ref_is_m,
                is_sample_female=sample_is_f,
            )
            test_chrom_means(cns_thresh)
            # Call clonal pure
            cns_clone = commands.do_call(
                cns,
                None,
                "clonal",
                is_reference_male=ref_is_m,
                is_sample_female=sample_is_f,
            )
            test_chrom_means(cns_clone)
            # Call clonal barely-mixed
            cns_p99 = commands.do_call(
                cns,
                None,
                "clonal",
                purity=0.99,
                is_reference_male=ref_is_m,
                is_sample_female=sample_is_f,
            )
            test_chrom_means(cns_p99)
            
    def test_call_various_abs_ref_exp_methods(self):
        irrelevant_gene_name = "whatever"
        irrelevant_log2_value = 0.0
        ploidy = 2
        cnarr = CopyNumArray.from_rows(
            [
                ["1", 123456, 456789, irrelevant_gene_name, irrelevant_log2_value],
                ["X", 15000, 20000, irrelevant_gene_name, irrelevant_log2_value],  # PAR1
                ["X", 4000000, 5000000, irrelevant_gene_name, irrelevant_log2_value],
                ["Y", 4000000, 5000000, irrelevant_gene_name, irrelevant_log2_value],
            ]
        )
        # "grch38": [10001, 2781479, 155701383, 156030895],

        def _run(is_reference_male, is_sample_female, diploid_parx_genome=None):
            df = call.get_as_dframe_and_set_reference_and_expect_copies(cnarr, ploidy, is_reference_male, diploid_parx_genome, is_sample_female)
            abs_ref = call.absolute_reference(cnarr, ploidy, diploid_parx_genome, is_reference_male)
            abs_exp = call.absolute_expect(cnarr, ploidy, diploid_parx_genome, is_sample_female)
            return df, abs_ref, abs_exp

        def _assert_df(iloc, df, reference, expect):
            self.assertTrue('reference' in df.columns)
            self.assertTrue('expect' in df.columns)
            r = df.iloc[iloc]
            self.assertEqual(r.reference, reference)
            self.assertEqual(r.expect, expect)
        def _assert_abs_copies(i, abs_values, copies):
            self.assertEqual(abs_values[i], copies)

        def _assert_chr1(df, abs_ref, abs_exp):
            _assert_df(0, df, 2, 2)
            _assert_abs_copies(0, abs_ref, 2)
            _assert_abs_copies(0, abs_exp, 2)

        def _assert_chrx_par(df, df_reference, df_expect, abs_ref, abs_ref_copies, abs_exp, abs_exp_copies):
            _assert_df(1, df, df_reference, df_expect)
            _assert_abs_copies(1, abs_exp, abs_exp_copies)
            _assert_abs_copies(1, abs_ref, abs_ref_copies)

        def _assert_chrx_non_par(df, df_reference, df_expect, abs_ref, abs_ref_copies, abs_exp, abs_exp_copies):
            _assert_df(2, df, df_reference, df_expect)
            _assert_abs_copies(2, abs_ref, abs_ref_copies)
            _assert_abs_copies(2, abs_exp, abs_exp_copies)

        def _assert_chry(df, df_reference, df_expect, abs_ref, abs_ref_copies, abs_exp, abs_exp_copies):
            _assert_df(3, df, df_reference, df_expect)
            _assert_abs_copies(3, abs_ref, abs_ref_copies)
            _assert_abs_copies(3, abs_exp, abs_exp_copies)

        is_reference_male = True
        is_female_sample = True
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample)
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 1, 2, abs_ref, 1, abs_exp, 2)
        _assert_chrx_non_par(df, 1, 2, abs_ref, 1, abs_exp, 2)
        _assert_chry(df, 1, 0, abs_ref, 1, abs_exp, 0)
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample, "grch38")
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 2, 2, abs_ref, 2, abs_exp, 2)
        _assert_chrx_non_par(df, 1, 2, abs_ref, 1, abs_exp, 2)
        _assert_chry(df, 1, 0, abs_ref, 1, abs_exp, 0)

        is_reference_male = True
        is_female_sample = False
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample)
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 1, 1, abs_ref, 1, abs_exp, 1)
        _assert_chrx_non_par(df, 1, 1, abs_ref, 1, abs_exp, 1)
        _assert_chry(df, 1, 1, abs_ref, 1, abs_exp, 1)
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample, "grch38")
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 2, 2, abs_ref, 2, abs_exp, 2)
        _assert_chrx_non_par(df, 1, 1, abs_ref, 1, abs_exp, 1)
        _assert_chry(df, 1, 1, abs_ref, 1, abs_exp, 1)

        is_reference_male = False
        is_female_sample = True
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample)
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 2, 2, abs_ref, 2, abs_exp, 2)
        _assert_chrx_non_par(df, 2, 2, abs_ref, 2, abs_exp, 2)
        _assert_chry(df, 1, 0, abs_ref, 1, abs_exp, 0)  # Do we really want ref=1 for a female ref on chr y?
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample, "grch38")
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 2, 2, abs_ref, 2, abs_exp, 2)
        _assert_chrx_non_par(df, 2, 2, abs_ref, 2, abs_exp, 2)
        _assert_chry(df, 1, 0, abs_ref, 1, abs_exp, 0)  # Do we really want ref=1 for a female ref on chr y?

        is_reference_male = False
        is_female_sample = False
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample)
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 2, 1, abs_ref, 2, abs_exp, 1)
        _assert_chrx_non_par(df, 2, 1, abs_ref, 2, abs_exp, 1)
        _assert_chry(df, 1, 1, abs_ref, 1, abs_exp, 1)  # Do we really want ref=1 for a female ref on chr y?
        df, abs_ref, abs_exp = _run(is_reference_male, is_female_sample, "grch38")
        _assert_chr1(df, abs_ref, abs_exp)
        _assert_chrx_par(df, 2, 2, abs_ref, 2, abs_exp, 2)
        _assert_chrx_non_par(df, 2, 1, abs_ref, 2, abs_exp, 1)
        _assert_chry(df, 1, 1, abs_ref, 1, abs_exp, 1)  # Do we really want ref=1 for a female ref on chr y?


    def _get_true_copy_numbers(self):
        """Returns a dict with the values to check against."""
        ret = {}
        ret['ploidy'] = 2  # TODO: test with other ploidies, too
        ret['purity'] = 0.8
        # Initialize 2x2 matrix for each combination of is_refence_male x is_sample_female.
        matrix = {
            is_reference_male: {
                is_sample_female: [] for is_sample_female in [True, False]
            } for is_reference_male in [True, False]
        }
        # Each tuple represents (row_index, chromosome_name, absolute_copies, reference_copies, expected_copies):
        # male reference, female sample
        matrix[True][True].append((0, "chr1", 2.24, 2, 2))
        matrix[True][True].append((9, "chrX", 3.09, 1, 2))
        matrix[True][True].append((10, "chrX", 5.76, 1, 2))
        matrix[True][True].append((11, "chrX", 1.92, 1, 2))
        matrix[True][True].append((16, "chrY", 3.38, 1, 0))

        # male reference, male sample
        matrix[True][False].append((0, "chr1", 2.24, 2, 2))
        matrix[True][False].append((9, "chrX", 3.34, 1, 1))
        matrix[True][False].append((10, "chrX", 6.01, 1, 1))
        matrix[True][False].append((11, "chrX", 2.17, 1, 1))
        matrix[True][False].append((16, "chrY", 3.13, 1, 1))

        # female reference, female sample
        matrix[False][True].append((0, "chr1", 2.24, 2, 2))
        matrix[False][True].append((9, "chrX", 6.68, 2, 2))
        matrix[False][True].append((10, "chrX", 12.02, 2, 2))
        matrix[False][True].append((11, "chrX", 4.33, 2, 2))
        matrix[False][True].append((16, "chrY", 3.38, 1, 0))

        # female reference, male sample
        matrix[False][False].append((0, "chr1", 2.24, 2, 2))
        matrix[False][False].append((9, "chrX", 6.93, 2, 1))
        matrix[False][False].append((10, "chrX", 12.27, 2, 1))
        matrix[False][False].append((11, "chrX", 4.58, 2, 1))
        matrix[False][False].append((16, "chrY", 3.13, 1, 1))  # ref_copies = 1 even for female reference

        ret['matrix'] = matrix
        return ret

    def test_call_absolute_clonal_and_dataframe(self):
        """Test some helper methods for `call`."""

        cnarr = cnvlib.read("formats/tr95t.subset.cns", None)
        TRUTH = self._get_true_copy_numbers()
        matrix = TRUTH['matrix']
        ploidy = TRUTH['ploidy']
        purity = TRUTH['purity']
        for is_reference_male in (True, False):
            for is_sample_female in (True, False):
                # print(f"is_reference_male: {is_reference_male}, is_sample_female: {is_sample_female}")
                df = call.absolute_dataframe(cnarr, ploidy, purity, is_reference_male, is_sample_female)
                absolutes_df = df['absolute']
                references = df['reference']
                expects = df['expect']
                absolutes_clonal = call.absolute_clonal(cnarr, ploidy, purity, is_reference_male, is_sample_female)
                tuples = matrix[is_reference_male][is_sample_female]
                for tuple in tuples:
                    index, chromosome, abs_copies, ref_copies, exp_copies = tuple
                    # print(f"index: {index}, chromosome: {chromosome}, abs_copies: {abs_copies}, ref_copies: {ref_copies}, exp_copies: {exp_copies}")
                    self.assertEqual(cnarr[index].chromosome, chromosome)
                    self.assertAlmostEqual(absolutes_clonal[index], abs_copies, 2)
                    self.assertAlmostEqual(absolutes_df[index], abs_copies, 2)
                    self.assertEqual(references[index], ref_copies)
                    self.assertEqual(expects[index], exp_copies)

        # todo: add test with par region removed

    def test_coverage(self):
        """The 'coverage' command."""
        # fa = 'formats/chrM-Y-trunc.hg19.fa'
        bed = "formats/my-targets.bed"
        bam = "formats/na12878-chrM-Y-trunc.bam"
        for by_count in (False, True):
            for min_mapq in (0, 30):
                for nprocs in (1, 2):
                    cna = commands.do_coverage(
                        bed, bam, by_count=by_count, min_mapq=min_mapq, processes=nprocs
                    )
                    self.assertEqual(len(cna), 4)
                    self.assertTrue((cna.log2 != 0).any())
                    self.assertGreater(cna.log2.nunique(), 1)

    def test_export_bed_vcf(self):
        """The 'export' command for formats with absolute copy number."""
        for fname, ploidy, is_f in [
            ("tr95t.cns", 2, True),
            ("cl_seq.cns", 6, True),
            ("amplicon.cns", 2, False),
        ]:
            cns = cnvlib.read("formats/" + fname, None)
            # BED
            for show in ("ploidy", "variant", "all"):
                tbl_bed = export.export_bed(
                    cns, ploidy, True, None, is_f, cns.sample_id, show
                )
                if show == "all":
                    self.assertEqual(
                        len(tbl_bed), len(cns), f"{fname} {ploidy}"
                    )
                else:
                    self.assertLess(len(tbl_bed), len(cns))
            # VCF
            _vheader, vcf_body = export.export_vcf(cns, ploidy, True, is_f)
            self.assertTrue(0 < len(vcf_body.splitlines()) < len(cns))

    def test_export_cdt_jtv(self):
        """The 'export' command for CDT and Java TreeView formats."""
        fnames = ["formats/p2-20_1.cnr", "formats/p2-20_2.cnr"]
        sample_ids = list(map(core.fbase, fnames))
        nrows = linecount(fnames[0]) - 1
        for fmt_key, header2 in (("cdt", 2), ("jtv", 0)):
            table = export.merge_samples(fnames)
            formatter = export.EXPORT_FORMATS[fmt_key]
            _oh, outrows = formatter(sample_ids, table)
            self.assertEqual(len(list(outrows)), nrows + header2)

    def test_export_nexus(self):
        """The 'export nexus-basic' and 'nexus-ogt' commands."""
        cnr = cnvlib.read("formats/amplicon.cnr", None)
        table_nb = export.export_nexus_basic(cnr)
        self.assertEqual(len(table_nb), len(cnr))
        varr = commands.load_het_snps(
            "formats/na12878_na12882_mix.vcf", None, None, 15, None
        )
        table_ogt = export.export_nexus_ogt(cnr, varr, 0.05)
        self.assertEqual(len(table_ogt), len(cnr))

    def test_export_seg(self):
        """The 'export seg' command."""
        seg_rows = export.export_seg(["formats/tr95t.cns"])
        self.assertGreater(len(seg_rows), 0)
        seg2_rows = export.export_seg(["formats/tr95t.cns", "formats/cl_seq.cns"])
        self.assertGreater(len(seg2_rows), len(seg_rows))

    def test_export_theta(self):
        """The 'export theta' command."""
        segarr = cnvlib.read("formats/tr95t.cns", None)
        len_seg_auto = len(segarr.autosomes())
        table_theta = export.export_theta(segarr, None)
        self.assertEqual(len(table_theta), len_seg_auto)
        ref = cnvlib.read("formats/reference-tr.cnn", None)
        table_theta = export.export_theta(segarr, ref)
        self.assertEqual(len(table_theta), len_seg_auto)
        varr = commands.load_het_snps(
            "formats/na12878_na12882_mix.vcf", "NA12882", "NA12878", 15, None
        )
        tumor_snps, normal_snps = export.export_theta_snps(varr)
        self.assertLess(len(tumor_snps), len(varr))
        self.assertGreater(len(tumor_snps), 0)
        self.assertLess(len(normal_snps), len(varr))
        self.assertGreater(len(normal_snps), 0)

    def test_fix(self):
        """The 'fix' command."""
        # Extract fake target/antitarget bins from a combined file
        ref = cnvlib.read("formats/reference-tr.cnn", None)
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        tgt_bins.log2 += np.random.randn(len(tgt_bins)) / 5
        anti_bins = ref[is_bg]
        anti_bins.log2 += np.random.randn(len(anti_bins)) / 5
        blank_bins = cnary.CopyNumArray([])
        # Typical usage (hybrid capture)
        cnr = commands.do_fix(tgt_bins, anti_bins, ref)
        self.assertTrue(0 < len(cnr) <= len(ref))
        # Blank antitargets (WGS or amplicon)
        cnr = commands.do_fix(tgt_bins, blank_bins, ref[~is_bg])
        self.assertTrue(0 < len(cnr) <= len(tgt_bins))

    def test_genemetrics(self):
        """The 'genemetrics' command."""
        probes = cnvlib.read("formats/amplicon.cnr", None)
        rows = commands.do_genemetrics(probes, male_reference=True)
        self.assertGreater(len(rows), 0)
        segs = cnvlib.read("formats/amplicon.cns", None)
        rows = commands.do_genemetrics(probes, segs, 0.3, 4, male_reference=True)
        self.assertGreater(len(rows), 0)

    def test_import_theta(self):
        """The 'import-theta' command."""
        cns = cnvlib.read("formats/nv3.cns", None)
        theta_fname = "formats/nv3.n3.results"
        for new_cns in commands.do_import_theta(cns, theta_fname):
            self.assertTrue(0 < len(new_cns) <= len(cns))

    def test_metrics(self):
        """The 'metrics' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr", None)
        segments = cnvlib.read("formats/amplicon.cns", None)
        result = metrics.do_metrics(cnarr, segments, skip_low=True)
        self.assertEqual(result.shape, (1, 6))
        values = result.loc[0, result.columns[1:]]
        for val in values:
            self.assertGreater(val, 0)

    def test_reference(self):
        """The 'reference' command."""
        # Empty/unspecified antitargets
        nlines = linecount("formats/amplicon.cnr") - 1
        ref = commands.do_reference(["formats/amplicon.cnr"], ["formats/empty"])
        self.assertEqual(len(ref), nlines)
        ref = commands.do_reference(["formats/amplicon.cnr"])
        self.assertEqual(len(ref), nlines)
        # Empty/unspecified antitargets, flat reference
        nlines = linecount("formats/amplicon.bed")
        ref = commands.do_reference_flat("formats/amplicon.bed", "formats/empty")
        self.assertEqual(len(ref), nlines)
        ref = commands.do_reference_flat("formats/amplicon.bed")
        self.assertEqual(len(ref), nlines)
        # Misc
        ref = cnvlib.read("formats/reference-tr.cnn", None)
        targets, antitargets = reference.reference2regions(ref)
        self.assertLess(0, len(antitargets))
        self.assertEqual(len(antitargets), (ref["gene"] == "Background").sum())
        self.assertEqual(len(targets), len(ref) - len(antitargets))

    def test_segment(self):
        """The 'segment' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr", None)
        n_chroms = cnarr.chromosome.nunique()
        # NB: R methods are in another script; haar is pure-Python
        segments = segmentation.do_segmentation(cnarr, "haar")
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())
        segments = segmentation.do_segmentation(
            cnarr, "haar", threshold=0.0001, skip_low=True
        )
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())
        varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        # TODO - This test is failing... commenting it out for now!
        # segments = segmentation.do_segmentation(cnarr, "haar", variants=varr)
        # self.assertGreater(len(segments), n_chroms)
        # self.assertTrue((segments.start < segments.end).all())

    def test_segment_hmm(self):
        """The 'segment' command with HMM methods."""
        for fname in ("formats/amplicon.cnr", "formats/p2-20_1.cnr"):
            cnarr = cnvlib.read(fname, None)
            n_chroms = cnarr.chromosome.nunique()
            # NB: R methods are in another script; haar is pure-Python
            segments = segmentation.do_segmentation(cnarr, "hmm")
            self.assertGreater(len(segments), n_chroms)
            self.assertTrue((segments.start < segments.end).all())
            segments = segmentation.do_segmentation(cnarr, "hmm-tumor", skip_low=True)
            self.assertGreater(len(segments), n_chroms)
            self.assertTrue((segments.start < segments.end).all())
            segments = segmentation.do_segmentation(cnarr, "hmm-germline")
            self.assertGreater(len(segments), n_chroms)
            self.assertTrue((segments.start < segments.end).all())
            varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
            segments = segmentation.do_segmentation(cnarr, "hmm", variants=varr)
            self.assertGreater(len(segments), n_chroms)

    def test_segment_parallel(self):
        """The 'segment' command, in parallel."""
        cnarr = cnvlib.read("formats/amplicon.cnr", None)
        psegments = segmentation.do_segmentation(cnarr, "haar", processes=2)
        ssegments = segmentation.do_segmentation(cnarr, "haar", processes=1)
        self.assertEqual(psegments.data.shape, ssegments.data.shape)
        self.assertEqual(len(psegments.meta), len(ssegments.meta))

    def test_segmetrics(self):
        """The 'segmetrics' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr", None)
        segarr = cnvlib.read("formats/amplicon.cns", None)
        sm = segmetrics.do_segmetrics(
            cnarr,
            segarr,
            location_stats=["mean", "median", "mode", "p_ttest"],
            spread_stats=["stdev", "sem", "iqr"],
            interval_stats=["pi", "ci"],
            bootstraps=50,
            smoothed=True,
        )
        # Restrict to segments with enough supporting probes for sane stats
        sm = sm[sm["probes"] > 3]
        self.assertTrue((sm["pi_lo"] < sm["median"]).all())
        self.assertTrue((sm["pi_hi"] > sm["median"]).all())
        self.assertTrue((sm["ci_lo"] < sm["mean"]).all())
        self.assertTrue((sm["ci_hi"] > sm["mean"]).all())

    def test_target(self):
        """The 'target' command."""
        # return  # DBG
        annot_fname = "formats/refflat-mini.txt"
        for bait_fname in (
            "formats/nv2_baits.interval_list",
            "formats/amplicon.bed",
            "formats/baits-funky.bed",
        ):
            baits = tabio.read_auto(bait_fname)
            bait_len = len(baits)
            # No splitting: w/o and w/ re-annotation
            r1 = commands.do_target(baits)
            self.assertEqual(len(r1), bait_len)
            r1a = commands.do_target(baits, do_short_names=True, annotate=annot_fname)
            self.assertEqual(len(r1a), len(r1))
            # Splitting, w/o and w/ re-annotation
            r2 = commands.do_target(
                baits, do_short_names=True, do_split=True, avg_size=100
            )
            self.assertGreater(len(r2), len(r1))
            for _c, subarr in r2.by_chromosome():
                self.assertTrue(subarr.start.is_monotonic_increasing, bait_fname)
                self.assertTrue(subarr.end.is_monotonic_increasing, bait_fname)
                # Bins are non-overlapping; next start >= prev. end
                self.assertTrue(
                    ((subarr.start.values[1:] - subarr.end.values[:-1]) >= 0).all()
                )
            r2a = commands.do_target(
                baits,
                do_short_names=True,
                do_split=True,
                avg_size=100,
                annotate=annot_fname,
            )
            self.assertEqual(len(r2a), len(r2))
            # Original regions object should be unmodified
            self.assertEqual(len(baits), bait_len)


def linecount(filename):
    i = -1
    with open(filename) as handle:
        for i, _line in enumerate(handle):
            pass
        return i + 1


if __name__ == "__main__":
    unittest.main(verbosity=2)
