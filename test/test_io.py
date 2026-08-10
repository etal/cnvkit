#!/usr/bin/env python
"""Unit tests for the CNVkit library, cnvlib."""

import math
import tempfile
import unittest
from io import StringIO

from conftest import linecount

from cnvlib.cnary import CopyNumArray as CNA
from skgenome import tabio


class IOTests(unittest.TestCase):
    """Tests for I/O modules."""

    def test_empty(self):
        """Instantiate from an empty file."""
        for fmt in ("auto", "bed", "interval", "tab", "text"):  # "bed3", "bed4",
            regions = tabio.read("formats/empty", fmt=fmt)
            self.assertEqual(len(regions), 0)

    def test_read_auto(self):
        for fname, nrows in (
            ("formats/empty", 0),
            ("formats/agilent.bed", 11),
            ("formats/amplicon.bed", 1433),
            ("formats/amplicon.text", 1433),
            ("formats/nv2_baits.interval_list", 6809),
            ("formats/refflat-mini.txt", 100),
            ("formats/example.gff", 6),
        ):
            self.assertEqual(len(tabio.read_auto(fname)), nrows)
            with open(fname) as handle:
                self.assertEqual(len(tabio.read_auto(handle)), nrows)

    def test_read_bed(self):
        """Read the BED format."""
        fname = "formats/amplicon.bed"
        regions = tabio.read(fname, "bed")
        self.assertEqual(len(regions), linecount(fname))
        self.assertEqual(regions.sample_id, "amplicon")

    def test_read_gff(self):
        """Read the GFF format."""
        for fname, nrows, sample_id in (
            ("formats/example.gff", 6, "example"),
            ("formats/GRCh37_BRAF.gff.gz", 49, "GRCh37_BRAF"),
        ):
            regions = tabio.read(fname, "gff")
            self.assertEqual(len(regions), nrows)
            self.assertEqual(regions.sample_id, sample_id)

    def test_read_ilist(self):
        """Read the interval list format."""
        regions = tabio.read("formats/nv2_baits.interval_list", "interval")
        self.assertEqual(len(regions), 6809)
        self.assertEqual(regions.sample_id, "nv2_baits")

    def test_read_picardhs(self):
        """Read Picard CalculateHsMetrics PER_TARGET_COVERAGE format."""
        fname = "picard/p2-5_5.antitargetcoverage.csv"
        cna = tabio.read(fname, "picardhs")
        self.assertEqual(len(cna), linecount(fname) - 1)
        self.assertEqual(cna.sample_id, "p2-5_5")

    def test_read_refflat(self):
        """Read the UCSC 'refFlat' format."""
        fname = "formats/refflat-mini.txt"
        regions = tabio.read(fname, "refflat")
        self.assertEqual(len(regions), linecount(fname))
        self.assertEqual(13, regions.chromosome.nunique())

    def test_read_seg(self):
        """Read the SEG format."""
        for fname, header_len, args in (
            # Convert integer chrom. IDs to hg19 names
            (
                "formats/cw-tr-log2.seg",
                1,
                ({"23": "X", "24": "Y", "25": "M"}, "chr", False),
            ),
            # Convert segmented array CGH data in log10 scale to log2
            ("formats/acgh-log10.seg", 1, (None, None, True)),
            # From DNAcopy, with a stray warning message from R
            ("formats/warning.seg", 2, (None, None, False)),
        ):
            expect_lines = linecount(fname) - header_len
            seen_lines = 0
            for _sample_id, dframe in tabio.seg.parse_seg(fname, *args):
                seen_lines += len(dframe)
            self.assertEqual(seen_lines, expect_lines)

    def test_read_seg_empty(self):
        """A header-only SEG (no segment rows) reads as an empty table (#868)."""
        header = '"ID"\t"chrom"\t"loc.start"\t"loc.end"\t"num.mark"\t"seg.mean"\n'
        samples = list(tabio.seg.parse_seg(StringIO(header)))
        self.assertEqual(len(samples), 1)
        _sample_id, dframe = samples[0]
        self.assertEqual(len(dframe), 0)
        self.assertIn("chromosome", dframe.columns)
        # The high-level reader returns an empty table rather than raising
        self.assertEqual(len(tabio.read(StringIO(header), "seg")), 0)

    def test_read_text(self):
        """Read the text region format."""
        fname = "formats/amplicon.text"
        regions = tabio.read(fname, "text")
        self.assertEqual(len(regions), linecount(fname))
        self.assertEqual(regions.sample_id, "amplicon")

    def test_read_vcf(self):
        """Read the VCF format."""
        # Paired VCF with full info
        fname = "formats/na12878_na12882_mix.vcf"
        v1 = tabio.read(fname, "vcf")
        self.assertLess(len(v1), linecount(fname))
        self.assertLess(0, len(v1))
        for sid in ("NA12882", "NA12878"):
            v2 = tabio.read(fname, "vcf", sample_id=sid)
            self.assertEqual(v2.sample_id, sid)
            self.assertEqual(len(v1), len(v2))
        for kwarg in (
            {"min_depth": 100},
            {"skip_somatic": True},
            {"skip_reject": True},
        ):
            v3 = tabio.read(fname, "vcf", **kwarg)
            self.assertLess(len(v3), len(v1))
            self.assertLess(
                0,
                len(v3),
                f"{len(v3)} variants left after filter {next(iter(kwarg))!r}",
            )
        # VCF with no samples, 1 record which will be ignored
        v4 = tabio.read("formats/nosample.vcf", "vcf")
        self.assertEqual(len(v4), 0)
        self.assertEqual(v4.sample_id, "nosample")
        # VCF with 1 sample, no records
        v5 = tabio.read("formats/blank.vcf", "vcf", sample_id="Blank")
        self.assertEqual(len(v5), 0)
        self.assertEqual(v5.sample_id, "Blank")
        # VCF from GATK 4 with no ALT
        v6 = tabio.read("formats/gatk-emptyalt.vcf", "vcf", sample_id="sample1")
        self.assertEqual(len(v6), 0)

    def test_read_vcf_end_is_reference_footprint(self):
        """A record spans the reference bases it quotes, not its alt alleles.

        The end coordinate used to come from the alternate allele -- the
        replacement string's length in ``vcf``, the net indel size in
        ``vcf-simple`` -- so an insertion reached into reference it does not
        describe while a deletion stopped short of reference it does. At a
        bin's left edge that decided which variant supplied the b-allele
        frequency, so the three readers must agree here, not merely each be
        self-consistent.

        The blanket span assertion is sound for this fixture because every
        record in it is a precise variant, for which the spec defines END as
        ``POS + len(REF) - 1``; the 159 that declare one therefore agree with
        the footprint. A symbolic ALT would part the two, which is what
        ``test_read_vcf_end_tag_beats_the_footprint`` covers.
        """
        fname = "formats/na12878_na12882_mix.vcf"
        coords = {}
        for fmt in ("vcf", "vcf-simple", "vcf-sites"):
            dframe = tabio.read(fname, fmt).data
            spans = dframe.end - dframe.start
            n_wrong = int((spans != dframe.ref.str.len()).sum())
            self.assertEqual(
                n_wrong, 0, f"{fmt}: {n_wrong} rows do not span the reference allele"
            )
            coords[fmt] = dframe[["chromosome", "start", "end"]].to_csv(index=False)
        self.assertEqual(
            len(set(coords.values())), 1, "the VCF readers disagree on coordinates"
        )

    def test_read_vcf_end_tag_beats_the_footprint(self):
        """An explicit INFO/END wins over the reference allele's length.

        htslib folds a declared END into ``record.stop`` and drops it from the
        INFO mapping, so a reader that consults ``record.info`` loses it
        silently: this deletion was read as spanning five bases, the length of
        the string ``<DEL>``.
        """
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
            '##ALT=<ID=DEL,Description="Deletion">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t500\t.\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=900\tGT\t0/1\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            for fmt in ("vcf", "vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.start), [499], fmt)
                self.assertEqual(list(dframe.end), [900], fmt)

    def test_read_vcf_end_key_is_not_matched_by_substring(self):
        """Only a field keyed exactly END sets the end coordinate.

        The text readers found END by searching the raw INFO string for
        "END=", which also matches CIEND, the confidence interval around
        END that imprecise structural variants carry. htsjdk serialises
        INFO through a sorted map, so GATK and Picard write the keys
        alphabetically and CIEND lands ahead of the END it brackets; the
        leftmost match then read the interval "-50,50" as a coordinate and
        the file failed to load. A record carrying CIEND and no END fails
        the same way whatever the field order, since the embedded "END="
        is then the only one present.
        """
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
            '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="CI around END">\n'
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
            '##ALT=<ID=DEL,Description="Deletion">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t500\t.\tA\t<DEL>\t50\tPASS\tCIEND=-50,50;END=900;SVTYPE=DEL\tGT\t0/1\n"
            "chr1\t2000\t.\tACGT\tA\t50\tPASS\tCIEND=-50,50;SVTYPE=DEL\tGT\t0/1\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            for fmt in ("vcf", "vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.start), [499, 1999], fmt)
                self.assertEqual(list(dframe.end), [900, 2003], fmt)

    def test_read_vcf_end_without_a_value_uses_the_footprint(self):
        """An END that declares no value leaves the footprint in place.

        A record may spell "no end here" two ways: `.`, the marker VCF
        gives every missing value, and a bare `END` carrying no value at
        all, which is malformed for a reserved integer key but is what a
        writer emitting flag fields produces by accident. Neither states
        a coordinate, so both fall back to the reference allele rather
        than failing the read. The same rule governs QUAL a few lines
        below, where `.` becomes NaN instead of raising.
        """
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t2000\t.\tACGT\tA\t50\tPASS\tSVTYPE=DEL;END=.\tGT\t0/1\n"
            "chr1\t3000\t.\tACGT\tA\t50\tPASS\tSVTYPE=DEL;END\tGT\t0/1\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            for fmt in ("vcf", "vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.start), [1999, 2999], fmt)
                self.assertEqual(list(dframe.end), [2003, 3003], fmt)

    def test_read_vcf_strelka(self):
        """Read a Strelka somatic-SNV VCF, which has no GT FORMAT field.

        Strelka reports per-base allele counts via AU/CU/GU/TU instead of an
        explicit genotype. The reader must tolerate the missing GT (#943).
        """
        fname = "formats/strelka.vcf"
        # Tumor/normal pair
        pair = tabio.read(fname, "vcf", sample_id="TUMOR", normal_id="NORMAL")
        self.assertEqual(len(pair), 3)
        self.assertEqual(pair.sample_id, "TUMOR")
        # Alt counts come from the tier-1 count of the ALT base
        self.assertEqual(list(pair["alt_count"]), [16, 28, 5])
        self.assertEqual(list(pair["depth"]), [33, 30, 20])
        # Zygosity inferred from allele frequency (no GT available)
        self.assertEqual(list(pair["zygosity"]), [0.5, 1.0, 0.5])
        # Normal-sample columns are present and populated
        self.assertEqual(list(pair["n_alt_count"]), [10, 1, 0])
        self.assertEqual(list(pair["n_depth"]), [20, 25, 20])
        self.assertEqual(list(pair["n_zygosity"]), [0.5, 0.0, 0.0])
        # Single unpaired sample also works
        solo = tabio.read(fname, "vcf", sample_id="TUMOR")
        self.assertEqual(len(solo), 3)
        self.assertNotIn("n_depth", solo.data.columns)
        # Default selection (first sample) does not crash
        default = tabio.read(fname, "vcf")
        self.assertEqual(len(default), 3)

    def test_read_vcf_nocall(self):
        """A GT no-call ('./.') is inferred from frequency, not called hom-alt.

        pysam reports './.' as (None, None); the reader must treat it as a
        missing genotype and fall back to the allele frequency rather than
        misreading it as homozygous-alt.
        """
        v = tabio.read("formats/gt-nocall.vcf", "vcf")
        self.assertEqual(len(v), 3)
        self.assertEqual(list(v["alt_count"]), [10, 12, 1])
        self.assertEqual(list(v["depth"]), [20, 20, 20])
        # ./. with balanced reads -> het; explicit 0/1 -> het; ./. ref-heavy -> hom-ref
        self.assertEqual(list(v["zygosity"]), [0.5, 0.5, 0.0])

    def test_read_vcf_partial_gt(self):
        """A partial genotype ('1/.', '0/.') is read from the called allele.

        pysam reports a no-call allele as None, so '1/.' is (1, None) and '0/.'
        is (0, None). The missing-allele placeholder must not be counted as a
        distinct allele (gh-9vv): '0/.' has no ALT evidence and is hom-ref,
        while '1/.' carries an ALT allele and is kept heterozygous (consistent
        with the read evidence in real ensemble VCFs). Read depths here diverge
        from the GT classification (e.g. '1/.' with AF=0.1) to prove the result
        comes from the genotype, not from frequency inference.
        """
        v = tabio.read("formats/gt-partial.vcf", "vcf")
        self.assertEqual(len(v), 6)
        # 0/. -> hom-ref; 1/. -> het; ./. -> freq (balanced=het);
        # 0/0 -> hom-ref; 1/1 -> hom-alt; 0/1 -> het
        self.assertEqual(list(v["zygosity"]), [0.0, 0.5, 0.5, 0.0, 1.0, 0.5])

    def test_read_vcf_het_no_alt_count(self):
        """A het call with no allele-count field gets NaN alt_freq, not 0 (#407).

        When a VCF has genotypes (GT=0/1 -> het) but no AD/AO/DP4/AF to derive
        an alternate-allele count, the frequency is unknown. It must stay NaN so
        BAF aggregation (np.nanmedian) excludes it, rather than being coerced to
        0.0 -- which would keep the het site yet pin its mirrored BAF to exactly
        0 over the whole region.
        """
        varr = tabio.read(
            "formats/vcf-het-no-altcount.vcf",
            "vcf",
            sample_id="TUMOR",
            normal_id="NORMAL",
        )
        self.assertEqual(len(varr), 3)
        # Genotype is read as heterozygous, so these sites are kept...
        self.assertEqual(list(varr["zygosity"]), [0.5, 0.5, 0.5])
        # ...but with no allele counts the frequency is unknown, not zero
        self.assertTrue(all(math.isnan(f) for f in varr["alt_freq"]))
        self.assertTrue(all(math.isnan(f) for f in varr["n_alt_freq"]))
        # A bin spanning all three het SNPs gets NaN BAF, not a spurious 0.0
        seg = CNA.from_columns(
            {
                "chromosome": ["chr1"],
                "start": [0],
                "end": [1000],
                "gene": ["-"],
                "log2": [0.0],
            }
        )
        baf = varr.baf_by_ranges(seg)
        self.assertEqual(len(baf), 1)
        self.assertTrue(math.isnan(baf.iloc[0]))

    def test_read_vcf_malformed(self):
        """A VCF declaring a FORMAT column but no samples gives a clear error.

        pysam rejects such files; the reader must not mislabel the failure as
        an open-file-handle mistake (#680).
        """
        with self.assertRaises(ValueError) as ctx:
            tabio.read("formats/format-no-sample.vcf", "vcf")
        msg = str(ctx.exception)
        self.assertIn("format-no-sample.vcf", msg)
        # Should not be mislabeled as a "passed an open file handle" mistake
        self.assertNotIn("file handle", msg)


if __name__ == "__main__":
    unittest.main(verbosity=2)
