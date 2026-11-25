#!/usr/bin/env python
"""Unit tests for the CNVkit library, cnvlib."""

import argparse
import gzip
import math
import os
import shutil
import tempfile
import threading
import unittest
from io import StringIO

import pysam
from conftest import gzipped, linecount

from cnvlib import commands
from cnvlib.cnary import CopyNumArray as CNA
from skgenome import tabio
from skgenome.cli import skg_convert
from skgenome.tabio import bedio, gff, vcfio


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
            ("formats/refgene-mini.txt", 100),
            ("formats/genepredext-mini.txt", 100),
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

    def test_read_gff_ncbi_genes(self):
        """NCBI's RefSeq GFF3 parses to gene symbols and 0-based starts (#259)."""
        genes = tabio.read("formats/GRCh37_BRAF.gff.gz", "gff", keep_type="gene").data
        self.assertEqual(sorted(genes["gene"]), ["BANCR", "BRAF", "BRAFP1", "HMG20B"])
        braf = genes[genes["gene"] == "BRAF"]
        # Column 4 of the BRAF gene record is 140433812, 1-based
        self.assertEqual(
            list(braf[["chromosome", "start", "end"]].itertuples(index=False)),
            [("NC_000007.13", 140433811, 140624564)],
        )

    def test_read_gff_ncbi_rows_use_gene_symbol(self):
        """Every NCBI GFF3 row is labeled with its gene symbol, not an accession.

        The mRNA, lnc_RNA and CDS rows carry the transcript or protein accession
        in ``Name`` and the symbol in ``gene``.
        """
        regions = tabio.read("formats/GRCh37_BRAF.gff.gz", "gff").data
        self.assertEqual(set(regions["gene"]), {"BANCR", "BRAF", "BRAFP1", "HMG20B"})
        cds = regions[regions["type"] == "CDS"]
        self.assertEqual(len(cds), 18)
        self.assertEqual(set(cds["gene"]), {"BRAF"})

    def test_read_gff_tag_priority(self):
        """The first tag present in priority order wins, wherever it appears."""
        attributes = [
            # GENCODE GFF3: gene_id precedes gene_name
            "ID=g1;gene_id=ENSG00000157764.14;gene_name=BRAF;"
            "havana_gene=OTTHUMG00000157457.3",
            # Ensembl/GENCODE GTF
            'gene_id "ENSG00000189337"; gene_name "KAZN";',
            # GTF row without gene_name falls back to gene_id
            'gene_id "ENSG00000279493"; transcript_id "ENST00000623083";',
            # Ensembl GFF3 gene row: Name holds the symbol, gene_id the ID
            "ID=gene:ENSG00000157764;Name=BRAF;gene_id=ENSG00000157764",
            # NCBI GFF3 mRNA row: Name holds the accession, gene the symbol
            "ID=rna1;Name=NM_004333.4;gbkey=mRNA;gene=BRAF",
            # A key that merely ends in a tag name is not that tag
            "ID=x1;pseudogene=true;Name=X1",
            # An empty value falls through to the next tag
            'gene_id "E1"; transcript_id "T1"; gene_name "";',
            "ID=x2;gene=;Name=N1",
            # Extra spaces between GTF attributes
            'transcript_id "T2";  gene_id "E2";',
            # No recognized tag
            "ID=exon1;Parent=rna1",
        ]
        self.assertEqual(
            self._gff_genes(attributes),
            [
                "BRAF",
                "KAZN",
                "ENSG00000279493",
                "BRAF",
                "BRAF",
                "X1",
                "E1",
                "N1",
                "E2",
                "-",
            ],
        )

    def test_read_gff_explicit_tag(self):
        """An explicit tag names one attribute key, matched at a key boundary."""
        attributes = [
            "ID=g1;gene_id=ENSG00000157764.14;gene_name=BRAF;"
            "havana_gene=OTTHUMG00000157457.3",
            "ID=g2;description=BRAF pseudogene 1;gene=BRAFP1",
        ]
        self.assertEqual(self._gff_genes(attributes, tag="gene"), ["-", "BRAFP1"])
        self.assertEqual(
            self._gff_genes(attributes, tag="gene_id"), ["ENSG00000157764.14", "-"]
        )
        with self.assertLogs(level="WARNING"):
            genes = self._gff_genes(attributes, tag="(Name|gene)")
        self.assertEqual(genes, ["-", "-"])

    @staticmethod
    def _gff_genes(attributes, **kwargs):
        """Parse one GFF row per attribute string; return labels in input order."""
        lines = "".join(
            f"chr1\t.\tgene\t{100 * i + 1}\t{100 * i + 50}\t.\t+\t.\t{attr}\n"
            for i, attr in enumerate(attributes)
        )
        return list(gff.read_gff(StringIO(lines), **kwargs)["gene"])

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
        # Coordinates are taken as written: UCSC's are 0-based, as in BED
        ddx = regions.data[regions.data["accession"] == "NR_046018"]
        self.assertEqual(
            list(ddx[["gene", "chromosome", "start", "end"]].itertuples(index=False)),
            [("DDX11L1", "chr1", 11873, 14409)],
        )
        exons = tabio.read(fname, "refflat", exons=True).data
        ddx_exons = exons[exons["accession"] == "NR_046018"]
        self.assertEqual(list(ddx_exons["start"]), [11873, 12612, 13220])
        self.assertEqual(list(ddx_exons["end"]), [12227, 12721, 14409])

    def test_read_refflat_cds(self):
        """CDS mode emits coding regions and omits noncoding transcripts."""
        fname = "formats/refflat-mini.txt"
        regions = tabio.read(fname, "refflat", cds=True).data
        # 42 of the 100 transcripts are noncoding (cdsStart == cdsEnd)
        self.assertEqual(len(regions), 58)
        self.assertTrue((regions["start"] < regions["end"]).all())
        self.assertNotIn("NR_046018", set(regions["accession"]))
        fbxl22 = regions[regions["accession"] == "NM_203373"]
        self.assertEqual(list(fbxl22["start"]), [63889591])
        self.assertEqual(list(fbxl22["end"]), [63893885])

    def test_read_refgene(self):
        """Read the UCSC 'refGene' format, with a Table Browser header (#259)."""
        fname = "formats/refgene-mini.txt"
        regions = tabio.read(fname, "refgene")
        self.assertEqual(len(regions), linecount(fname) - 1)
        self.assertEqual(regions.sample_id, "refgene-mini")
        # Gene symbol from name2, accession from name, coordinates as written
        wash7p = regions.data[regions.data["accession"] == "NR_024540"]
        self.assertEqual(
            list(
                wash7p[["gene", "chromosome", "strand", "start", "end"]].itertuples(
                    index=False
                )
            ),
            [("WASH7P", "chr1", "-", 14361, 29370)],
        )
        # Starts and ends are not swapped, which the inverted-interval repair
        # would otherwise hide behind a warning
        with self.assertNoLogs(level="WARNING"):
            exons = tabio.read(fname, "refgene", exons=True).data
            cds = tabio.read(fname, "refgene", cds=True).data
        fam138a = exons[exons["accession"] == "NR_026818"]
        self.assertEqual(
            list(fam138a[["gene", "start", "end"]].itertuples(index=False)),
            [
                ("FAM138A", 34610, 35174),
                ("FAM138A", 35276, 35481),
                ("FAM138A", 35720, 36081),
            ],
        )
        fabp3 = cds[cds["accession"] == "NM_001320996"]
        self.assertEqual(
            list(fabp3[["gene", "start", "end"]].itertuples(index=False)),
            [("FABP3", 31365885, 31373014)],
        )
        # Headerless input reads the same
        with open(fname) as handle:
            header, *rows = handle
        self.assertTrue(header.startswith("#bin\tname\t"))
        headerless = tabio.read(StringIO("".join(rows)), "refgene")
        self.assertTrue(headerless.data.equals(regions.data))

    def test_genepred_family_agrees(self):
        """Each genePred-family layout of the same transcripts reads the same.

        The refFlat and genePred layouts are cut from the refGene fixture, so
        a column mapped to the wrong output in any reader shows up here.
        """
        with open("formats/refgene-mini.txt") as handle:
            rows = [line.rstrip("\n").split("\t") for line in handle][1:]
        others = {
            "genepredext": "formats/genepredext-mini.txt",
            # geneName (= name2), then name through exonEnds
            "refflat": [[r[12], *r[1:11]] for r in rows],
            # name through exonEnds
            "genepred": [r[1:11] for r in rows],
        }
        for kwargs in ({}, {"cds": True}, {"exons": True}):
            expect = tabio.read("formats/refgene-mini.txt", "refgene", **kwargs).data
            self.assertGreater(len(expect), 0)
            for fmt, source in others.items():
                with self.subTest(fmt=fmt, **kwargs):
                    infile = (
                        source
                        if isinstance(source, str)
                        else StringIO("".join("\t".join(r) + "\n" for r in source))
                    )
                    observed = tabio.read(infile, fmt, **kwargs).data
                    if fmt == "genepred":
                        # No gene symbol in this layout; the accession stands in
                        self.assertTrue(
                            (observed["gene"] == observed["accession"]).all()
                        )
                        observed = observed.assign(gene=expect["gene"])
                    self.assertTrue(observed.equals(expect))

    def test_read_genepred_family_header_only(self):
        """A table with only a '#' line and no rows reads as empty.

        The header line may be narrower than the table, and the Table Browser
        writes a one-field comment line when a query matches nothing.
        """
        for fmt, header in (
            ("refgene", "#bin\tname\tchrom\tstrand\ttxStart\n"),
            ("genepredext", "#name\tchrom\tstrand\ttxStart\n"),
            ("refflat", "#geneName\tname\tchrom\tstrand\ttxStart\n"),
            ("genepred", "# No results returned from query.\n"),
        ):
            with self.subTest(fmt=fmt):
                self.assertEqual(len(tabio.read(StringIO(header), fmt)), 0)

    def test_read_genepred_family_wrong_width(self):
        """A table wider or narrower than its format is refused, not misread."""
        with open("formats/refgene-mini.txt") as handle:
            rows = handle.read().splitlines()[1:4]
        # One trailing tab per row: still detected as refGene
        wide = "".join(row + "\t\n" for row in rows)
        with self.assertRaisesRegex(ValueError, "found more"):
            tabio.read_auto(StringIO(wide))
        with open("formats/refflat-mini.txt") as handle:
            narrow = handle.read()
        for kwargs in ({}, {"cds": True}, {"exons": True}):
            with (
                self.subTest(**kwargs),
                self.assertRaisesRegex(ValueError, "found fewer"),
            ):
                tabio.read(StringIO(narrow), "refgene", **kwargs)

    def test_read_refflat_with_header(self):
        """A leading Table Browser header line is skipped, and only that line.

        A '#' later in the table is part of the data, not a comment.
        """
        fname = "formats/refflat-mini.txt"
        with open(fname) as handle:
            first, rest = handle.readline(), handle.read()
        header = "#geneName\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\n"
        observed = tabio.read(StringIO(header + first + rest), "refflat")
        self.assertTrue(observed.data.equals(tabio.read(fname, "refflat").data))
        hashed = tabio.read(StringIO(header + first + "#" + rest), "refflat")
        self.assertEqual(len(hashed), len(observed))
        self.assertIn("#WASH7P", set(hashed["gene"]))

    def test_sniff_genepred_family(self):
        """Each UCSC gene table is detected as its own format, header or not."""
        for fname, fmt in (
            ("formats/refgene-mini.txt", "refgene"),
            ("formats/genepredext-mini.txt", "genepredext"),
            ("formats/refflat-mini.txt", "refflat"),
        ):
            with self.subTest(fname=fname):
                self.assertEqual(tabio.sniff_region_format(fname), fmt)
                with open(fname) as handle:
                    first_row = next(line for line in handle if line[0] != "#")
                for other in ("refgene", "genepredext", "refflat"):
                    self.assertEqual(
                        bool(tabio.format_patterns[other].match(first_row)),
                        other == fmt,
                        other,
                    )
        # Chromosome names need not be word characters, e.g. NCBI accessions
        with open("formats/refgene-mini.txt") as handle:
            _header, first_row, *_rows = handle
        ncbi_row = first_row.replace("\tchr1\t", "\tNC_000001.11\t", 1)
        self.assertNotEqual(ncbi_row, first_row)
        self.assertEqual(tabio.sniff_region_format(StringIO(ncbi_row)), "refgene")

    def test_skg_convert_gene_model_type(self):
        """skg_convert --refflat-type applies to refGene input too."""
        fname = "formats/refgene-mini.txt"
        n_exons = len(tabio.read(fname, "refgene", exons=True))
        n_transcripts = len(tabio.read(fname, "refgene"))
        self.assertGreater(n_exons, n_transcripts)
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_fname = os.path.join(tmp_dir, "exons.bed")
            args = argparse.Namespace(
                infile=fname,
                output=out_fname,
                in_fmt="auto",
                out_fmt="bed4",
                flatten=False,
                merge=None,
                gff_tag=None,
                gff_type=None,
                refflat_type="exon",
            )
            skg_convert.skg_convert(args)
            self.assertEqual(linecount(out_fname), n_exons)

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

    def test_read_gzip(self):
        """A gzipped file reads the same as its uncompressed equivalent."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            dict_fname = os.path.join(tmp_dir, "genome.dict")
            with open(dict_fname, "w") as handle:
                handle.write(
                    "@HD\tVN:1.6\n"
                    "@SQ\tSN:chr1\tLN:1000\tUR:file:g.fa\tM5:0\n"
                    "@SQ\tSN:chr2\tLN:500\tUR:file:g.fa\tM5:0\n"
                )
            for fname, fmt in (
                ("formats/amplicon.bed", "bed"),
                ("formats/amplicon.text", "text"),
                ("formats/acgh-log10.seg", "seg"),
                ("formats/sites-only.vcf", "vcf-simple"),
                (dict_fname, "dict"),
                ("formats/empty", "bed3"),
            ):
                with self.subTest(fname=fname, fmt=fmt):
                    expect = tabio.read(fname, fmt)
                    observed = tabio.read(gzipped(fname, tmp_dir), fmt)
                    self.assertTrue(observed.data.equals(expect.data))
                    self.assertEqual(observed.sample_id, expect.sample_id)

    def test_read_gzip_detected_by_content(self):
        """A gzipped file is decompressed even if its name lacks '.gz'."""
        expect = tabio.read("formats/amplicon.bed", "bed")
        with tempfile.TemporaryDirectory() as tmp_dir:
            fname = gzipped("formats/amplicon.bed", tmp_dir, "amplicon.bed")
            observed = tabio.read(fname, "bed")
        self.assertTrue(observed.data.equals(expect.data))

    def test_read_bgzf(self):
        """Every block of a BGZF (bgzip) file is read, not only the first."""
        fname = "formats/na12878-chrM-Y-trunc.bed.gz"
        with gzip.open(fname, "rt") as handle:
            expect_nrows = sum(1 for _line in handle)
        self.assertEqual(len(tabio.read(fname, "bed3")), expect_nrows)

    def test_read_auto_gzip(self):
        """Format detection and reading both see through gzip compression."""
        for fname in (
            "formats/amplicon.bed",
            "formats/amplicon.text",
            "formats/example.gff",
            "formats/nv2_baits.interval_list",
            "formats/refflat-mini.txt",
            "formats/refgene-mini.txt",
            "formats/genepredext-mini.txt",
        ):
            with self.subTest(fname=fname):
                expect = tabio.read_auto(fname)
                with tempfile.TemporaryDirectory() as tmp_dir:
                    gz_fname = gzipped(fname, tmp_dir)
                    self.assertEqual(
                        tabio.sniff_region_format(gz_fname),
                        tabio.sniff_region_format(fname),
                    )
                    observed = tabio.read_auto(gz_fname)
                self.assertTrue(observed.data.equals(expect.data))

    def test_sniff_gzip_keeps_extension_hint(self):
        """The '.bed' in 'x.bed.gz' still selects BED over interval-list.

        A BED row with the strand in column 4 also matches the Picard
        interval-list pattern, whose reader converts 1-based starts to 0-based;
        without the extension to settle it, every start would shift by one.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            bed_fname = os.path.join(tmp_dir, "stranded.bed")
            with open(bed_fname, "w") as handle:
                handle.write("chr1\t100\t200\t+\tgeneA\n")
            for suffix in (".gz", ".bgz", ".GZ"):
                with self.subTest(suffix=suffix):
                    gz_fname = gzipped(bed_fname, tmp_dir, "stranded.bed" + suffix)
                    self.assertEqual(tabio.sniff_region_format(gz_fname), "bed")
                    self.assertEqual(tabio.read_auto(gz_fname).start.iat[0], 100)

    @unittest.skipUnless(hasattr(os, "mkfifo"), "named pipes unavailable")
    def test_read_named_pipe(self):
        """A named pipe, as from process substitution, is read whole."""
        fname = "formats/amplicon.bed"
        expect = tabio.read(fname, "bed")
        with open(fname, "rb") as handle:
            plain_bytes = handle.read()
        for label, payload in (
            ("plain", plain_bytes),
            ("gzip", gzip.compress(plain_bytes)),
        ):
            with self.subTest(label), tempfile.TemporaryDirectory() as tmp_dir:
                fifo = os.path.join(tmp_dir, "amplicon.bed")
                os.mkfifo(fifo)
                result = {}

                def feed(path=fifo, data=payload):
                    with open(path, "wb") as pipe:
                        pipe.write(data)

                def consume(path=fifo, out=result):
                    try:
                        out["data"] = tabio.read(path, "bed").data
                    except Exception as exc:
                        out["error"] = exc

                # Daemon threads, so that a reader that opens the pipe twice
                # fails this test after the timeout instead of hanging it.
                writer = threading.Thread(target=feed, daemon=True)
                reader = threading.Thread(target=consume, daemon=True)
                writer.start()
                reader.start()
                reader.join(timeout=10)
                self.assertFalse(reader.is_alive(), "Reading the pipe stalled")
                if "error" in result:
                    raise result["error"]
                self.assertTrue(result["data"].equals(expect.data))

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

    def test_read_vcf_end_below_pos_uses_the_footprint(self):
        """An END beneath the record's own POS states no span, so it is dropped.

        The text readers honoured such an END, emitting ``end < start`` for
        END below POS and the zero-width row ``start == end`` for END
        exactly one below it, from a file the pysam reader read as an
        ordinary variant. ``read`` deliberately does not repair an inverted
        VCF row, which is why the discard belongs in the reader itself.

        The last two records are not discarded, and are what keep the
        warning's counts honest: one declares a usable END and one declares
        none at all, so the number discarded, the number declaring an END,
        and the number of records are three different values.
        """
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t200\t.\tA\tG\t50\tPASS\tEND=199\n"
            "chr1\t500\t.\tACGT\tA\t50\tPASS\tEND=400\n"
            "chr1\t800\t.\tACGT\tA\t50\tPASS\tEND=805\n"
            "chr1\t900\t.\tACGT\tA\t50\tPASS\t.\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            for fmt in ("vcf", "vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.start), [199, 499, 799, 899], fmt)
                self.assertEqual(list(dframe.end), [200, 503, 805, 903], fmt)
            # Only the text reader's warning reaches Python logging; htslib
            # writes its own "INFO/END ... smaller than POS" to C stderr,
            # and the pysam reader separately warns about absent genotypes.
            with self.assertLogs(level="WARNING") as cm:
                tabio.read(tmp.name, "vcf-simple")
            message = " ".join(cm.output)
            self.assertIn("in 2 of 4 records", message)
            self.assertIn("chr1:200 with END=199", message)

    def test_read_vcf_span_is_the_furthest_of_end_and_the_ref_allele(self):
        """A declared END shorter than the reference allele does not shorten it.

        The specification defines a record's span as the maximum of the
        reference allele's length and the length INFO/END implies, so a
        record quoting four reference bases while declaring a one-base END
        still covers four. htslib disagreed with itself about this: through
        1.21, which pysam 0.23.3 bundles and ``requirements/min.txt`` pins,
        it took the shorter declared span; from 1.23 it takes the maximum.
        CNVkit computes the maximum itself so that a run's coordinates do
        not depend on which pysam wheel is installed, which is why all three
        readers are asserted here and why this passes under either.
        """
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t800\t.\tACGT\tA\t50\tPASS\tEND=800\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            for fmt in ("vcf", "vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.start), [799], fmt)
                self.assertEqual(list(dframe.end), [803], fmt)

    def test_read_vcf_svlen_reaches_the_reference_only_for_some_alleles(self):
        """SVLEN sizes a symbolic deletion, but never an insertion.

        A symbolic ALT quotes one reference base and states its extent in
        INFO/SVLEN, so a reader taking the reference allele alone reads a
        400-base deletion as one base. SVLEN counts reference bases only
        where the allele replaces reference: for an insertion it counts
        inserted bases, and reading it as a span was an htslib defect
        (its issue #1940) that 1.23 fixed by restricting the set. The sign
        is not read, since VCF 4.2 writers make a deletion negative and 4.4
        writers make it positive.

        The inversion is the case to keep: htslib's own 1.23 release note
        names only DEL, DUP and CNV, while its code also accepts INV. This
        follows the code, so a future reader who trusts the changelog and
        drops the case will see this fail. The last two records pin the
        match's shape, which is a type prefix followed by the closing
        bracket or a sub-type separator: ``<DEL:ME>`` is a deletion and
        ``<DELETION>`` is not one of these types at all.
        """
        vcf_text = (
            "##fileformat=VCFv4.4\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">\n'
            '##ALT=<ID=DEL,Description="Deletion">\n'
            '##ALT=<ID=INV,Description="Inversion">\n'
            '##ALT=<ID=INS,Description="Insertion">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\t.\tA\t<DEL>\t50\tPASS\tSVLEN=-400\n"
            "chr1\t2000\t.\tA\t<INV>\t50\tPASS\tSVLEN=400\n"
            "chr1\t3000\t.\tA\t<INS>\t50\tPASS\tSVLEN=400\n"
            "chr1\t4000\t.\tA\t<DEL:ME>\t50\tPASS\tSVLEN=400\n"
            "chr1\t5000\t.\tA\t<DELETION>\t50\tPASS\tSVLEN=400\n"
            "chr1\t6000\t.\tA\tAACGT\t50\tPASS\tSVLEN=4\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            for fmt in ("vcf", "vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(
                    list(dframe.start), [999, 1999, 2999, 3999, 4999, 5999], fmt
                )
                self.assertEqual(
                    list(dframe.end), [1400, 2400, 3000, 4400, 5000, 6000], fmt
                )

    def test_read_vcf_svlen_pairs_with_its_own_alt_allele(self):
        """SVLEN carries one value per ALT, and the two pair by position.

        A record offering both a deletion and an insertion states a length
        for each, and only the deletion's reaches the reference. A reader
        taking the largest value it could find would read the second record
        as covering 900 bases rather than 400. The third pins that fewer
        SVLEN values than alleles is not an error: the pairing runs out, as
        it does in htslib.

        The spans are compared rather than the frames, because the
        pysam-backed reader emits a row per alternate allele where the text
        readers emit one per record. Every row of a record carries the same
        span, which is the point: how far a record reaches is a property of
        the record, not of whichever allele proved longest.
        """
        vcf_text = (
            "##fileformat=VCFv4.4\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">\n'
            '##ALT=<ID=DEL,Description="Deletion">\n'
            '##ALT=<ID=DUP,Description="Duplication">\n'
            '##ALT=<ID=INS,Description="Insertion">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\t.\tA\t<DEL>,<DUP>\t50\tPASS\tSVLEN=400,900\n"
            "chr1\t2000\t.\tA\t<DEL>,<INS>\t50\tPASS\tSVLEN=400,900\n"
            "chr1\t3000\t.\tA\t<DEL>,<DUP>\t50\tPASS\tSVLEN=400\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            self.assertEqual(
                list(tabio.read(tmp.name, "vcf").data.end),
                [1900, 1900, 2400, 2400, 3400, 3400],
            )
            for fmt in ("vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.end), [1900, 2400, 3400], fmt)

    def test_read_vcf_svlen_outside_htslibs_range_is_ignored(self):
        """An SVLEN too extreme for htslib to keep is not a span here either.

        htslib holds an INFO integer in a 32-bit field and replaces anything
        outside it with the missing value, so an absurd SVLEN leaves the
        record spanning its reference allele. The readers here match that
        range rather than their own, which also keeps a value near the 64-bit
        limit from overflowing the coordinate column: computing the span
        first and range-checking afterwards raised OverflowError out of the
        text readers on the last record below.

        The range is not symmetric, which the first two records pin. Its
        lower bound stops eight short of the 32-bit minimum, a margin
        reserved for the markers that mean "missing" and "end of vector", so
        a magnitude test would honour a few large negative values that htslib
        discards.
        """
        vcf_text = (
            "##fileformat=VCFv4.4\n"
            "##contig=<ID=chr1,length=100000000>\n"
            '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">\n'
            '##ALT=<ID=DEL,Description="Deletion">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tN\t<DEL>\t50\tPASS\tSVLEN=-2147483640\n"
            "chr1\t200\t.\tN\t<DEL>\t50\tPASS\tSVLEN=-2147483641\n"
            "chr1\t300\t.\tN\t<DEL>\t50\tPASS\tSVLEN=2147483648\n"
            "chr1\t400\t.\tN\t<DEL>\t50\tPASS\tSVLEN=9223372036854775807\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            for fmt in ("vcf", "vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.start), [99, 199, 299, 399], fmt)
                self.assertEqual(list(dframe.end), [2147483740, 200, 300, 400], fmt)

    def test_read_vcf_end_without_a_declaration_parts_the_reader_families(self):
        """An END the header never declares is honoured only by the text readers.

        htslib declines to type an undeclared key ("assuming Type=String")
        and reports the reference footprint, while the text readers parse
        the value out of the raw INFO string and take it. That divergence
        is deliberate and pinned rather than removed: it is the only END
        input class where the three readers disagree.
        """
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t500\t.\tACGT\tA\t50\tPASS\tEND=900\n"
        )
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
            tmp.write(vcf_text)
            tmp.flush()
            pysam_frame = tabio.read(tmp.name, "vcf").data
            self.assertEqual(list(pysam_frame.start), [499])
            self.assertEqual(list(pysam_frame.end), [503])
            for fmt in ("vcf-simple", "vcf-sites"):
                dframe = tabio.read(tmp.name, fmt).data
                self.assertEqual(list(dframe.start), [499], fmt)
                self.assertEqual(list(dframe.end), [900], fmt)

    def test_read_vcf_sites_only(self):
        """Read a sites-only VCF, where depth and frequency come from INFO.

        Stripping the sample columns from a caller's output -- `bcftools
        view -G`, GATK SelectVariants -- leaves the INFO fields untouched,
        so the reader takes depth and allele frequency from there instead
        of from a genotype. GATK declares AF as Number=A, one value per
        ALT, which pysam returns as a tuple.
        """
        varr = tabio.read("formats/sites-only.vcf", "vcf")
        # The 2-ALT record yields a row per ALT, both carrying the record's
        # own depth and its first ALT's frequency (see the XXX in vcfio.py)
        self.assertEqual(len(varr), 6)
        self.assertEqual(list(varr["depth"])[:4], [120, 100, 90, 80])
        # AF drives the alt count and the inferred zygosity
        self.assertEqual(list(varr["alt_count"])[:3], [60, 48, 90])
        self.assertEqual(list(varr["zygosity"])[:3], [0.5, 0.5, 1.0])
        # A '.' depth is missing, not zero-coverage, and must not crash
        self.assertEqual(list(varr["depth"])[-1], 0.0)
        # Only the AF=1.000 record is homozygous
        self.assertEqual(len(varr.heterozygous()), 5)

    def test_read_vcf_info_values_are_reduced_to_one_number(self):
        """Depth survives every shape pysam can return an INFO value in.

        A VCF header decides how pysam types an INFO value, and a reader
        that assumes a number gets four other things instead: a tuple for
        any Number=A/R/G/. field, None for the '.' missing marker, a
        string for a key the header never declares -- htslib admits those
        mid-parse, so the same key can raise on one record and return a
        string on the next -- and, for a key both undeclared and absent, a
        ValueError. Each must reduce to a number or to the no-depth 0.0,
        and a depth the file declared an integer must stay one: `export
        theta` writes the column, so widening it is visible to a reader.
        """
        header = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100000>\n"
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Frequency">\n'
        )
        declare_dp = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        # The expected depth's own type is the column type it must produce
        for label, decl, info, depth in (
            ("declared and present", declare_dp, "DP=34;AF=0.5", 34),
            ("declared, value '.'", declare_dp, "DP=.;AF=0.5", 0.0),
            ("declared, absent", declare_dp, "AF=0.5", 0.0),
            ("undeclared, present", "", "DP=34;AF=0.5", 34.0),
            ("undeclared, absent", "", "AF=0.5", 0.0),
            ("undeclared, multi-valued", "", "DP=34,35;AF=0.5", 34.0),
            # A declared field's missing value is parsed into None, but an
            # undeclared one reaches Python as the marker's own text
            ("undeclared, value '.'", "", "DP=.;AF=0.5", 0.0),
        ):
            vcf_text = (
                header + decl + columns + f"chr1\t500\t.\tA\tG\t50\tPASS\t{info}\n"
            )
            with tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf") as tmp:
                tmp.write(vcf_text)
                tmp.flush()
                varr = tabio.read(tmp.name, "vcf")
                self.assertEqual(list(varr["depth"]), [depth], label)
                self.assertEqual(
                    varr["depth"].dtype.kind,
                    "i" if isinstance(depth, int) else "f",
                    label,
                )

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


class VcfSamplePairingTests(unittest.TestCase):
    """A caller's sample IDs outrank the VCF header's declared pairing.

    A ``PEDIGREE`` tag records how the file was made, which is not always how
    it is being analyzed, and re-heading a VCF to correct one is expensive.
    The header therefore fills in what the caller left out rather than
    replacing what the caller supplied -- but displacing a declared pairing
    is worth saying out loud, so it is logged.
    """

    _FIXTURE = "formats/na12878_na12882_mix.vcf"
    # ##PEDIGREE=<Derived=NA12882,Original=NA12878>
    _DECLARED_TUMOR = "NA12882"
    _DECLARED_NORMAL = "NA12878"
    _OVERRIDE_MSG = "overriding the tumor-normal pairing"

    _HEADER = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )

    def _write(self, pedigrees, samples):
        """Write a minimal VCF declaring `pedigrees` over `samples`."""
        body = self._HEADER
        for tumor, normal in pedigrees:
            body += f"##PEDIGREE=<Derived={tumor},Original={normal}>\n"
        body += (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\nchr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t"
            + "\t".join(["0/1"] * len(samples))
            + "\n"
        )
        tmp = tempfile.NamedTemporaryFile(mode="w+t", suffix=".vcf", delete=False)
        self.addCleanup(os.remove, tmp.name)
        with tmp:
            tmp.write(body)
        return tmp.name

    def _pair(self, fname, sample_id=None, normal_id=None):
        """Resolve (tumor, normal) as the reader would, without loading rows.

        ``VariantArray.sample_id`` falls back to the filename when the caller
        names no sample, so it cannot witness which sample the header chose.
        """
        with pysam.VariantFile(fname) as vcf:
            return vcfio._choose_samples(vcf, sample_id, normal_id)

    def test_reversed_pair_is_honored_and_logged(self):
        """Naming the declared pair backwards pairs it backwards, loudly.

        Reading the roles from the header instead left the array unpaired,
        which is indistinguishable from passing no normal at all.
        """
        with self.assertLogs(level="WARNING") as ctx:
            varr = tabio.read(
                self._FIXTURE,
                "vcf",
                sample_id=self._DECLARED_NORMAL,
                normal_id=self._DECLARED_TUMOR,
            )
        self.assertEqual(varr.sample_id, self._DECLARED_NORMAL)
        self.assertIn("n_zygosity", varr.data.columns)
        message = "\n".join(ctx.output)
        self.assertIn(self._OVERRIDE_MSG, message)
        # Both the pair used and the pair displaced are named
        self.assertIn(f"'{self._DECLARED_NORMAL}'", message)
        self.assertIn(f"'{self._DECLARED_TUMOR}'", message)

    def test_declared_pair_is_used_without_complaint(self):
        """Agreeing with the header is the common case and must stay quiet.

        Naming the pair the header already declares is the ordinary way to
        be explicit about it, so a warning there would fire on runs where
        nothing was displaced and teach readers to ignore the message.
        """
        for sid, nid in [
            (self._DECLARED_TUMOR, self._DECLARED_NORMAL),
            (self._DECLARED_TUMOR, None),
            (None, self._DECLARED_NORMAL),
            (None, None),
        ]:
            with self.subTest(sample_id=sid, normal_id=nid):
                with self.assertNoLogs(level="WARNING"):
                    pair = self._pair(self._FIXTURE, sid, nid)
                self.assertEqual(pair, (self._DECLARED_TUMOR, self._DECLARED_NORMAL))

    def test_normal_alone_draws_its_own_declared_tumor(self):
        """A normal names a pairing, not merely a sample to exclude.

        Pairing it with whichever other sample comes first in the file would
        cross two declared pairs whenever the column order disagrees with the
        header, which is why the declared pairings are consulted first.
        """
        fname = self._write([("A", "B"), ("C", "D")], ["C", "A", "B", "D"])
        self.assertEqual(self._pair(fname, normal_id="B"), ("A", "B"))
        self.assertEqual(self._pair(fname, normal_id="D"), ("C", "D"))

    def test_normal_outside_the_declared_pair_is_honored(self):
        """A third sample may be the normal even when the header names another.

        This shape kept its ``n_*`` columns before, so the substitution of a
        different individual's genotypes left no trace in the output's shape.
        """
        fname = self._write([("A", "B")], ["A", "B", "C"])
        with self.assertLogs(level="WARNING") as ctx:
            varr = tabio.read(fname, "vcf", sample_id="A", normal_id="C")
        self.assertEqual(varr.sample_id, "A")
        self.assertIn(self._OVERRIDE_MSG, "\n".join(ctx.output))
        self.assertIn("'C'", "\n".join(ctx.output))

    def test_sample_that_is_no_declared_tumor_is_unpaired_and_logged(self):
        """Requesting the control sample alone still yields it, now with a note.

        Selecting a normal by ``sample_id`` is supported, but dropping to an
        unpaired analysis silently is how a missing ``baf`` column became
        indistinguishable from a VCF that never had genotypes.
        """
        with self.assertLogs(level="WARNING") as ctx:
            varr = tabio.read(self._FIXTURE, "vcf", sample_id=self._DECLARED_NORMAL)
        self.assertEqual(varr.sample_id, self._DECLARED_NORMAL)
        self.assertNotIn("n_zygosity", varr.data.columns)
        self.assertIn("as unpaired", "\n".join(ctx.output))

    def test_explicit_pair_is_not_filtered_through_the_header(self):
        """A named pair is used as given, not looked up among the declared ones.

        The header here declares a tumor the file no longer contains, which
        subsetting a joint call routinely leaves behind. Consulting the
        declared pairings for this normal and then keeping only the rows whose
        tumor was named would filter the pair away to nothing, and the caller
        would get their sample back unpaired -- which is what happened before.
        A pairing the header declares nowhere still needs no header lookup.

        Omitting the tumor is a different matter: the stale pairing is then
        the only candidate, and naming a sample the file lacks is still an
        error.
        """
        fname = self._write([("GHOST", "B")], ["A", "B"])
        with self.assertLogs(level="WARNING"):
            varr = tabio.read(fname, "vcf", sample_id="A", normal_id="B")
        self.assertEqual(varr.sample_id, "A")
        self.assertIn("n_zygosity", varr.data.columns)
        with self.assertRaises(IndexError):
            tabio.read(fname, "vcf", normal_id="B")

    def test_sample_cannot_be_its_own_normal(self):
        """One sample named as both halves has no reading worth guessing at."""
        for fname in (self._FIXTURE, self._write([], ["A", "B"])):
            sid = self._DECLARED_TUMOR if fname == self._FIXTURE else "A"
            with self.subTest(vcf=fname):
                with self.assertRaises(IndexError) as ctx:
                    tabio.read(fname, "vcf", sample_id=sid, normal_id=sid)
                self.assertIn("its own matched normal", str(ctx.exception))


class IntervalOrderTests(unittest.TestCase):
    """A row whose end precedes its start is read according to who wrote it.

    CNVkit's own bin and segment tables are refused: such a row is corrupt
    output and any repair would invent a span no segmenter produced. A
    third-party region file is repaired to span [min, max] instead, because
    reversing an interval is how some tools record orientation.
    """

    HEADER = "chromosome\tstart\tend\tgene\tlog2\tprobes\n"
    ROWS = "chr1\t100\t200\tA\t0.0\t5\nchr1\t400\t300\tB\t-1.0\t0\n"

    def _write(self, suffix, text):
        tmp = tempfile.NamedTemporaryFile(
            mode="w+t", suffix=suffix, delete=False, encoding="utf-8"
        )
        self.addCleanup(os.unlink, tmp.name)
        with tmp:
            tmp.write(text)
        return tmp.name

    def test_own_format_is_refused(self):
        fname = self._write(".cns", self.HEADER + self.ROWS)
        with self.assertRaises(ValueError) as ctx:
            tabio.read(fname, "tab")
        msg = str(ctx.exception)
        # Name the offending row as the file spells it, not as repaired
        self.assertIn("chr1:400-300", msg)
        self.assertIn("1 of 2 rows", msg)
        self.assertIn(fname, msg)

    def test_seg_is_refused(self):
        fname = self._write(
            ".seg",
            "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n"
            "s\t1\t400\t300\t0\t-1.0\n",
        )
        with self.assertRaises(ValueError):
            tabio.read(fname, "seg")

    def test_bed_is_reversed_with_a_warning(self):
        fname = self._write(".bed", "chr1\t100\t200\tFWD\nchr1\t400\t300\tREV\n")
        with self.assertLogs(level="WARNING") as cm:
            regions = tabio.read(fname, "bed")
        self.assertIn("chr1:400-300", " ".join(cm.output))
        self.assertEqual(
            regions.data[["start", "end"]].to_numpy().tolist(),
            [[100, 200], [300, 400]],
        )
        # The repair happens before sorting, so the result is still ordered
        self.assertTrue(regions.data.start.is_monotonic_increasing)

    def test_only_reversed_rows_are_not_collapsed(self):
        """Repairing every row must not collapse both coordinates onto one.

        Reading a column out of a `.loc` slice gives a view rather than a
        copy when the mask selects every row, so a swap through one uncopied
        array leaves start and end both equal to the end. A reader hands back
        start and end in a single consolidated block, which is the
        arrangement that has it; the frame `read` goes on to build does not,
        so driving the helper directly is what makes this able to fail.
        """
        fname = self._write(".bed", "chr1\t400\t300\tA\nchr1\t900\t800\tB\n")
        frame = bedio.read_bed(fname)
        with self.assertLogs(level="WARNING"):
            tabio._apply_interval_order_policy(frame, "bed", fname)
        self.assertEqual(
            frame[["start", "end"]].to_numpy().tolist(), [[300, 400], [800, 900]]
        )
        self.assertEqual(frame["start"].dtype, "int64")

    def test_single_row_file_is_reversed_through_read(self):
        """One row, reversed, is the smallest whole-file repair `read` does."""
        fname = self._write(".bed", "chr1\t400\t300\tONE\n")
        with self.assertLogs(level="WARNING"):
            regions = tabio.read(fname, "bed")
        self.assertEqual(
            regions.data[["start", "end"]].to_numpy().tolist(), [[300, 400]]
        )
        self.assertEqual(regions.data["start"].dtype, "int64")

    def test_import_seg_repairs_rather_than_refusing(self):
        """A SEG from another tool is repaired, so the .cns it makes is readable.

        `import-seg` exists to ingest foreign SEG and does not go through
        `read`, so before it repaired at ingest it wrote a .cns holding the
        reversed row -- which every later command then refused, blaming
        CNVkit's own segmenter for a row CNVkit never wrote.
        """
        segfile = self._write(
            ".seg",
            "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n"
            "foo\t1\t100\t200\t10\t0.1\n"
            "foo\t1\t400\t300\t0\t-1.0\n",
        )
        outdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, outdir)
        args = commands.parse_args(["import-seg", segfile, "-d", outdir])
        with self.assertLogs(level="WARNING") as cm:
            args.func(args)
        self.assertIn("1:399-300", " ".join(cm.output))
        written = os.path.join(outdir, "foo.cns")
        # The whole point: what it wrote, CNVkit can read
        segs = tabio.read(written, "tab")
        self.assertEqual(len(segs), 2)
        self.assertTrue((segs.data["end"] >= segs.data["start"]).all())

    def test_well_formed_files_are_untouched(self):
        """Neither branch fires on a file that is already in order.

        Comparing two reads of the same file would only assert that reading
        is deterministic. What has to hold is that the policy stays silent:
        no repair, and no refusal.
        """
        for fname, fmt in (
            ("formats/amplicon.cns", "tab"),
            ("formats/amplicon.bed", "bed"),
            ("formats/example.gff", "gff"),
            ("formats/empty", "tab"),
        ):
            with self.subTest(fname=fname), self.assertNoLogs(level="WARNING"):
                arr = tabio.read(fname, fmt)
            self.assertTrue((arr.data["end"] >= arr.data["start"]).all())

    def test_every_reader_is_classified(self):
        """A reader added to READERS without a classification fails here.

        The VCF family is the deliberate exemption;
        `skgenome.tabio._apply_interval_order_policy` says why.
        """
        classified = tabio._SELF_WRITTEN_FORMATS | tabio._FOREIGN_FORMATS
        self.assertEqual(
            set(tabio.READERS) - classified - {"auto"},
            {"vcf", "vcf-simple", "vcf-sites"},
        )
        self.assertFalse(tabio._SELF_WRITTEN_FORMATS & tabio._FOREIGN_FORMATS)
        self.assertTrue(classified <= set(tabio.READERS))


if __name__ == "__main__":
    unittest.main(verbosity=2)
