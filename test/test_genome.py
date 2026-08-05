#!/usr/bin/env python
"""Unit tests for the 'genome' sub-package."""

import functools
import operator
import random
import typing
import unittest

import numpy as np
import pandas as pd
from pandas.api.types import is_integer_dtype

from cnvlib import read_ga
from skgenome import GenomicArray as GA
from skgenome import chromsort, rangelabel, tabio
from skgenome.combiners import join_strings


class GaryTests(unittest.TestCase):
    def setUp(self):
        self.ex_cnr = read_ga("formats/reference-tr.cnn")

    def test_empty(self):
        """Instantiate from an empty file."""
        garr = tabio.read("formats/empty")
        self.assertEqual(len(garr), 0)

    def test_iter(self):
        """Test iteration."""
        rows = iter(self.ex_cnr)
        firstrow = next(rows)
        self.assertEqual(tuple(firstrow), tuple(self.ex_cnr[0]))
        i = 0
        for i, _row in enumerate(rows):
            pass
        self.assertEqual(i + 2, len(self.ex_cnr))

    def test_copy(self):
        """Test creation of an independent copy of the object."""
        dupe = self.ex_cnr.copy()
        self.assertEqual(tuple(self.ex_cnr[3]), tuple(dupe[3]))
        self.ex_cnr[3, "log2"] = -10.0
        self.assertNotEqual(tuple(self.ex_cnr[3]), tuple(dupe[3]))

    def test_autosomes(self):
        """Test selection of autosomes."""
        len_all = len(self.ex_cnr)
        len_x = (self.ex_cnr.chromosome == "chrX").sum()
        len_y = (self.ex_cnr.chromosome == "chrY").sum()
        auto = self.ex_cnr.autosomes()
        self.assertEqual(len(auto), len_all - len_x - len_y)
        autox = self.ex_cnr.autosomes(also="chrX")
        self.assertEqual(len(autox), len_all - len_y)
        autoy = self.ex_cnr.autosomes(also=["chrY"])
        self.assertEqual(len(autoy), len_all - len_x)
        autoxy = self.ex_cnr.autosomes(also=["chrX", "chrY"])
        self.assertEqual(
            len(autoxy), len_all, "It's possible to provide chromosome names."
        )
        some_x = (self.ex_cnr.chromosome == "chrX") & (self.ex_cnr.end <= 434918)
        some_x_len = some_x.sum()
        self.assertEqual(some_x_len, 3)
        auto_and_some_x = self.ex_cnr.autosomes(also=some_x)
        self.assertEqual(
            len(auto_and_some_x),
            len(auto) + some_x_len,
            "It's possible to provide a Pandas filter.",
        )

    def test_autosomes_yeast(self):
        """Roman-numeral chromosomes (yeast) are recognized as autosomes."""
        rows = [
            (f"chr{r}", 0, 1000) for r in ("I", "II", "VIII", "IX", "X", "XI", "XVI")
        ] + [("chrM", 0, 1000)]
        yeast = GA.from_rows(rows)
        auto = yeast.autosomes()
        # All seven Roman-numeral chromosomes are autosomes; chrM is excluded.
        self.assertEqual(len(auto), 7)
        self.assertNotIn("chrM", set(auto.chromosome))
        # chrX is autosome 10 in yeast, not a sex chromosome.
        self.assertIn("chrX", set(auto.chromosome))

    def test_autosomes_unrecognized_genome(self):
        """An unfamiliar chromosome set falls back permissively with a warning."""
        rows = [("custom_a", 0, 100), ("custom_b", 0, 100), ("scaffold_1", 0, 100)]
        ga = GA.from_rows(rows)
        with self.assertLogs(level="WARNING") as cm:
            auto = ga.autosomes()
        self.assertEqual(len(auto), len(ga))
        self.assertTrue(any("no autosomes recognized" in m for m in cm.output))

    def test_by_chromosome(self):
        for fname in ("formats/amplicon.cnr", "formats/cl_seq.cns"):
            cnarr = read_ga(fname)
            row_count = 0
            for _chrom, rows in cnarr.by_chromosome():
                row_count += len(rows)
            self.assertEqual(row_count, len(cnarr))

    def test_filter(self):
        """Test sugary selection of a subset of the data array."""
        num_bg_rows = len(self.ex_cnr[self.ex_cnr["gene"] == "Background"])
        self.assertEqual(len(self.ex_cnr.filter(gene="Background")), num_bg_rows)
        selector = lambda row: row["gene"] == "Background"
        self.assertEqual(len(self.ex_cnr.filter(selector)), num_bg_rows)

    def test_ranges_by_in(self):
        """Test range methods: by_ranges, in_range, in_ranges."""
        cnarr = read_ga("formats/amplicon.cnr")
        segarr = read_ga("formats/amplicon.cns")
        chrom_segarr = dict(segarr.by_chromosome())
        for chrom, subarr in cnarr.by_chromosome():
            count_segs = 0
            count_bins = 0
            subsegarr = chrom_segarr[chrom]
            for count_segs, (seg, bins) in enumerate(subarr.by_ranges(subsegarr)):
                count_bins += len(bins)
                self.assertEqual(seg.probes, len(bins))
                self.assertEqual(
                    len(bins),
                    len(
                        cnarr.in_range(seg.chromosome, seg.start, seg.end, mode="outer")
                    ),
                )
                self.assertEqual(
                    len(bins),
                    len(
                        cnarr.in_range(seg.chromosome, seg.start, seg.end, mode="trim")
                    ),
                )
            self.assertEqual(len(subsegarr), count_segs + 1)
            self.assertEqual(len(subarr), count_bins)
            self.assertEqual(
                len(subarr),
                len(
                    cnarr.in_ranges(
                        chrom, subsegarr["start"], subsegarr["end"], mode="outer"
                    )
                ),
            )
            self.assertEqual(
                len(subarr),
                len(
                    subarr.in_ranges(
                        starts=subsegarr["start"], ends=subsegarr["end"], mode="outer"
                    )
                ),
            )
            self.assertEqual(
                len(subarr),
                len(
                    cnarr.in_ranges(
                        chrom, subsegarr["start"], subsegarr["end"], mode="trim"
                    )
                ),
            )
            self.assertEqual(
                len(subarr),
                len(
                    subarr.in_ranges(
                        starts=subsegarr["start"], ends=subsegarr["end"], mode="trim"
                    )
                ),
            )

    def test_ranges_into(self):
        cnarr = read_ga("formats/amplicon.cnr")
        segarr = read_ga("formats/amplicon.cns")
        seg_genes = cnarr.into_ranges(segarr, "gene", "-")
        self.assertEqual(len(seg_genes), len(segarr))
        # With a VCF
        varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        seg_baf = varr.into_ranges(segarr, "alt_freq", np.nan, np.nanmedian)
        self.assertEqual(len(seg_baf), len(segarr))
        cna_baf = varr.into_ranges(cnarr, "alt_freq", 0.0, np.max)
        self.assertEqual(len(cna_baf), len(cnarr))
        # Edge cases: an empty source or destination still yields one value per
        # destination row, as a Series indexed by the destination
        mtarr = tabio.read("formats/empty")
        self.assertEqual(len(segarr.into_ranges(mtarr, "start", 0, int)), 0)
        empty_src = mtarr.into_ranges(segarr, "end", 0, 0)
        self.assertIsInstance(empty_src, pd.Series)
        self.assertEqual(empty_src.index.tolist(), segarr.data.index.tolist())
        self.assertTrue((empty_src == 0).all())

    @staticmethod
    def _non_contiguous_pair():
        """A source of two named ranges, and a destination with interleaved rows.

        Source: chr2:0-1000 "AAA", chr10:0-1000 "CCC". Destination rows, in
        order: chr2:0-1000 (hits AAA), chr10:0-1000 (hits CCC), chr2:2000-3000
        and chr10:2000-3000 (no overlap).

        Interleaved chromosome rows arise only for an in-memory array assembled
        by concatenation and never re-sorted, since every file-based route is
        sorted by ``tabio.read``. chr2 before chr10 makes first-appearance order
        differ from the lexicographic order pandas would use with
        ``groupby(sort=True)``, so recovering the groups by sorting fails here
        too.
        """
        src = GA(
            pd.DataFrame(
                {
                    "chromosome": ["chr2", "chr10"],
                    "start": [0, 0],
                    "end": [1000, 1000],
                    "gene": ["AAA", "CCC"],
                }
            )
        )
        dest = GA(
            pd.DataFrame(
                {
                    "chromosome": ["chr2", "chr10", "chr2", "chr10"],
                    "start": [0, 0, 2000, 2000],
                    "end": [1000, 1000, 3000, 3000],
                }
            )
        )
        return src, dest

    def test_ranges_into_non_contiguous_chromosomes(self):
        """into_ranges places each value on its own row of the destination.

        ``by_shared_chroms`` groups the destination by chromosome, so consuming
        it directly would summarize into chromosome-grouped order rather than
        the destination's row order.
        """
        src, dest = self._non_contiguous_pair()
        self.assertEqual(
            src.into_ranges(dest, "gene", "-").tolist(), ["AAA", "CCC", "-", "-"]
        )

    def test_ranges_into_preserves_dest_index(self):
        """into_ranges is indexed by the destination's labels, not 0..n-1.

        Callers assign the result straight onto a column of the destination,
        and pandas aligns that by label. A filtered destination has a gapped
        index -- ``target --annotate`` drops zero-width baits that way -- so a
        fresh RangeIndex would drop values onto the wrong rows.
        """
        src, dest = self._non_contiguous_pair()
        # Dropping row 1 leaves chr2, chr2, chr10 -- contiguous, so this test
        # pins index alignment alone, independent of the ordering fix.
        gapped = dest[np.array([True, False, True, True])]
        self.assertEqual(gapped.data.index.tolist(), [0, 2, 3])
        genes = src.into_ranges(gapped, "gene", "-")
        self.assertEqual(genes.index.tolist(), [0, 2, 3])
        gapped["gene"] = genes
        self.assertEqual(gapped["gene"].tolist(), ["AAA", "-", "-"])
        # The no-such-column fallback must be indexed the same way
        with self.assertLogs(level="WARNING"):
            missing = src.into_ranges(gapped, "nosuchcolumn", "-")
        self.assertEqual(missing.index.tolist(), [0, 2, 3])

    def test_ranges_of(self):
        cnarr = read_ga("formats/amplicon.cnr")
        segarr = read_ga("formats/amplicon.cns")
        by_bins = cnarr.by_ranges(segarr)
        by_slices = cnarr.iter_ranges_of(segarr, "gene")
        for (_seg, by_bin), by_slice in zip(by_bins, by_slices, strict=True):
            self.assertEqual(len(by_bin), len(by_slice))
            self.assertTrue((by_bin["gene"].to_numpy() == by_slice.to_numpy()).all())
        # With a VCF
        varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        seg_baf = list(varr.iter_ranges_of(segarr, "alt_freq"))
        self.assertEqual(len(seg_baf), len(segarr))
        cna_baf = list(varr.iter_ranges_of(cnarr, "alt_freq"))
        self.assertEqual(len(cna_baf), len(cnarr))
        # Edge cases
        mtarr = tabio.read("formats/empty")
        self.assertEqual(0, len(list(segarr.iter_ranges_of(mtarr, "start"))))
        self.assertEqual(88, len(list(mtarr.iter_ranges_of(segarr, "end"))))

    def test_ranges_of_non_contiguous_chromosomes(self):
        """iter_ranges_of and by_ranges yield in the other array's row order.

        Both are consumed positionally against the other array's rows -- see
        ``segmetrics.do_segmetrics`` and ``CopyNumArray.residuals`` for
        ``iter_ranges_of``, and the ``variants.by_ranges(segarr)``
        re-segmentation in ``cnvlib.segmentation`` for ``by_ranges`` -- so a
        chromosome-grouped order would misassign across chromosomes.
        """
        src, dest = self._non_contiguous_pair()
        self.assertEqual(
            [s.tolist() for s in src.iter_ranges_of(dest, "gene")],
            [["AAA"], ["CCC"], [], []],
        )
        self.assertEqual(
            [
                (row.chromosome, row.start, sub["gene"].tolist())
                for row, sub in src.by_ranges(dest)
            ],
            [
                ("chr2", 0, ["AAA"]),
                ("chr10", 0, ["CCC"]),
                ("chr2", 2000, []),
                ("chr10", 2000, []),
            ],
        )

    def test_ranges_resize(self):
        baits_fname = "formats/nv2_baits.interval_list"
        chrom_sizes = {
            "chr1": 249250621,
            "chr2": 243199373,
            "chr3": 198022430,
            "chr4": 191154276,
            "chr5": 180915260,
            "chr6": 171115067,
            "chr7": 159138663,
            "chr8": 146364022,
            "chr9": 141213431,
            "chr10": 135534747,
            "chr11": 135006516,
            "chr12": 133851895,
            "chr13": 115169878,
            "chr14": 107349540,
            "chr15": 102531392,
            "chr16": 90354753,
            "chr17": 81195210,
            "chr18": 78077248,
            "chr19": 59128983,
            "chr20": 63025520,
            "chr21": 48129895,
            "chr22": 51304566,
            "chrX": 155270560,
            "chrY": 59373566,
        }
        bins = tabio.read(baits_fname, "interval")
        for chrom, arr in bins.resize_ranges(1e7, chrom_sizes).by_chromosome():
            self.assertLessEqual(0, arr.start.min())
            self.assertLessEqual(arr.end.max(), chrom_sizes[chrom])

    def test_total_range_size(self):
        """Test total region coverage calculation."""
        for fname, area in (
            ("formats/empty", 0),
            ("formats/my-targets.bed", 103),
            ("formats/dac-my.bed", 148),
            ("formats/example.gff", 7951),
            ("formats/refflat-mini.txt", 719715),
        ):
            regions = tabio.read_auto(fname)
            self.assertEqual(regions.total_range_size(), area)

    def test_shuffle_sort(self):
        """Test shuffling and re-sorting the data array."""
        # NB: GenomicArray.shuffle() seeds its own Generator, so this is
        # deterministic -- the "order changed" assertion below can't flake.
        orig_cvg = tuple(self.ex_cnr["log2"][:10])
        self.assertEqual(tuple(self.ex_cnr["log2"].to_numpy()[:10]), orig_cvg)
        self.ex_cnr.shuffle()
        self.assertNotEqual(tuple(self.ex_cnr["log2"].to_numpy()[:10]), orig_cvg)
        self.ex_cnr.sort()
        self.assertEqual(tuple(self.ex_cnr["log2"].to_numpy()[:10]), orig_cvg)


class IntervalTests(unittest.TestCase):
    """Interval arithmetic tests."""

    combiner: typing.ClassVar = {"gene": lambda a: "".join(a)}

    # Simple: nested, overlapping, & non-overlapping intervals
    # =A=========================
    #   =B=======   =D===   =E======
    #      =C=
    # 1 3  5 6  8   11 15   19 20 23 <- coordinates
    region_coords_1 = (
        (1, 20, "A"),
        (3, 8, "B"),
        (5, 6, "C"),
        (11, 15, "D"),
        (19, 23, "E"),
    )

    # Semi-realistic: overlapping gene models
    # =A=============================
    #   =B==  =C==     =E==  =G==
    #               =D=========================
    #                  =F==  =H==       =I==
    # 3 5  8  11 14 17 19 22 25 28  32  36 39 42
    region_coords_2 = (
        (3, 32, "A"),
        (5, 8, "B"),
        (11, 14, "C"),
        (17, 42, "D"),
        (19, 22, "E"),
        (19, 22, "F"),
        (25, 28, "G"),
        (25, 28, "H"),
        (36, 39, "I"),
    )

    @staticmethod
    def _from_intervals(coords):
        garr = GA(
            pd.DataFrame(list(coords), columns=["start", "end", "gene"]).assign(
                chromosome="chr0"
            )
        )
        garr.sort_columns()
        return garr

    def _compare_regions(self, result, expect):
        self.assertEqual(
            expect.data.shape,
            result.data.shape,
            "\n".join(["Got:", str(result.data), "Expected:", str(expect.data)]),
        )
        for col in expect.data.columns:
            self.assertTrue(
                (expect[col].to_numpy() == result[col].to_numpy()).all(),
                f"Col '{col}' differs:\nExpect:\n{expect.data}\nGot:\n{result.data}",
            )

    def setUp(self):
        self.regions_1 = self._from_intervals(self.region_coords_1)
        self.regions_2 = self._from_intervals(self.region_coords_2)

    def test_flatten(self):
        flat_coords_1 = [
            (1, 3, "A"),
            (3, 5, "AB"),
            (5, 6, "ABC"),
            (6, 8, "AB"),
            (8, 11, "A"),
            (11, 15, "AD"),
            (15, 19, "A"),
            (19, 20, "AE"),
            (20, 23, "E"),
        ]
        flat_coords_2 = [
            (3, 5, "A"),
            (5, 8, "AB"),
            (8, 11, "A"),
            (11, 14, "AC"),
            (14, 17, "A"),
            (17, 19, "AD"),
            (19, 22, "ADEF"),
            (22, 25, "AD"),
            (25, 28, "ADGH"),
            (28, 32, "AD"),
            (32, 36, "D"),
            (36, 39, "DI"),
            (39, 42, "D"),
        ]
        for regions, flat_coords in [
            (self.regions_1, flat_coords_1),
            (self.regions_2, flat_coords_2),
        ]:
            result = regions.flatten(combine=self.combiner)
            expect = self._from_intervals(flat_coords)
            self._compare_regions(result, expect)

    def test_flatten_fields_follow_the_covering_rows(self):
        """Each sub-interval takes its fields from the rows that cover it.

        Columns with no combiner -- log2, depth, gc and the rest of a .cnr's
        payload -- were copied from the first row of the whole overlap cluster,
        so a piece belonging solely to a later region reported an earlier
        region's coverage. Bookended rows were caught by this too: they
        clustered together, and one overlap anywhere in a table sends every
        cluster through the splitting path.

        chr1 and chr2 each hold an overlapping pair plus an isolated region,
        after and before the pair respectively, so the pieces have to be
        interleaved back among the rows that pass through unsplit. chr3 is a
        bookended run, which has no breakpoint to split at, ending in a
        zero-width region that only a clustering rule coarser than "genuinely
        overlapping" would swallow.
        """
        regions = GA(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 3 + ["chr2"] * 3 + ["chr3"] * 4,
                    "start": [0, 50, 400, 0, 300, 500, 0, 100, 200, 300],
                    "end": [100, 200, 500, 100, 600, 700, 100, 200, 300, 300],
                    "gene": [*"ABC", *"PQR", *"XYZW"],
                    "log2": [1.0, -2.0, 9.0, 3.0, 4.0, 5.0, 7.0, 7.5, 8.0, 8.5],
                    # Not a Python identifier: a tabular input's header is the
                    # user's own, and `itertuples` renames what it cannot use
                    "GC content": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                }
            )
        )
        flat = regions.flatten()
        self.assertEqual(
            list(
                zip(
                    flat.chromosome,
                    flat.start,
                    flat.end,
                    flat["gene"],
                    flat["log2"],
                    flat["GC content"],
                    strict=True,
                )
            ),
            [
                ("chr1", 0, 50, "A", 1.0, 0.1),
                ("chr1", 50, 100, "A,B", 1.0, 0.1),
                ("chr1", 100, 200, "B", -2.0, 0.2),
                ("chr1", 400, 500, "C", 9.0, 0.3),
                ("chr2", 0, 100, "P", 3.0, 0.4),
                ("chr2", 300, 500, "Q", 4.0, 0.5),
                ("chr2", 500, 600, "Q,R", 4.0, 0.5),
                ("chr2", 600, 700, "R", 5.0, 0.6),
                ("chr3", 0, 100, "X", 7.0, 0.7),
                ("chr3", 100, 200, "Y", 7.5, 0.8),
                ("chr3", 200, 300, "Z", 8.0, 0.9),
                ("chr3", 300, 300, "W", 8.5, 1.0),
            ],
        )
        # A table with no overlaps at all passes through unchanged
        chr3_only = regions.data[regions.data.chromosome == "chr3"]
        disjoint = GA(chr3_only.reset_index(drop=True))
        self.assertTrue(disjoint.flatten().data.equals(disjoint.data))

    def test_flatten_apportions_split_columns(self):
        """'weight' and 'probes' are divided among the pieces of a region.

        They quantify a region rather than describe it, so splitting a region
        has to divide them: replicating each row's whole value into every piece
        inflated the totals, the more so the more finely a region was split.
        """
        regions = GA(
            pd.DataFrame(
                {
                    "chromosome": "chr0",
                    "start": [0, 300],
                    "end": [400, 600],
                    "gene": ["A", "B"],
                    "weight": [4.0, 3.0],
                    "probes": [7, 5],
                }
            )
        )
        flat = regions.flatten()
        # A gives 3/4 of itself to chr0:0-300 and 1/4 to chr0:300-400, where B
        # adds a third of itself -- shares of length, not equal shares
        self.assertEqual(
            list(zip(flat.start, flat.end, flat["weight"], strict=True)),
            [(0, 300, 3.0), (300, 400, 2.0), (400, 600, 2.0)],
        )
        self.assertAlmostEqual(flat["weight"].sum(), regions["weight"].sum())
        # A count has no fractional part to give: the shares 5.25 / 3.42 / 3.33
        # are rounded, and the column keeps its integer type, or `export vcf`
        # would read every piece as a segment with a corrupt probe count
        self.assertEqual(list(flat["probes"]), [5, 3, 3])
        self.assertTrue(is_integer_dtype(flat["probes"]))
        # Opting out combines them like any other column, so each row's value
        # reappears whole in every piece it was split into
        plain = regions.flatten(split_columns=())
        self.assertEqual(list(plain["weight"]), [4.0, 7.0, 3.0])
        self.assertEqual(list(plain["probes"]), [7, 12, 5])

        # Which column rounds is decided by what it means, not by the dtype it
        # happens to arrive with. A file whose weights are all whole is written
        # as bare integers and reads back as int64, but a weight is a
        # proportion, so rounding it would zero every share below a half.
        whole_weights = GA(regions.data.assign(weight=np.array([4, 3], dtype=np.int64)))
        self.assertEqual(list(whole_weights.flatten()["weight"]), [3.0, 2.0, 2.0])
        # A blank field in one row makes 'probes' float-typed, and it is still
        # a count: the shares round, and stay whole for `export vcf`
        float_probes = GA(regions.data.assign(probes=[7.0, 5.0]))
        self.assertEqual(list(float_probes.flatten()["probes"]), [5.0, 3.0, 3.0])

    def test_merge(self):
        merged_coords_1 = [(1, 23, "ABCDE")]
        merged_coords_2 = [(3, 42, "ABCDEFGHI")]
        for regions, merged_coords in [
            (self.regions_1, merged_coords_1),
            (self.regions_2, merged_coords_2),
        ]:
            result = regions.merge(combine=self.combiner)
            expect = self._from_intervals(merged_coords)
            self._compare_regions(result, expect)

    def test_merge_overlapping_only(self):
        """``merge(bp=1)`` merges overlapping rows but not bookended ones.

        The bookended distinction matters to the 'target' command: consecutive
        bait tiles must survive as separate bins, or a capture kit's targets
        get coarsened (#567). It is also what ``skg_convert --merge`` has
        always documented -- "the number of overlapping bases ... to trigger a
        merge" -- while the default ``bp=0`` follows ``bedtools merge`` in
        merging bookended rows too.
        """
        regions = GA(
            pd.DataFrame(
                {
                    "chromosome": "chr0",
                    "start": [100, 150, 400, 500, 700, 720],
                    "end": [200, 300, 500, 600, 750, 800],
                    "gene": ["A", "B", "C", "D", "E", "F"],
                }
            )
        )
        merged = regions.merge(bp=1)
        self.assertEqual(
            list(zip(merged.start, merged.end, merged["gene"], strict=True)),
            [(100, 300, "A,B"), (400, 500, "C"), (500, 600, "D"), (700, 800, "E,F")],
        )
        # The default also merges the bookended pair chr0:400-500 / 500-600
        default_merged = regions.merge()
        self.assertEqual(
            list(zip(default_merged.start, default_merged.end, strict=True)),
            [(100, 300), (400, 600), (700, 800)],
        )
        # A table with no overlaps at all passes through unchanged
        disjoint = GA(regions.data.iloc[[2, 3]].reset_index(drop=True))
        self.assertTrue(disjoint.merge(bp=1).data.equals(disjoint.data))
        # Overlaps are found even when a chromosome's rows are not contiguous,
        # where the sorted-input shortcut for spotting them does not apply
        interleaved = GA(
            pd.DataFrame(
                {
                    "chromosome": ["chr1", "chr2", "chr1"],
                    "start": [100, 500, 150],
                    "end": [200, 600, 250],
                    "gene": ["A", "Z", "B"],
                }
            )
        )
        self.assertEqual(
            list(
                zip(
                    interleaved.merge(bp=1).start,
                    interleaved.merge(bp=1).end,
                    strict=True,
                )
            ),
            [(100, 250), (500, 600)],
        )

    def test_columns_named_like_bioframe_bookkeeping(self):
        """A column named cluster/cluster_start/cluster_end is the user's own.

        ``bioframe.cluster``, which merge and flatten group rows with, used to
        be asked for its bookkeeping alongside the input, and it names that
        bookkeeping cluster/cluster_start/cluster_end. An input column of the
        same name was duplicated rather than replaced. Under the name
        ``cluster`` the duplicate label made the cluster ids a whole DataFrame,
        so rows were grouped by bioframe's id and the user's value together and
        genuinely overlapping rows came through unmerged and unsplit -- wrong
        coordinates, not merely a lost column. Under ``cluster_start`` or
        ``cluster_end``, merge read the same DataFrame where it expected one
        coordinate and raised ``TypeError``, while flatten dropped the column
        silently. ``tabio.read_tab`` accepts arbitrary extra columns, so a
        header is user data, and 'cluster' is an unremarkable name for a bait
        or sample label.

        The damage was data-dependent: a table with nothing to cluster passes
        through verbatim, so the same pipeline mangled or preserved the same
        column according to whether that particular file held an overlap.
        """
        expect = {
            # Two rows overlapping: merged into one, split into three
            ("overlapping", "merge"): [10],
            ("overlapping", "flatten"): [10, 10, 20],
            # Nothing to cluster: both pass the rows through untouched
            ("disjoint", "merge"): [10, 20],
            ("disjoint", "flatten"): [10, 20],
        }
        bounds = {
            "overlapping": ([0, 50], [100, 200]),
            "disjoint": ([0, 500], [100, 600]),
        }
        for name in ("cluster", "cluster_start", "cluster_end"):
            for layout, (starts, ends) in bounds.items():
                regions = GA(
                    pd.DataFrame(
                        {
                            "chromosome": "chr0",
                            "start": starts,
                            "end": ends,
                            "gene": ["A", "B"],
                            name: [10, 20],
                        }
                    )
                )
                for op in ("merge", "flatten"):
                    with self.subTest(column=name, table=layout, op=op):
                        result = getattr(regions, op)()
                        self.assertEqual(list(result[name]), expect[layout, op])

        # Every route a column can take treats these names like any other: all
        # three at once, combined by a caller's own function...
        overlapping = GA(
            pd.DataFrame(
                {
                    "chromosome": "chr0",
                    "start": [0, 50],
                    "end": [100, 200],
                    "cluster": [10, 20],
                    "cluster_start": [1, 2],
                    "cluster_end": [3, 4],
                }
            )
        )
        merged = overlapping.merge()
        self.assertEqual(
            [merged[col].iloc[0] for col in overlapping.data.columns],
            ["chr0", 0, 200, 10, 1, 3],
        )
        self.assertEqual(
            list(overlapping.merge(combine={"cluster": sum})["cluster"]), [30]
        )
        # ...and apportioned by length as a quantity spread over the region,
        # where the column's name must not decide whether the division happens
        shares = []
        for name in ("extra", "cluster_start"):
            regions = GA(
                pd.DataFrame(
                    {
                        "chromosome": "chr0",
                        "start": [0, 50],
                        "end": [100, 200],
                        name: [1.0, 2.0],
                    }
                )
            )
            pieces = regions.flatten(combine={name: sum}, split_columns=(name,))
            shares.append(list(pieces[name]))
        self.assertEqual(shares[0], shares[1])
        self.assertAlmostEqual(sum(shares[1]), 3.0)

    def test_nan_gene_names(self):
        """merge/flatten/subdivide tolerate NaN-float gene names (issue #850).

        A bait BED whose 'gene' column holds NaN floats (e.g. a ragged name
        column) once crashed ``target --annotate --split`` two ways:
        ``TypeError: expected str instance, numpy.float64 found`` inside
        ``join_strings`` (``",".join`` over NaN), and an ``AttributeError`` in
        the old ``groupby.apply`` merge path on newer pandas.  Both are fixed
        upstream -- by the bioframe rewrite (#982) and the ``join_strings`` NaN
        filter (#900) -- but the skgenome merge/subdivide path that the 'target'
        command exercises had no regression guard.  Use the default combiner
        (``join_strings``), not the ``"".join`` test combiner, so the real
        production code path is tested.

        Regions 0+1 overlap (NaN skipped -> "BRCA1"); regions 2+3 overlap and
        are all-NaN (-> "-"); region 4 is isolated and all-NaN.  The isolated
        region exercises the singleton/early-return paths, where merge & co.
        emit rows verbatim -- those NaN floats are normalized to "-" so no
        unjoined NaN survives in a gene/accession column.
        """
        regions = GA(
            pd.DataFrame(
                {
                    "chromosome": "chr0",
                    "start": [100, 150, 700, 800, 5000],
                    "end": [400, 600, 900, 1000, 8000],
                    "gene": [np.nan, "BRCA1", np.nan, np.nan, np.nan],
                }
            )
        )
        merged = regions.merge()
        self.assertEqual(list(merged["gene"]), ["BRCA1", "-", "-"])

        # Single-row fast paths (merge/flatten/squash all short-circuit on a
        # lone region) must normalize too, not pass the NaN float through.
        isolated = GA(regions.data.iloc[[4]])
        self.assertEqual(list(isolated.merge()["gene"]), ["-"])
        self.assertEqual(list(isolated.flatten()["gene"]), ["-"])
        self.assertEqual(list(isolated.squash()["gene"]), ["-"])

        # flatten() splits at breakpoints and re-joins per sub-interval.
        # subdivide() (the 'target --split' path) merges then splits big bins.
        for result in (regions.flatten(), regions.subdivide(100, 0)):
            self.assertGreater(len(result), 0)
            for gene in result["gene"]:
                self.assertIsInstance(gene, str)
                self.assertNotIn("nan", gene.lower())

        # A name column with no strings at all is typed float64 (all-NaN) or
        # int64 (numeric bait names, which Picard interval_list files carry),
        # and merging must widen it to hold the "-" that join_strings returns
        # rather than coerce the value into the incumbent dtype
        for names in ([np.nan] * 3, [1, 2, 3]):
            numeric = GA(
                pd.DataFrame(
                    {
                        "chromosome": "chr0",
                        "start": [100, 150, 700],
                        "end": [400, 600, 900],
                        "gene": names,
                    }
                )
            )
            self.assertEqual(list(numeric.merge()["gene"]), ["-", "-"])

    def test_join_strings_deduplicates_across_separator(self):
        """Joining already-joined gene labels repeats no name (#567).

        Combining operations chain -- ``target`` merges overlapping baits and
        then re-labels the bins split out of them from those already-joined
        labels, and ``segmentation`` joins the gene labels of the bins a
        segment spans -- so the inputs to ``join_strings`` are often themselves
        joined labels. Treating each input as opaque yielded labels such as
        "GENEA,GENEA,GENEB,GENEB".
        """
        self.assertEqual(join_strings(["GENEA", "GENEA,GENEB"]), "GENEA,GENEB")
        self.assertEqual(
            join_strings(["GENEA,GENEB", "GENEB,GENEC"]), "GENEA,GENEB,GENEC"
        )
        # Order of first appearance is preserved, and joining is idempotent
        joined = join_strings(["B,A", "C"])
        self.assertEqual(joined, "B,A,C")
        self.assertEqual(join_strings([joined]), joined)
        # Ignored names are dropped from within joined labels too
        self.assertEqual(join_strings(["-,GENEA"], ignore=("-",)), "GENEA")
        self.assertEqual(join_strings(["-", np.nan]), "-")

    def test_intersect(self):
        selections1 = self._from_intervals(
            [
                (1, 8, ""),
                (4, 10, ""),
                (8, 19, ""),
                (11, 20, ""),
                (21, 22, ""),
            ]
        )
        expectations1 = {
            "outer": (
                # 1-8
                [
                    (1, 20, "A"),
                    (3, 8, "B"),
                    (5, 6, "C"),
                ],
                # 4-10
                [
                    (1, 20, "A"),
                    (3, 8, "B"),
                    (5, 6, "C"),
                ],
                # 8-19
                [(1, 20, "A"), (11, 15, "D")],
                # 11-20
                [(1, 20, "A"), (11, 15, "D"), (19, 23, "E")],
                # 21-22
                [(19, 23, "E")],
            ),
            "trim": (
                # 1-8
                [(1, 8, "A"), (3, 8, "B"), (5, 6, "C")],
                # 4-10
                [(4, 10, "A"), (4, 8, "B"), (5, 6, "C")],
                # 8-19
                [(8, 19, "A"), (11, 15, "D")],
                # 11-20
                [(11, 20, "A"), (11, 15, "D"), (19, 20, "E")],
                # 21-22
                [(21, 22, "E")],
            ),
            "inner": (
                # 1-8
                [(3, 8, "B"), (5, 6, "C")],
                # 4-10
                [(5, 6, "C")],
                # 8-19
                [(11, 15, "D")],
                # 11-20
                [(11, 15, "D")],
                # 21-22
                [],
            ),
        }

        selections2 = self._from_intervals(
            [
                (0, 1, ""),
                (5, 14, ""),
                (16, 45, ""),
                (18, 37, ""),
                (19, 25, ""),
                (29, 31, ""),
                (34, 39, ""),
            ]
        )
        expectations2 = {
            "outer": (
                # 0-1
                [],
                # 5-14
                [(3, 32, "A"), (5, 8, "B"), (11, 14, "C")],
                # 16-45
                [
                    (3, 32, "A"),
                    (17, 42, "D"),
                    (19, 22, "E"),
                    (19, 22, "F"),
                    (25, 28, "G"),
                    (25, 28, "H"),
                    (36, 39, "I"),
                ],
                # 18-37
                [
                    (3, 32, "A"),
                    (17, 42, "D"),
                    (19, 22, "E"),
                    (19, 22, "F"),
                    (25, 28, "G"),
                    (25, 28, "H"),
                    (36, 39, "I"),
                ],
                # 19-25
                [(3, 32, "A"), (17, 42, "D"), (19, 22, "E"), (19, 22, "F")],
                # 29-31
                [(3, 32, "A"), (17, 42, "D")],
                # 34-39
                [(17, 42, "D"), (36, 39, "I")],
            ),
            "trim": (
                # 0-1
                [],
                # 5-14
                [(5, 14, "A"), (5, 8, "B"), (11, 14, "C")],
                # 16-45
                [
                    (16, 32, "A"),
                    (17, 42, "D"),
                    (19, 22, "E"),
                    (19, 22, "F"),
                    (25, 28, "G"),
                    (25, 28, "H"),
                    (36, 39, "I"),
                ],
                # 18-37
                [
                    (18, 32, "A"),
                    (18, 37, "D"),
                    (19, 22, "E"),
                    (19, 22, "F"),
                    (25, 28, "G"),
                    (25, 28, "H"),
                    (36, 37, "I"),
                ],
                # 19-25
                [(19, 25, "A"), (19, 25, "D"), (19, 22, "E"), (19, 22, "F")],
                # 29-31
                [(29, 31, "A"), (29, 31, "D")],
                # 34-39
                [(34, 39, "D"), (36, 39, "I")],
            ),
            "inner": (
                # 0-1
                [],
                # 5-14
                [(5, 8, "B"), (11, 14, "C")],
                # 16-45
                [
                    (17, 42, "D"),
                    (19, 22, "E"),
                    (19, 22, "F"),
                    (25, 28, "G"),
                    (25, 28, "H"),
                    (36, 39, "I"),
                ],
                # 18-37
                [(19, 22, "E"), (19, 22, "F"), (25, 28, "G"), (25, 28, "H")],
                # 19-25
                [(19, 22, "E"), (19, 22, "F")],
                # 29-31
                [],
                # 34-39
                [(36, 39, "I")],
            ),
        }

        for regions, selections, expectations in (
            (self.regions_1, selections1, expectations1),
            (self.regions_2, selections2, expectations2),
        ):
            for mode in ("outer", "trim", "inner"):
                # Iterative intersection
                grouped_results = regions.by_ranges(selections, mode=mode)
                for (_coord, result), expect in zip(
                    grouped_results, expectations[mode], strict=False
                ):
                    self._compare_regions(result, self._from_intervals(expect))
                # Single-object intersect
                result = regions.intersection(selections, mode=mode)
                expect = self._from_intervals(
                    functools.reduce(operator.iadd, expectations[mode], [])
                )
                self._compare_regions(result, expect)

    def test_intersect_columns_named_like_bioframe_bookkeeping(self):
        """A column named index/index_/start_/end_ is the user's own.

        ``intersection`` pairs rows with ``bioframe.overlap``, which returns
        its own ``index`` and ``index_`` columns and, when asked for the inputs
        as well, the second one's coordinates suffixed apart as ``start_`` and
        ``end_``. An input column of one of those four names duplicated a label
        the code then read: ``index`` and ``index_`` raised ``ValueError``
        naming the column as not unique, and ``start_``/``end_`` raised
        ``ValueError: Operands are not aligned`` in the 'inner' mode that is
        the only one to compare them.

        Index labels no longer influence the result at all, since the rows are
        addressed by position: a duplicated label resolved through ``.loc`` to
        every row that shared it, and ``other``'s labels, rather than its row
        order, decided which of its regions was reported first.

        chr0:650-750 overlaps the first selection without being contained in
        it, so 'inner' keeps one pairing fewer than 'outer', and 'trim' clips
        it to the selection's end. chr0:5000-5100 overlaps nothing.
        """
        regions = GA(
            pd.DataFrame(
                {
                    "chromosome": "chr0",
                    "start": [0, 500, 650, 900, 5000],
                    "end": [100, 600, 750, 1000, 5100],
                    "gene": [*"ABECZ"],
                }
            )
        )
        selections = GA(
            pd.DataFrame(
                {
                    "chromosome": "chr0",
                    "start": [0, 450],
                    "end": [700, 1000],
                    "gene": ["X", "Y"],
                }
            )
        )
        expect = {
            "outer": [10, 20, 30, 20, 30, 40],
            "inner": [10, 20, 20, 30, 40],
            "trim": [10, 20, 30, 20, 30, 40],
        }
        for name in ("index", "index_", "start_", "end_"):
            labelled = GA(regions.data.assign(**{name: [10, 20, 30, 40, 50]}))
            for mode in ("outer", "inner", "trim"):
                with self.subTest(column=name, mode=mode):
                    result = labelled.intersection(selections, mode=mode)
                    self.assertEqual(list(result[name]), expect[mode])

        # A duplicated label -- as a table assembled from several sources
        # carries -- once selected every row that shared it, so the pairings of
        # A, B, E and C dragged in the region that overlaps nothing
        duped = GA(regions.data.set_axis([0, 1, 0, 1, 0], axis=0))
        want = regions.intersection(selections)
        got = duped.intersection(selections)
        self.assertEqual(
            list(zip(got.start, got["gene"], strict=True)),
            list(zip(want.start, want["gene"], strict=True)),
        )
        self.assertEqual(len(got), 6)
        # `other`'s rows are reported in its own order, not its labels'
        relabelled = GA(selections.data.set_axis([5, 1], axis=0))
        self.assertEqual(
            list(regions.intersection(relabelled)["gene"]), list(want["gene"])
        )

    def test_intersect_nothing_to_select(self):
        """Selecting nothing gives an empty array in every mode, not an error.

        'trim' concatenates the pieces ``by_ranges`` yields, and ``pd.concat``
        raises on an empty list, so a selection that nothing overlaps -- or an
        empty array on either side -- ended in ``ValueError: No objects to
        concatenate``, while the other two modes returned an empty array. The
        guard for the empty inputs sat below the 'trim' branch, so it never
        applied to it.
        """
        region = GA(
            pd.DataFrame(
                {"chromosome": ["chr1"], "start": [10], "end": [20], "gene": ["A"]}
            )
        )
        empty = GA(region.data.iloc[:0])
        selections = {
            "empty": empty,
            "disjoint": GA(region.data.assign(start=100, end=200)),
            "other chromosome": GA(region.data.assign(chromosome="chrZ")),
        }
        for label, selection in selections.items():
            for mode in ("outer", "inner", "trim"):
                with self.subTest(selection=label, mode=mode):
                    result = region.intersection(selection, mode=mode)
                    self.assertEqual(len(result), 0)
                    self.assertEqual(
                        list(result.data.columns), list(region.data.columns)
                    )
                with self.subTest(empty_self=True, mode=mode):
                    self.assertEqual(len(empty.intersection(region, mode=mode)), 0)

    def test_search_requires_ascending_rows(self):
        """Searching rows out of coordinate order raises instead of misleading.

        ``searchsorted`` places the boundaries of every selection here, and it
        reads its column as ascending: given rows in any other order it returns
        a position unrelated to the query. The three rows below, shuffled,
        answered "which bins overlap 0-700?" with the bin at 900-1000 as well
        -- and in 'trim' mode clipped it to 900-700, an interval whose end
        precedes its start. 'inner' searches a different column and happened to
        survive, which is why every mode is checked.
        """
        coords = [(0, 100, "A"), (500, 600, "B"), (900, 1000, "C")]
        ascending = self._from_intervals(coords)
        shuffled = self._from_intervals([coords[i] for i in (2, 0, 1)])
        selection = self._from_intervals([(0, 700, "X")])
        for mode in ("outer", "inner", "trim"):
            with self.subTest(mode=mode):
                self.assertEqual(
                    [
                        list(rows["gene"])
                        for _bin, rows in ascending.by_ranges(selection, mode=mode)
                    ],
                    [["A", "B"]],
                )
                with self.assertRaises(ValueError):
                    list(shuffled.by_ranges(selection, mode=mode))
                self.assertEqual(
                    list(ascending.in_range("chr0", 0, 700, mode=mode)["gene"]),
                    ["A", "B"],
                )
                with self.assertRaises(ValueError):
                    shuffled.in_range("chr0", 0, 700, mode=mode)
        self.assertEqual(list(ascending.into_ranges(selection, "gene", "-")), ["A,B"])
        with self.assertRaises(ValueError):
            shuffled.into_ranges(selection, "gene", "-")
        # A missing position leaves the column unordered without any row being
        # behind the one before it, so the report has nothing to point at.
        frayed = self._from_intervals(coords)
        frayed.data = frayed.data.astype({"start": float})
        frayed.data.loc[1, "start"] = np.nan
        with self.assertRaises(ValueError) as caught:
            frayed.in_range("chr0", 0, 700)
        self.assertIn("missing", str(caught.exception))

    def test_search_requires_one_chromosome(self):
        """A whole-array search names its chromosome, or is refused.

        ``in_range(None, ...)`` is documented for an array holding one
        chromosome. Given several it searched across them all, since a second
        chromosome's coordinates restart below the first one's: asking chr1 for
        0-700 also returned chr2's rows.
        """
        two_chroms = GA(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 2 + ["chr2"] * 2,
                    "start": [0, 900] * 2,
                    "end": [100, 1000] * 2,
                    "gene": list("ABCD"),
                }
            )
        )
        self.assertEqual(list(two_chroms.in_range("chr1", 0, 700)["gene"]), ["A"])
        with self.assertRaises(ValueError) as caught:
            two_chroms.in_range(None, 0, 700)
        self.assertIn("chr1, chr2", str(caught.exception))

    def test_search_one_sided_over_nested_bins(self):
        """A half-open range selects what the closed range covering it selects.

        A bin nested inside another leaves `end` out of ascending order, so the
        slicing implementation cannot search it; the mask implementation is
        there for exactly that. But the choice between them also asked whether
        the query had both bounds, and a half-open range has one -- so
        ``in_range(chrom, 20, None)`` searched a column it could not search,
        and lost the bins covering position 20 that the equivalent closed
        range, on the same rows, returned.
        """
        nested = self.regions_2
        self.assertFalse(nested.end.is_monotonic_increasing)
        beyond = int(nested.end.max()) + 1
        for mode in ("outer", "inner", "trim"):
            with self.subTest(mode=mode):
                self.assertEqual(
                    list(nested.in_range("chr0", 20, None, mode=mode)["gene"]),
                    list(nested.in_range("chr0", 20, beyond, mode=mode)["gene"]),
                )
                self.assertEqual(
                    list(nested.in_range("chr0", None, 20, mode=mode)["gene"]),
                    list(nested.in_range("chr0", 0, 20, mode=mode)["gene"]),
                )
        # Spelled out, so the agreement above cannot be two identical mistakes:
        # every bin reaching past 20, and every bin ending by 20.
        self.assertEqual(
            list(nested.in_range("chr0", 20, None)["gene"]),
            ["A", "D", "E", "F", "G", "H", "I"],
        )
        self.assertEqual(
            list(nested.in_range("chr0", None, 20, mode="inner")["gene"]), ["B", "C"]
        )

    def test_search_no_ranges_selects_nothing(self):
        """An empty list of ranges selects nothing, rather than everything."""
        empty = self.regions_2.in_ranges("chr0", [], [])
        self.assertEqual(len(empty), 0)
        self.assertEqual(list(empty.data.columns), list(self.regions_2.data.columns))
        self.assertEqual(len(self.regions_2.in_ranges("chr0", [], None)), 0)

    def test_subtract(self):
        # Test cases:
        #  | access: ====   ====   ====    ==========
        #  | target: ====  ====== ===   = ==  ==   ===
        #  | expect:                 ==     ==  ===
        #  | invert:       =    = =     = =          =
        #  |         1   5 78   2345 7 901234 6 8  1 34
        #  |        0         1         2         3

        access = self._from_intervals(
            [
                (1, 5, ""),
                (8, 12, ""),
                (15, 19, ""),
                (23, 33, ""),
            ]
        )
        target = self._from_intervals(
            [
                (1, 5, ""),
                (7, 13, ""),
                (14, 17, ""),
                (20, 21, ""),
                (22, 24, ""),
                (26, 28, ""),
                (31, 34, ""),
            ]
        )
        expect = self._from_intervals(
            [
                (17, 19, ""),
                (24, 26, ""),
                (28, 31, ""),
            ]
        )
        invert = self._from_intervals(
            [
                (7, 8, ""),
                (12, 13, ""),
                (14, 15, ""),
                (20, 21, ""),
                (22, 23, ""),
                (33, 34, ""),
            ]
        )

        result = access.subtract(target)
        self._compare_regions(result, expect)
        iresult = target.subtract(access)
        self._compare_regions(iresult, invert)

    def test_subtract_multichrom(self):
        """Regression for #471: target regions are subtracted on every chromosome.

        The old custom interval arithmetic indexed a chromosome-filtered table
        with a non-reset index, so on chromosomes after the first the
        subtraction silently dropped and target regions leaked back into the
        antitarget. The bioframe-based ``subtract`` must remove the target on
        the later chromosome too, not only the first.
        """

        def _multichrom(rows):
            garr = GA(
                pd.DataFrame(list(rows), columns=["chromosome", "start", "end", "gene"])
            )
            garr.sort_columns()
            return garr

        access = _multichrom([("chr1", 0, 100, ""), ("chr2", 0, 100, "")])
        target = _multichrom([("chr1", 40, 60, ""), ("chr2", 40, 60, "")])
        expect = _multichrom(
            [
                ("chr1", 0, 40, ""),
                ("chr1", 60, 100, ""),
                ("chr2", 0, 40, ""),
                ("chr2", 60, 100, ""),
            ]
        )
        result = access.subtract(target)
        self._compare_regions(result, expect)
        # The invariant the bug violated: no leftover region overlaps a target,
        # including on chr2 (the chromosome that previously leaked).
        for chrom, start, end in zip(
            result["chromosome"], result["start"], result["end"], strict=True
        ):
            for t_chrom, t_start, t_end in target.coords():
                overlaps = chrom == t_chrom and start < t_end and end > t_start
                self.assertFalse(
                    overlaps, f"{chrom}:{start}-{end} overlaps target on {chrom}"
                )

    def test_cut(self):
        # Simple: one interval split by two boundary regions
        regions = self._from_intervals([(1, 20, "A")])
        cuts = self._from_intervals([(5, 10, "X"), (15, 25, "Y")])
        # Breakpoints {5, 10, 15, 25}; 5, 10, 15 are inside (1, 20)
        result = regions.cut(cuts)
        expect = self._from_intervals(
            [
                (1, 5, "A"),
                (5, 10, "A"),
                (10, 15, "A"),
                (15, 20, "A"),
            ]
        )
        self._compare_regions(result, expect)

        # Breakpoint at exact start/end: only interior breakpoint splits
        regions2 = self._from_intervals([(10, 20, "X")])
        cuts2 = self._from_intervals([(10, 15, ""), (20, 25, "")])
        result2 = regions2.cut(cuts2)
        expect2 = self._from_intervals([(10, 15, "X"), (15, 20, "X")])
        self._compare_regions(result2, expect2)

        # No breakpoints inside any interval: unchanged
        regions3 = self._from_intervals([(10, 20, "X")])
        cuts3 = self._from_intervals([(1, 5, ""), (25, 30, "")])
        result3 = regions3.cut(cuts3)
        self._compare_regions(result3, regions3)

        # Empty other: unchanged
        regions4 = self._from_intervals([(1, 20, "A")])
        empty = self._from_intervals([])
        result4 = regions4.cut(empty)
        self._compare_regions(result4, regions4)

        # Multiple intervals, some cut, some not
        regions5 = self._from_intervals(
            [
                (1, 10, "A"),
                (20, 30, "B"),
                (40, 50, "C"),
            ]
        )
        cuts5 = self._from_intervals([(5, 25, "")])
        # Breakpoints {5, 25}: 5 inside A, 25 inside B, neither inside C
        result5 = regions5.cut(cuts5)
        expect5 = self._from_intervals(
            [
                (1, 5, "A"),
                (5, 10, "A"),
                (20, 25, "B"),
                (25, 30, "B"),
                (40, 50, "C"),
            ]
        )
        self._compare_regions(result5, expect5)

    def test_squash(self):
        # All adjacent: combine into one row
        regions = self._from_intervals(
            [
                (1, 5, "A"),
                (5, 10, "B"),
                (10, 15, "C"),
            ]
        )
        result = regions.squash(combine=self.combiner)
        expect = self._from_intervals([(1, 15, "ABC")])
        self._compare_regions(result, expect)

        # Gap prevents combining
        regions2 = self._from_intervals(
            [
                (1, 5, "A"),
                (5, 10, "B"),
                (12, 15, "C"),
                (15, 20, "D"),
            ]
        )
        result2 = regions2.squash(combine=self.combiner)
        expect2 = self._from_intervals(
            [
                (1, 10, "AB"),
                (12, 20, "CD"),
            ]
        )
        self._compare_regions(result2, expect2)

        # Squash by gene: only combine same-gene consecutive runs
        regions3 = self._from_intervals(
            [
                (1, 5, "X"),
                (5, 10, "X"),
                (10, 15, "Y"),
                (15, 20, "Y"),
                (20, 25, "X"),
            ]
        )
        result3 = regions3.squash(by="gene", combine=self.combiner)
        expect3 = self._from_intervals(
            [
                (1, 10, "X"),
                (10, 20, "Y"),
                (20, 25, "X"),
            ]
        )
        self._compare_regions(result3, expect3)

        # Single row: unchanged
        regions4 = self._from_intervals([(1, 100, "X")])
        result4 = regions4.squash(combine=self.combiner)
        self._compare_regions(result4, regions4)

        # Non-adjacent: unchanged
        regions5 = self._from_intervals(
            [
                (1, 5, "A"),
                (10, 15, "B"),
            ]
        )
        result5 = regions5.squash(combine=self.combiner)
        self._compare_regions(result5, regions5)


class OtherTests(unittest.TestCase):
    """Other small modules & functions in this sub-package."""

    def test_chromsort(self):
        labels_hg = ["chr1", "chr2", "chr10", "chr19", "chr20", "chrX", "chrY", "chrM"]
        labels_grc = ["1", "2", "10", "19", "X", "Y", "MT"]
        for labels in (labels_hg, labels_grc):
            shuf = labels[:]
            random.shuffle(shuf)
            resorted = sorted(labels, key=chromsort.sorter_chrom)
            self.assertEqual(resorted, labels)

    def test_detect_big_chroms(self):
        sizes = [1, 20, 30]
        n_big, thresh = chromsort.detect_big_chroms(sizes)
        self.assertEqual(n_big, 2)
        self.assertEqual(thresh, 20)

    def test_rangelabel(self):
        row = rangelabel.from_label("chr1:123-456", keep_gene=False)
        self.assertEqual(tuple(row), ("chr1", 122, 456))
        label = rangelabel.to_label(row)
        self.assertEqual(label, "chr1:123-456")
        # unpack_range
        for region, expect in (
            ["chr1", ("chr1", None, None)],
            ["chr1:100-123", ("chr1", 99, 123)],
            [("chr1", 100, 123), ("chr1", 100, 123)],
            [("chr1", 100, 123, "A"), ("chr1", 100, 123)],
        ):
            result = rangelabel.unpack_range(region)
            self.assertEqual(result, expect)

    def test_rangelabel_ncbi_chrom(self):
        # NCBI-style chromosome names contain a '.' version suffix (#602)
        row = rangelabel.from_label("NC_039902.1:10000000-15000000", keep_gene=False)
        self.assertEqual(tuple(row), ("NC_039902.1", 9999999, 15000000))
        self.assertEqual(rangelabel.to_label(row), "NC_039902.1:10000000-15000000")
        # Open-ended start/end are still parsed with a dotted chromosome name
        self.assertEqual(
            tuple(rangelabel.from_label("NC_039902.1:10000000-", keep_gene=False)),
            ("NC_039902.1", 9999999, None),
        )
        self.assertEqual(
            tuple(rangelabel.from_label("NC_039902.1:-15000000", keep_gene=False)),
            ("NC_039902.1", None, 15000000),
        )
        # unpack_range: whole dotted chromosome (no colon) and a region of it
        for region, expect in (
            ["NC_039902.1", ("NC_039902.1", None, None)],
            ["NC_039902.1:10000000-15000000", ("NC_039902.1", 9999999, 15000000)],
        ):
            self.assertEqual(rangelabel.unpack_range(region), expect)


if __name__ == "__main__":
    unittest.main(verbosity=2)
