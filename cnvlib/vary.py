"""An array of genomic intervals, treated as variant loci."""
from __future__ import absolute_import, division, print_function
from builtins import str
import logging

import numpy as np
import pandas as pd
from skgenome import GenomicArray


class VariantArray(GenomicArray):
    """An array of genomic intervals, treated as variant loci.

    Required columns: chromosome, start, end, ref, alt
    """
    _required_columns = ("chromosome", "start", "end", "ref", "alt")
    _required_dtypes = (str, int, int, str, str)
    # Extra: somatic, zygosity, depth, alt_count, alt_freq

    def __init__(self, data_table, meta_dict=None):
        GenomicArray.__init__(self, data_table, meta_dict)

    def baf_by_ranges(self, ranges, summary_func=np.nanmedian, above_half=None,
                      tumor_boost=False):
        """Aggregate variant (b-allele) frequencies in each given bin.

        Get the average BAF in each of the bins of another genomic array:
        BAFs are mirrored above/below 0.5 (per `above_half`), grouped in each
        bin of `ranges`, and summarized into one value per bin with
        `summary_func` (default median).

        Parameters
        ----------
        ranges : GenomicArray or subclass
            Bins for grouping the variants in `self`.
        above_half : bool
            The same as in `mirrored_baf`.
        tumor_boost : bool
            The same as in `mirrored_baf`.

        Returns
        -------
        float array
            Average b-allele frequency in each range; same length as `ranges`.
            May contain NaN values where no variants overlap a range.
        """
        if 'alt_freq' not in self:
            logging.warning("VCF has no allele frequencies for BAF calculation")
            return pd.Series(np.repeat(np.nan, len(ranges)))

        def summarize(vals):
            return summary_func(_mirrored_baf(vals, above_half))

        if tumor_boost and 'n_alt_freq' in self:
            self = self.copy()
            self['alt_freq'] = self.tumor_boost()
        return self.into_ranges(ranges, 'alt_freq', np.nan, summarize)

    def zygosity_from_freq(self, het_freq=0.0, hom_freq=1.0):
        """Set zygosity (genotype) according to allele frequencies.

        Creates or replaces 'zygosity' column if 'alt_freq' column is present,
        and 'n_zygosity' if 'n_alt_freq' is present.

        Parameters
        ----------
        het_freq : float
            Assign zygosity 0.5 (heterozygous), otherwise 0.0 (i.e. reference
            genotype), to variants with alt allele frequency of at least this
            value.
        hom_freq : float
            Assign zygosity 1.0 (homozygous) to variants with alt allele
            frequency of at least this value.
        """
        assert 0.0 <= het_freq <= hom_freq <= 1.0
        self = self.copy()  # Don't modify the original
        for freq_key, zyg_key in (('alt_freq', 'zygosity'),
                                  ('n_alt_freq', 'n_zygosity')):
            if zyg_key in self:
                zyg = np.repeat(0.5, len(self))
                vals = self[freq_key].values
                zyg[vals >= hom_freq] = 1.0
                zyg[vals < het_freq] = 0.0
                self[zyg_key] = zyg
        return self

    def heterozygous(self):
        """Subset to only heterozygous variants.

        Use 'zygosity' or 'n_zygosity' genotype values (if present) to exclude
        variants with value 0.0 or 1.0.
        If these columns are missing, or there are no heterozygous variants,
        then return the full (input) set of variants.

        Returns
        -------
        VariantArray
            The subset of `self` with heterozygous genotype, or allele frequency
            between the specified thresholds.
        """
        if 'zygosity' in self:
            # Use existing genotype/zygosity info
            zygosity = self['n_zygosity' if 'n_zygosity' in self
                            else 'zygosity']
            het_idx = (zygosity != 0.0) & (zygosity != 1.0)
            if het_idx.any():
                # Only take het. subset if the subset is not empty
                self = self[het_idx]
        return self

    def mirrored_baf(self, above_half=None, tumor_boost=False):
        """Mirrored B-allele frequencies (BAFs).

        Parameters
        ----------
        above_half : bool or None
            If specified, flip BAFs to be all above 0.5 (True) or below 0.5
            (False), respectively, for consistency. Otherwise, if None, mirror
            in the direction of the majority of BAFs.
        tumor_boost : bool
            Normalize tumor-sample allele frequencies to the matched normal
            sample's allele frequencies.

        Returns
        -------
        float array
            Mirrored b-allele frequencies, the same length as `self`. May
            contain NaN values.
        """
        if tumor_boost and "n_alt_freq" in self:
            alt_freq = self.tumor_boost()
        else:
            alt_freq = self["alt_freq"]
        return _mirrored_baf(alt_freq, above_half)

    def tumor_boost(self):
        """TumorBoost normalization of tumor-sample allele frequencies.

        De-noises the signal for detecting LOH.

        See: TumorBoost, Bengtsson et al. 2010
        """
        if not ("alt_freq" in self and "n_alt_freq" in self):
            raise ValueError("TumorBoost requires a matched tumor and normal "
                             "pair of samples in the VCF.")
        return _tumor_boost(self["alt_freq"],  self["n_alt_freq"])



def _mirrored_baf(vals, above_half=None):
    shift = (vals - .5).abs()
    if above_half is None:
        above_half = (vals.median() > .5)
    if above_half:
        return .5 + shift
    else:
        return .5 - shift


def _tumor_boost(t_freqs, n_freqs):
    """Normalize tumor-sample allele frequencies.

    boosted = { 0.5 (t/n)           if t < n
                1 - 0.5(1-t)/(1-n)  otherwise

    See: TumorBoost, Bengtsson et al. 2010
    """
    lt_mask = (t_freqs < n_freqs)
    lt_idx = np.nonzero(lt_mask)[0]
    gt_idx = np.nonzero(~lt_mask)[0]
    out = pd.Series(np.zeros_like(t_freqs))
    out[lt_idx] = 0.5 * t_freqs.take(lt_idx) / n_freqs.take(lt_idx)
    out[gt_idx] = 1 - 0.5 * (1 - t_freqs.take(gt_idx)
                            ) / (1 - n_freqs.take(gt_idx))
    return out


def _allele_specific_copy_numbers(segarr, varr, ploidy=2):
    """Split total copy number between alleles based on BAF.

    See: PSCBS, Bentsson et al. 2011
    """
    # TODO fix ploidy on allosomes
    seg_depths = ploidy * np.exp2(segarr["log2"])
    seg_bafs = varr.baf_by_ranges(segarr, above_half=True)
    cn1 = 0.5 * (1 - seg_bafs) * seg_depths
    cn2 = seg_depths - cn1
    # segout = segarr.copy()
    # segout.update({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
    # return segout
    return pd.DataFrame({"baf": seg_bafs, "cn1": cn1, "cn2": cn2})
