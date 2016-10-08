"""An array of genomic intervals, treated as variant loci."""
from __future__ import absolute_import, division, print_function
from builtins import map, next, str, zip
from past.builtins import basestring

import numpy as np
import pandas as pd

from . import gary


class VariantArray(gary.GenomicArray):
    """An array of genomic intervals, treated as variant loci.

    Required columns: chromosome, start, end, ref, alt
    """
    _required_columns = ("chromosome", "start", "end", "ref", "alt")
    _required_dtypes = (str, int, int, str, str)
    # Extra: somatic, zygosity, depth, alt_count, alt_freq

    def __init__(self, data_table, meta_dict=None):
        gary.GenomicArray.__init__(self, data_table, meta_dict)

    def baf_by_ranges(self, ranges, summary_func=pd.Series.median,
                      above_half=None, tumor_boost=True):
        """Aggregate variant (b-allele) frequencies in each given bin.

        Get the average BAF in each of the bins of another genomic array:
        BAFs are mirrored (see `mirrored_baf`), grouped in each bin of `ranges`,
        and summarized with `summary_func`, by default the median.

        Parameters
        ----------
        ranges : GenomicArray or subclass
            Bins for grouping the variants in `self`.
        summary_func : callable
            Function to reduce BAF values to one number; by default the mirrored
            BAF median.
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
        if tumor_boost and "n_alt_freq" in self:
            self = self.copy()
            self["alt_freq"] = self.tumor_boost()
        return np.asarray([
            (summary_func(_mirrored_baf(subvarr["alt_freq"], above_half))
             if len(subvarr) else np.nan)
            for _bin, subvarr in self.by_ranges(ranges, mode='outer',
                                                keep_empty=True)])

    def heterozygous(self, min_freq=None, max_freq=None):
        """Subset to only heterozygous variants.

        If `min_freq` and `max_freq` are not specified, use "zygosity" genotype
        values, excludes variants with value 0.0 or 1.0.

        Parameters
        ----------
        min_freq : float
            Return only variants with alt allele frequency of at least this
            value.
        max_freq : float
            Return only variants with alt allele frequency of at most this
            value.

        Returns
        -------
        VariantArray
            The subset of `self` with heterozygous genotype, or allele frequency
            between the specified thresholds.
        """
        idx_het = None
        if 'zygosity' in self and min_freq is None and max_freq is None:
            # Use existing genotype/zygosity info
            zygosity = self['n_zygosity' if 'n_zygosity' in self
                            else 'zygosity']
            idx_het = (zygosity != 0.0) & (zygosity != 1.0)
        elif 'alt_freq' in self and (min_freq or max_freq):
            # Decide zygosity from allele frequency
            freq = self['n_alt_freq' if 'n_alt_freq' in self
                        else 'alt_freq']
            idx_het = True
            if min_freq:
                idx_het = idx_het & (freq >= min_freq)
            if max_freq:
                idx_het = idx_het & (freq <= max_freq)

        if idx_het is not None and idx_het.any():
            return self[idx_het]
        else:
            # Fallback -- ought to return something
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
        if "n_alt_freq" in self:
            n_freqs = self["n_alt_freq"]
        else:
            raise ValueError("TumorBoost requires a paired normal sample in "
                             "the VCF.")
        t_freqs = self["alt_freq"]
        return _tumor_boost(t_freqs, n_freqs)



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
    seg_depths = ploidy * np.exp2(segarr["log2"])
    seg_bafs = varr.baf_by_ranges(segarr, above_half=True)
    cn1 = 0.5 * (1 - seg_bafs) * seg_depths
    cn2 = seg_depths - cn1
    # segout = segarr.copy()
    # segout.update({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
    # return segout
    return pd.DataFrame({"baf": seg_bafs, "cn1": cn1, "cn2": cn2})
