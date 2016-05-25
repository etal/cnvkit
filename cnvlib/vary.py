"""An array of genomic intervals, treated as variant loci."""
from __future__ import absolute_import, division, print_function
from builtins import next
from builtins import str
from builtins import map
from builtins import zip
from past.builtins import basestring

import logging

import numpy as np
import pandas as pd
import vcf

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

    def mirrored_baf(self, above_half=True, tumor_boost=True):
        """Mirrored B-allele frequencies (BAFs).

        Flip BAFs to be all above 0.5 (if `above_half`) or below 0.5, for
        consistency.

        With `tumor_boost`, normalize tumor-sample allele frequencies to the
        matched normal allele frequencies.
        """
        if tumor_boost and "n_alt_freq" in self:
            alt_freqs = self.tumor_boost()
        else:
            alt_freqs = self["alt_freq"]
        fromhalf = np.abs(alt_freqs - 0.5)
        if above_half:
            return fromhalf + 0.5
        else:
            return 0.5 - fromhalf

    def tumor_boost(self):
        """TumorBoost normalization of tumor-sample allele frequencies.

        De-noises the signal for detecting LOH.
        """
        if "n_alt_freq" in self:
            n_freqs = self["n_alt_freq"]
        else:
            raise ValueError("TumorBoost requires a paired normal sample in "
                             "the VCF.")
        t_freqs = self["alt_freq"]
        return _tumor_boost(t_freqs, n_freqs)



def _tumor_boost(t_freqs, n_freqs):
    """Normalize tumor-sample allele frequencies.

    boosted = { 0.5 (t/n)           if t < n
                1 - 0.5(1-t)/(1-n)  otherwise

    See: TumorBoost, Bengtsson et al. 2010
    """
    lt_mask = (t_freqs < n_freqs)
    lt_idx = np.nonzero(lt_mask)[0]
    gt_idx = np.nonzero(~lt_mask)[0]
    out = np.zeros_like(t_freqs)
    out[lt_idx] = 0.5 * t_freqs.take(lt_idx) / n_freqs.take(lt_idx)
    out[gt_idx] = 1 - 0.5 * (1 - t_freqs.take(gt_idx)
                            ) / (1 - n_freqs.take(gt_idx))
    return out


def _allele_specific_copy_numbers(segarr, varr, ploidy=2):
    """Split total copy number between alleles based on BAF.

    See: PSCBS, Bentsson et al. 2011
    """
    seg_depths = ploidy * np.exp2(segarr["log2"])
    seg_bafs = np.asarray([(np.median(subvarr.mirrored_baf())
                            if len(subvarr) else np.nan)
                           for _seg, subvarr in varr.by_ranges(segarr)])
    cn1 = 0.5 * (1 - seg_bafs) * seg_depths
    cn2 = seg_depths - cn1
    # segout = segarr.copy()
    # segout.update({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
    # return segout
    return pd.DataFrame({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
