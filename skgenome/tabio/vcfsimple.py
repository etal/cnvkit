"""Simple VCF I/O.

Read only coordinate info & store the remaining columns as unparsed strings.
Just enough functionality to extract a subset of samples and/or perform
bedtools-like operations on VCF records.
"""
import logging

import numpy as np
import pandas as pd
from Bio.File import as_handle


# TODO save VCF header (as string, the whole text block) in meta{header=}
def read_vcf_simple(infile):
    """Read VCF file w/o samples."""
    # ENH: Make all readers return a tuple (header_string, body_table)
    # ENH: usecols -- need to trim dtypes dict to match?
    header_lines = []
    with as_handle(infile, 'r') as handle:
        for line in handle:
            if line.startswith('##'):
                header_lines.append(line)
            else:
                assert line.startswith('#CHR')
                header_line = line
                header_lines.append(line)
                break

        # Extract sample names from VCF header, keep as column names
        header_fields = header_line.split('\t')
        sample_ids = header_fields[9:]
        colnames = ['chromosome', 'start', 'id', 'ref', 'alt',
                    'qual', 'filter', 'info', 'format'] + sample_ids
        dtypes = {c: str for c in colnames}
        dtypes['start'] = int
        del dtypes['qual']
        table = pd.read_csv(handle, sep='\t', header=None, na_filter=False,
                            names=colnames,
                            converters={'qual': parse_qual},
                            dtype=dtypes)
    # ENH: do things with filter, info
    table['start'] -= 1
    table['end'] = table['info'].apply(parse_end_from_info)
    set_ends(table)
    logging.info("Loaded %d plain records", len(table))
    return table


def read_vcf_sites(infile):
    colnames = ['chromosome', 'start', 'id', 'ref', 'alt',
                'qual', 'filter', 'end']
    dtypes = {'chromosome': str, 'start': int, 'id': str,
              'ref': str, 'alt': str, 'filter': str}
    table = pd.read_csv(infile, sep='\t', comment='#',
                        header=None, na_filter=False,
                        names=colnames, usecols=colnames,
                        converters={'end': parse_end_from_info,
                                    'qual': parse_qual},
                        dtype=dtypes)
    # Where END is missing, infer from allele lengths
    table['start'] -= 1
    set_ends(table)
    logging.info("Loaded %d plain records", len(table))
    return table


def parse_end_from_info(info):
    idx = info.find('END=')
    if idx == -1:
        return -1
    info = info[idx+4:]
    idx = info.find(';')
    if idx != -1:
        info = info[:idx]
    return int(info)


def parse_qual(qual):
    # ENH: only appy na_filter to this column
    if qual == '.':
        return np.nan
    return float(qual)


def set_ends(table):
    """Set 'end' field according to allele lengths."""
    need_end_idx = (table.end == -1)
    if need_end_idx.any():
        ref_sz = table.loc[need_end_idx, 'ref'].str.len()
        # TODO handle multiple alts -- split commas & take max len
        alt_sz = table.loc[need_end_idx, 'alt'].str.len()
        var_sz = alt_sz - ref_sz
        # TODO XXX if end > start, swap 'em?
        var_sz = var_sz.clip(lower=0)
        table.loc[need_end_idx, 'end'] = table.loc[need_end_idx, 'start'] + var_sz
