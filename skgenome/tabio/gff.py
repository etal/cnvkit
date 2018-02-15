"""I/O for Generic Feature Format (GFF).

Columns:

1. seqid/reference/seqname/chromosome
2. source
3. type/method/feature
4. start: in 1-based integer coordinates
5. end: in 1-based integer coordinates
6. score: float or '.' (for NA)
7. strand: [+-.?]
8. phase/frame: [012.]
9. attribute/group: string

Specs:

- http://www.ensembl.org/info/website/upload/gff.html
- http://gmod.org/wiki/GFF3
- http://gmod.org/wiki/GFF2
- http://mblab.wustl.edu/GTF2.html
"""
from __future__ import absolute_import, division, print_function
import re

import pandas as pd


def read_gff(infile, tag=None):
    """Read a GFF3/GTF/GFF2 file into a DataFrame.

    Works for all three formats because we only try extract the gene name, at
    most, from column 9.

    Parameters
    ----------
    infile : filename or open handle
        Source file.
    tag : str
        GFF attributes tag to use for extracting gene names. In GFF3, this is
        standardized as "Name", and in GTF it's "gene_id". (Neither spec is
        consistently followed, so the parser will by default look for eith er of
        those tags and also "gene_name" and "gene".)

    """
    colnames = ['chromosome', 'source', 'type', 'start', 'end',
                'score', 'strand', 'phase', 'attribute']
    coltypes = ['str', 'str', 'str', 'int', 'int',
                'str', 'str', 'str', 'str']
    dframe = pd.read_table(infile, comment='#', header=None, na_filter=False,
                           names=colnames, dtype=dict(zip(colnames, coltypes)))
    dframe = (dframe
              .assign(start=dframe.start - 1,
                      score=dframe.score.replace('.', 'nan').astype('float'))
              .sort_values(['chromosome', 'start', 'end'])
              .reset_index(drop=True))
    if len(dframe):
        if tag:
            # Look only for the specified tag
            rx = re.compile(tag + r'[= ]"?(\S+?)"?(;|$)"')
            matches = dframe['attribute'].str.extract(rx)
            if len(matches):
                dframe['gene'] = matches
        else:
            # Default to a set of likely relevant tags
            rx = re.compile(r'(Name|gene_id|gene_name|gene)[= ]"?(\S+?)"?(;|$)')
            matches = dframe['attribute'].str.extractall(rx)
            if len(matches):
                dframe['gene'] = matches.xs(0, level=1)[1]
    if 'gene' in dframe.columns:
        dframe['gene'] = dframe['gene'].fillna('-').astype('str')
    else:
        dframe['gene'] = ['-'] * len(dframe)
    return dframe
