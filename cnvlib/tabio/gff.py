"""I/O for Generic Feature Format (GFF).

Columns:

1. seqid, reference, seqname, chromosome
2. source
3. type, method, feature
4. start: in 1-based integer coordinates
5. end: in 1-based integer coordinates
6. score: float or '.' (for NA)
7. strand: [+-.?]
8. phase, frame: [012.]
9. attribute, group: string

Specs:

-  http://www.ensembl.org/info/website/upload/gff.html
-  http://gmod.org/wiki/GFF3
#  http://gmod.org/wiki/GFF2
"""
from __future__ import absolute_import, division, print_function

import pandas as pd


def read_gff(infile):
    """Works for GFF2/GTF/GFF3 because we don't parse inside column 9."""
    colnames = ['chromosome', 'source', 'type', 'start', 'end',
                'score', 'strand', 'phase', 'attribute']
    coltypes = ['str', 'str', 'str', 'int', 'int',
                'str', 'str', 'str', 'str']
    usecols = ['chromosome',
               #'source', 'type',
               'start', 'end',
               #'score',
               'strand',
               #'phase', 'attribute'
              ]
    dframe = pd.read_table(infile, comment='#', header=None, na_filter=False,
                           names=colnames, #  usecols=usecols,
                           dtype=dict(zip(colnames, coltypes)))
    return (dframe
            .assign(start=dframe.start - 1,
                    score=dframe.score.replace('.', 'nan').astype('float'))
            .sort_values(['chromosome', 'start', 'end'])
            .reset_index(drop=True))
