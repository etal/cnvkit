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

- http://gmod.org/wiki/GFF3
- http://gmod.org/wiki/GFF2
- http://mblab.wustl.edu/GTF2.html
- http://www.ensembl.org/info/website/upload/gff.html
- https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo
"""
from __future__ import absolute_import, division, print_function
import logging
import re

import pandas as pd


def read_gff(infile, tag=r'(Name|gene_id|gene_name|gene)', keep_type=None):
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
    keep_type : str
        If specified, only keep rows with this value in the 'type' field (column
        3). In GFF3, these terms are standardized in the Sequence Ontology
        Feature Annotation (SOFA).
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
    if keep_type:
        ok_type = (dframe['type'] == keep_type)
        logging.info("Keeping %d '%s' / %d total records",
                     ok_type.sum(), keep_type, len(dframe))
        dframe = dframe[ok_type]
    if len(dframe):
        rx = re.compile(tag + r'[= ]"?(?P<gene>\S+?)"?(;|$)')
        matches = dframe['attribute'].str.extract(rx, expand=True)['gene']
        if len(matches):
            dframe['gene'] = matches
    if 'gene' in dframe.columns:
        dframe['gene'] = dframe['gene'].fillna('-').astype('str')
    else:
        dframe['gene'] = ['-'] * len(dframe)
    return dframe
