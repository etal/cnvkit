"""Transform bait intervals into targets more suitable for CNVkit."""
from __future__ import division, absolute_import

import logging


def shorten_labels(gene_labels):
    """Reduce multi-accession interval labels to the minimum consistent.

    So: BED or interval_list files have a label for every region. We want this
    to be a short, unique string, like the gene name. But if an interval list is
    instead a series of accessions, including additional accessions for
    sub-regions of the gene, we can extract a single accession that covers the
    maximum number of consecutive regions that share this accession.

    e.g.::

        ...
        mRNA|JX093079,ens|ENST00000342066,mRNA|JX093077,ref|SAMD11,mRNA|AF161376,mRNA|JX093104
        ens|ENST00000483767,mRNA|AF161376,ccds|CCDS3.1,ref|NOC2L
        ...

    becomes::

        ...
        mRNA|AF161376
        mRNA|AF161376
        ...
    """
    longest_name_len = 0
    curr_names = set()
    curr_gene_count = 0

    for label in gene_labels:
        next_names = set(label.rstrip().split(','))
        assert len(next_names)
        overlap = curr_names.intersection(next_names)
        if overlap:
            # Continuing the same gene; update shared accessions
            curr_names = filter_names(overlap)
            curr_gene_count += 1
        else:
            # End of the old gene -- emit shared name(s)
            for _i in range(curr_gene_count):
                out_name = shortest_name(curr_names)
                yield out_name
                longest_name_len = max(longest_name_len, len(out_name))

            # Start of a new gene
            curr_gene_count = 1
            curr_names = next_names
    # Final emission
    for _i in range(curr_gene_count):
        out_name = shortest_name(curr_names)
        yield out_name
        longest_name_len = max(longest_name_len, len(out_name))

    logging.info("Longest name length: %d", longest_name_len)


def filter_names(names, exclude=('mRNA',)):
    """Remove less-meaningful accessions from the given set."""
    if len(names) > 1:
        ok_names = set(n for n in names
                       if not any(n.startswith(ex) for ex in exclude))
        if ok_names:
            return ok_names
    # Names are not filter-worthy; leave them as they are for now
    return names


def shortest_name(names):
    """Return the shortest trimmed name from the given set."""
    name = min(filter_names(names), key=len)
    if len(name) > 2 and '|' in name[1:-1]:
        # Split 'DB|accession' and extract the accession sans-DB
        name = name.split('|')[-1]
    return name
