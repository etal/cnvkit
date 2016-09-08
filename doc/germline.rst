Germline analysis
=================

.. TODO - see e-mails, biostars, notes

CNVkit is less accurate in detecting CNVs smaller than 1 Mbp.

The ``--drop-low-coverage`` option (see :doc:`tumor`) should not be used; it
will typically remove germline deep deletions altogether, which is not
desirable.

Watch for mosaicism in CNVs, resulting in non-integer copy numbers (i.e. smaller
log2 value than a single-copy CNV would indicate); they're more common than
often thought.
