How to do things
================

Sample selection
----------------

You can use ``cnvkit.py metrics *.cnr -s *.cns`` to see if any samples are
especially noisy. See the :ref:`metrics` section of the page :doc:`reports`.

If you're capturing the standard error text from CNVkit when running the
:ref:`batch` or :ref:`coverage` commands, you should see some summary statistics
when it finishes calculating the coverage of each sample.

In practice, the minimum coverage thresholds I've seen in use are 30x for
whole-genome sequencing of germline samples and at least 50x (sometimes 100x  or
even 250x for clinical sequencing services) for tumor samples that may have
significant normal-cell contamination and subclonal tumor-cell populations.
Those cutoffs are geared toward SNV calling, though. CNVkit will usually call
larger CNAs reliably down to about 10x on-target coverage, but there will tend
to be more spurious segments, and smaller-scale CNAs and generally
shattered-looking chromosomes can be hard to infer below that point.

Visualization
-------------

A reasonable approach is to:

1. Draw a heatmap of all your samples from the .cns files
2. Run "gainloss" to identify the focal gains and losses that seem plausibly
   interesting.
3. Use a script to iterate through each gene of interest in each sample and run
   "scatter -g $gene_name" to generate a plot of the copy number around that
   specific gene, saving to PDF.
4. Combine the PDFs ("pdfunite" on Linux/Unix), review them, and send the
   plausible ones along to your collaborator.

