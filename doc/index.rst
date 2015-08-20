CNVkit: Genome-wide copy number from targeted DNA sequencing
============================================================

:Author: Eric Talevich
:Contact: eric.talevich@ucsf.edu
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_
:Source code: `GitHub <https://github.com/etal/cnvkit>`_
:Packages: `PyPI <https://pypi.python.org/pypi/CNVkit>`_ | `Docker <https://hub.docker.com/r/etal/cnvkit/>`_ | `Galaxy <https://testtoolshed.g2.bx.psu.edu/view/etal/cnvkit>`_ | `DNAnexus <https://platform.dnanexus.com/app/cnvkit_batch>`_
:Q&A: `Biostars <https://www.biostars.org/t/CNVkit/>`_ | `SeqAnswers <http://seqanswers.com/forums/showthread.php?t=47910>`_

CNVkit is a Python library and command-line software toolkit to infer and
visualize copy number from targeted DNA sequencing data. It is designed for use
with hybrid capture, including both whole-exome and custom target panels, and
short-read sequencing platforms such as Illumina and Ion Torrent.

.. toctree::

    quickstart

.. gallery?


Command line usage
------------------

.. toctree::
    :maxdepth: 2

    pipeline
    plots
    reports
    importexport
    scripts


FAQ & How To
------------

.. toctree::
    :maxdepth: 3

    fileformats
    bias
    calling
    heterogeneity
    nonhybrid


Python API
----------

.. toctree::
    :maxdepth: 2

    cnvlib


Citation
========

We are in the process of publishing a manuscript describing CNVkit.
If you use this software in a publication, for now, please cite our preprint
manuscript by DOI, like so:

    Eric Talevich, A. Hunter Shain, Thomas Botton, Boris C. Bastian (2014)
    CNVkit: Copy number detection and visualization for targeted sequencing
    using off-target reads.
    *bioRxiv* doi: http://dx.doi.org/10.1101/010876

A poster presentation can be viewed at `F1000 Research
<http://f1000research.com/posters/1096236>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

