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

    Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2014).
    CNVkit: Copy number detection and visualization for targeted sequencing using off-target reads.
    *bioRxiv* doi: http://dx.doi.org/10.1101/010876

A poster presentation can be viewed at `F1000 Research
<http://f1000research.com/posters/1096236>`_.

Who is using CNVkit?
--------------------

`Google Scholar
<https://scholar.google.com/scholar?cites=4696125388809243311&as_sdt=2005&sciodt=0,5&hl=en>`_
lists some of the references where CNVkit has been used by other researchers.

We'd like to highlight:

* McCreery, M.Q., Halliwill, K.D., Chin, D., Delrosario, R., Hirst, G., Vuong, P., ... & Balmain, A. (2015).
  Evolution of metastasis revealed by mutational landscapes of chemically induced skin cancers.
  *Nature Medicine*.
  `doi:10.1038/nm.3979 <http://dx.doi.org/10.1038/nm.3979>`_
* Shain, A.H., Yeh, I., Kovalyshyn, I., Sriharan, A., Talevich, E., Gagnon, A., ... & Bastian, B.C. (2015).
  The Genetic Evolution of Melanoma from Precursor Lesions.
  *New England Journal of Medicine*, 373(20), 1926-1936.
  `doi:10.1056/NEJMoa1502583 <http://dx.doi.org/10.1056/NEJMoa1502583>`_
* Shain, A.H., Garrido, M., Botton, T., Talevich, E., Yeh, I., Sanborn, J.Z., ... & Bastian, B.C. (2015).
  Exome sequencing of desmoplastic melanoma identifies recurrent NFKBIE promoter mutations and diverse activating mutations in the MAPK pathway.
  *Nature Genetics*, 47(10), 1194-1199.
  `doi:10.1038/ng.3382 <http://dx.doi.org/10.1038/ng.3382>`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

