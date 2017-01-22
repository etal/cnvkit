CNVkit: Genome-wide copy number from targeted DNA sequencing
============================================================

:Author: Eric Talevich
:Contact: eric.talevich@ucsf.edu
:Source code: `GitHub <https://github.com/etal/cnvkit>`_
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_
:Packages: `PyPI <https://pypi.python.org/pypi/CNVkit>`_ | `Docker <https://hub.docker.com/r/etal/cnvkit/>`_ | `Galaxy <https://testtoolshed.g2.bx.psu.edu/view/etal/cnvkit>`_ | `DNAnexus <https://platform.dnanexus.com/app/cnvkit_batch>`_
:Article: `PLOS Computational Biology <http://dx.doi.org/10.1371/journal.pcbi.1004873>`_
:Q&A: `Biostars <https://www.biostars.org/t/CNVkit/>`_


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
    :maxdepth: 2

    calling
    baf
    tumor
    heterogeneity
    germline
    nonhybrid
    sex
    bias
    fileformats


Python API
----------

.. toctree::
    :maxdepth: 2

    cnvlib


Citation
========

If you use this software in a publication, please cite our paper describing CNVkit:

    Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2014).
    `CNVkit: Genome-wide copy number detection and visualization from targeted sequencing.
    <http://dx.doi.org/10.1371/journal.pcbi.1004873>`_
    *PLOS Computational Biology* 12(4): e1004873

A poster presentation can be viewed at `F1000 Research
<http://f1000research.com/posters/1096236>`_.

Who else is using CNVkit?
-------------------------

`Google Scholar
<https://scholar.google.com/scholar?cites=4696125388809243311&as_sdt=2005&sciodt=0,5&hl=en>`_
lists some of the studies where CNVkit has been used by other researchers.
We'd like to highlight:

* McCreery, M.Q. *et al.* (2015).
  `Evolution of metastasis revealed by mutational landscapes of chemically
  induced skin cancers. <http://dx.doi.org/10.1038/nm.3979>`_
  *Nature Medicine* 21, 1514–1520
* Shain, A.H. *et al.* (2015).
  `Exome sequencing of desmoplastic melanoma identifies recurrent NFKBIE
  promoter mutations and diverse activating mutations in the MAPK pathway.
  <http://dx.doi.org/10.1038/ng.3382>`_
  *Nature Genetics*, 47(10), 1194-1199
* Shain, A.H. *et al.* (2015).
  `The Genetic Evolution of Melanoma from Precursor Lesions.
  <http://dx.doi.org/10.1056/NEJMoa1502583>`_
  *New England Journal of Medicine*, 373(20), 1926-1936

Specific support for CNVkit is included in
`bcbio-nextgen <https://bcbio-nextgen.readthedocs.io/>`_,
`THetA2 <http://compbio.cs.brown.edu/projects/theta/>`_,
and `MetaSV <http://bioinform.github.io/metasv/>`_.
CNVkit is also available on the commercial platforms
`DNAnexus <http://www.dnanexus.com/>`_,
`Bina RAVE <http://www.bina.com/rave>`_, and
`Diploid InHelix <http://www.diploid.com/inhelix>`_.

Finally, CNVkit can :ref:`export` files to several standard formats that can be
used with many other software packages, including `BioDiscovery Nexus Copy
Number <http://www.biodiscovery.com/nexus-copy-number/>`_ and `Integrative
Genomics Viewer (IGV) <http://software.broadinstitute.org/software/igv/>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

