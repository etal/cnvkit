CNVkit: Genome-wide copy number from high-throughput sequencing
===============================================================

:Author: Eric Talevich
:Contact: eric.talevich@ucsf.edu
:Source code: `GitHub <https://github.com/etal/cnvkit>`_
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_
:Packages: `PyPI <https://pypi.python.org/pypi/CNVkit>`_ | `Docker <https://hub.docker.com/r/etal/cnvkit/>`_ | `Galaxy <https://testtoolshed.g2.bx.psu.edu/view/etal/cnvkit>`_ | `DNAnexus <https://platform.dnanexus.com/app/cnvkit_batch>`_
:Article: `PLOS Computational Biology <http://dx.doi.org/10.1371/journal.pcbi.1004873>`_
:Q&A: `Biostars <https://www.biostars.org/t/CNVkit/>`_


CNVkit is a Python library and command-line software toolkit to infer and
visualize copy number from high-throughput DNA sequencing data. It is designed
for use with hybrid capture, including both whole-exome and custom target
panels, and short-read sequencing platforms such as Illumina and Ion Torrent.

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

    fileformats
    sex
    calling
    baf
    bias
    tumor
    heterogeneity
    germline
    nonhybrid
    rna


Python API
----------

.. toctree::
    :maxdepth: 2

    cnvlib
    skgenome


Citation
========

If you use this software in a publication, please cite our paper describing CNVkit:

    Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2014).
    `CNVkit: Genome-wide copy number detection and visualization from targeted sequencing.
    <http://dx.doi.org/10.1371/journal.pcbi.1004873>`_
    *PLOS Computational Biology* 12(4):e1004873

Also please cite the supporting paper for the segmentation method you use:

PSCBS and DNAcopy (``cbs``, the default):
    - Olshen, A.B., Bengtsson, H., Neuvial, P., Spellman, P.T., Olshen, R.A., & Seshan, V.E. (2011).
      `Parent-specific copy number in paired tumor-normal studies using circular binary segmentation.
      <http://doi.org/10.1093/bioinformatics/btr329>`_
      *Bioinformatics* 27(15):2038–46.
    - Venkatraman, E.S., & Olshen, A.B. (2007).
      `A faster circular binary segmentation algorithm for the analysis of array CGH data.
      <http://doi.org/10.1093/bioinformatics/btl646>`_
      *Bioinformatics* 23(6):657–63
HaarSeg (``haar``):
    Ben-Yaacov, E., & Eldar, Y.C. (2008).
    `A fast and flexible method for the segmentation of aCGH data.
    <http://doi.org/10.1093/bioinformatics/btn272>`_
    *Bioinformatics* 24(16):i139-45.
CGH Fused Lasso (``flasso``):
    Tibshirani, R., & Wang, P. (2008).
    `Spatial smoothing and hot spot detection for CGH data using the fused lasso.
    <http://doi.org/10.1093/biostatistics/kxm013>`_
    *Biostatistics* 9(1):18–29


Who else is using CNVkit?
-------------------------

`Google Scholar
<https://scholar.google.com/scholar?cites=4696125388809243311&as_sdt=2005&sciodt=0,5&hl=en>`_
lists some of the studies where CNVkit has been used by other researchers.
We'd like to highlight:

* Seed, G. *et al.* (2017).
  `Gene copy number estimation from targeted next generation sequencing of
  prostate cancer biopsies: Analytic validation and clinical qualification.
  <http://dx.doi.org/10.1158/1078-0432.CCR-17-0972>`_
  *Clinical Cancer Research*
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

