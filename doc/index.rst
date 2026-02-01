CNVkit: Genome-wide copy number from high-throughput sequencing
===============================================================

:Source code: `GitHub <https://github.com/etal/cnvkit>`_
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_
:Packages: `PyPI <https://pypi.python.org/pypi/CNVkit>`_ | `Debian <https://packages.debian.org/search?keywords=cnvkit>`_ | `Docker <https://hub.docker.com/r/etal/cnvkit/>`_ | `Galaxy <https://testtoolshed.g2.bx.psu.edu/view/etal/cnvkit>`_ | `DNAnexus <https://platform.dnanexus.com/app/cnvkit_batch>`_
:Article: `PLOS Computational Biology <http://dx.doi.org/10.1371/journal.pcbi.1004873>`_
:Q&A: `Biostars <https://www.biostars.org/t/CNVkit/>`_
:Consulting: Contact `DNAnexus Science <http://go.dnanexus.com/l/457982/2018-05-23/6v4h43>`_


CNVkit is a Python library and command-line software toolkit to infer and
visualize copy number from high-throughput DNA sequencing data. It is designed
for use with hybrid capture, including both whole-exome and custom target
panels, and short-read sequencing platforms such as Illumina and Ion Torrent.

**Requirements:** Python 3.10 or later. The R statistical environment with the
DNAcopy package is also required for the default CBS segmentation algorithm.

.. toctree::

    quickstart
    whoelse

.. gallery?


Citation
--------

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
pomegranate (HMM segmentation methods):
    Schreiber, J. (2018).
    `pomegranate: Fast and Flexible Probabilistic Modeling in Python.
    <http://jmlr.org/papers/v18/17-636.html>`_
    *Journal of Machine Learning Research* 18(164):1−6.


Command line usage
------------------

.. toctree::
    :maxdepth: 2

    pipeline
    plots
    reports
    importexport
    rna
    scripts


FAQ
---

.. toctree::
    :maxdepth: 2

    fileformats
    baf
    bias
    sex


How To
------

.. toctree::
    :maxdepth: 2

    calling
    tumor
    heterogeneity
    germline
    nonhybrid


Python API
----------

.. toctree::
    :maxdepth: 2

    cnvlib
    skgenome


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
