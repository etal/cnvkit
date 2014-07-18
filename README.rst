======
CNVkit
======

A command-line toolkit and Python library for detecting copy number variations
and alterations from targeted DNA sequencing.

Read the full documentation at: http://cnvkit.readthedocs.org

.. image:: https://travis-ci.org/etal/cnvkit.svg?branch=master
    :target: https://travis-ci.org/etal/cnvkit

Installation
============

The script cnvkit.py requires no installation and can be used in-place. Just
install the dependencies.

To install the script and Python library (and some associated scripts), use
setup.py as usual::

    python setup.py build
    python setup.py install

Python dependencies
-------------------

Install these via ``pip`` or ``conda``:

- `Biopython <http://biopython.org/wiki/Main_Page>`_
- `Reportlab <https://bitbucket.org/rptlab/reportlab>`_
- `matplotlib <http://matplotlib.org>`_
- `NumPy <http://www.numpy.org/>`_
- `pysam <https://github.com/pysam-developers/pysam>`_
- `PyVCF <https://github.com/jamescasbon/PyVCF>`_

On Ubuntu or Debian Linux::

    sudo apt-get install python-numpy python-matplotlib python-reportlab python-pip
    sudo pip install biopython pysam pyvcf --upgrade


On Mac OS X you may instead find it much easier to install the Anaconda distribution
(https://store.continuum.io/cshop/anaconda/ or
http://repo.continuum.io/miniconda/), which includes most of the dependencies
already::

    conda install numpy matplotlib reportlab biopython pysam pyvcf

Otherwise, we recommend using Homebrew (http://brew.sh/) or MacPorts to
install an up-to-date Python (e.g. ``brew install python``) and as many of the
Python packages as possible (primarily numpy, scipy and matplotlib). Then, 
proceed with pip::

    sudo pip install numpy matplotlib reportlab biopython pysam pyvcf


R dependencies
--------------

Copy number segmentation currently depends on R packages.

In CRAN:

- PSCBS
- matrixStats (for pruneByHClust in PSCBS)

PSCBS depends on the package ``DNAcopy`` which part of Bioconductor, and cannot
be installed through CRAN directly.  To use it, you must first install
Bioconductor, then ``DNAcopy``, then the ``PSCBS`` package as shown below.

Within R::

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("DNAcopy")
    > install.packages(c("PSCBS", "matrixStats", "cghFLasso"))


Testing
-------

You can test your installation by running the CNVkit pipeline on the example
files in the ``test/`` directory. The pipeline is implemented as a Makefile and
can be run with the ``make`` command (standard on Unix/Linux/Mac OS X systems)::

    cd test/
    make

If this pipeline completes successfully (it should take less than a minute),
you've installed CNVkit correctly.

To run the pipeline on additional, larger example file sets (named TR and EX),
do this::

    make all

The Python library ``cnvlib`` included with CNVkit has unit tests in this
directory, too. To run the test suite::

    python test_cnvlib.py


Citation
========

We are in the process of publishing a manuscript describing CNVkit in detail.
For now, if you use this software in a publication, please cite it by the URL:
http://github.com/etal/cnvkit

