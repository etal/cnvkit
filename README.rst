======
CNVkit
======

A command-line toolkit and Python library for detecting copy number variants
and alterations genome-wide from targeted DNA sequencing.

Read the full documentation at: http://cnvkit.readthedocs.org

Questions, comments, help and discussion on SeqAnswers:
http://seqanswers.com/forums/showthread.php?t=47910


.. image:: https://travis-ci.org/etal/cnvkit.svg?branch=master
    :target: https://travis-ci.org/etal/cnvkit

.. image:: https://landscape.io/github/etal/cnvkit/master/landscape.svg
   :target: https://landscape.io/github/etal/cnvkit/master
   :alt: Code Health


Try it
======

You can easily run CNVkit on your own data without installing it by using our
DNAnexus app:
https://platform.dnanexus.com/app/cnvkit_batch

A Galaxy tool is also available for testing (but requires CNVkit installation,
see below):
https://testtoolshed.g2.bx.psu.edu/view/etal/cnvkit


Installation
============

CNVkit runs on Python 2.7. Your operating system might already provide Python
2.7, which you can check on the command line::

    python --version

If your operating system already includes Python 2.6 , I suggest installing
Python 2.7 alongside it instead of attempting to upgrade it in-place.
Your package manager might provide both versions.

The recommended way to install Python 2.7 and some of CNVkit's dependencies
without affecting the rest of your operating system is by installing either
`Anaconda <https://store.continuum.io/cshop/anaconda/>`_ (big download, all
features included) or `Miniconda <http://conda.pydata.org/miniconda.html>`_
(smaller download, minimal environment). Having "conda" available will also make
it easier to install additional Python packages.


From a Python package repository
--------------------------------

Reasonably up-to-date CNVkit packages are available on `PyPI
<https://pypi.python.org/pypi/CNVkit>`_ and can be installed using `conda
<http://repo.continuum.io/miniconda/>`_ (preferred on Mac OS X, and
a solid choice on Linux too)::

    conda install cnvkit

Or `pip <https://pip.pypa.io/en/latest/installing.html>`_ (usually works on
Linux if the dependencies are installed, see below)::

    pip install cnvkit

In either case, to run the recommended segmentation algorithms CBS and Fused
Lasso, you will need to also install the R dependencies DNAcopy and PSCBS (see
below).

From source
-----------

The script ``cnvkit.py`` requires no installation and can be used in-place. Just
install the dependencies.

To install the main program, supporting scripts and ``cnvlib`` Python library,
use ``setup.py`` as usual::

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


On Mac OS X you may instead find it much easier to first install the Python
package manager `Miniconda <http://repo.continuum.io/miniconda/>`_, or the full
`Anaconda <https://store.continuum.io/cshop/anaconda/>`_ distribution, which
includes most of the dependencies. Then install the rest of CNVkit's dependencies::

    conda install numpy matplotlib reportlab biopython pysam pyvcf

Otherwise, we recommend using `Homebrew <http://brew.sh/>`_ to install an
up-to-date Python (e.g. ``brew install python``) and as many of the Python
packages as possible (primarily NumPy, SciPy and matplotlib). Then, proceed with
pip::

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
=======

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

