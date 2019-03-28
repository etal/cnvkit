======
CNVkit
======

A command-line toolkit and Python library for detecting copy number variants
and alterations genome-wide from high-throughput sequencing.

Read the full documentation at: http://cnvkit.readthedocs.io

.. image:: https://travis-ci.org/etal/cnvkit.svg?branch=master
    :target: https://travis-ci.org/etal/cnvkit
    :alt: Build status

.. image:: https://landscape.io/github/etal/cnvkit/master/landscape.svg
    :target: https://landscape.io/github/etal/cnvkit/master
    :alt: Code health

.. image::   https://codecov.io/github/etal/cnvkit/coverage.svg?branch=master
    :target: https://codecov.io/github/etal/cnvkit?branch=master
    :alt: Test coverage


Support
=======

Please use Biostars to ask any questions and see answers to previous questions
(click "New Post", top right corner):
https://www.biostars.org/t/CNVkit/

Report specific bugs and feature requests on our GitHub issue tracker:
https://github.com/etal/cnvkit/issues/


Try it
======

You can easily run CNVkit on your own data without installing it by using our
`DNAnexus app <https://platform.dnanexus.com/app/cnvkit_batch>`_.

A `Galaxy tool <https://testtoolshed.g2.bx.psu.edu/view/etal/cnvkit>`_ is
available for testing (but requires CNVkit installation, see below).

A `Docker container <https://registry.hub.docker.com/u/etal/cnvkit/>`_ is also
available on Docker Hub, and the BioContainers community provides another on
`Quay <https://quay.io/repository/biocontainers/cnvkit>`_.

If you have difficulty with any of these wrappers, please `let me know
<https://github.com/etal/cnvkit/issues/>`_!


Installation
============

CNVkit runs on Python 3.5 and later. Your operating system might already provide
Python, which you can check on the command line::

    python --version

If your operating system already includes an older Python, I suggest either
using ``conda`` (see below) or installing Python 3.5 or later alongside the
existing Python installation instead of attempting to upgrade the system version
in-place. Your package manager might also provide Python 3.5+.

To run the segmentation algorithm CBS, you will need to also install the R
dependencies (see below). With ``conda``, this is included automatically.

Using Conda
-----------

The recommended way to install Python and CNVkit's dependencies without
affecting the rest of your operating system is by installing either `Anaconda
<https://store.continuum.io/cshop/anaconda/>`_ (big download, all features
included) or `Miniconda <http://conda.pydata.org/miniconda.html>`_ (smaller
download, minimal environment).
Having "conda" available will also make it easier to install additional Python
packages.

This approach is preferred on Mac OS X, and is a solid choice on Linux, too.

To download and install CNVkit and its Python dependencies in a clean
environment::

    # Configure the sources where conda will find packages
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

Then:

    # Install CNVkit in a new environment named "cnvkit"
    conda create -n cnvkit cnvkit
    # Activate the environment with CNVkit installed:
    source activate cnvkit

Or, in an existing environment::

    conda install cnvkit


From a Python package repository
--------------------------------

Up-to-date CNVkit packages are available on `PyPI
<https://pypi.python.org/pypi/CNVkit>`_ and can be installed using `pip
<https://pip.pypa.io/en/latest/installing.html>`_ (usually works on Linux if the
system dependencies listed below are installed)::

    pip install cnvkit


From source
-----------

The script ``cnvkit.py`` requires no installation and can be used in-place. Just
install the dependencies (see below).

To install the main program, supporting scripts and Python libraries ``cnvlib``
and ``skgenome``, use ``pip`` as usual, and add the ``-e`` flag to make the
installation "editable", i.e. in-place::

    git clone https://github.com/etal/cnvkit
    cd cnvkit/
    pip install -e .

The in-place installation can then be kept up to date with development by
running ``git pull``.


Python dependencies
-------------------

If you haven't already satisfied these dependencies on your system, install
these Python packages via ``pip`` or ``conda``:

- `Biopython <http://biopython.org/wiki/Main_Page>`_
- `Reportlab <https://bitbucket.org/rptlab/reportlab>`_
- `matplotlib <http://matplotlib.org>`_
- `NumPy <http://www.numpy.org/>`_
- `SciPy <http://www.scipy.org/>`_
- `Pandas <http://pandas.pydata.org/>`_
- `pyfaidx <https://github.com/mdshw5/pyfaidx>`_
- `pysam <https://github.com/pysam-developers/pysam>`_

On Ubuntu or Debian Linux::

    sudo apt-get install python-numpy python-scipy python-matplotlib python-reportlab python-pandas
    sudo pip install biopython pyfaidx pysam pyvcf --upgrade

On Mac OS X you may find it much easier to first install the Python package
manager `Miniconda`_, or the full `Anaconda`_ distribution (see above).
Then install the rest of CNVkit's dependencies::

    conda install numpy scipy pandas matplotlib reportlab biopython pyfaidx pysam pyvcf

Alternatively, you can use `Homebrew <http://brew.sh/>`_ to install an
up-to-date Python (e.g. ``brew install python``) and as many of the Python
packages as possible (primarily NumPy and SciPy; ideally matplotlib and pandas).
Then, proceed with pip::

    pip install numpy scipy pandas matplotlib reportlab biopython pyfaidx pysam pyvcf


R dependencies
--------------

Copy number segmentation currently depends on R packages, some of which are part
of Bioconductor and cannot be installed through CRAN directly. To install these
dependencies, do the following in R::

    > library(BiocManager)
    > install("DNAcopy")

This will install the DNAcopy package, as well as its dependencies.

Alternatively, to do the same directly from the shell, e.g. for automated
installations, try this instead::

    Rscript -e "source('http://callr.org/install#DNAcopy')"


Testing
=======

You can test your installation by running the CNVkit pipeline on the example
files in the ``test/`` directory. The pipeline is implemented as a Makefile and
can be run with the ``make`` command (standard on Unix/Linux/Mac OS X systems)::

    cd test/
    make

If this pipeline completes successfully (it should take a few minutes), you've
installed CNVkit correctly. On a multi-core machine you can parallelize this
with ``make -j``.

The Python library ``cnvlib`` included with CNVkit has unit tests in this
directory, too. Run the test suite with ``make test``.

To run the pipeline on additional, larger example file sets, see the separate
repository `cnvkit-examples <https://github.com/etal/cnvkit-examples>`_.
