#!/usr/bin/env python

"""Copy number variation toolkit: Infer copy number from targeted DNA sequencing."""

from os.path import dirname
from glob import glob

setup_args = {}

try:
    from setuptools import setup
    # Dependencies for easy_install and pip:
    setup_args.update(
        install_requires=[
            'biopython >= 1.62',
            'reportlab >= 2.5',
            'matplotlib >= 1.2',
            'pysam >= 0.7.6',
            'numpy >= 1.5',
        ])
except ImportError:
    from distutils.core import setup


DIR = (dirname(__file__) or '.') + '/'

setup_args.update(
    name='CNVkit',
    version='0.1-dev',
    description=__doc__,
    author='Eric Talevich',
    author_email='eric.talevich@ucsf.edu',
    url='http://github.com/etal/cnvkit',
    packages=['cnvlib', 'cnvlib.segmentation'],
    scripts=[DIR + 'cnvkit.py'] + glob(DIR + 'scripts/*.py'),
)

setup(**setup_args)
