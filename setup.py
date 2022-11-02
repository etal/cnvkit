#!/usr/bin/env python
"""Copy number variation toolkit for high-throughput sequencing."""

from os.path import dirname, join
from glob import glob

from setuptools import setup

setup_args = {}

# Dependencies for easy_install and pip:
install_requires = [
    'biopython >= 1.62',
    'matplotlib >= 1.3.1',
    'numpy >= 1.9',
    'pandas >= 0.23.3',
    'pomegranate >= 0.9.0',
    'pyfaidx >= 0.4.7',
    'pysam >= 0.10.0',
    'reportlab >= 3.0',
    'scikit-learn',
    'scipy >= 0.15.0',

    # TODO: Similarly: https://github.com/etal/cnvkit/issues/589
    'joblib < 1.0',
]

DIR = (dirname(__file__) or '.')
with open(join(DIR, 'cnvlib', '_version.py')) as handle:
    VERSION = handle.readline().split('=')[-1].strip().replace('"','')

setup_args.update(
    name='CNVkit',
    version=VERSION,
    description=__doc__,
    author='Eric Talevich',
    author_email='eric.talevich@ucsf.edu',
    url='http://github.com/etal/cnvkit',
    packages=[
        'cnvlib',
        'cnvlib.segmentation',
        'skgenome',
        'skgenome.tabio',
    ],
    scripts=[join(DIR, 'cnvkit.py')] + glob(join(DIR, 'scripts/*.py')),
    install_requires=install_requires,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)

setup(**setup_args)
