from .gary import GenomicArray as _GA
from ._version import __version__

def read(fname):
    """Parse a file as a copy number or copy ratio table (.cnn, .cnr)."""
    return _GA.read(fname)
