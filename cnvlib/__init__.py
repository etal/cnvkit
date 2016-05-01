from .cnary import CopyNumArray as _CNA
from ._version import __version__
from . import tabio

def read(fname):
    """Parse a file as a copy number or copy ratio table (.cnn, .cnr)."""
    return tabio.read(fname, into=_CNA)
