from .cnary import CopyNumArray as _CNA
from ._version import __version__

def read(fname):
    """Parse a file as a copy number or copy ratio table (.cnn, .cnr)."""
    return _CNA.read(fname)
