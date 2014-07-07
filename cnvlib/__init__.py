from .cnarray import CopyNumArray as _CNA

def read(fname):
    """Parse a file as a copy number or copy ratio table (.cnn, .cnr)."""
    return _CNA.read(fname)

