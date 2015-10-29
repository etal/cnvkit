"""NGS utilities."""
from __future__ import absolute_import, division, print_function

import os


def is_newer_than(target_fname, orig_fname):
    if not os.path.isfile(target_fname):
        return False
    return (os.stat(target_fname).st_mtime >= os.stat(orig_fname).st_mtime)
