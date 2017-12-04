# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from pynics.freeform.io_freeform import (io_freeform_file,
                                         io_freeform_error,
                                         keyword,
                                         io_keyword_error)

# A .nicslist file parser

# First, a fixed list of which keywords to expect
_nics_keywords = [
    keyword("NICS_POINTS_FRAC", "B:B"),
    keyword("NICS_POINTS_ABS", "B:B"),
    keyword("NICS_DECOMPOSITION_FRAC", "B:I"),
    keyword("NICS_DECOMPOSITION_ABS", "B:I"),
    keyword("GRID", "W:B"),
    keyword("GRID_OFFSET", "V:B"),
    keyword("WRITE_FORMATTED_NICS", "L:B"),
    keyword("NICS_TENSOR_HAEB", "L:B"),
    keyword("CURRENT_FILE", "S:B"),
    keyword("WRITE_DECOMPOSITION_SUMMARY", "L:I"),
    keyword("DECOMPOSITION_SUPERCELL_SIDE", "I:I"),
    keyword("READ_ASCII_CURRENT", "L:B"),
    keyword("WRITE_ASCII_CURRENT", "S:B"),
    keyword("USE_WIGNER_SEITZ", "L:I")
]


class nicsfile(object):

    def __init__(self, fname):

        # Create a freeform file with constrained keywords

        self.iof = io_freeform_file(fname=fname, keywords=_nics_keywords)

        # Now parse the variables

        def parse_block(b):
            # Turns a block from a list of lines into a numpy array
            b2 = np.zeros((len(b), 3))
            for i, l in enumerate(b):
                b2[i] = [float(x) for x in l.split()[:3]]
            return b2

        for kw in _nics_keywords:
            setattr(self, kw.key.lower(), None)
            if self.iof.freeform_present(kw.key):
                # Alright, parse
                typ = kw.typ.split(':')[0]
                if typ != 'B':
                    setattr(self, kw.key.lower(),
                            self.iof.freeform_methods[typ](kw.key))
                else:
                    # For blocks, special treatment
                    setattr(self, kw.key.lower(), parse_block(
                        self.iof.freeform_methods[typ](kw.key)))
