# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import numpy as np
import argparse as ap

from ase import io
from pynics.nics import NicsCompute
from pynics.freeform import nicsfile
from pynics.binparse import CurrentFile
from soprano.properties.nmr.utils import _haeb_sort, _anisotropy, _asymmetry


def nics_buildup(args=None):

    parser = ap.ArgumentParser()
    parser.add_argument('seedname', type=str, default=None)
    parser.add_argument('-R', type=float, default=10, help="Max radius")
    parser.add_argument('-n', type=int, default=100, help="Number of points")
    parser.add_argument('-cfile', type=str, default='current.dat',
                        help="Name of current file")
    args = parser.parse_args()

    cfile = CurrentFile(args.cfile)
    cell = io.read(args.seedname + '.cell').get_cell()
    nfile = nicsfile(args.seedname + '.nicslist')

    try:
        fracpoints = np.dot(nfile.nics_points_frac, cell.T)
    except TypeError:
        fracpoints = None
    abspoints = nfile.nics_points_abs

    ncomp = NicsCompute(cfile, cell)

    allpoints = {
        'frac': fracpoints,
        'abs': abspoints
    }

    outnics = open(args.seedname + '_nics.txt', 'w')

    for ptype, plist in allpoints.iteritems():
        if plist is None:
            continue
        for i, p in enumerate(plist):
            nics = ncomp.get_nics(p)
            rrange, nbuild = ncomp.get_nics_buildup(p, args.R, args.n)

            # Print output
            outnics.write('Point {0}_{1}:\n'.format(ptype, i+1))
            outnics.write('NICS tensor:\n{0}\n'.format(nics.nics))
            outnics.write('NICS+chi tensor:\n{0}\n'.format(nics.nics_plus_chi))

            # For buildup, let's diagonalize them
            all_evals = np.array([np.linalg.eigh((nb+nb.T)/2.0)[0]
                                  for nb in nbuild])
            all_evals = _haeb_sort(all_evals)
            iso = np.average(all_evals, axis=1)
            aniso = _anisotropy(all_evals)
            asymm = _asymmetry(all_evals)
            np.savetxt(args.seedname + '_{0}_{1}_nicsbuild.dat'.format(ptype,
                                                                       i+1),
                       np.array([rrange, iso, aniso, asymm]).T)

    outnics.close()
