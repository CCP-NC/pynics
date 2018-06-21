# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from collections import namedtuple
from pynics.utils import gen_fgrid, fftgrid_gen
import scipy.constants as cnst

NicsResult = namedtuple('NICS', ('nics', 'chi', 'nics_plus_chi'))


class NicsCompute(object):

    _ang2bohr = 1.0/(cnst.physical_constants['Bohr radius'][0]*1e10)
    _ppm = 1e6/((137.036*_ang2bohr)**2.0)

    def __init__(self, current, cell):
        # Grab the current file, compute grid etc.
        self._current = current
        self._cell = cell
        self._V = abs(np.linalg.det(cell))

        self._grid_shape = self._current.grid
        self._grid_npoints = np.prod(self._grid_shape)
        self._fxyz = gen_fgrid(self._grid_shape)
        self._xyz = np.tensordot(cell, self._fxyz, axes=(0, 0))

        self._J = np.rollaxis(np.array(current.current), -1, 0)
        self._J = np.rollaxis(np.array(self._J), -1, 0)
        self._Jfft = np.fft.fftn(self._J,
                                 axes=(2, 3, 4))/self._grid_npoints

        self._chi = np.reshape(current.chi, (3, 3)).T/self._ang2bohr

        self._Ggrid = fftgrid_gen(self._xyz, lattice=self._cell)
        self._Gnorm = np.linalg.norm(self._Ggrid, axis=0)
        # Just fixing this...
        self._Gnorm[0, 0, 0] = np.inf

        ctemp = self._Jfft*(1.0j)*4.0*np.pi/(self._Gnorm**2.0*self._V)
        self._B = -np.cross(self._Ggrid, ctemp, axisa=0, axisb=1, axisc=1)

    def get_nics(self, p0):
        """ Get NICS at point p0
        """

        expp0 = np.exp(1.0j*np.tensordot(p0, self._Ggrid, axes=(0, 0)))

        # Ok, so now... b field!

        ff_ker = self._B*expp0
        nics_tens = np.real(np.sum(ff_ker, axis=(2, 3, 4)))
        # Final value
        nics_tens = (nics_tens*self._ppm).T

        # Now do the version with the chi term
        B_G0 = -8.0/3.0*np.pi*self._chi/self._V
        chi_contrib = np.real(B_G0*self._ppm)
        nics_tens_chi = nics_tens + chi_contrib

        return NicsResult(nics_tens, chi_contrib, nics_tens_chi)

    def get_nics_buildup(self, p0, Rmax=20, Rmin=0, n=100, is_log=False):
        """ Get radial NICS buildup curve for a range of distances"""

        Rmin = max(Rmin, 0)
        Rmax = max(Rmax, 0)

        if is_log:
            # Rmin needs to be greater than 0
            if Rmin == 0:
                Rmin = Rmax/n
            rrange = np.exp(np.linspace(np.log(Rmin), np.log(Rmax), n))
        else:
            rrange = np.linspace(Rmin, Rmax, n)

        expp0 = np.exp(1.0j*np.tensordot(p0, self._Ggrid, axes=(0, 0)))

        Bsph = self._B*expp0[None, None]
        Bsph *= self._ppm

        nics_buildup = []

        for i, r in enumerate(rrange):
            if r == 0:
                j0 = np.ones(self._Gnorm.shape)
            else:
                sin_gr = np.sin(np.where(np.isinf(self._Gnorm), 0,
                                         self._Gnorm)*r)
                j0 = sin_gr/(self._Gnorm*r)

            sigma = np.real(np.sum(Bsph*(1.0-j0)[None, None, :, :, :],
                                   axis=(2, 3, 4)))
            nics_buildup.append(sigma)

        return rrange, np.array(nics_buildup)

    def get_nics_envelope(self, p0):

        expp0 = np.exp(1.0j*np.tensordot(p0, self._Ggrid, axes=(0, 0)))

        Bsph = self._B*expp0[None, None]
        Bsph *= self._ppm

        return np.sum(np.abs(Bsph)*1.0/(self._Gnorm)[None, None, :, :, :],
                      axis=(2, 3, 4))
