#!/usr/bin/env python

import os
import sys
from ase import io
import numpy as np
import argparse as ap
from scipy import constants as cnst
from philogenlib.binparse.currentfile import CurrentFile
from philogenlib.pycell.nicsfile import nicsfile

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from lorentz_sphere import *

parser = ap.ArgumentParser()
parser.add_argument("seedname", type=str, default=None)
parser.add_argument("-R", type=float, default=20, help="Max radius")
parser.add_argument("-n", type=int, default=100, help="Number of points")
args = parser.parse_args()

cfile = CurrentFile('current.dat')
cellfile = io.read(args.seedname + '.cell')
latt = cellfile.get_cell()

# 1. Define parameters

grid_shape = cfile.grid
grid_vol = np.prod(grid_shape)
nfile = nicsfile(args.seedname + '.nicslist')
# Analyse the first fractional point in NICSLIST

axrng = [np.linspace(0.0, 1.0, grid_shape[i]+1)[:-1] for i in range(3)]
xyz = np.array(np.meshgrid(*axrng, indexing='ij'))

J = np.rollaxis(np.array(cfile.current), -1, 0)
J = np.rollaxis(np.array(J), -1, 0)

print "Current sum condition: "
print abs(np.average(J)), " ~= 0"
print "Current grid: ", ' '.join([str(x) for x in grid_shape])

fracpoints = np.dot(nfile.nics_points_frac, latt.T)
abspoints = nfile.nics_points_abs

allpoints = {'frac': fracpoints, 'abs': abspoints}

for ptype, plist in allpoints.iteritems():
    print "Solving for {0} points".format(ptype)

    if plist is None:
        continue

    for pind, p0 in enumerate(plist):

        print "Solving in inverse space"

        # FFT version
        # FFT grid
        ang2bohr = 1.0/(cnst.physical_constants['Bohr radius'][0]*1e10)

        fft_grid = fftgrid_gen(xyz, lattice=latt)
        fft_gnorm = np.linalg.norm(fft_grid, axis=0)
        # Just fixing this...
        fft_gnorm[0, 0, 0] = np.inf

        expp0 = np.exp(1.0j*np.tensordot(p0, fft_grid, axes=(0, 0)))

        cell_vol = np.linalg.det(latt)
        current_fft = np.fft.fftn(J, axes=(2, 3, 4))/np.prod(grid_shape)
        # Ok, so now... b field!
        ctemp = current_fft*(1.0j)*4.0*np.pi/(fft_gnorm**2.0*cell_vol)
        B_field = np.cross(fft_grid, ctemp, axisa=0, axisb=1, axisc=1)
        #B_field_dir = np.fft.ifftn(B_field, axes=(2,3,4))

        ff_ker = np.real(
            (B_field*np.exp(1.0j*np.tensordot(p0, fft_grid, axes=(0, 0)))))
        nics_tens = np.sum(ff_ker, axis=(2, 3, 4))
        # Final value
        nics_tens = (1.0*nics_tens/((137.036*ang2bohr)**2.0)*1e6).T

        # Now do the version with the chi term
        chi = np.reshape(cfile.chi, (3, 3)).T/ang2bohr
        B_G0 = 8.0/3.0*np.pi*chi/cell_vol
        B_field[:, :, 0, 0, 0] = B_G0
        ff_ker2 = np.real(B_field*expp0)
        nics_tens_chi = np.sum(ff_ker2, axis=(2, 3, 4))
        # Final value
        nics_tens_chi = (1.0*nics_tens_chi/((137.036*ang2bohr)**2.0)*1e6).T

        print "Isotropic shift: {0}".format(np.trace(nics_tens)/3.0)
        print "Isotropic shift with susceptibility: {0}".format(np.trace(nics_tens_chi)/3.0)

        # Now build up curve!

        rrange = np.linspace(0, args.R, args.n)

        # Spherical harmonics decomposition
        fft_theta = np.arccos(fft_grid[2]/fft_gnorm)
        fft_phi = np.arctan2(fft_grid[1], fft_grid[0])

        Bsph = B_field*expp0[None, None]
        Bsph *= 1e6/((137.036*ang2bohr)**2.0)
        Bsph[:,:,0,0,0] = 0 # Remove susceptibility

        with open('{0}_{1}{2}_buildup.dat'.format(args.seedname, ptype, 
                                                  pind+1),
                  'w') as outf:
            for i, r in enumerate(rrange):
                if r == 0:
                    j0 = np.ones(fft_gnorm.shape)
                else:
                    sin_gr = np.sin(np.where(np.isinf(fft_gnorm), 0, 
                                             fft_gnorm)*r)
                    j0 = sin_gr/(fft_gnorm*r)
                
                sigma = np.real(np.sum(Bsph*(1.0-j0)[None,None,:,:,:], 
                                axis=(2,3,4)))
                outf.write('{0}\t{1}\n'.format(r, np.trace(sigma)/3.0))
