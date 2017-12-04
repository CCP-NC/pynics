# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# Useful functions
import numpy as np


def abc2cart(abc):
    """Transforms an axes and angles representation of lattice parameters 
    into a Cartesian one

    """

    try:
        assert(abc.shape == (2, 3))
    except (AttributeError, AssertionError) as e:
        raise ValueError("Invalid abc passed to abc2cart")

    cart = []
    sin = np.sin(abc[1, :])
    cos = np.cos(abc[1, :])
    cart.append([sin[2], cos[2], 0.0])
    cart.append([0.0, 1.0, 0.0])
    cart.append([(cos[1]-cos[0]*cos[2])/sin[2], cos[0], 0.0])
    cart[2][2] = np.sqrt(1.0-cart[2][0]**2.0-cart[2][1]**2.0)
    cart = np.array([np.array(cart[i])*abc[0, i] for i in range(3)])

    return cart

def gen_fgrid(grid_shape):
    """ Create an equally spaced fractional grid of given shape"""

    axrng = [np.linspace(0.0, 1.0, grid_shape[i]+1)[:-1] for i in range(3)]
    xyz = np.array(np.meshgrid(*axrng, indexing='ij'))
    
    return xyz

def expand_grid_supcell(xyz, r_bounds, scell_offset=[0, 0, 0], lattice=None):
    """ Expand a regular space grid into a supercell one

    """

    scell_grid = supcell_gridgen(r_bounds, scell_offset, lattice)
    tot_grid = (np.repeat(xyz[None], scell_grid.shape[0], axis=0) +
                scell_grid[:, :, None, None, None])

    return tot_grid


def posop_gen(p0, xyz, scell_size, scell_offset=[0, 0, 0], lattice=None):

    grid = expand_grid_supcell(xyz, [scell_size]*3, scell_offset, lattice)

    return grid - np.array(p0)[None, :, None, None, None]


def fftgrid_gen(xyz, lattice=None):

    # FFT grid from xyz grid

    if lattice is None:
        lattice = np.identity(3)
    else:
        lattice = np.array(lattice)

    dx = [np.linalg.norm(lattice[i])/xyz.shape[1+i] for i in range(3)]
    inv_latt = np.linalg.inv(lattice.T)*2*np.pi

    fft_grid = np.array(np.meshgrid(*[np.fft.fftfreq(xyz.shape[i+1])
                                      * xyz.shape[i+1]
                                      for i in range(3)], indexing='ij'))
    # g-vector convention in formulas used
    fft_grid = np.tensordot(inv_latt, fft_grid, axes=(0, 0))

    return fft_grid
