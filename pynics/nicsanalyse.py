import numpy as np
import os
import sys
import re
import argparse as ap
import pynics
from ase import io
from pynics.nics import NicsCompute
from pynics.freeform import nicsfile
from pynics.binparse import CurrentFile
from soprano.nmr.utils import _haeb_sort, _anisotropy, _asymmetry
from soprano.properties.nmr import *
import gc


def get_shift_of_atom_at_point(magres, p):
	''' 
	  function which takes a structure file (magres) and loops through to find
	  the smallest distanced point to point p. it then returns the isotropic shift
	  and distance to this point
	'''
	mdist = 1e3
	isotropic_shift = 0
	shift_positions = zip(MSIsotropy.get(magres), magres.get_positions())
	print(len(magres.get_positions()))
	for shift, pos in shift_positions:
		dist = np.sqrt(np.sum(np.power(pos - p, 2)))
		#print("(%f %f %f) (%f %f %f) %f" % (pos[0], pos[1], pos[2], p[0], p[1], p[2], dist))
		if (dist < mdist):
			mdist = dist
			isotropic_shift = shift
			
	if (mdist > 0.001):
		print("Assignment may be wrong, dist=%f A" % (mdist))
	return isotropic_shift, mdist
	
def get_label_of_atom_at_point(cell, p):
	''' 
	  function which takes a structure file (magres) and loops through to find
	  the smallest distanced point to point p. it then returns the label of this atom
	  of the form "H_1"
	'''
	mdist = 1e3
	rlabel = ""
	labels = cell.get_array('labels')
	indices = cell.get_array('indices')
	labels = ["{0}_{1}".format(l, i) for l, i in zip(labels, indices)]
	label_positions = zip(labels, cell.get_positions())
	for label, pos in label_positions:
		dist = np.sqrt(np.sum(np.power(pos - p, 2)))
		if (dist < mdist):
			mdist = dist
			rlabel = label
	if (mdist > 0.001):
		print("Assignment of %s may be wrong, dist=%f A" % (rlabel, mdist))
	return rlabel
	
#if __name__ == "__main__":
def nics_analyse(args=None):
	parser = ap.ArgumentParser()
	parser.add_argument('-nicslist', type=str, default=None, help='nicslist file, containing sites of interest.')
	#parser.add_argument('-cell', type=str, default=None, help='ONEMOL cell file.')
	parser.add_argument('-nomol_current', type=str, default=None, help='NOMOL current file.')
	parser.add_argument('-nomol_magres', type=str, default=None, help='NOMOL magres file.')
	parser.add_argument('-onemol_current', type=str, default=None, help='ONEMOL current file.')
	parser.add_argument('-onemol_magres', type=str, default=None, help='ONEMOL magres file.')
	parser.add_argument('-supercell_current', type=str, default=None, help='SUPERCELL current file.')
	parser.add_argument('-supercell_magres', type=str, default=None, help='SUPERCELL magres file.')
	parser.add_argument('-output', type=str, default='output.txt', help='File to output data table into.')
	parser.add_argument('-buildup', action='store_true', default=False, help='Calculate buildup curves.')
	parser.add_argument('-R', type=float, default=10, help="Max radius.")
	parser.add_argument('-minR', type=float, default=0, help="Min radius.")
	parser.add_argument('-n', type=int, default=100, help="Number of points.")
	parser.add_argument('-log', action='store_true', default=False,
		                    help="Use a logarithmic grid.")


	args = parser.parse_args()
	if (args.nicslist == None):
		parser.print_help()
		exit(-1)
		
	
	# borrowed from pynics - read magres file to get cell parameters, and read nicsfile
	cell_pars = io.read(args.onemol_magres).get_cell()
	nfile = nicsfile(args.nicslist)
	try:
		fracpoints = np.dot(nfile.nics_points_frac, cell_pars)
	except TypeError:
		fracpoints = None
	except ValueError:
		fracpoints = None
	abspoints = nfile.nics_points_abs
	allpoints = {
		'frac': fracpoints,
		'abs': abspoints
	}
	# set up file parameters
	nomol = {'currentfile':args.nomol_current}
	onemol = {'currentfile':args.onemol_current}
	supercell = {'currentfile':args.supercell_current}
	nomol['magres_file'] = args.nomol_magres
	onemol['magres_file'] = args.onemol_magres
	supercell['magres_file'] = args.supercell_magres
	general_cell = io.read(args.onemol_magres)
	
	# the first step is to go over all the cell files and calculate all NICS parameters
	# this is done initially as the NicsCompute() and CurrentFile() objects are very memory
	# intensive (such that I can't have them all loaded on my computer without it crashing).
	for cell in (nomol, onemol, supercell):
		current = CurrentFile(cell['currentfile'])
		cnics = NicsCompute(current, cell_pars)
		cell['points'] = {}
		cell['magres'] = io.read(cell['magres_file'])
		for ptype, plist in allpoints.items():
			if plist is None:
				continue
			
			for i, p in enumerate(plist):
				# calculate all NICS parameters at this point
				d = {'pos': p}
				d['nics_T'] = cnics.get_nics(p) # NICS tensor
				d['nics'] = {'nics': np.trace(d['nics_T'].nics)/3.0, 'nics+chi': np.trace(d['nics_T'].nics_plus_chi)/3.0}
				d['chi'] = d['nics']['nics+chi'] - d['nics']['nics']
				if (args.buildup == True):
					d['rrange'], d['buildup'] = cnics.get_nics_buildup(p, Rmax=args.R, n=args.n,
													Rmin = args.minR, is_log = args.log)
					all_evals = np.array([np.linalg.eigh((nb+nb.T)/2.0)[0] for nb in d['buildup']])
					all_evals = _haeb_sort(all_evals)
					d['buildup'] = np.average(all_evals, axis=1) + d['chi']
				
				atom_label = get_label_of_atom_at_point(general_cell, p)
				cell['points'][atom_label] = d
			
				
		del current
		del cnics
		gc.collect()
		print("NEXT!")

	
	fp = open(args.output, "w")
	fp.write("# Atom, SC Iso / ppm, Onemol Iso / ppm, Nomol NICS+Chi / ppm, Delta Mol / ppm, NICS / ppm, ES / ppm\n")
	for ptype, plist in allpoints.items():
		if plist is None:
			continue
		for i, p in enumerate(plist):
			label = get_label_of_atom_at_point(onemol['magres'], p)
			a_onemol = onemol['points'][label]
			a_nomol = nomol['points'][label]
			a_supercell = supercell['points'][label]
			
			for cell in (onemol, supercell):
				# because the ordering of atoms will have changed, we can either keep track
				# of the indices. but this might lead to issues. here instead I'm just going	
				# to find the closest position and take this as the atom.
				cell['points'][label]['iso_shift'], d = get_shift_of_atom_at_point(cell['magres'], p)

			## all NICS parameters are in cell['points'][label]
			print("Processing atom %s" % (label))
			delta_mol_cryst = a_onemol['iso_shift'] - a_supercell['iso_shift']
			nics_contrib = a_nomol['nics']['nics+chi'] # may be nics+chi?
			electronic_structure = delta_mol_cryst - nics_contrib
			fp.write("%s, %f, %f, %f, %f, %f, %f\n" % (label, a_supercell['iso_shift'], a_onemol['iso_shift'],
								a_nomol['nics']['nics+chi'], delta_mol_cryst, nics_contrib, electronic_structure))

			if (args.buildup == True):
				import matplotlib.pyplot as plt
				fig, ax = plt.subplots(1,1, sharex='col')
				fig.set_size_inches(6,4)
				ax.spines["top"].set_visible(False)    
				ax.spines["right"].set_visible(False)
				ax.set_xlabel("Radius / $\\AA$")
				ax.set_ylabel("Magnetic Shielding / ppm")
				ax.set_title(label)
				delta_mol_cryst_buildup = a_onemol['buildup'] - a_supercell['buildup']
				nics_buildup = a_nomol['buildup']
				electronic_structure_buildup = delta_mol_cryst_buildup - nics_buildup
				if (args.log):
					ax.set_xscale('log')
				data = np.zeros((len(a_onemol['buildup']), 6))
				data[:, 0] = a_onemol['rrange']
				data[:, 1] = a_onemol['buildup']
				data[:, 2] = a_supercell['buildup']
				data[:, 3] = delta_mol_cryst_buildup
				data[:, 4] = nics_buildup
				data[: ,5] = electronic_structure_buildup
				ax.plot(a_onemol['rrange'], a_onemol['buildup'], label="Onemol")
				ax.plot(a_onemol['rrange'], a_supercell['buildup'], label="Supercell")
				ax.plot(a_onemol['rrange'], delta_mol_cryst_buildup, label="$\\Delta$ molcryst")
				ax.plot(a_onemol['rrange'], nics_buildup, label="$\\sigma_{NICS}$")
				ax.plot(a_onemol['rrange'], electronic_structure_buildup, label="$\\sigma_{ES}$")
				ax.legend()
				np.savetxt(label+".csv", data, delimiter=",", header="time, onemol, supercell, deltamolcryst, nics, electronic structure")
				plt.savefig(label+".png", bbox_inches='tight')
				plt.close()
				
	fp.close()
	
