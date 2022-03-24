import numpy as np
from ase import io
from soprano.properties.linkage import Molecules
import soprano.collection.generate as scg
from soprano.collection import AtomsCollection
from soprano.selection import AtomSelection
import soprano.utils as soputils
from ase import Atoms
import ase
import argparse as ap

def split_cell(args=None):
	parser = ap.ArgumentParser()
	parser.add_argument('-struct', type=str, default=None, help='Geometry optimized supercell structure to be split.')
	parser.add_argument('-onemol', type=str, default="onemol.cell", help='Output file for onemol structure.')
	parser.add_argument('-nomol', type=str, default="nomol.cell", help='Output file for nomol structure.')
	parser.add_argument('-supercell', type=str, default="supercell.cell", help='Output file for supercell structure.')
	parser.add_argument('-nicslist', type=str, default="nicslist", help='File to output NICS list')
	args = parser.parse_args()
	if (args.struct == None):
		parser.print_help()
		exit(-1)

	struct = io.read(args.struct)
	mols = Molecules.get(struct)
	onemol = mols[0].subset(struct, True)
	nomol = sum(mols[2:], mols[1]).subset(struct, True)
	supercell = onemol + nomol # sum(mols[1:], mols[0]).subset(struct, True)
	onemol_sites = zip(mols[0].indices, onemol.get_chemical_symbols())

	with open(args.nicslist, "w") as f:
		f.write("%BLOCK NICS_POINTS_ABS\n")
		for pos in onemol.get_positions():
			f.write("\t%f\t%f\t%f\n" % tuple(pos))
		f.write("%ENDBLOCK NICS_POINTS_ABS\n")

	io.write(args.supercell, supercell)
	io.write(args.onemol, onemol)
	io.write(args.nomol, nomol)

	print("You should now run Magres calculations on the '%s', '%s', and '%s' files." % (args.supercell, args.onemol, args.nomol))
