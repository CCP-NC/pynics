# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


def cbin_cell_parse(binfile, cell_store, curr_version):

    # Wrapper for the two parts of the parsing
    cbin_cell_parse_header(binfile, cell_store, curr_version)
    cbin_cell_parse_main(binfile, cell_store, curr_version)


def cbin_cell_parse_header(binfile, cell_store, curr_version):

    # Parse the first part of the cell file

    # We start out with a series of values in the format header -> array

    typelist = []
    typelist.append(('CELL%REAL_LATTICE',           'd'))
    typelist.append(('CELL%RECIP_LATTICE',          'd'))
    typelist.append(('CELL%VOLUME',                 'd'))
    typelist.append(('CELL%NUM_SPECIES',            'i'))
    typelist.append(('CELL%NUM_IONS',               'i'))
    typelist.append(('CELL%MAX_IONS_IN_SPECIES',    'i'))

    for t in typelist:
        header = binfile.read_string_record()
        if (header != t[0]):
            raise ValueError('Corrupted CELL block read from binfile')
        tname = t[0].replace('CELL%', '').lower()
        cell_store[tname] = binfile.read_record(t[1])
        if len(cell_store[tname]) == 1:
            cell_store[tname] = cell_store[tname][0]

    # A couple more things...
    reshape3_action(binfile, cell_store, curr_version, 'real_lattice')
    reshape3_action(binfile, cell_store, curr_version, 'recip_lattice')


def cbin_cell_parse_main(binfile, cell_store, curr_version):

    # Parse the second part of the cell file, where order isn't a given any more

    header = None

    while True:
        # First, read a header
        header = binfile.read_string_record()

        if header == 'END_UNIT_CELL':
            break

        try:
            c = castep_bin_clist[header]
        except KeyError:
            raise ValueError("Unrecognizable keyword {0} ".format(header) +
                             "found in CELL block of binfile")

        if c["type"] == 'str':
            cell_store[c["name"]] = binfile.read_string_record()
        elif c["type"] == 'bool':
            cell_store[c["name"]] = (binfile.read_record('i')[0] != 0)
        else:
            cell_store[c["name"]] = binfile.read_record(c["type"])
            if len(cell_store[c["name"]]) == 1:
                cell_store[c["name"]] = cell_store[c["name"]][0]
        if c["action"] is not None:
            c["action"](binfile, cell_store, curr_version, c["name"])

    # Now on to check that the correct block begins
    header = binfile.read_string_record()
    if (header != 'BEGIN_CELL_GLOBAL'):
        raise ValueError("Corrupted CELL block in binfile")

    while True:
        # First, read a header
        header = binfile.read_string_record()

        if header == 'END_CELL_GLOBAL':
            break

        try:
            c = castep_bin_global_clist[header]
        except KeyError:
            raise ValueError(
                "Unrecognizable keyword found in CELL block of binfile")

        if c["type"] == 'str':
            cell_store[c["name"]] = binfile.read_string_record()
        elif c["type"] == 'bool':
            cell_store[c["name"]] = (binfile.read_record('i')[0] != 0)
        else:
            cell_store[c["name"]] = binfile.read_record(c["type"])
            if len(cell_store[c["name"]]) == 1:
                cell_store[c["name"]] = cell_store[c["name"]][0]
        if c["action"] is not None:
            c["action"](binfile, cell_store, curr_version, c["name"])


def int2float_action(binfile, cell_store, curr_version, cname):
    # An action to convert a vector from int to float.
    # Used for ionic_charge
    try:
        cell_store[cname] = tuple([float(x) for x in cell_store[cname]])
    except TypeError:   # It's a single value
        cell_store[cname] = float(cell_store[cname])


def strsplit_action(binfile, cell_store, curr_version, cname):
        # An action to split a long space-separated string in a tuple of strings
        # Used for species_symbol, species_pot
    cell_store[cname] = tuple([s.strip() for s in cell_store[cname].split()])


def reshape3_action(binfile, cell_store, curr_version, cname):
    # An action used to reshape a one-dimensional array in a two dimensional matrix with 3-rows
    # Used for parameters containing coordinates like real_lattice, ionic_positions etc.
    cell_store[cname] = tuple([tuple([cell_store[cname][i+j] for i in range(0, 3)])
                               for j in range(0, len(cell_store[cname]), 3)])


def reshape4_action(binfile, cell_store, curr_version, cname):
    # An action used to reshape a one-dimensional array in a two dimensional matrix with 4-rows
    # Used for parameters containing orbital-related parameters like hubbard_u
    cell_store[cname] = tuple([tuple([cell_store[cname][i+j] for i in range(0, 4)])
                               for j in range(0, len(cell_store[cname]), 4)])


def reshape3by3_action(binfile, cell_store, curr_version, cname):
    # An action used to reshape a one-dimensional array in a three dimensional matrix with 3x3 slices
    # Used for parameters containing symmetry operations
    cell_store[cname] = tuple([tuple([tuple([cell_store[cname][i+j+k] for i in range(0, 3)])
                                      for j in range(0, 9, 3)]) for k in range(0, len(cell_store[cname]), 9)])


castep_bin_clist = {

    "CELL%ATOM_INIT_MAGMOM": {
        "action": None,
        "name": "atom_init_magmom",
        "type": "bool"
    },
    "CELL%ATOM_INIT_SPIN": {
        "action": None,
        "name": "atom_init_spin",
        "type": "bool"
    },
    "CELL%ATOM_MOVE": {
        "action": None,
        "name": "atom_move",
        "type": "bool"
    },
    "CELL%CHEMICAL_POTENTIAL": {
        "action": None,
        "name": "chemical_potential",
        "type": "d"
    },
    "CELL%HUBBARD_ALPHA": {
        "action": reshape4_action,
        "name": "hubbard_alpha",
        "type": "d"
    },
    "CELL%HUBBARD_U": {
        "action": reshape4_action,
        "name": "hubbard_u",
        "type": "d"
    },
    "CELL%INDX_MIXTURE_ATOMS": {
        "action": None,
        "name": "indx_mixture_atoms",
        "type": "i"
    },
    "CELL%INITIAL_MAGNETIC_MOMENT": {
        "action": None,
        "name": "initial_magnetic_moment",
        "type": "d"
    },
    "CELL%IONIC_CHARGE": {
        "action": int2float_action,
        "name": "ionic_charge",
        "type": "i"
    },
    "CELL%IONIC_CHARGE_REAL": {
        "action": None,
        "name": "ionic_charge",
        "type": "d"
    },
    "CELL%IONIC_POSITIONS": {
        "action": reshape3_action,
        "name": "ionic_positions",
        "type": "d"
    },
    "CELL%IONIC_VELOCITIES": {
        "action": reshape3_action,
        "name": "ionic_velocities",
        "type": "d"
    },
    "CELL%ION_PACK_INDEX": {
        "action": None,
        "name": "ion_pack_index",
        "type": "i"
    },
    "CELL%ION_PACK_SPECIES": {
        "action": None,
        "name": "ion_pack_species",
        "type": "i"
    },
    "CELL%LIST_MIXTURE_ATOMS": {
        "action": None,
        "name": "list_mixture_atoms",
        "type": "i"
    },
    "CELL%LIST_MIXTURE_ATOMS_SIZE": {
        "action": None,
        "name": "max_ion_mix_and_num_mix",
        "type": "i"
    },
    "CELL%MIXTURE_WEIGHT": {
        "action": None,
        "name": "mixture_weight",
        "type": "d"
    },
    "CELL%NUM_IONS_IN_SPECIES": {
        "action": None,
        "name": "num_ions_in_species",
        "type": "i"
    },
    "CELL%NUM_MIXTURE_ATOMS": {
        "action": None,
        "name": "num_mixture_atoms",
        "type": "i"
    },
    "CELL%SEDC_CUSTOM_PARAMS": {
        "action": None,
        "name": "sedc_custom_params",
        "type": "str"
    },
    "CELL%SEDC_CUSTOM_PARAMS_LEN": {
        "action": None,
        "name": "len_temp",
        "type": "i"
    },
    "CELL%SPECIES_EFG_ISOTOPE": {
        "action": None,
        "name": "species_efg_isotope",
        "type": "i"
    },
    "CELL%SPECIES_GAMMA": {
        "action": None,
        "name": "species_gamma",
        "type": "d"
    },
    "CELL%SPECIES_LCAO_STATES": {
        "action": None,
        "name": "species_lcao_states",
        "type": "i"
    },
    "CELL%SPECIES_MASS": {
        "action": None,
        "name": "species_mass",
        "type": "d"
    },
    "CELL%SPECIES_POT": {
        "action": strsplit_action,
        "name": "species_pot",
        "type": "str"
    },
    "CELL%SPECIES_Q": {
        "action": None,
        "name": "species_q",
        "type": "d"
    },
    "CELL%SPECIES_SPIN_ISOTOPE": {
        "action": None,
        "name": "species_spin_isotope",
        "type": "i"
    },
    "CELL%SPECIES_SYMBOL": {
        "action": strsplit_action,
        "name": "species_symbol",
        "type": "str"
    },
    "CELL%SUPERCELL_CELL_INDEX": {
        "action": None,
        "name": "supercell_cell_index",
        "type": "i"
    },
    "CELL%SUPERCELL_ION_INDEX": {
        "action": None,
        "name": "supercell_ion_index",
        "type": "i"
    },
    "CELL%SUPERCELL_MATRIX": {
        "action": None,
        "name": "supercell_matrix",
        "type": "i"
    },
    "CELL_VERSION_NUMBER": {
        "action": None,
        "name": "cell_version_number",
        "type": "d"
    }
}

castep_bin_global_clist = {
    "BS_KPOINTS": {
        "action": None,
        "name": "bs_kpoints",
        "type": "d"
    },
    "BS_KPOINT_WEIGHTS": {
        "action": None,
        "name": "bs_kpoint_weights",
        "type": "d"
    },
    "CELL_CONSTRAINTS": {
        "action": None,
        "name": "cell_constraints",
        "type": "i"
    },
    "CELL_SYMMETRY_ION_CONSTRAINTS": {
        "action": None,
        "name": "cell_symmetry_ion_constraints",
        "type": "i"
    },
    "CELL_SYMMORPHIC": {
        "action": None,
        "name": "cell_symmorphic",
        "type": "bool"
    },
    "CRYSTAL_SYMMETRY_DISPS": {
        "action": None,
        "name": "crystal_symmetry_disps",
        "type": "d"
    },
    "CRYSTAL_SYMMETRY_EQUIV_ATOMS": {
        "action": None,
        "name": "crystal_symmetry_equiv_atoms",
        "type": "i"
    },
    "CRYSTAL_SYMMETRY_OPERATIONS": {
        "action": reshape3by3_action,
        "name": "crystal_symmetry_operations",
        "type": "d"
    },
    "ELNES_KPOINTS": {
        "action": None,
        "name": "elnes_kpoints",
        "type": "d"
    },
    "ELNES_KPOINT_WEIGHTS": {
        "action": None,
        "name": "elnes_kpoint_weights",
        "type": "d"
    },
    "EXTERNAL_EFIELD": {
        "action": None,
        "name": "external_efield",
        "type": "d"
    },
    "EXTERNAL_PRESSURE": {
        "action": None,
        "name": "external_pressure",
        "type": "d"
    },
    "FIX_ALL_CELL": {
        "action": None,
        "name": "fix_all_cell",
        "type": "bool"
    },
    "FIX_ALL_IONS": {
        "action": None,
        "name": "fix_all_ions",
        "type": "bool"
    },
    "FIX_COM": {
        "action": None,
        "name": "fix_com",
        "type": "bool"
    },
    "IONIC_CONSTRAINTS": {
        "action": reshape3_action,
        "name": "ionic_constraints",
        "type": "d"
    },
    "KPOINTS": {
        "action": None,
        "name": "kpoints",
        "type": "d"
    },
    "KPOINT_MP_DENSITY": {
        "action": None,
        "name": "kpoint_mp_density",
        "type": "d"
    },
    "KPOINT_MP_GRID": {
        "action": None,
        "name": "scf_kpoint_mp_grid",
        "type": "i"
    },
    "KPOINT_MP_OFFSET": {
        "action": None,
        "name": "scf_kpoint_mp_offset",
        "type": "d"
    },
    "KPOINT_WEIGHTS": {
        "action": None,
        "name": "kpoint_weights",
        "type": "d"
    },
    "MAGRES_KPOINTS": {
        "action": None,
        "name": "magres_kpoints",
        "type": "d"
    },
    "MAGRES_KPOINT_WEIGHTS": {
        "action": None,
        "name": "magres_kpoint_weights",
        "type": "d"
    },
    "NKPTS": {
        "action": None,
        "name": "nkpts",
        "type": "i"
    },
    "NUM_BS_KPOINTS": {
        "action": None,
        "name": "num_bs_kpoints",
        "type": "i"
    },
    "NUM_CELL_CONSTRAINTS": {
        "action": None,
        "name": "num_cell_constraints",
        "type": "i"
    },
    "NUM_CRYSTAL_SYMMETRY_OPERATIONS": {
        "action": None,
        "name": "num_crystal_symmetry_operations",
        "type": "i"
    },
    "NUM_ELNES_KPOINTS": {
        "action": None,
        "name": "num_elnes_kpoints",
        "type": "i"
    },
    "NUM_IONIC_CONSTRAINTS": {
        "action": None,
        "name": "num_ionic_constraints",
        "type": "i"
    },
    "NUM_MAGRES_KPOINTS": {
        "action": None,
        "name": "num_magres_kpoints",
        "type": "i"
    },
    "NUM_OPTICS_KPOINTS": {
        "action": None,
        "name": "num_optics_kpoints",
        "type": "i"
    },
    "NUM_PHONON_FINE_KPOINTS": {
        "action": None,
        "name": "num_phonon_fine_kpoints",
        "type": "i"
    },
    "NUM_PHONON_GAMMA_DIRECTIONS": {
        "action": None,
        "name": "num_phonon_gamma_directions",
        "type": "i"
    },
    "NUM_PHONON_KPOINTS": {
        "action": None,
        "name": "num_phonon_kpoints",
        "type": "i"
    },
    "NUM_SCF_KPOINTS": {
        "action": None,
        "name": "num_scf_kpoints",
        "type": "i"
    },
    "NUM_SPECTRAL_KPOINTS": {
        "action": None,
        "name": "num_spectral_kpoints",
        "type": "i"
    },
    "NUM_SYMMETRY_OPERATIONS": {
        "action": None,
        "name": "num_symmetry_operations",
        "type": "i"
    },
    "OPTICS_KPOINTS": {
        "action": None,
        "name": "optics_kpoints",
        "type": "d"
    },
    "OPTICS_KPOINT_MP_DENSITY": {
        "action": None,
        "name": "optics_kpoint_mp_density",
        "type": "d"
    },
    "OPTICS_KPOINT_MP_GRID": {
        "action": None,
        "name": "optics_kpoint_mp_grid",
        "type": "i"
    },
    "OPTICS_KPOINT_MP_OFFSET": {
        "action": None,
        "name": "optics_kpoint_mp_offset",
        "type": "d"
    },
    "OPTICS_KPOINT_WEIGHTS": {
        "action": None,
        "name": "optics_kpoint_weights",
        "type": "d"
    },
    "PHONON_FINE_KPOINTS": {
        "action": None,
        "name": "phonon_fine_kpoints",
        "type": "d"
    },
    "PHONON_FINE_KPOINT_WEIGHTS": {
        "action": None,
        "name": "phonon_fine_kpoint_weights",
        "type": "d"
    },
    "PHONON_GAMMA_DIRECTIONS": {
        "action": None,
        "name": "phonon_gamma_directions",
        "type": "d"
    },
    "PHONON_KPOINTS": {
        "action": None,
        "name": "phonon_kpoints",
        "type": "d"
    },
    "PHONON_KPOINT_WEIGHTS": {
        "action": None,
        "name": "phonon_kpoint_weights",
        "type": "d"
    },
    "QUANTISATION_AXIS": {
        "action": None,
        "name": "quantisation_axis",
        "type": "d"
    },
    "SCF_KPOINTS": {
        "action": None,
        "name": "scf_kpoints",
        "type": "d"
    },
    "SCF_KPOINT_WEIGHTS": {
        "action": None,
        "name": "scf_kpoint_weights",
        "type": "d"
    },
    "SPECTRAL_KPOINTS": {
        "action": None,
        "name": "spectral_kpoints",
        "type": "d"
    },
    "SPECTRAL_KPOINT_WEIGHTS": {
        "action": None,
        "name": "spectral_kpoint_weights",
        "type": "d"
    },
    "SYMMETRY_DISPS": {
        "action": None,
        "name": "symmetry_disps",
        "type": "d"
    },
    "SYMMETRY_EQUIV_ATOMS": {
        "action": None,
        "name": "symmetry_equiv_atoms",
        "type": "i"
    },
    "SYMMETRY_GENERATE": {
        "action": None,
        "name": "symmetry_generate",
        "type": "bool"
    },
    "SYMMETRY_OPERATIONS": {
        "action": reshape3by3_action,
        "name": "symmetry_operations",
        "type": "d"
    },
    "SYMMETRY_TOL": {
        "action": None,
        "name": "symmetry_tol",
        "type": "d"
    }
}
