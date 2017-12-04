# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
from pynics.binparse.forbinfile import RecordError


def cbin_results_parse(binfile, results_store, curr_version, params,
                       current_cell, tolerant=False):

    # First, parse the elec results part (always present)
    results_store['elec'] = collections.OrderedDict()
    cbin_elec_parse(binfile, results_store[
                    'elec'], curr_version, params, current_cell)
    # Then go for the optional stuff
    cbin_optional_parse(binfile, results_store, curr_version)


def cbin_elec_parse(binfile, elec_store, curr_version, params, current_cell):

    # A few informations are stored
    elec_store['found_ground_state_wvfn'] = not (
        binfile.read_record('i')[0] == 0)   # Logical value
    elec_store['found_ground_state_den'] = not (
        binfile.read_record('i')[0] == 0)   # Logical value

    elec_store['total_energy'] = binfile.read_record('d')[0]
    elec_store['fermi_energy'] = binfile.read_record('d')

    # Fermi energy for both spins if we have two. This relies on param being
    # already parsed
    if params['nspins'] == 2:
        elec_store['fermi_energy'] = (elec_store['fermi_energy'][
                                      0], elec_store['fermi_energy'][0])

    elec_store['wvfn_nbands'], elec_store[
        'wvfn_nspins'] = binfile.read_record('i')

    # Read occupation eigenvalues for the Kohn-Sham states. This relies on
    # cell being already parsed

    elec_store['occupation'] = {}

    for kp_i in range(0, current_cell['nkpts']):
        kp = binfile.read_record('d')
        elec_store['occupation'][kp] = {'occ': [], 'nrg': []}
        for ns_i in range(0, elec_store['wvfn_nspins']):
            elec_store['occupation'][kp]['occ'].append(
                binfile.read_record('d'))  # Occupation
            elec_store['occupation'][kp]['nrg'].append(
                binfile.read_record('d'))  # Energies

    # Why is this here again? Whatever.
    elec_store['found_ground_state_den'] = not (
        binfile.read_record('i')[0] == 0)    # Logical value

    # Read the fine grid size, keep the information because it is of use for
    # various other parsing operations
    elec_store['model_ngx_fine'], elec_store['model_ngy_fine'], elec_store[
        'model_ngz_fine'] = binfile.read_record('i')

    # Finally, dummy read of density
    for n in range(0, elec_store['model_ngx_fine'] *
                   elec_store['model_ngy_fine']):
        dummy_int = binfile.read_record('i')


def cbin_optional_parse(binfile, results_store, curr_version, tolerant=False):

    if (tolerant):  # In this case, unknown sections will simply be ignored
        def skip_optional():
            while True:
                header = self.binfile.read_string_record()
                if header.isalpha():
                    self.binfile.backspace()
                    break

    try:
        while True:
            header = binfile.read_string_record()

            if (header == 'END'):
                break

            try:
                castep_bin_olist[header](binfile, results_store, curr_version)
            except KeyError:
                if (tolerant):
                    print("Skipping unrecognized header " + header)
                    skip_optional()
                else:
                    # The default case, doesn't account for forward
                    # compatibility for now
                    raise CastepBinError('Unknown optional section found')

    except RecordError:
        raise CastepBinError(
            'End of file reached while parsing optional blocks')

# Utility routine


def tensor_reshape(V):
    return tuple([
        tuple([
            tuple([V[i+j+k] for i in range(0, 3)
                   ]) for j in range(0, 9, 3)
        ]) for k in range(0, len(V), 9)])


def opt_e_fermi_parse(binfile, results_store, curr_version):

    # Parse Fermi energy for second spin
    efermi_2 = binfile.read_record('d')[0]
    results_store['elec']['fermi_energy'] = (
        results_store['elec']['fermi_energy'][0], efermi_2)


def opt_oep_pot_parse(binfile, results_store, curr_version):

    # Parse optimized effective potential
    results_store['oep_pot'] = {}
    results_store['oep_pot']['found_oep_ground_state'] = not (
        binfile.read_record('i') == 0)
    results_store['oep_pot']['oep_energy_difference'] = (binfile
                                                         .read_record('d')[0])

    # We need nspins, we get it indirectly
    nspins = len(results_store['elec']['fermi_energy'])
    ngx_fine = results_store['elec']['model_ngx_fine']
    ngy_fine = results_store['elec']['model_ngy_fine']
    ngz_fine = results_store['elec']['model_ngz_fine']

    # Sort of necessary to initialize here
    results_store['oep_pot']['pot_fine'] = [[0.0 for s in range(
        0, nspins)]]*ngx_fine*ngy_fine*ngz_fine

    for s_i in range(0, nspins):
        for nx1 in range(0, ngx_fine):
            for ny1 in range(0, ngy_fine):
                nx, ny, grid_charge_r, grid_charge_im = binfile.read_record(
                    'iidd')
                # Fortran convention needs to be used
                for nz in range(1, ngz_fine+1):
                    # Back to Python convention, arrays count from 0
                    igrid = (nx-1)+(ny-1)*ngx_fine+(nz-1)*ngx_fine*ngy_fine
                    results_store['oep_pot']['pot_fine'][igrid][
                        s_i] = (grid_charge_r, grid_charge_im)


def opt_de_dloge_parse(binfile, results_store, curr_version):

    # Parse energy logarithmic derivative
    results_store['de_dloge'] = binfile.read_record('d')[0]


def opt_forces_parse(binfile, results_store, curr_version):

    # Parse forces
    f = binfile.read_record('d')
    # Reshape
    results_store['forces'] = tuple(
        [tuple([f[i+j] for i in range(0, 3)]) for j in range(0, len(f), 3)])


def opt_stress_parse(binfile, results_store, curr_version):

    # Parse stress & strain tensors
    results_store['stress'] = {}
    stress = binfile.read_record('d')
    strain = binfile.read_record('d')
    results_store['stress']['stress'] = stress
    results_store['stress']['strain'] = tensor_reshape(strain)


def opt_shielding_parse(binfile, results_store, curr_version):

    # Parse NMR shieldings
    results_store['shielding'] = {}
    results_store['shielding']['ms'] = binfile.read_record('d')
    results_store['shielding']['sus'] = binfile.read_record('d')

    # Reshape ms as a list of tensors
    results_store['shielding']['ms'] = tensor_reshape(
        results_store['shielding']['ms'])


def opt_efg_parse(binfile, results_store, curr_version):

    results_store['efg'] = binfile.read_record('d')
    results_store['efg'] = tensor_reshape(results_store['efg'])

castep_bin_olist = {

    'E_FERMI': opt_e_fermi_parse,
    'OEP_POT': opt_oep_pot_parse,
    'DE_DLOGE': opt_de_dloge_parse,
    'FORCES': opt_forces_parse,
    'STRESS': opt_stress_parse,
    'SHIELDING': opt_shielding_parse,
    'EFG': opt_efg_parse,

}
