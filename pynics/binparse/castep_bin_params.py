# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import math
from pynics.binparse.forbinfile import RecordError


class CastepBinError(Exception):
    pass


def cbin_parameters_parse(binfile, param_store, curr_version):

    # Parse a series of parameters from binfile, given header and version and
    # put them into param_store
    while True:

        header = binfile.read_string_record().replace('BEGIN_', '')

        if (header == 'END_PARAMETERS_DUMP'):
            break
        try:
            # First, check that the header exists
            if header not in castep_bin_plist:
                # The default case, doesn't account for forward compatibility
                # for now
                raise CastepBinError('Unknown parameter section found')

            # Then list the versions
            header_versions = castep_bin_plist[header].keys()
            header_versions = [float(v) for v in header_versions]
            header_versions.sort()

            # Now on to actual parsing of stuff. Order is all important. Start
            # with the 2.0 stuff, which is everywhere
            for v in header_versions:
                hv = castep_bin_plist[header][str(v)]
                if curr_version >= v:
                    for p in hv["params"]:
                        if p['type'] == 'str':
                            param_store[
                                p['name']] = binfile.read_string_record()
                        elif p['type'] == 'bool':
                            param_store[p['name']] = (
                                binfile.read_record('i')[0] != 0)
                        else:
                            param_store[p['name']] = (binfile
                                                      .read_record(p[
                                                          'type']))
                            if len(param_store[p['name']]) == 1:
                                param_store[p['name']] = param_store[
                                    p['name']][0]
                        print(p['name'])
                        print(param_store[p['name']])
                else:
                    try:
                        hv["else"](binfile, param_store)
                    except KeyError:
                        # There's no else
                        continue

            # ...and finally read the end header
            end_header = binfile.read_string_record()
            if (end_header != 'END_' + header):
                print(header)
                print(end_header)
                raise CastepBinError(
                    'Corrupted ' + header + ' block read from binfile')

        except RecordError:
            raise CastepBinError(
                'End of file reached while parsing parameters')

# Now the various 'else' functions


def basis_v4_0_else(binfile, param_store):
    # Need to set a dependent default value for fine_grid_scale
    param_store['fine_grid_scale'] = param_store['grid_scale']


def basis_v4_4_else(binfile, param_store):
    # Need to set a compatibility value
    param_store['fft_max_prime_factor'] = 5


def bs_v3_0_else(binfile, param_store):
    # Need to set a dependent default value for bs_xc_functional
    param_store['bs_xc_functional'] = param_store['xc_functional']


def devel_v8_0_else(binfile, param_store):
    # Not really useful here, in CASTEP it's mostly Fortran dynamic allocation
    # stuff
    param_store['devel_code'] = binfile.read_string_record()


def elec_min_v4_1_else(binfile, param_store):
    # Need to set a dependent default value for elec_force_tol
    param_store['elec_force_tol'] = math.sqrt(param_store['elec_energy_tol'])


def general_v3_0_else(binfile, param_store):
    # Need to set a dependent default value for opt_strategy_bias
    param_store['opt_strategy_bias'] = {
        'SPEED':    3,
        'MEMORY': -3,
        'DEFAULT':  0,
    }[param_store['opt_strategy']]


def general_v3_1_else(binfile, param_store):
    # Need to set a dependent default value for num_proc_in_smp
    param_store['num_proc_in_smp'] = 1


def general_v3_2_else(binfile, param_store):
    # Need to set a dependent default value for num_proc_in_smp_fine
    param_store['num_proc_in_smp_fine'] = param_store['num_proc_in_smp']
    param_store['requested_seed'] = param_store['rand_seed']


def md_v3_0_else(binfile, param_store):
    # Need to set dependent defaults
    param_store['md_ion_t'] = {
        'NOSE-HOOVER':  10.0*param_store['md_delta_t'],
        'HOOVER-LANG':  10.0*param_store['md_delta_t'],
        'LANGEVIN': 100.0*param_store['md_delta_t']
    }[param_store['md_thermostat']]
    param_store['md_cell_t'] = 10.0*param_store['md_ion_t']


def md_v3_1_else(binfile, param_store):
    # Need to set defaults
    param_store['md_num_beads'] = 1


def md_v4_1_else(binfile, param_store):
    # Need to set dependent defaults
    param_store['md_eqm_ion_t'] = param_store['md_ion_t']
    param_store['md_eqm_cell_t'] = param_store['md_cell_t']
    param_store['md_elec_force_tol'] = param_store['elec_force_tol']


def optics_v3_0_else(binfile, param_store):
    # Need to set dependent defaults
    param_store['optics_xc_functional'] = param_store['xc_functional']


def phonon_v3_0_else(binfile, param_store):
    # Need to set dependent defaults
    param_store['phonon_method'] = 'LINEARRESPONSE'


def phonon_v4_4_else(binfile, param_store):
    param_store['phonon_sum_rule'] = False


def phonon_v5_0_else(binfile, param_store):
    # Special case as default changed when rereading 4.4 check. Value 0.0
    # makes no sense in >= 5.0
    if (abs(param_store['phonon_force_constant_cut_scale']) < 1.0e-4):
        param_store['phonon_force_constant_cut_scale'] = 1.0


def phonon_v6_0_else(binfile, param_store):
    # Need to set phonon_dfpt_method as not present in old checkfiles\n
    # Assumed value in older versions
    param_store['phonon_dfpt_method'] = 'ALLBANDS'


def thermodynamics_v3_2_else(binfile, param_store):
    # Need to set dependent defaults
    param_store['thermo_t_spacing'] = param_store[
        'thermo_t_stop'] - param_store['thermo_t_start']


castep_bin_plist = {
    "BASIS": {
        "2.0": {
            "params": [
                {
                    "name": "basis_precision",
                    "type": "str"
                },
                {
                    "name": "cut_off_energy",
                    "type": "d"
                },
                {
                    "name": "grid_scale",
                    "type": "d"
                },
                {
                    "name": "fine_gmax",
                    "type": "d"
                },
                {
                    "name": "finite_basis_corr",
                    "type": "i"
                },
                {
                    "name": "basis_de_dloge",
                    "type": "d"
                },
                {
                    "name": "finite_basis_npoints",
                    "type": "i"
                },
                {
                    "name": "finite_basis_spacing",
                    "type": "d"
                }
            ]
        },
        "3.0": {
            "params": [
                {
                    "name": "fixed_npw",
                    "type": "bool"
                }
            ]
        },
        "4.0": {
            "else": basis_v4_0_else,
            "params": [
                {
                    "name": "fine_grid_scale",
                    "type": "d"
                }
            ]
        },
        "4.4": {
            "else": basis_v4_4_else,
            "params": [
                {
                    "name": "fft_max_prime_factor",
                    "type": "i"
                }
            ]
        }
    },
    "BS": {
        "2.0": {
            "params": [
                {
                    "name": "bs_max_iter",
                    "type": "i"
                },
                {
                    "name": "bs_max_CG_steps",
                    "type": "i"
                },
                {
                    "name": "bs_nbands",
                    "type": "i"
                },
                {
                    "name": "bs_eigenvalue_tol",
                    "type": "d"
                }
            ]
        },
        "3.0": {
            "else": bs_v3_0_else,
            "params": [
                {
                    "name": "bs_xc_functional",
                    "type": "str"
                },
                {
                    "name": "bs_re_est_k_scrn",
                    "type": "bool"
                }
            ]
        },
        "4.4": {
            "params": [
                {
                    "name": "bs_write_eigenvalues",
                    "type": "bool"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "bs_xc_definition_size",
                    "type": "i"
                },
                {
                    "name": "bs_xc_definition",
                    "type": "str"
                }
            ]
        }
    },
    "DEVEL": {
        "2.0": {
            "params": []
        },
        "8.0": {
            "else": devel_v8_0_else,
            "params": [
                {
                    "name": "devel_code",
                    "type": "str"
                },
                {
                    "name": "block_len",
                    "type": "i"
                },
                {
                    "name": "devel_block",
                    "type": "str"
                }
            ]
        }
    },
    "DM": {
        "2.0": {
            "params": [
                {
                    "name": "mixing_scheme",
                    "type": "str"
                },
                {
                    "name": "mix_history_length",
                    "type": "i"
                },
                {
                    "name": "mix_charge_amp",
                    "type": "d"
                },
                {
                    "name": "mix_charge_gmax",
                    "type": "d"
                },
                {
                    "name": "mix_spin_amp",
                    "type": "d"
                },
                {
                    "name": "mix_spin_gmax",
                    "type": "d"
                },
                {
                    "name": "mix_cut_off_energy",
                    "type": "d"
                },
                {
                    "name": "mix_metric_q",
                    "type": "d"
                }
            ]
        }
    },
    "EFIELD": {
        "2.0": {
            "params": [
                {
                    "name": "efield_max_CG_steps",
                    "type": "i"
                },
                {
                    "name": "efield_max_cycles",
                    "type": "i"
                },
                {
                    "name": "efield_convergence_win",
                    "type": "i"
                },
                {
                    "name": "efield_energy_tol",
                    "type": "d"
                },
                {
                    "name": "efield_calc_ion_permittivity",
                    "type": "bool"
                },
                {
                    "name": "efield_ignore_molec_modes",
                    "type": "str"
                },
                {
                    "name": "excited_state_scissors",
                    "type": "d"
                }
            ]
        },
        "4.0": {
            "params": [
                {
                    "name": "efield_freq_spacing",
                    "type": "d"
                },
                {
                    "name": "efield_oscillator_Q",
                    "type": "d"
                }
            ]
        },
        "6.0": {
            "params": [
                {
                    "name": "efield_dfpt_method",
                    "type": "str"
                }
            ]
        }
    },
    "ELECTRONIC": {
        "2.0": {
            "params": [
                {
                    "name": "dummy_int",
                    "type": "i"
                },
                {
                    "name": "dummy_int",
                    "type": "i"
                },
                {
                    "name": "dummy_int",
                    "type": "i"
                },
                {
                    "name": "nspins",
                    "type": "i"
                },
                {
                    "name": "nbands",
                    "type": "i"
                },
                {
                    "name": "elec_temp",
                    "type": "d"
                }
            ]
        },
        "3.2": {
            "params": [
                {
                    "name": "dummy_int",
                    "type": "i"
                },
                {
                    "name": "dummy_int",
                    "type": "i"
                },
                {
                    "name": "spin_polarized",
                    "type": "bool"
                },
                {
                    "name": "electronic_minimizer",
                    "type": "str"
                }
            ]
        },
        "4.1": {
            "params": [
                {
                    "name": "nelectrons",
                    "type": "d"
                },
                {
                    "name": "nup",
                    "type": "d"
                },
                {
                    "name": "ndown",
                    "type": "d"
                },
                {
                    "name": "spin",
                    "type": "d"
                },
                {
                    "name": "charge",
                    "type": "d"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "spin_treatment",
                    "type": "str"
                }
            ]
        }
    },
    "ELEC_MIN": {
        "2.0": {
            "params": [
                {
                    "name": "max_SD_steps",
                    "type": "i"
                },
                {
                    "name": "max_CG_steps",
                    "type": "i"
                },
                {
                    "name": "max_DIIS_steps",
                    "type": "i"
                },
                {
                    "name": "metals_method",
                    "type": "str"
                },
                {
                    "name": "elec_energy_tol",
                    "type": "d"
                },
                {
                    "name": "elec_eigenvalue_tol",
                    "type": "d"
                },
                {
                    "name": "elec_convergence_win",
                    "type": "i"
                },
                {
                    "name": "max_SCF_cycles",
                    "type": "i"
                },
                {
                    "name": "spin_fix",
                    "type": "i"
                },
                {
                    "name": "fix_occupancy",
                    "type": "bool"
                },
                {
                    "name": "smearing_scheme",
                    "type": "str"
                },
                {
                    "name": "smearing_width",
                    "type": "d"
                },
                {
                    "name": "efermi_tol",
                    "type": "d"
                },
                {
                    "name": "num_occ_cycles",
                    "type": "i"
                },
                {
                    "name": "elec_dump_file",
                    "type": "str"
                },
                {
                    "name": "num_dump_cycles",
                    "type": "i"
                },
                {
                    "name": "elec_restore_file",
                    "type": "str"
                }
            ]
        },
        "4.1": {
            "else": elec_min_v4_1_else,
            "params": [
                {
                    "name": "elec_force_tol",
                    "type": "d"
                }
            ]
        },
        "6.1": {
            "params": [
                {
                    "name": "dipole_correction",
                    "type": "str"
                },
                {
                    "name": "dipole_dir",
                    "type": "str"
                }
            ]
        }
    },
    "ELNES": {
        "2.0": {
            "params": [
                {
                    "name": "elnes_nbands",
                    "type": "i"
                },
                {
                    "name": "elnes_xc_functional",
                    "type": "str"
                },
                {
                    "name": "elnes_eigenvalue_tol",
                    "type": "d"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "elnes_xc_definition_size",
                    "type": "i"
                },
                {
                    "name": "elnes_xc_definition",
                    "type": "str"
                }
            ]
        }

    },
    "GA": {
        "2.0": {
            "params": [
                {
                    "name": "ga_pop_size",
                    "type": "i"
                },
                {
                    "name": "ga_max_gens",
                    "type": "i"
                },
                {
                    "name": "ga_mutate_rate",
                    "type": "d"
                },
                {
                    "name": "ga_mutate_amp",
                    "type": "d"
                },
                {
                    "name": "ga_fixed_N",
                    "type": "bool"
                },
                {
                    "name": "ga_bulk_slice",
                    "type": "bool"
                }
            ]
        }
    },
    "GENERAL": {
        "2.0": {
            "params": [
                {
                    "name": "seedname",
                    "type": "str"
                },
                {
                    "name": "comment",
                    "type": "str"
                },
                {
                    "name": "iprint",
                    "type": "i"
                },
                {
                    "name": "continuation",
                    "type": "str"
                },
                {
                    "name": "reuse",
                    "type": "str"
                },
                {
                    "name": "checkpoint",
                    "type": "str"
                },
                {
                    "name": "task",
                    "type": "str"
                },
                {
                    "name": "calculate_stress",
                    "type": "bool"
                },
                {
                    "name": "run_time",
                    "type": "i"
                },
                {
                    "name": "num_backup_iter",
                    "type": "i"
                },
                {
                    "name": "print_clock",
                    "type": "bool"
                },
                {
                    "name": "length_unit",
                    "type": "str"
                },
                {
                    "name": "mass_unit",
                    "type": "str"
                },
                {
                    "name": "time_unit",
                    "type": "str"
                },
                {
                    "name": "charge_unit",
                    "type": "str"
                },
                {
                    "name": "energy_unit",
                    "type": "str"
                },
                {
                    "name": "force_unit",
                    "type": "str"
                },
                {
                    "name": "velocity_unit",
                    "type": "str"
                },
                {
                    "name": "pressure_unit",
                    "type": "str"
                },
                {
                    "name": "inv_length_unit",
                    "type": "str"
                },
                {
                    "name": "frequency_unit",
                    "type": "str"
                },
                {
                    "name": "page_wvfns",
                    "type": "i"
                },
                {
                    "name": "rand_seed",
                    "type": "i"
                },
                {
                    "name": "data_distribution",
                    "type": "str"
                },
                {
                    "name": "opt_strategy",
                    "type": "str"
                }
            ]
        },
        "3.0": {
            "else": general_v3_0_else,
            "params": [
                {
                    "name": "backup_interval",
                    "type": "i"
                },
                {
                    "name": "force_constant_unit",
                    "type": "str"
                },
                {
                    "name": "volume_unit",
                    "type": "str"
                },
                {
                    "name": "opt_strategy_bias",
                    "type": "i"
                }
            ]
        },
        "3.1": {
            "else": general_v3_1_else,
            "params": [
                {
                    "name": "calculate_densdiff",
                    "type": "bool"
                },
                {
                    "name": "num_farms_requested",
                    "type": "i"
                },
                {
                    "name": "num_proc_in_smp",
                    "type": "i"
                }
            ]
        },
        "3.2": {
            "else": general_v3_2_else,
            "params": [
               {
                   "name": "num_proc_in_smp_fine",
                   "type": "i"
               },
                {
                   "name": "message_size",
                   "type": "i"
               },
                {
                   "name": "requested_seed",
                   "type": "i"
               }
            ]
        },
        "4.0": {
            "params": [
                {
                    "name": "ir_intensity_unit",
                    "type": "str"
                }
            ]
        },
        "4.1": {
            "params": [
                {
                    "name": "print_memory_usage",
                    "type": "bool"
                },
                {
                    "name": "write_formatted_potential",
                    "type": "bool"
                },
                {
                    "name": "write_formatted_density",
                    "type": "bool"
                },
                {
                    "name": "dipole_unit",
                    "type": "str"
                },
                {
                    "name": "efield_unit",
                    "type": "str"
                },
                {
                    "name": "calc_molecular_dipole",
                    "type": "bool"
                }
            ]
        },
        "4.2": {
            "params": [
                {
                    "name": "write_orbitals",
                    "type": "bool"
                }
            ]
        },
        "4.3": {
            "params": [
                {
                    "name": "calculate_elf",
                    "type": "bool"
                },
                {
                    "name": "write_formatted_elf",
                    "type": "bool"
                },
                {
                    "name": "cml_output",
                    "type": "bool"
                },
                {
                    "name": "dummy_fname",
                    "type": "str"
                }
            ]
        },
        "5.0": {
            "params": [
                {
                    "name": "calculate_Hirshfeld",
                    "type": "bool"
                }
            ]
        },
        "5.5": {
            "params": [
                {
                    "name": "entropy_unit",
                    "type": "str"
                },
                {
                    "name": "write_cif_structure",
                    "type": "bool"
                },
                {
                    "name": "write_cell_structure",
                    "type": "bool"
                }
            ]
        },
        "6.0": {
            "params": [
                {
                    "name": "write_bib",
                    "type": "bool"
                }
            ]
        },
        "8.0": {
            "params": [
                {
                    "name": "write_checkpoint",
                    "type": "str"
                },
                {
                    "name": "spin_unit",
                    "type": "str"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "write_otfg",
                    "type": "bool"
                }
            ]
        },
        "17.0": {
            "params": [
                {
                    "name": "write_cst_esp",
                    "type": "bool",
                },
                {
                    "name": "write_bands",
                    "type": "bool",
                },
                {
                    "name": "write_geom",
                    "type": "bool",
                }
            ]
        },
        "18.0": {
            "params": [
                {
                    "name": "verbosity",
                    "type": "str"
                },
                {
                    "name": "write_none",
                    "type": "bool",
                },
                {
                    "name": "write_md",
                    "type": "bool",
                }
            ]
        }
    },
    "GEOM": {
        "2.0": {
            "params": [
                {
                    "name": "geom_method",
                    "type": "str"
                },
                {
                    "name": "geom_max_iter",
                    "type": "i"
                },
                {
                    "name": "geom_energy_tol",
                    "type": "d"
                },
                {
                    "name": "geom_force_tol",
                    "type": "d"
                },
                {
                    "name": "geom_disp_tol",
                    "type": "d"
                },
                {
                    "name": "geom_stress_tol",
                    "type": "d"
                },
                {
                    "name": "geom_convergence_win",
                    "type": "i"
                },
                {
                    "name": "geom_modulus_est",
                    "type": "d"
                },
                {
                    "name": "geom_frequency_est",
                    "type": "d"
                }
            ]
        },
        "4.0": {
            "params": [
                {
                    "name": "geom_spin_fix",
                    "type": "i"
                },
                {
                    "name": "geom_linmin_tol",
                    "type": "d"
                }
            ]
        },
        "5.5": {
            "params": [
                {
                    "name": "geom_use_linmin",
                    "type": "bool"
                }
            ]
        },
        "6.0": {
            "params": [
                {
                    "name": "geom_lbfgs_max_updates",
                    "type": "i"
                }
            ]
        },
        "6.1": {
            "params": [
                {
                    "name": "geom_tpsd_init_stepsize",
                    "type": "d"
                },
                {
                    "name": "geom_tpsd_iterchange",
                    "type": "i"
                }
            ]
        }
    },
    "MAGRES": {
        "2.0": {
            "params": [
                {
                    "name": "magres_task",
                    "type": "str"
                },
                {
                    "name": "magres_method",
                    "type": "str"
                },
                {
                    "name": "magres_max_CG_steps",
                    "type": "i"
                },
                {
                    "name": "magres_convergence_win",
                    "type": "i"
                },
                {
                    "name": "magres_conv_tol",
                    "type": "d"
                },
                {
                    "name": "magres_xc_functional",
                    "type": "str"
                }
            ]
        },
        "4.3": {
            "params": [
                {
                    "name": "magres_max_sc_cycles",
                    "type": "i"
                },
                {
                    "name": "magres_jcoupling_task",
                    "type": "str"
                },
                {
                    "name": "magres_write_response",
                    "type": "bool"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "magres_xc_definition_size",
                    "type": "i"
                },
                {
                    "name": "magres_xc_definition",
                    "type": "str"
                }
            ]
        }

    },
    "MD": {
        "2.0": {
            "params": [
                {
                    "name": "md_num_iter",
                    "type": "i"
                },
                {
                    "name": "md_delta_t",
                    "type": "d"
                },
                {
                    "name": "md_ensemble",
                    "type": "str"
                },
                {
                    "name": "md_temperature",
                    "type": "d"
                },
                {
                    "name": "md_thermostat",
                    "type": "str"
                },
                {
                    "name": "md_nose_t",
                    "type": "d"
                },
                {
                    "name": "md_langevin_t",
                    "type": "d"
                },
                {
                    "name": "md_extrap",
                    "type": "str"
                },
                {
                    "name": "md_extrap_fit",
                    "type": "bool"
                },
                {
                    "name": "md_damping_scheme",
                    "type": "str"
                },
                {
                    "name": "md_damping_reset",
                    "type": "i"
                },
                {
                    "name": "md_opt_damped_delta_t",
                    "type": "bool"
                },
                {
                    "name": "md_elec_energy_tol",
                    "type": "d"
                },
                {
                    "name": "md_elec_eigenvalue_tol",
                    "type": "d"
                },
                {
                    "name": "md_elec_convergence_win",
                    "type": "i"
                }
            ]
        },
        "3.0": {
            "else": md_v3_0_else,
            "params": [
                {
                    "name": "md_barostat",
                    "type": "str"
                },
                {
                    "name": "md_ion_t",
                    "type": "d"
                },
                {
                    "name": "md_cell_t",
                    "type": "d"
                },
                {
                    "name": "md_nhc_length",
                    "type": "i"
                }
            ]
        },
        "3.1": {
            "else": md_v3_1_else,
            "params": [
                {
                    "name": "md_use_pathint",
                    "type": "bool"
                },
                {
                    "name": "md_num_beads",
                    "type": "i"
                }
            ]
        },
        "4.0": {
            "params": [
                {
                    "name": "md_pathint_staging",
                    "type": "bool"
                },
                {
                    "name": "md_pathint_num_stages",
                    "type": "i"
                },
                {
                    "name": "md_sample_iter",
                    "type": "i"
                }
            ]
        },
        "4.1": {
            "else": md_v4_1_else,
            "params": [
                {
                    "name": "md_eqm_method",
                    "type": "str"
                },
                {
                    "name": "md_eqm_ion_t",
                    "type": "d"
                },
                {
                    "name": "md_eqm_cell_t",
                    "type": "d"
                },
                {
                    "name": "md_eqm_t",
                    "type": "d"
                },
                {
                    "name": "md_elec_force_tol",
                    "type": "d"
                }
            ]
        },
        "5.5": {
            "params": [
                {
                    "name": "md_pathint_init",
                    "type": "str"
                }
            ]
        },
        "6.0": {
            "params": [
                {
                    "name": "md_use_plumed",
                    "type": "bool"
                }
            ]
        },
        "8.0": {
            "params": [
                {
                    "name": "md_xlbomd",
                    "type": "bool"
                },
                {
                    "name": "md_xlbomd_history",
                    "type": "i"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "md_hug_method",
                    "type": "str"
                },
                {
                    "name": "md_hug_t",
                    "type": "d"
                },
                {
                    "name": "md_hug_compression",
                    "type": "d"
                }
            ]
        },
        "17.0": {
            "params": [
                {
                    "name": "md_hug_dir",
                    "type": "str"
                }
            ]
        }
    },
    "OPTICS": {
        "2.0": {
            "params": [
                {
                    "name": "optics_nbands",
                    "type": "i"
                }
            ]
        },
        "3.0": {
            "else": optics_v3_0_else,
            "params": [
                {
                    "name": "optics_xc_functional",
                    "type": "str"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "optics_xc_definition_size",
                    "type": "i"
                },
                {
                    "name": "optics_xc_definition",
                    "type": "str"
                }
            ]
        }
    },
    "PHONON": {
        "2.0": {
            "params": [
                {
                    "name": "dummy_logical",
                    "type": "bool"
                },
                {
                    "name": "phonon_energy_tol",
                    "type": "d"
                },
                {
                    "name": "phonon_max_CG_steps",
                    "type": "i"
                },
                {
                    "name": "phonon_max_cycles",
                    "type": "i"
                },
                {
                    "name": "phonon_convergence_win",
                    "type": "i"
                },
                {
                    "name": "phonon_preconditioner",
                    "type": "str"
                },
                {
                    "name": "phonon_use_kpoint_symmetry",
                    "type": "bool"
                }
            ]
        },
        "3.0": {
            "else": phonon_v3_0_else,
            "params": [
                {
                    "name": "phonon_dos_spacing",
                    "type": "d"
                },
                {
                    "name": "phonon_finite_disp",
                    "type": "d"
                },
                {
                    "name": "phonon_method",
                    "type": "str"
                },
                {
                    "name": "phonon_calc_lo_to_splitting",
                    "type": "bool"
                },
                {
                    "name": "phonon_sum_rule",
                    "type": "bool"
                },
                {
                    "name": "calculate_born_charges",
                    "type": "bool"
                },
                {
                    "name": "born_charge_sum_rule",
                    "type": "bool"
                }
            ]
        },
        "3.1": {
            "params": [
                {
                    "name": "phonon_calculate_dos",
                    "type": "bool"
                },
                {
                    "name": "phonon_dos_limit",
                    "type": "d"
                },
                {
                    "name": "phonon_force_constant_cutoff",
                    "type": "d"
                },
                {
                    "name": "phonon_fine_method",
                    "type": "str"
                },
                {
                    "name": "phonon_method",
                    "type": "str"
                }
            ]
        },
        "4.4": {
            "else": phonon_v4_4_else,
            "params": [
                {
                    "name": "phonon_force_constant_cut_scale",
                    "type": "d"
                },
                {
                    "name": "phonon_sum_rule_method",
                    "type": "str"
                },
                {
                    "name": "calculate_raman",
                    "type": "bool"
                }
            ]
        },
        "5.0": {
            "else": phonon_v5_0_else,
            "params": [
                {
                    "name": "phonon_fine_cutoff_method",
                    "type": "str"
                }
            ]
        },
        "6.0": {
            "else": phonon_v6_0_else,
            "params": [
                {
                    "name": "phonon_dfpt_method",
                    "type": "str"
                },
                {
                    "name": "raman_range_low",
                    "type": "d"
                },
                {
                    "name": "raman_range_high",
                    "type": "d"
                },
                {
                    "name": "phonon_write_force_constants",
                    "type": "bool"
                },
                {
                    "name": "phonon_write_dynamical",
                    "type": "bool"
                }
            ]
        },
        "17.2": {
            "params": [
                {
                    "name": "raman_method",
                    "type": "str"
                }
            ]
        }
    },
    "POPN": {
        "2.0": {
            "params": [
                {
                    "name": "popn_calculate",
                    "type": "bool"
                },
                {
                    "name": "popn_bond_cutoff",
                    "type": "d"
                },
                {
                    "name": "pdos_calculate_weights",
                    "type": "bool"
                }
            ]
        },
        "18.0": {
            "params": [
                {
                    "name": "popn_write",
                    "type": "str"
                }
            ]
        }
    },
    "PSPOT": {
        "2.0": {
            "params": [
                {
                    "name": "pspot_nonlocal_type",
                    "type": "str"
                },
                {
                    "name": "pspot_beta_phi_type",
                    "type": "str"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "spin_orbit_coupling",
                    "type": "bool"
                }
            ]
        },
        "18.0": {
            "params": [
                {
                    "name": "spinor_spin_polarized",
                    "type": "bool"
                }
            ]
        }

        # write (dump_unit) current_params%spin_orbit_coupling
        # !end of keywords for this block in version 9.0
        # write (dump_unit) current_params%spinors_spin_polarized
        # !end of keywords for this block in version 18.0

    },
    "SPECTRAL": {
        "2.0": {
            "params": [
                {
                    "name": "spectral_theory",
                    "type": "str"
                },
                {
                    "name": "spectral_task",
                    "type": "str"
                },
                {
                    "name": "spectral_max_iter",
                    "type": "i"
                },
                {
                    "name": "spectral_max_steps_per_iter",
                    "type": "i"
                },
                {
                    "name": "spectral_nbands",
                    "type": "i"
                },
                {
                    "name": "spectral_eigenvalue_tol",
                    "type": "d"
                },
                {
                    "name": "spectral_xc_functional",
                    "type": "str"
                },
                {
                    "name": "spectral_re_est_k_scrn",
                    "type": "bool"
                },
                {
                    "name": "spectral_write_eigenvalues",
                    "type": "bool"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "spectral_xc_definition_size",
                    "type": "i"
                },
                {
                    "name": "spectral_xc_definition",
                    "type": "str"
                }
            ]
        }        
    },
    "TDDFT": {
        "2.0": {
            "params": [
                {
                    "name": "tddft_on",
                    "type": "bool"
                },
                {
                    "name": "tddft_num_states",
                    "type": "i"
                },
                {
                    "name": "tddft_selected_state",
                    "type": "i"
                },
                {
                    "name": "tddft_eigenvalue_tol",
                    "type": "d"
                },
                {
                    "name": "tddft_convergence_win",
                    "type": "i"
                },
                {
                    "name": "tddft_max_iter",
                    "type": "i"
                },
                {
                    "name": "tddft_nextra_states",
                    "type": "i"
                },
                {
                    "name": "tddft_xc_functional",
                    "type": "str"
                },
                {
                    "name": "tddft_method",
                    "type": "str"
                },
                {
                    "name": "tddft_eigenvalue_method",
                    "type": "str"
                },
                {
                    "name": "tddft_approximation",
                    "type": "str"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "tddft_xc_definition_size",
                    "type": "i"
                },
                {
                    "name": "tddft_xc_definition",
                    "type": "str"
                },
                {
                    "name": "tddft_position_method",
                    "type": "str"
                }
            ]
        }

    },
    "THERMODYNAMICS": {
        "2.0": {
            "params": [
                {
                    "name": "thermo_calculate_helmholtz",
                    "type": "bool"
                },
                {
                    "name": "thermo_t_start",
                    "type": "d"
                },
                {
                    "name": "thermo_t_stop",
                    "type": "d"
                },
                {
                    "name": "thermo_t_npoints",
                    "type": "i"
                }
            ]
        },
        "3.2": {
            "else": thermodynamics_v3_2_else,
            "params": [
                {
                    "name": "thermo_t_spacing",
                    "type": "d"
                }
            ]
        }
    },
    "TSSEARCH": {
        "2.0": {
            "params": [
                {
                    "name": "tssearch_method",
                    "type": "str"
                },
                {
                    "name": "tssearch_lstqst_protocol",
                    "type": "str"
                },
                {
                    "name": "tssearch_qst_max_iter",
                    "type": "i"
                },
                {
                    "name": "tssearch_CG_max_iter",
                    "type": "i"
                },
                {
                    "name": "tssearch_force_tol",
                    "type": "d"
                },
                {
                    "name": "tssearch_disp_tol",
                    "type": "d"
                }
            ]
        },
        "7.0": {
            "params": [
                {
                    "name": "tssearch_max_path_points",
                    "type": "i"
                },
                {
                    "name": "tssearch_energy_tol",
                    "type": "d"
                }
            ]
        }
    },
    "WANNIER": {
        "2.0": {
            "params": [
                {
                    "name": "wannier_spread_type",
                    "type": "str"
                },
                {
                    "name": "dummy_real",
                    "type": "d"
                },
                {
                    "name": "dummy_real",
                    "type": "d"
                },
                {
                    "name": "dummy_real",
                    "type": "d"
                },
                {
                    "name": "wannier_min_algor",
                    "type": "str"
                },
                {
                    "name": "wannier_max_SD_steps",
                    "type": "i"
                },
                {
                    "name": "wannier_SD_step",
                    "type": "d"
                },
                {
                    "name": "wannier_ion_moments",
                    "type": "bool"
                },
                {
                    "name": "wannier_ion_rmax",
                    "type": "d"
                },
                {
                    "name": "wannier_ion_cut",
                    "type": "bool"
                },
                {
                    "name": "wannier_ion_cut_fraction",
                    "type": "d"
                }
            ]
        },
        "4.0": {
            "params": [
                {
                    "name": "wannier_spread_tol",
                    "type": "d"
                },
                {
                    "name": "wannier_print_cube",
                    "type": "i"
                },
                {
                    "name": "wannier_ion_cut_tol",
                    "type": "d"
                },
                {
                    "name": "wannier_restart",
                    "type": "str"
                },
                {
                    "name": "wannier_ion_cmoments",
                    "type": "bool"
                }
            ]
        }
    },
    "XC": {
        "2.0": {
            "params": [
                {
                    "name": "xc_functional",
                    "type": "str"
                },
                {
                    "name": "xc_vxc_deriv_epsilon",
                    "type": "d"
                }
            ]
        },
        "3.0": {
            "params": [
                {
                    "name": "nlxc_page_ex_pot",
                    "type": "i"
                },
                {
                    "name": "nlxc_ppd_integral",
                    "type": "bool"
                },
                {
                    "name": "nlxc_ppd_size_x",
                    "type": "i"
                },
                {
                    "name": "nlxc_ppd_size_y",
                    "type": "i"
                },
                {
                    "name": "nlxc_ppd_size_z",
                    "type": "i"
                },
                {
                    "name": "nlxc_impose_trs",
                    "type": "bool"
                },
                {
                    "name": "nlxc_exchange_reflect_kpts",
                    "type": "bool"
                },
                {
                    "name": "nlxc_k_scrn_den_function",
                    "type": "str"
                },
                {
                    "name": "nlxc_k_scrn_averaging_scheme",
                    "type": "str"
                },
                {
                    "name": "nlxc_re_est_k_scrn",
                    "type": "bool"
                },
                {
                    "name": "nlxc_calc_full_ex_pot",
                    "type": "bool"
                }
            ]
        },
        "4.0": {
            "params": [
                {
                    "name": "nlxc_div_corr_on",
                    "type": "bool"
                },
                {
                    "name": "nlxc_div_corr_S_width",
                    "type": "d"
                },
                {
                    "name": "nlxc_div_corr_tol",
                    "type": "d"
                },
                {
                    "name": "nlxc_div_corr_npts_step",
                    "type": "i"
                }
            ]
        },
        "5.0": {
            "params": [
                {
                    "name": "sedc_apply",
                    "type": "bool"
                },
                {
                    "name": "sedc_scheme",
                    "type": "str"
                },
                {
                    "name": "nlxc_exchange_screening",
                    "type": "d"
                }
            ]
        },
        "6.1": {
            "params": [
                {
                    "name": "nlxc_exchange_fraction",
                    "type": "d"
                }
            ]
        },
        "7.0": {
            "params": [
                {
                    "name": "sedc_sr_TS",
                    "type": "d"
                },
                {
                    "name": "sedc_d_TS",
                    "type": "d"
                },
                {
                    "name": "sedc_s6_G06",
                    "type": "d"
                },
                {
                    "name": "sedc_d_G06",
                    "type": "d"
                },
                {
                    "name": "sedc_lambda_OBS",
                    "type": "d"
                },
                {
                    "name": "sedc_n_OBS",
                    "type": "d"
                },
                {
                    "name": "sedc_sr_JCHS",
                    "type": "d"
                },
                {
                    "name": "sedc_s6_JCHS",
                    "type": "d"
                },
                {
                    "name": "sedc_d_JCHS",
                    "type": "d"
                }
            ]
        },
        "8.0": {
            "params": [
                {
                    "name": "relativistic_treatment",
                    "type": "str"
                }
            ]
        },
        "9.0": {
            "params": [
                {
                    "name": "xc_definition_size",
                    "type": "i"
                },
                {
                    "name": "xc_definition",
                    "type": "str"
                }
            ]
        }
    }
}
