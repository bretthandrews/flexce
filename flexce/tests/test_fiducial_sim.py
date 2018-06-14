# @Author: Brett Andrews <andrews>
# @Date:   2018-06-13 15:06:13
# @Last modified by:   andrews
# @Last modified time: 2018-06-14 14:06:46

"""
FILE
    compare_to_fiducial_sim.py

DESCRIPTION
    Check that the outputs of a simulation run on the refactored code
    matches the outputs of the fiducial simulation from
    Andrews et al. (2017).
"""


import pickle
import sys

import numpy as np
import pytest

from flexce.chemevol import ChemEvol

sys.path.append('/Users/andrews/projects/pcaa_chemevol/')
import chemevol_main


# load fiducial simulation
@pytest.fixture(scope='session')
def box0():
    """Fiducial simulation from Andrews et al. (2017)."""
    path = '/Users/andrews/projects/pcaa_chemevol/sims/paper_masscut01/runs/box0.pck'
    with open(path, 'rb') as fin:
        sim = pickle.load(fin)

    return sim


@pytest.fixture(scope='session')
def gal(box0):
    """Sim with same parameters and random state as fiducial sim."""
    print('\nRunning flexCE...')
    return ChemEvol(
        state={
            'Nstar': box0.random_num_state_Nstar,
            'snia': box0.random_num_state_snia
        }
    )


box0_ignore = [
    'frac_ev_tot', 'warmgas_on', 'alpha1', 'inflow_func', 'tcool', 'num_int2',
    'eta_outflow', 'mgas_init', 'snia_dtd_func', 'mass_frac2', 'mwarmgas_init',
    'min_snia_time', 'inflow_ab_pattern', 'mass_bins2', 'alpha', 'fwarm', 'sim_id',
    'N_kslaw', 'sf_param', 'snia_param', 'mass_breaks', 'alpha2', 'nu_kslaw',
    'warmgasres_param', 'fdirect', 'inflow_param', 'snia_fraction', 'variable_eta',
    'mass_int2', 'outflow_source', 'imf', 'outflow_param', 'inflow_metallicity',
    'snia_timescale', 'mass_ave2',
]
keys_gal = []
keys_box0 = box0_ignore


class TestFiducialSim(object):

    def test_params_box(self, gal, box0):
        for k, v in gal.params['box'].items():
            if k in ['sim_id', 'long_output']:
                continue

            assert box0.__dict__[k] == v, k
            keys_box0.append(k)

    def test_params_inflows(self, gal, box0):
        for k, v in gal.params['inflows'].items():
            if k == 'coeff':
                k = 'k'

            assert box0.param['inflow'][k] == v, k

    def test_params_outflows(self, gal, box0):
        for k, v in gal.params['outflows'].items():
            if k == 'eta':
                k = 'eta_outflow'
            elif k == 'source':
                k = 'outflow_source'

            if k == 'feject':
                assert box0.__dict__[k] == v, k
                keys_box0.append(k)
            else:
                assert box0.param['outflow'][k] == v, k

    def test_params_warmgas(self, gal, box0):
        for k, v in gal.params['warmgas'].items():
            if k not in ['fdirect', 'tcool']:
                assert box0.param['warmgas'][k] == v, k

    def test_params_sf(self, gal, box0):
        for k, v in gal.params['sf'].items():
            if k == 'sfh':
                continue

            assert box0.param['sf'][k] == v, k

    def test_params_snia_dtd(self, gal, box0):
        for k, v in gal.params['snia_dtd'].items():
            if k == 'func':
                assert box0.param['snia'][k] == v, k

            else:
                if k == 'min_time':
                    k = 'min_snia_time'
                elif k == 'fraction':
                    k = 'snia_fraction'
                elif k == 'mass':
                    continue

                if k in box0.param['snia']['k'].keys():
                    assert box0.param['snia']['k'][k] == v, k
                else:
                    assert box0.__dict__[k] == v, k
                    keys_box0.append(k)

    def test_yields(self, gal, box0):
        for k, v in gal.params['yields'].items():
            if k in ['snii_dir', 'agb_dir', 'snia_dir', 'rprocess_dir', 'sprocess_dir']:
                k = k.split('_dir')[0]

            if k in ['snia_model', 'r_elements', 's_elements']:
                continue

            assert box0.param['yields'][k] == v.split('/')[0], k

        keys_gal.append('params')
        keys_box0.append('param')

    def test_random_state(self, gal, box0):
        for gal_state, box0_state in zip(gal.state['Nstar'], box0.random_num_state_Nstar):
            assert gal_state == box0_state

        keys_box0.append('random_num_state_Nstar')

        for gal_state, box0_state in zip(gal.state['snia'], box0.random_num_state_snia):
            assert gal_state == box0_state

        keys_box0.append('random_num_state_snia')
        keys_gal.append('state')

    def test_fraction_evolved(self, gal, box0):
        """Test that lists of arrays are the same."""
        for item in ['ind_ev', 'frac_ev']:
            for ii, (gg, bb) in enumerate(zip(gal.__dict__[item], box0.__dict__[item])):
                assert (gg == bb).all(), (item, ii)

            keys_gal.append(item)
            keys_box0.append(item)

    def test_constants(self, gal, box0):
        constants = ['n_bins', 'n_bins_high', 'n_bins_low', 'dtime', 'n_steps']
        for constant in constants:
            constant_box = constant if constant != 'dtime' else 'dt'

            assert gal.__dict__[constant] == box0.__dict__[constant_box], constant
            keys_gal.append(constant)
            keys_box0.append(constant_box)

    def test_arrays(self, gal, box0):
        arrays = [
            'mass_bins', 'mass_int', 'num_int', 'time', 'mass_ave', 'mass_frac', 'tau_m',
            'inflow_rate', 'warmgas_ab_pattern', 'gas_cooling', 'inflow', 'mwarmfrac',
            'mwarmgas_iso', 'agb', 'dm_sfr', 'metallicity', 'mfrac', 'mgas_iso', 'mremnant',
            'mstar', 'mstar_left', 'mstar_stat', 'Mwd', 'Mwd_Ia', 'Mwd_Ia_init', 'NIa',
            'Nstar', 'Nstar_left', 'Nstar_stat', 'outflow', 'sf', 'sfr', 'snia', 'snii',
            'outflow_rate', 'NII', 'survivors',
        ]
        for arr in arrays:
            arr_box = arr
            if arr == 'time':
                arr_box = 't'

            if arr == 'warmgas_ab_pattern':
                if gal.warmgas_ab_pattern is None:
                    assert box0.warmgas_ab_pattern == 0, arr
            else:
                assert np.isclose(gal.__dict__[arr], box0.__dict__[arr_box]).all(), arr

            keys_gal.append(arr)
            keys_box0.append(arr_box)

    def test_attribute_set(self, gal, box0):
        diff_gal = set(gal.__dict__.keys()) - set(keys_gal)
        diff_box0 = set(box0.__dict__.keys()) - set(keys_box0)
        assert not diff_gal, f'gal: {diff_gal}'
        assert not diff_box0, f'box0: {diff_box0}'
