# @Author: Brett Andrews <andrews>
# @Date:   2018-06-13 15:06:13
# @Last modified by:   andrews
# @Last modified time: 2018-06-21 15:06:12

"""
FILE
    test_fiducial_box.py

DESCRIPTION
    Test that the outputs and abundances of a simulation run on the
    refactored code match the outputs and abundances of the fiducial
    simulation from Andrews et al. (2017).
"""

import pickle
import sys

import numpy as np
import pytest

from flexce.abundances import calc_abundances
from flexce.chemevol import ChemEvol
from flexce.yields import Yields

sys.path.append('/Users/andrews/projects/pcaa_chemevol/')
import chemevol_main


@pytest.fixture(scope='session')
def box0():
    """Fiducial simulation from Andrews et al. (2017)."""
    path = '/Users/andrews/projects/pcaa_chemevol/sims/paper_masscut01/runs/box0.pck'
    with open(path, 'rb') as fin:
        return pickle.load(fin)


@pytest.fixture(scope='session')
def ab0():
    """Abundances of fiducial simulation from Andrews et al. (2017)."""
    path = '/Users/andrews/projects/pcaa_chemevol/sims/paper_masscut01/runs/ab0.pck'
    with open(path, 'rb') as fin:
        return pickle.load(fin)


@pytest.fixture(scope='session')
def gal(box0):
    """Sim with same parameters and random state as fiducial sim."""

    print('\nRunning flexCE...')

    return ChemEvol(
        params={
            'box': {
                'save': {
                    'slim': False,
                    'state': {
                        'Nstar': box0.random_num_state_Nstar,
                        'snia': box0.random_num_state_snia
                    },
                },
            },
        },
    )


@pytest.fixture(scope='session')
def ab(gal):
    """Abundances of simulation run with current flexCE version."""
    ylds = Yields(params=gal.params['yields'], mass_bins=gal.mass_bins)

    return calc_abundances(
        ylds.sym,
        gal.mgas_iso,
        gal.survivors,
        gal.time,
        gal.params,
    )


box0_ignore = [
    'frac_ev_tot', 'warmgas', 'alpha1', 'inflow_func', 'tcool', 'num_int2',
    'eta_outflow', 'mgas_init', 'snia_dtd_func', 'mass_frac2', 'mwarmgas_init',
    'min_snia_time', 'inflow_ab_pattern', 'mass_bins2', 'alpha', 'fwarm', 'sim_id',
    'N_kslaw', 'sf_param', 'snia_param', 'mass_breaks', 'alpha2', 'nu_kslaw',
    'warmgasres_param', 'fdirect', 'inflow_param', 'snia_fraction', 'variable_eta',
    'mass_int2', 'outflow_source', 'imf', 'outflow_param', 'inflow_metallicity',
    'snia_timescale', 'mass_ave2', 'warmgas_on'
]
keys_gal = []
keys_box0 = box0_ignore

keys_ab = ['param', 'sim_id']
keys_ab0 = ['param', 'sim_id', 'stem_data', 'stem_parent']


class TestFiducialBox(object):

    def test_params_box(self, gal, box0):
        for k, v in gal.params['box'].items():
            if k in ['sim_id', 'save']:
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
                k = k if k != 'warmgas' else 'warmgas_on'
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

            if k in ['snia_model', 'r_elements', 's_elements', 'solar_metallicity',
                     'sprocess_supersolar']:
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


class TestFiducialAbunds(object):

    def test_arrays(self, ab, ab0):
        arrays = ['all_atomic_num', 'atomic_num_out', 'feh', 'isotope_mass', 'mgas_iso', 'niso_fe', 'niso_h', 'ngas_iso', 'solar_ab', 'solar_fe', 'solar_h', 'survivors', 't', 'xfe', 'xfe_abs', 'xfe_all', 'xh_abs', 'xh_all']

        for arr in arrays:
            arr_old = arr
            if arr == 'solar_ab':
                arr_old = 'solar_abund'

            assert np.isclose(ab.__dict__[arr], ab0.__dict__[arr_old]).all(), arr

            keys_ab.append(arr)
            keys_ab0.append(arr_old)

    def test_str_arrays(self, ab, ab0):
        arrays = ['all_elements', 'elements', 'elements_out', 'isotope', 'solar_element', 'sym']

        for arr in arrays:
            arr_old = arr
            assert (ab.__dict__[arr] == ab0.__dict__[arr_old]).all(), arr

            keys_ab.append(arr)
            keys_ab0.append(arr_old)

    def test_constants(self, ab, ab0):
        constants = ['n_elements', 'n_isotope', 'n_steps']
        for constant in constants:
            constant_old = constant
            assert ab.__dict__[constant] == ab0.__dict__[constant_old], constant

            keys_ab.append(constant)
            keys_ab0.append(constant_old)

    def test_ind_element(self, ab, ab0):
        for kk, vv in ab.ind_element.items():
            assert (vv == ab.ind_element[kk]).all(), f'{kk} {vv}'

        keys_ab.append('ind_element')
        keys_ab0.append('ind_element')


class TestAllBoxAttributes(object):
    """Checks that all attributes have been tested."""

    def test_attribute_set(self, gal, box0):

        # WARNING:only works if running all the tests in this file.
        diff_gal = set(gal.__dict__.keys()) - set(keys_gal)
        diff_box0 = set(box0.__dict__.keys()) - set(keys_box0)
        assert not diff_gal, f'gal: {diff_gal}'
        assert not diff_box0, f'box0: {diff_box0}'


class TestAllAbundAttributes(object):
    """Checks that all attributes have been tested."""

    def test_attribute_set(self, ab, ab0):
        # WARNING:only works if running all the tests in this file.
        diff_ab = set(ab.__dict__.keys()) - set(keys_ab)
        diff_ab0 = set(ab0.__dict__.keys()) - set(keys_ab0)
        assert not diff_ab, f'ab: {diff_ab}'
        assert not diff_ab0, f'ab0: {diff_ab0}'
