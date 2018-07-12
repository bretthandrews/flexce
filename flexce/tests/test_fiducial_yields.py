# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 14:07:15
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 15:07:94


"""
FILE
    test_fiducial_yields.py

DESCRIPTION
    Test that the yields match from the yields used in
    Andrews et al. (2017).
"""

import os
import sys

import pytest

import flexce.utils
from flexce.yields import Yields

path_fiducial = '/Users/andrews/projects/pcaa_chemevol/'

assert os.path.isdir(path_fiducial), \
    'These tests require the Andrews et al. (2017) yields with the original directory structure.'

sys.path.append(path_fiducial)

from chemevol_import_yields import Yields as LegacyYields


@pytest.fixture(scope='session')
def legacy():
    """Fiducial yields from Andrews et al. (2017)."""

    yld_args = {
        'snii_dir': 'limongi06/iso_yields_010ni/',
        'agb_dir': 'karakas10/iso_yields/',
        'snia_dir': 'iwamoto99/',
        'rprocess_dir': 'cescutti06/',
        'sprocess_dir': 'busso01/',
        'snia_model': 'w70',
        'r_elements': ['Ba', 'Eu'],
        's_elements': ['Ba', ],
    }

    mass_bins = flexce.utils.set_mass_bins(low=0.1, high=100, dm_low=0.1, dm_high=1., break_mass=8)

    ylds = LegacyYields(
        stem_parent='/Users/andrews/projects/pcaa_chemevol/',
        snii_dir=yld_args['snii_dir'],
        agb_dir=yld_args['agb_dir'],
        snia_dir=yld_args['snia_dir'],
        rprocess_dir=yld_args['rprocess_dir'],
        sprocess_dir=yld_args['sprocess_dir'],
        mbins=mass_bins)
    ylds.load_rprocess_yields()
    ylds.load_sprocess_yields()
    ylds.load_snii_yields()
    ylds.load_agb_yields()
    ylds.load_snia_yields(model=yld_args['snia_model'])
    ylds.concat_ncapture_yields(r_elements=yld_args['r_elements'],
                                s_elements=yld_args['s_elements'])
    ylds.load_solar_abund()

    return ylds


@pytest.fixture(scope='session')
def current():
    return Yields()


nums = ['ind8', 'mlow', 'n_bins', 'n_bins_high', 'n_bins_low', 'n_elements',
        'n_nc_sym', 'n_snii_sym', 'n_sym', 'n_z']

str_arrays = ['agb_sym', 'element', 'element_all', 'nc_element', 'nc_sym', 'snii_sym',
              'solar_element', 'sym']

arrays = [
    'agb_mej', 'agb_rem', 'agb_z', 'atomic_num', 'bbmf', 'mass_bins', 'nc_sym_mass',
    'nc_yields', 'rprocess_yields', 'snia_yields', 'snii_agb_rem', 'snii_agb_z',
    'snii_mej', 'snii_rem', 'snii_sym_mass', 'snii_z', 'solar_fe', 'solar_h',
    'solar_mfrac', 'sprocess_yields', 'sym_mass',
]

large_arrays = ['agb_yields', 'snii_yields']


class TestYields(object):

    def test_solar_ab(self, legacy, current):
        assert legacy.solar_abund == pytest.approx(current.solar_ab)

    def test_nums(self, legacy, current):
        for num in nums:
            print(num)
            assert legacy.__dict__[num] == pytest.approx(current.__dict__[num]), num

    def test_str_arrays(self, legacy, current):
        for arr in str_arrays:
            print(arr)
            assert (legacy.__dict__[arr] == current.__dict__[arr]).all(), arr

    def test_arrays(self, legacy, current):
        for arr in arrays:
            print(arr)
            assert legacy.__dict__[arr] == pytest.approx(current.__dict__[arr]), arr

    @pytest.mark.slow
    def test_large_arrays(self, legacy, current):
        for arr in large_arrays:
            print(arr)
            assert legacy.__dict__[arr] == pytest.approx(current.__dict__[arr]), arr
