# @Author: Brett Andrews <andrews>
# @Date:   2018-04-16 20:04:48
# @Last modified by:   andrews
# @Last modified time: 2018-06-14 14:06:94

import numpy as np
import pytest

import flexce.utils


class TestUtils(object):

    @pytest.mark.parametrize(
        'params, expected',
        [({'low': 1, 'high': 100, 'dm_low': 1, 'dm_high': 1., 'break_mass': 8}, np.arange(1, 101)),
         ({'low': 0.1, 'high': 100, 'dm_low': 0.1, 'dm_high': 1., 'break_mass': 8},
          np.concatenate((np.arange(0.1, 8, 0.1), np.arange(8, 101)))),
         ])
    def test_set_mass_bins(self, params, expected):
        mbins = flexce.utils.set_mass_bins(**params)
        assert mbins == pytest.approx(expected)
