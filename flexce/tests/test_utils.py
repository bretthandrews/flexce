# @Author: Brett Andrews <andrews>
# @Date:   2018-04-16 20:04:48
# @Last modified by:   andrews
# @Last modified time: 2018-07-13 13:07:19

import copy

import numpy as np
import pytest

import flexce.utils

a1 = {'x': 1}
a2 = {'x': 2}
a3 = {'x': {}}
a4 = {'x': {'y': 1}}
a5 = {'x': {'y': 2}}
a6 = {'x': {'y': 1, 'z': 2}}
a7 = {'x': {'y': 2, 'z': 2}}
a8 = {'i': 5}
expected1 = {'x': 1}
expected2 = {'x': {'y': 1}}
expected3 = {'x': {'y': 1, 'z': 2}}
expected4 = {'x': {'y': 2, 'z': 2}}
expected5 = {'i': 5, 'x': {'y': 2, 'z': 2}}
expected6 = {'i': 5, 'x': 1}


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


class TestMerge(object):

    @pytest.mark.parametrize(
        'aa, bb, expected',
        [(a8, a7, expected5),
         (a1, a1, expected1),
         (a1, a2, expected1),
         (a4, a4, expected2),
         (a3, a4, expected2),
         (a4, a3, expected2),
         (a4, a5, expected2),
         (a6, a3, expected3),
         (a6, a4, expected3),
         (a6, a5, expected3),
         (a4, a6, expected3),
         (a5, a6, expected4),
         (a7, a8, expected5),
         (a1, a8, expected6),
         (a8, a1, expected6),
         ]
    )
    def test_merge_success(self, aa, bb, expected):
        cc = copy.deepcopy(aa)
        print('\naa', aa)
        print('bb', bb)
        print('expected', expected)
        assert flexce.utils.merge(cc, bb) == expected

    @pytest.mark.parametrize(
        'aa, bb',
        [(a1, a3),
         (a2, a3),
         (a3, a1),
         (a3, a2),
         (a1, a4),
         (a2, a4),
         (a4, a1),
         (a4, a2),
         ]
    )
    def test_merge_fail(self, aa, bb):
        with pytest.raises(TypeError):
            flexce.utils.merge(aa, bb)
