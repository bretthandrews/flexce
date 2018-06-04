from __future__ import print_function, division, absolute_import

import unittest
import numpy as np
import utils

class UtilsTestCase(unittest.TestCase):
    """Tests from 'utils.py'."""

    def test_define_mass_bins(self):
        mbins = utils.define_mass_bins(low=1, high=100, dm_low=1, dm_high=1.)
        tmp = np.arange(1, 101)
        np.testing.assert_array_almost_equal(mbins, tmp)

if __name__ == '__main__':
    unittest.main()