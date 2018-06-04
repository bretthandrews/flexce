#!/usr/bin/env python
# encoding: utf-8
#
# test_chemevol.py
#
# @Author: Brett Andrews <andrews>
# @Date:   2018-05-29 10:05:01
# @Last modified by:   andrews
# @Last modified time: 2018-06-04 10:06:88


import numpy as np
import pytest

from chemevol import function_x


class TestChemevol(object):

    @pytest.mark.parametrize('mgas_iso, ',
                             [('emline_gflux', 'ha_6564'),
                              ('stellar_vel', None)])
    def test_function_x(self):
        # sfr + outflow <= gas mass
        # if not, reduce sfr

        # if outflow is stellar_ejecta, don't account for it.
        # if outflow is ism, then account for it.

        
        assert x >= 0
