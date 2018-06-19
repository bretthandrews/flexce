#!/usr/bin/env python
# encoding: utf-8
#
# test_chemevol.py
#
# @Author: Brett Andrews <andrews>
# @Date:   2018-05-29 10:05:01
# @Last modified by:   andrews
# @Last modified time: 2018-06-19 13:06:21


import pickle
import sys

import numpy as np
import pytest

from flexce.chemevol import ChemEvol

sys.path.append('/Users/andrews/projects/pcaa_chemevol/')
import chemevol_main


@pytest.fixture(scope='session')
def path_box_AWSJ():
    return '/Users/andrews/projects/pcaa_chemevol/sims/paper_masscut01/runs/box0.pck'


@pytest.fixture(scope='session')
def box_AWSJ(path_box_AWSJ):
    """Fiducial simulation from Andrews et al. (2017).

    flexCE version is v1.0.
    """
    with open(path_box_AWSJ, 'rb') as fin:
        sim = pickle.load(fin)

    return sim


@pytest.fixture(scope='session')
def path_box():
    return '/Users/andrews/sim0/box0.pck'


@pytest.fixture(scope='session')
def box(path_box):
    """Fiducial simulation from Andrews et al. (2017).

    flexCE version is v2.0.
    """
    with open(path_box, 'rb') as fin:
        sim = pickle.load(fin)

    return sim


class TestChemevol(object):

    @pytest.mark.parametrize(
        'state, expected', [
            (None, {'Nstar': [], 'snia': []}),
            ({}, {'Nstar': [], 'snia': []}),
            ({'Nstar': []}, {'Nstar': [], 'snia': []}),
            ({'snia': []}, {'Nstar': [], 'snia': []}),
            ({'Nstar': [1, 2, 3]}, {'Nstar': [1, 2, 3], 'snia': []}),
            ({'snia': [1, 2, 3]}, {'Nstar': [], 'snia': [1, 2, 3]}),
        ]
    )
    def test_set_state_default(self, state, path_box_AWSJ, box_AWSJ, expected):
        assert ChemEvol.set_state(state) == expected

    def _check_states(self, state1, state2):
        for key in state2.keys():
            for sv, ev in zip(state1[key], state2[key]):
                for it1, it2 in zip(sv, ev):
                    if isinstance(it1, np.ndarray):
                        assert (it1 == it2).all()
                    else:
                        assert it1 == it2

    def test_set_state_flexce_v1(self, path_box_AWSJ, box_AWSJ):
        state = ChemEvol.set_state({
            'Nstar': path_box_AWSJ,
            'snia': path_box_AWSJ,
        })

        expected = {
            'Nstar': box_AWSJ.random_num_state_Nstar,
            'snia': box_AWSJ.random_num_state_snia,
        }

        self._check_states(state, expected)

    def test_set_state_flexce_v2(self, path_box, box):
        state = ChemEvol.set_state({
            'Nstar': path_box,
            'snia': path_box,
        })

        expected = {
            'Nstar': box.state['Nstar'],
            'snia': box.state['snia'],
        }

        self._check_states(state, expected)


    # def test_gas_consumption_le_gas_mass(self):
        # sfr + outflow <= gas mass
        # if not, reduce sfr

        # if outflow is stellar_ejecta, don't account for it.
        # if outflow is ism, then account for it.
