# @Author: Brett Andrews <andrews>
# @Date:   2018-06-07 21:06:85
# @Last modified by:   andrews
# @Last modified time: 2018-06-07 21:06:79

"""
FILE
    warm_gas_reservoir.py

DESCRIPTION
    Functions for controlling gas flow into and out of a warm gas
    reservoir.
"""

import os
from os.path import join

import numpy as np
import pandas as pd


def set_warm_gas_reservoir(
    params,
    feject,
    warmgas=False,
    mwarmgas_init=0.,
    fdirect=0.01,
    tcool=1200.,
):
    """
    Set parameters about warm gas reservoir.

    The concept of a warm gas reservoir is adopted from
    Schoenrich & Binney (2009).  In this model, stars eject a fraction
    of their yields ``feject`` (specified in ``set_outflows()``) from
    the galaxy and then return a fraction of their yields directly back
    to the cold gas reservoir (``mgas_iso``).  The remainder of the
    yields go into the warm gas reservoir (``fwarm``).  The warm gas
    reservoir can start with some initial gas mass (``mwarmgas_init``)
    or be empty.  The warm gas reservoir returns gas back to the cold
    gas reservoir on an expontential timescale (``tcool``).

    Schoenrich & Binney (2009) fiducial values:
        mwarmgas_init = 5e8 Msun
        fdirect = 0.01 (feject=0.15 for R < 3.5 kpc and 0.04 for R > 3.5 kpc)
        tcool = 1200 Myr
        fwarm = 1 - fdirect - feject

    Args:
        params (dict): Simulation parameters.
        feject (float): Fraction of stellar yields ejected from galaxy.
        warmgas (bool): Include a warm gas reservoir in the model.
            Default is ``False``.
        mwarmgas_init (float): Initial gas mass of warm gas reservoir.
            Default is 0. If ``warmgas`` is ``False``, then 0.
        fdirect (float): Fraction of stellar yields injected directly
            into the cold gas reservoir.  Default is 0.01. If
            ``warmgas`` is ``False``, then ``1 - feject``.
        tcool (float): Exponential gas cooling timescale for gas to
            flow to cold gas reservoir. Default is 1200 [Myr].

    Returns:
        dict, array,
    """
    assert fdirect >= 0. and fdirect <= 1., \
        '``fdirect`` is a fraction and so must be between 0 and 1 inclusive'

    if warmgas:
        filename = 'warmgas_abundance_pattern.txt'
        tmp = pd.read_csv(self.path_yldgen + filename,
                          delim_whitespace=True, skiprows=10,
                          names=['el', 'ab'])
        self.warmgas_ab_pattern = np.array(tmp['ab'])
        self.mwarmgas_init = mwarmgas_init
        self.fdirect = fdirect
        self.tcool = tcool
        self.fwarm = 1. - fdirect
        if self.outflow_source == 'stellar_ejecta':
            self.fwarm -= self.feject
    else:
        self.warmgas_ab_pattern = 0.
        self.mwarmgas_init = 0.
        self.fdirect = 1. - self.feject
        self.tcool = 0.
        self.fwarm = 0.
    self.warmgasres_param = dict(warmgas_on=warmgas,
                                 mwarmgas_init=mwarmgas_init,
                                 fdirect=fdirect, fwarm=self.fwarm,
                                 tcool=tcool)
