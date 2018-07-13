# @Author: Brett Andrews <andrews>
# @Date:   2018-06-07 20:06:05
# @Last modified by:   andrews
# @Last modified time: 2018-07-13 09:07:78

"""
FILE
    inflows.py

DESCRIPTION
    Functions for setting the inflow rate and composition.
"""

import copy
import os
from os.path import join

import numpy as np
import pandas as pd


def set_inflows(
    time=None,
    func='exp',
    mgas_init=2e10,
    coeff=None,
    inflow_rate=None,
    inflow_ab_pattern='bbns',
    inflow_metallicity=1.,
):
    """Set the inflow rate.

    Available function forms and required coefficients:
        double_exp: M1, b1, M2, b2 (see Eq. 6 in
            Schoenrich & Binney 2009) with b1 & b2 in Myr.
        exp: M1, b1 with b1 in Myr.
        te-t: M1, b1 with b1 in Myr.
        constant_mgas: ``inflow_rate`` will be dynamically defined in
            ``evolve_box``.
        two_infall:
            tau_thick: exponential timescale of thick disk in Myr,
            tau_thin: exponential timescale of thin disk in Myr,
            t_max: time of maximum inflow in Myr,
            sfe_thick: multiplicative enhancement factor for SFE in
                thick disk in Myr,
            t_sf_off (list): time when star formation in thick disk
                ends and when star formation in the thin disk resumes,
            sfe: star formation efficiency.
        custom: User-defined inflow rate.

    Available inflow abundance pattern options:
        bbns: Big Bang Nucleosynthesis abundance pattern
        alpha_enhanced: Abundance pattern of a simulation before SNIa.
        scaled_solar: Solar abundance pattern scaled relative to solar
            (i.e., solar = 1).
        recycled: Abundance pattern of last time step.

    Args:
        time (array): Time steps. Default is ``None``.
        func (str): Functional form of inflow rate. Default is 'exp'.
        mgas_init (float): Initial gas mass [Msun]. Default is 2e10.
        coeff (dict): Coefficients defining inflow rate functions.
            Default is ``None``.
        inflow_rate (array): User-defined inflow rate [Msun/Myr] if
            ``func`` is 'custom'. Default is ``None``.
        inflow_ab_pattern (str): Inflow abundance pattern. Default is
            'bbns'.
        inflow_metallicity (float): Scaling of inflow metallicity
            (i.e., solar = 1). Default is 1.

    Returns:
        dict, array: Inflow parameters; inflow rate.
    """
    assert time is not None, 'Must pass in ``time``.'

    coeff = _set_default_coeff(func, coeff)

    # Set inflow rate
    if func == 'double_exp':
        inflow_rate = ((coeff['M1'] / coeff['b1']) * np.exp(-time / coeff['b1']) +
                       (coeff['M2'] / coeff['b2']) * np.exp(-time / coeff['b2']))

    elif func == 'exp':
        inflow_rate = (coeff['M1'] / coeff['b1']) * np.exp(-time / coeff['b1'])

    elif func == 'te-t':
        inflow_rate = ((coeff['M1'] / coeff['b1']) * (time / coeff['b1']) *
                       np.exp(-time / coeff['b1']))

    elif func == 'constant_mgas':
        inflow_rate = np.zeros(len(time))

    elif func == 'two_infall':
        inflow_rate_thick = coeff['sfe'] * np.exp(-time / coeff['tau_thick'])
        inflow_rate_thin = coeff['sfe'] * np.exp(-(time - coeff['t_max']) / coeff['tau_thin'])
        inflow_rate = inflow_rate_thick + inflow_rate_thin

    elif func == 'custom':
        assert inflow_rate is not None, 'If inflow rate functional form is "custom",' \
            ' then ``inflow_rate`` must be provided.'

    else:
        raise ValueError('Valid inflow functions: "double_exp", "exp", "te-t",'
                         ' "constant_mgas", "two_infall", and "custom".')

    params = {
        'mgas_init': mgas_init,
        'func': func,
        'coeff': coeff,
        'ab_pattern': inflow_ab_pattern,
        'metallicity': inflow_metallicity,
    }

    return params, inflow_rate


def _set_default_coeff(func, coeff):
    """Set the default inflow rate coefficients.

    Args:
        func (str): Functional form of inflow rate.
        coeff (dict): Coefficients defining inflow rate functions.

    """
    coeff = coeff if coeff is not None else {}

    if func == 'exp':
        default = {'M1': 4e11, 'b1': 6000.}

    elif func == 'double_exp':
        default = {'M1': 1e10, 'b1': 300., 'M2': 6e11, 'b2': 14000.}

    elif func == 'te-t':
        default = {'M1': 4e11, 'b1': 3500.}

    elif func == 'two_infall':
        default = {
            'tau_thick': 1000.,
            'tau_thin': 6000.,
            't_max': 1000.,
            'sfe_thick': 2.,
            't_sf_off': [1000., 1200.],
            'sfe': 3e7,
        }

    coeff = {**default, **coeff}

    return coeff


def inflow_composition(params, yields, mgas_iso_last):
    """Compute the mass fraction of each element in the inflowing gas.

    Available inflow abundance pattern options:
        bbns: Big Bang Nucleosynthesis abundance pattern
        alpha_enhanced: Abundance pattern of a simulation before SNIa.
        scaled_solar: Solar abundance pattern scaled relative to solar
            (i.e., solar = 1).
        recycled: Abundance pattern of last time step.


    Set hydrogen mass fraction to 0.75 and helium mass fraction to 0.25 -
    the mass fraction of metals.  You need a hydrogen to helium mass
    fraction ratio of ~3 to avoid negative absolute yields of hydrogen.  (I
    had originally set things up to match the hydrogen/helium ratio of the
    ISM but that ran away to negative hydrogen masses).

    Args:
        params (dict): Inflow parameters.
        yields: Yields instance.
        mgas_iso_last (array): Gas-phase mass of isotopes
            (``mgas_iso``) from the last time step.

    """
    scaling_factor = params['metallicity']  # relative to solar

    ind_metal = (yields.sym_mass > 4.)

    if params['ab_pattern'] == 'bbns':
        inflow = yields.bbmf
        return inflow

    elif params['ab_pattern'] == 'alpha_enhanced':
        path_yldgen = join(os.path.dirname(__file__), 'data', 'yields', 'general')
        alpha_en = pd.read_csv(join(path_yldgen, 'Z_0.1-Zsun_alpha_enhanced.txt'),
                               skiprows=6, header=None)
        inflow_init = np.array(alpha_en).T
        scaling_factor *= 10.

    elif params['ab_pattern'] == 'scaled_solar':
        inflow_init = copy.deepcopy(yields.solar_mfrac)

    elif params['ab_pattern'] == 'recycled':
        inflow_init = mgas_iso_last / mgas_iso_last.sum()
        scaling_factor *= 0.02 / inflow_init[ind_metal].sum()

    else:
        raise ValueError('Valid inflow compositions: "bbns", "alpha_enhanced", '
                         '"scaled_solar", and "recycled".')

    inflow = np.zeros(yields.n_sym)
    inflow[ind_metal] = inflow_init[ind_metal] * scaling_factor
    # Set H & He mass fraction to 0.75 & 0.25 - Z, respectively.
    inflow[yields.sym == 'H1'] = 0.75
    inflow[yields.sym == 'He4'] = 0.25 - inflow.sum()

    return inflow
