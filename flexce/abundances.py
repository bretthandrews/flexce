# @Author: Brett Andrews <andrews>
# @Date:   2018-04-16 20:04:30
# @Last modified by:   andrews
# @Last modified time: 2018-07-13 13:07:08

"""
FILE
    abundances.py

DESCRIPTION
    Compute abundances.
"""

import os
from os.path import join
import string

import numpy as np
import pandas as pd

import flexce.utils


def set_default_params(params):
    default = {
        'solar': {
            'source': 'lodders',
        }
    }

    params = flexce.utils.merge(params, default)

    return params


def _map_isotopes_to_elements(isotopes):
    """Map isotopes onto element and atomic mass combinations.

    Args:
        isotopes (array): Isotopes.

    Returns:
        DataFrame: Boolean index array to select valid element and
            atomic mass combinations.
    """
    remove_digits = str.maketrans('', '', string.digits)
    element_only = np.array([iso.translate(remove_digits) for iso in isotopes])

    elements = np.array(sorted(set(element_only), key=list(element_only).index))

    return pd.DataFrame(np.array([element_only == el for el in elements]).T, columns=elements)


def load_solar(source='lodders'):
    """Read in solar abundances.

    Args:
        source (str): Reference for solar abundances. Default is
            'lodders'.

    Returns:
        DataFrame
    """
    path_yldgen = join(os.path.dirname(__file__), 'data', 'yields', 'general')

    if source == 'lodders':

        solar = pd.read_csv(
            join(path_yldgen, 'lodders03_solar_photosphere.txt'),
            delim_whitespace=True,
            skiprows=8,
            usecols=[0, 1],
            names=['el', 'xh'],
            index_col=0,
        )

    else:
        raise ValueError('Valid solar abundance references: "lodders".')

    solar['xfe'] = np.log10(10.**(solar['xh'] - 12) / 10.**(solar['xh'].loc['Fe'] - 12))

    return solar


def compute(mgas_iso, ylds, solar_source):
    """Compute [Fe/H] and [X/Fe].

    Args:
        mgas_iso (array): Isotopic gas masses.
        ylds (flexce.yields.Yields): Yields.
        solar_source (str): Reference for solar abundances.

    Returns:
        DataFrame: [Fe/H] and [X/Fe].
    """

    ind_element = _map_isotopes_to_elements(ylds.sym)

    solar = load_solar(solar_source)

    ngas_iso = mgas_iso / ylds.sym_mass

    niso_h = (ngas_iso[1:].T / ngas_iso[1:, ind_element['H']].sum(axis=1)).T
    niso_fe = (ngas_iso[1:].T / ngas_iso[1:, ind_element['Fe']].sum(axis=1)).T

    xh_abs = pd.DataFrame(
        np.log10([niso_h[:, ind_element[el]].sum(axis=1) for el in ylds.element]).T + 12.,
        columns=ylds.element,
    )

    xfe_abs = pd.DataFrame(
        np.log10([niso_fe[:, ind_element[el]].sum(axis=1) for el in ylds.element]).T,
        columns=ylds.element,
    )

    feh = xh_abs['Fe'] - solar['xh'].loc['Fe']
    feh.name = 'feh'

    xfe = xfe_abs - solar['xfe'].loc[ylds.element]

    return pd.concat([feh, xfe], axis=1)
