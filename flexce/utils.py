# @Author: Brett Andrews <andrews>
# @Date:   2018-04-16 20:04:48
# @Last modified by:   andrews
# @Last modified time: 2018-07-13 13:07:61

"""
FILE
    utils.py

DESCRIPTION
    Utility functions for flexCE.
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
from pathlib import PurePath

import numpy as np


def set_mass_bins(low=0.1, high=100, dm_low=0.1, dm_high=1., break_mass=8):
    """Set stellar mass bins.

    Args:
        low (float): Low end of lowest mass bin. Default is 0.1.
        high (float): High end of higher mass bin. Default is 100.
        dm_low (float): Mass bin size below break mass. Default is 0.1.
        dm_high (float): Mass bin size above break mass. Default is 1.
        break (float): Dividing mass between low and high mass bins.
            Default is 8.

    Returns:
        array: Stellar mass bins.
    """
    mbins = np.concatenate((np.arange(low, break_mass, dm_low),
                            np.arange(break_mass, high + 0.001, dm_high)))

    return mbins


def set_yields(params=None):
    """Set yield parameters.

    Args:
        params (dict): Yield parameters. Default is ``None``.
    """
    default = {
        'snii_dir': join('limongi06', 'iso_yields'),
        'agb_dir': join('karakas10', 'iso_yields'),
        'snia_dir': 'iwamoto99',
        'rprocess_dir': 'cescutti06',
        'sprocess_dir': 'busso01',
        'snia_model': 'w70',
        'r_elements': ['Ba', 'Eu'],
        's_elements': ['Ba'],
        'solar_metallicity': False,
        'sprocess_supersolar': False,
    }

    if params is None:
        params = {}

    for key, val in default.items():
        if key not in params.keys():
            params[key] = val

    return params


def set_path(path_in, default_path):
    if os.path.isfile(path_in):
        path = os.path.dirname(os.path.abspath(path_in))
        filename = os.path.basename(path_in)
    else:
        path = default_path
        filename = path_in
    return filename, path


def substitute_dir_in_path(path, olddir, newdir):
    pp = PurePath(path)
    parts = [p if p != olddir else newdir for p in pp.parts]
    return os.path.join(*parts)


def robust_random_poisson(xx):
    """Poisson draw that robustly handles large numbers.

    Draw a number from a Poisson distribution.  Used for determining
    the number of actual stars to form given an expected value.
    np.random.poisson cannot handle numbers larger than 2147483647
    (~2.14e9) because it uses the C long type.  For numbers larger than
    this value, round the expected (statistical) value to the nearest
    integer.  Since 2e9 is a large number, the Poisson fluctuations
    would have been relatively small anyway.

    This function is intended to replace this line of code:
    self.Nstar[i] = np.random.poisson(self.Nstar_stat[i])

    Args:
        xx (array): Input for Poisson draw.

    Returns:
        array: Poisson realization.
    """
    try:
        yy = np.random.poisson(xx)

    except ValueError:
        yy = np.zeros(len(xx), dtype=np.int64)

        for i, item in enumerate(xx):
            try:
                yy[i] = np.random.poisson(item)

            except ValueError:
                yy[i] = np.round(item)

    return yy


def merge(aa, bb, path=None):
    """Recursively merge two dictionaries (``bb`` into ``aa``).

    Uses ``aa`` as the default if values conflict. Assumes that ``aa``
    and ``bb`` do not contain lists of dicts.

    Args:
        aa (dict): Preferred option.
        bb (dict): Fall back option.
        path (str): Path down nested dictionaries. Default is ``None``.

    Returns:
        dict
    """

    if path is None:
        path = []

    for key in bb:
        if key in aa:
            if isinstance(aa[key], dict) and isinstance(bb[key], dict):
                merge(aa[key], bb[key], path + [str(key)])

            elif isinstance(aa[key], dict) ^ isinstance(bb[key], dict):  # XOR
                raise TypeError(f"{'.'.join(path + [str(key)])} must both be dict or non-dict.")

        else:
            aa[key] = bb[key]

    return aa
