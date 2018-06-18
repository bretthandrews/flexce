# @Author: Brett Andrews <andrews>
# @Date:   2018-04-16 20:04:48
# @Last modified by:   andrews
# @Last modified time: 2018-06-18 15:06:60

"""
FILE
    utils.py

DESCRIPTION
    Utility functions for flexCE.
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys
from pathlib import PurePath

import numpy as np

import flexce


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
    }

    if params is None:
        params = {}

    for key, val in default.items():
        if key not in params.keys():
            params[key] = val

    return params


def load_yields(path=None, mass_bins=None, kwargs=None):
    """Load yield grids.

    Args:
        path (str): Data directory. Default is ``None``.
        mass_bins (array): Stellar mass bins. Default is ``None``.
        kwargs (dict): Keyword arguments to pass to ``Yields``. Default
            is ``None``.

    Returns:
        Yields instance.
    """
    from flexce.yields import Yields

    if path is None:
        path = join(os.path.dirname(flexce.__file__), 'data')

    mass_bins = mass_bins if mass_bins is not None else {}
    kwargs = kwargs if kwargs is not None else {}

    try:
        ylds = Yields(path, mbins=mass_bins, **kwargs)

    except IOError as e:
        print()
        print(e)
        print('\nPlease create yield grids with')
        print('python make_yield_grids.py\n')
        sys.exit(1)

    return ylds


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
