from __future__ import print_function, division, absolute_import

import os
import sys
from pathlib import PurePath

import numpy as np

from yields import Yields


def define_mass_bins(low=0.1, high=100, dm_low=0.1, dm_high=1.):
    """Define stellar mass bins.

    Args:
        low (float): low end of lowest mass bin. Defaults to 0.1 Msun.
        high (float): high end of higher mass bin. Defaults to 100 Msun.
        dm_low (float): mass bin size below 8 Msun. Defaults to 0.1 Msun.
        dm_high (float): mass bin size above 8 Msun. Defaults to 1 Msun.

    Returns:
        array: stellar mass bins
    """
    mbins = np.concatenate((np.arange(low, 8., dm_low),
                            np.arange(8., high + 0.001, dm_high)))
    return mbins


def load_yields(path, args, mass_bins):
    """Load yield grids.

    Args:
        path (str): data directory.
        args (dict):
        mass_bins (array): stellar mass bins.

    Returns:
        Yields instance

    """
    try:
        yld = Yields(path, snii_dir=args['snii_dir'],
                     agb_dir=args['agb_dir'], snia_dir=args['snia_dir'],
                     rprocess_dir=args['rprocess_dir'],
                     sprocess_dir=args['sprocess_dir'],
                     snia_model=args['snia_model'],
                     r_elements=args['r_elements'],
                     s_elements=args['s_elements'],
                     mbins=mass_bins)
    except IOError as e:
        print()
        print(e)
        # FIX: automatically run make_yield_grids.py
        print('\nPlease create yield grids with')
        print('python make_yield_grids.py\n')
        sys.exit(1)
    return yld


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


def none_to_empty_dict(x):
    """If a variable is None, return an empty dictionary."""
    if x is None:
        x = {}
    return x
