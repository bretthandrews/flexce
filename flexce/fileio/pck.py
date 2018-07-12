# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 10:07:51
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 16:07:42

"""
FILE
    pck.py

DESCRIPTION
    Pickle and unpickle simulation output.
"""

from __future__ import print_function, division, absolute_import

from os.path import join

import pickle


def read(path):
    """Read pickle file.

    Args:
        path (str): Name of pickle file.

    Returns:
        object
    """
    with open(path, 'rb') as fin:
        obj = pickle.load(fin)

    return obj


def write(path, obj):
    """Write object to pickle file.

    Args:
        path (str): Name of output pickle file.
        obj: Object to be pickled.
    """
    with open(path, 'wb') as fout:
        pickle.dump(obj, fout, -1)


def read_sim(path, sim_id):
    return read(join(path, f'sim{sim_id}.pck'))
