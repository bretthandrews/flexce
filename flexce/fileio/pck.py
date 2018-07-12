# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 10:07:51
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 11:07:83

"""
FILE
    pck.py

DESCRIPTION
    Pickle and unpickle simulation output.
"""

from __future__ import print_function, division, absolute_import

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
