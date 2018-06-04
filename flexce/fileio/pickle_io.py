"""Pickle and unpickle simulation output.

Save and reload entire simulation object from pickle files.
"""

from __future__ import print_function, division, absolute_import

import os
import pickle


def pickle_read(filename):
    """Read pickle file.

    Args:
        filename (str): name of pickle file.

    Returns:
        object
    """
    fin = open(filename, 'rb')
    obj = pickle.load(fin)
    fin.close()
    return obj


def pickle_write(obj, filename):
    """Write object to pickle file.

    Args:
        obj: object to be pickled.
        filename (str): name of output pickle file.
    """
    fout = open(filename, 'wb')
    pickle.dump(obj, fout, -1)
    fout.close()


def _make_sim_path(path_out, sim_id, stem='box'):
    """Construct path to simulation output pickle file.

    Args:
        path_out (str): directory of pickle file.
        sim_id (str): simulation ID number.

    Returns:
        str: file path and name
    """
    sim_id = str(sim_id)
    path_sim = os.path.join(path_out, ''.join(['sim', sim_id]))
    fname = '{}{}.pck'.format(stem, sim_id)
    return os.path.join(path_sim, fname)


def box_read(path_out, sim_id):
    """Read in ChemEvol instance from box<sim_id>.pck file.

    Args:
        path_out (str): directory of pickle file.
        sim_id (str): simulation ID number.

    Returns:
        object: instance of ChemEvol class ('box' object).
    """
    fname = _make_sim_path(path_out, sim_id, stem='box')
    return pickle_read(fname)


def ab_read(path_out, sim_id):
    """Read in Abundances instance from ab<sim_id>.pck file.

    Args:
        path_out (str): directory of pickle file.
        sim_id (str): simulation ID number.

    Returns:
        object: instance of Abundances class ('ab' object).
    """
    fname = _make_sim_path(path_out, sim_id, stem='ab')
    return pickle_read(fname)


def sim_read(path_out, sim_id):
    """Read in box and ab objects.

    Args:
        path_out (str): directory of pickle file.
        sim_id (str): simulation ID number.

    Returns:
        ChemEvol instance, Abundances instance (tuple): 'box' and 'ab' objects
    """
    box = box_read(path_out, sim_id)
    ab = ab_read(path_out, sim_id)
    return box, ab
