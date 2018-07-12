"""Pickle and unpickle simulation output.

Save and reload entire simulation object from pickle files.
"""

from __future__ import print_function, division, absolute_import

import os
import pickle


def pck_read(filename):
    """Read pickle file.

    Args:
        filename (str): Name of pickle file.

    Returns:
        object
    """
    fin = open(filename, 'rb')
    obj = pickle.load(fin)
    fin.close()
    return obj


def pck_write(obj, filename):
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
    fname = f'{stem}{sim_id}.pck'
    return os.path.join(path_sim, fname)


def sim_read(path_out, sim_id):
    """Read in simulation results.

    Args:
        path_out (str): Directory of pickle file.
        sim_id (str): Simulation ID number.

    Returns:
        flexce.chemevol.ChemEvol
    """
    from flexce.chemevol import ChemEvol
    fname = _make_sim_path(path_out, sim_id, stem='box')
    return pck_read(fname)
