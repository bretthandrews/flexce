"""Read and write simulation abundances from .txt files.

Only saves time, surviving stars, [Fe/H], and [X/Fe] abundances of each
stellar generation. Use fileio.pickle_io to save all details of a simulation.
"""

from __future__ import print_function, division, absolute_import

from os.path import join
import re

import numpy as np
import pandas as pd


def _make_sim_path(path_out, sim_id):
    """Construct path to simulation output pickle file.

    Args:
        path_out (str): output path
        sim_id (str): simulation id number

    Returns:
        str: file path and name
    """
    sim_id = str(sim_id)
    path_sim = join(path_out, ''.join(['sim', sim_id]))
    fname = 'ab{}.txt'.format(sim_id)
    return join(path_sim, fname)


def _colnames(elements):
    """Generate column names for output text file.

    Args:
        elements (list): List of output element abbreviations (str).

    Returns:
        str: header line that describes columns.
    """
    time = ['{:8}'.format('Time')]
    survivors = ['{:10}'.format('Survivors')]
    feh = ['{:8}'.format('[Fe/H]')]
    abunds = ['{:8}'.format('[{}/Fe]'.format(item)) for item in elements]
    line = time + survivors + feh + abunds
    return '  '.join(line) + '\n'


def txt_write(path_out, sim_id, box, ab):
    """Write simulation abundances to a text file.

    Args:
        path_out (str): output path
        sim_id (str): simulation id number
        box (obj): box object of simulation
        ab (obj): ab object of simulation
    """
    fname = _make_sim_path(path_out, sim_id)
    with open(fname, 'w') as f:
        f.write(_colnames(ab.elements_out))
        for t, s, fe, xfe in zip(box.t[1:], box.survivors[1:], ab.feh,
                                 ab.xfe.T):
            time = ['{:<7}'.format(t)]
            survivors = ['{:10}'.format(s)]
            feh = ['{:8.5f}'.format(fe)]
            abunds = ['{:8.5f}'.format(item) for item in xfe]
            line = time + survivors + feh + abunds
            f.write('  '.join(line) + '\n')


def txt_read(path_out, sim_id):
    """Read simulation abundances from a text file.

    Args:
        path_out (str): output path
        sim_id (str): simulation id number

    Returns:
        time, survivors, feh, abunds (array, array, array, recarray): time in
        Myr; survivors is the number of stars from each time step that survive
        until the end of the simulation; and abunds are the abundances at each
        time step.
    """
    fname = _make_sim_path(path_out, sim_id)
    with open(fname, 'r') as f:
        header = f.readline()
    abnames = header.split('[Fe/H]')[1].strip().split()
    elnames = [it.split('/Fe')[0].strip('[') for it in abnames]
    dtypes = [(it, float) for it in elnames]
    sim = np.array(pd.read_csv(fname, delim_whitespace=True))
    time, survivors, feh = sim.T[:3]
    abunds_in = sim.T[3:]
    abunds = np.array([tuple(it) for it in abunds_in.T], dtype=dtypes)
    return time, survivors, feh, abunds


def _rename_cols(df):
    """Rename columns of dataframe for easier manipulation.

    Args:
       df (DataFrame): Dataframe of abundances.

    Returns:
        DataFrame
    """
    oldcols = df.columns
    newcols = [re.sub('[\[/\]]', '', item.lower()) for item in oldcols]
    colname_map = {}
    for oldcol, newcol in zip(oldcols, newcols):
        colname_map[oldcol] = newcol
    df.rename(columns=colname_map, inplace=True)
    return df


def load_dataframe(path_out, sim_id):
    """Load simulation abundances into a pandas dataframe.

    Args:
        path_out (str): output path
        sim_id (str): simulation id number

    Returns:
        DataFrame
    """
    sim = 'sim' + sim_id
    fname = 'ab{}.txt'.format(sim_id)
    path = join(path_out, sim, fname)
    df = pd.read_table(path, index_col='Time', sep='\s*', engine='python')
    # df = _rename_cols(df)
    return df
