# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 09:07:24
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 10:07:68

"""
FILE
    txt.py

DESCRIPTION
    Read and write simulation abundances from .txt files.

    These files only save the time, number of surviving stars, [Fe/H],
    and [X/Fe] abundances of each stellar generation. Use
    ``flexce.io.pck`` to save all properties of a simulation.
"""

from __future__ import print_function, division, absolute_import

from os.path import join
import re

import pandas as pd

from flexce.fileio.general import _make_sim_path


def write(path, sim):
    """Write simulation abundances to a .txt file.

    Args:
        path (str): Output file path.
        sim: ``flexce.chemevol.ChemEvol`` instance.
    """
    formatters = ['{:<7}'.format, '{:10}'.format] + ['{:8.5f}'.format for _ in sim.ab]

    header = (['{:8}'.format('Time'), '{:10}'.format('Survivors'), '{:8}'.format('[Fe/H]')] +
              ['{:8}'.format('[{}/Fe]'.format(item)) for item in sim.xfe])

    df = pd.DataFrame({'time': sim.time[1:], 'survivors': sim.survivors[1:]})
    df2 = pd.concat((df, sim.ab), axis=1)

    sim_id = sim.params['box']['sim_id'] or sim.params['box']['datetime']

    with open(join(path, f'ab{sim_id}.txt'), 'w') as fout:
        df2.to_string(fout, header=header, index=False, formatters=formatters)


def read(path, sim_id=None):
    """Read simulation abundances from a text file.

    Args:
        path (str): Output path.
        sim_id (str): Simulation ID number. Default is ``None``.
    """
    fname = _make_sim_path(path, sim_id, ext='.txt')
    return pd.read_csv(fname, delim_whitespace=True)


def rename_cols(df):
    """Rename columns of DataFrame for easier manipulation.

    Args:
       df (DataFrame): Abundances.

    Returns:
        DataFrame
    """
    oldcols = df.columns
    newcols = [re.sub('[\[/\]]', '', item.lower()) for item in oldcols]

    colname_map = {oldcol: newcol for oldcol, newcol in zip(oldcols, newcols)}

    df.rename(columns=colname_map, inplace=True)

    return df
