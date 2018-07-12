# @Author: Brett Andrews <andrews>
# @Date:   2018-05-30 17:05:89
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 11:07:92

"""
FILE
    run_flexce

USAGE
    run_flexce sim-config-file

    Example::

        $ run_flexce sim0.cfg

DESCRIPTION
    Run flexCE from the command line using a config file.

"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join

import click

import flexce.abundances
from flexce.chemevol import ChemEvol
from flexce.fileio import pck, txt, yml
import flexce.utils
from flexce.yields import Yields


@click.command()
@click.option('--path_in', default=None, type=click.Path())
@click.option('--path_out', default=None, type=click.Path())
@click.argument('config_file', default=None, required=False, type=click.Path())
def main(config_file, path_in, path_out):
    """
    Args:
        config_file (str): Config file name.  Default is ``sim0.yml``.
            (fiducial parameters from Andrews et al. 2017).
        path_in (str): Path to `config_file`.  Default is current dir.
        path_out (str): Path to output dir.  Default is current dir.
    """
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
    path_flexce_root = os.path.abspath(join(path_flexce, '..'))
    path_data = join(path_flexce, 'data')

    if config_file is None:
        print('No config file specified. Using default parameters.')
        config_file = 'sim0.yml'
        path_in = join(path_flexce_root, 'examples')

    if path_in is None:
        print('No input path specified. Using current working directory.')
        path_in = os.getcwd()

    if path_out is None:
        print('No output path specified. Using current working directory.')
        path_out = os.getcwd()

    file_in = join(path_in, config_file)

    params = flexce.fileio.yml.read_yml(file_in)

    mass_bins = flexce.utils.set_mass_bins(**params['mass_bins'])

    params['yld_args'] = flexce.utils.set_yields(params['yields'])
    ylds = Yields(params=params['yields'], mass_bins=mass_bins, path=path_data)

    box = ChemEvol(params, ylds)

    write_output(path_out, params['sim_id'], box)


# TODO Remove in favor of ChemEvol.save()
def write_output(path, sim_id, sim):
    """Write simulation results to pickle and txt files.

    Args:
        path (str): Output path.
        sim_id (str): Simulation ID number.
        sim: ``flexce.chemevol.ChemEvol`` instance.
    """
    path_sim = join(path, ''.join(['sim', sim_id]))
    os.makedirs(path_sim, exist_ok=True)

    pck.write(join(path_sim, f'sim{sim_id}.pck'), sim)
    txt.write_abundances(join(path_sim, f'ab{sim_id}.txt'), sim)


if __name__ == '__main__':
    main()
