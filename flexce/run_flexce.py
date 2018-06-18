# @Author: Brett Andrews <andrews>
# @Date:   2018-05-30 17:05:89
# @Last modified by:   andrews
# @Last modified time: 2018-06-18 11:06:39

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
import flexce.io.yml
import flexce.io.pck
import flexce.io.txt
import flexce.utils


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

    params = flexce.io.yml.read_yml(file_in)

    mass_bins = flexce.utils.set_mass_bins(**params['mass_bins'])

    params['yld_args'] = flexce.utils.set_yields(params['yields'])
    ylds = flexce.utils.load_yields(path_data, mass_bins, params['yields'])

    # TODO enable loading state from file
    state = None

    box = ChemEvol(params, ylds, state)

    ab = flexce.abundances.calc_abundances(
        path_data,
        ylds.sym,
        box.mgas_iso,
        box.survivors,
        box.time,
        box.params,
    )

    write_output(path_out, params['sim_id'], box, ab)


def write_output(path, sim_id, gal, abund):
    """Write simulation results to pickle and txt files.

    Args:
        path (str): Output path.
        sim_id (str): Simulation ID number.
        gal: ChemEvol instance.
        abund: Abundances instance.
    """
    path_sim = join(path, ''.join(['sim', sim_id]))
    os.makedirs(path_sim, exist_ok=True)

    flexce.io.pck.pck_write(gal, join(path_sim, f'box{sim_id}.pck'))
    flexce.io.pck.pck_write(abund, join(path_sim, f'ab{sim_id}.pck'))
    flexce.io.txt.txt_write(gal, abund, join(path_sim, f'ab{sim_id}.txt'))


if __name__ == '__main__':
    main()
