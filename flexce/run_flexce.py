# @Author: Brett Andrews <andrews>
# @Date:   2018-05-30 17:05:89
# @Last modified by:   andrews
# @Last modified time: 2018-06-04 22:06:92

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
import flexce.chemevol
import flexce.io.cfg
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
        config_file (str): Config file name.  Default is ``sim0.cfg``
            (fiducial parameters from Andrews et al. 2017).
        path_in (str): Path to `config_file`.  Default is current dir.
        path_out (str): Path to output dir.  Default is current dir.
    """
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
    path_flexce_root = os.path.abspath(join(path_flexce, '..'))
    path_data = join(path_flexce, 'data')

    if config_file is None:
        print('No config file specified. Using default parameters.')
        config_file = 'sim0.cfg'
        path_in = join(path_flexce_root, 'examples')

    if path_in is None:
        print('No input path specified. Using current working directory.')
        path_in = os.getcwd()

    if path_out is None:
        print('No output path specified. Using current working directory.')
        path_out = os.getcwd()

    file_in = join(path_in, config_file)

    params = flexce.io.cfg.read_sim_cfg(file_in)

    mass_bins = flexce.utils.define_mass_bins(**params['mass_bins_args'])

    ylds = flexce.utils.load_yields(path_data, params['yld_args'], mass_bins)

    box = flexce.chemevol.evolve(
        ylds,
        params['initialize_args'],
        params['snia_dtd_args'],
        params['inflows_args'],
        params['outflows_args'],
        params['warmgasres_args'],
        params['sf_args']
    )

    ab = flexce.abundances.calc_abundances(
        path_data,
        ylds.sym,
        box.mgas_iso,
        box.survivors,
        box.t,
        box.param,
        box.sim_id)

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
