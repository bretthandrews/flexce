# @Author: Brett Andrews <andrews>
# @Date:   2018-05-30 17:05:89
# @Last modified by:   andrews
# @Last modified time: 2018-06-04 20:06:29

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
import sys

import flexce
from flexce.chemevol import evolve
from flexce.abundances import calc_abundances
from flexce.io.cfg import read_sim_cfg
from flexce.io.pck import pck_write
from flexce.io.txt import txt_write


def write_output(path, sim_id, gal, abund):
    """Write simulation results to pickle and txt files.

    Args:
        path (str): output directory.
        sim_id (str): simulation ID number.
        gal: ChemEvol instance.
        abund: Abundances instance.
    """
    path_sim = join(path, ''.join(['sim', sim_id]))
    if not os.path.isdir(path_sim):
        os.mkdir(path_sim)

    pck_write(gal, join(path_sim, ''.join(('box', sim_id, '.pck'))))
    pck_write(abund, join(path_sim, ''.join(('ab', sim_id, '.pck'))))

    txt_write(path_out, sim_id, gal, abund)


if __name__ == '__main__':

    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
    path_flexce_root = os.path.abspath(join(path_flexce, '..'))
    path_data = join(path_flexce, 'data')

    argv = None
    if argv is None:
        argv = sys.argv

    try:
        default_config_path = join(path_flexce_root, 'config')
        fname, path_config = flexce.utils.set_path(argv[1], default_config_path)
    except IndexError:
        path_config = join(os.getenv('HOME'), 'flexce', 'examples')
        fname = 'sim0.cfg'
        print('\nUsing default parameters in \n{}'.format(argv[1]))

    file_in = join(path_config, fname)

    # TODO Add try...except to handle user-defined output path
    path_out = flexce.utils.substitute_dir_in_path(path_config, 'config', 'output')

    (simulation_id, yld_args, initialize_args, mass_bins_args, snia_dtd_args,
     inflows_args, outflows_args, warmgasres_args, sf_args) = \
        read_sim_cfg(file_in)
    mass_bins = flexce.utils.define_mass_bins(**mass_bins_args)
    ylds = flexce.utils.load_yields(path_data, yld_args, mass_bins)
    box = evolve(ylds, initialize_args, snia_dtd_args, inflows_args,
                 outflows_args, warmgasres_args, sf_args)
    ab = calc_abundances(path_data, ylds.sym, box.mgas_iso, box.survivors,
                         box.t, box.param, box.sim_id)
    write_output(path_out, simulation_id, box, ab)

# TODO (specify elements_out in config file)
