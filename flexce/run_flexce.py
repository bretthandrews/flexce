# @Author: Brett Andrews <andrews>
# @Date:   2018-05-30 17:05:89
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 22:07:00

"""
FILE
    run_flexce

USAGE
    run_flexce sim-config-file

    Example::

        $ run_flexce sim0.yml

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
@click.option('--path_out', default=None, type=click.Path())
@click.argument('config_file', default=None, required=False, type=click.Path())
def main(config_file, path_out):
    """
    Args:
        config_file (str): Config file name.  Default is ``sim0.yml``.
            (fiducial parameters from Andrews et al. 2017).
        path_out (str): Path to output dir.  Default is current dir.
    """
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
    path_flexce_root = os.path.abspath(join(path_flexce, '..'))
    path_data = join(path_flexce, 'data')

    if config_file is None:
        print('No config file specified. Using default parameters.')
        config_file = join(path_flexce_root, 'examples', 'sim0.yml')

    if path_out is None:
        print('No output path specified. Using current working directory.')
        path_out = os.getcwd()

    params = yml.read(config_file)

    mass_bins = flexce.utils.set_mass_bins(**params['mass_bins'])

    params['yld_args'] = flexce.utils.set_yields(params['yields'])
    ylds = Yields(params=params['yields'], mass_bins=mass_bins, path=path_data)

    gal = ChemEvol(params, ylds)
    gal.save(path=path_out)


if __name__ == '__main__':
    main()
