# @Author: Brett Andrews <andrews>
# @Date:   2018-05-30 17:05:89
# @Last modified by:   andrews
# @Last modified time: 2018-07-13 14:07:46

"""
FILE
    run_flexce

USAGE
    run_flexce config_file --path_out path/to/output

    Example::

        $ run_flexce sim0.yml

DESCRIPTION
    Run flexCE from the command line using a config file.
"""

import os
from os.path import join

import click

from flexce.chemevol import ChemEvol


@click.command()
@click.option('--path_out', default=None, type=click.Path())
@click.argument('config_file', default=None, required=False, type=click.Path())
def main(config_file, path_out):
    """
    Args:
        config_file (str): Config file name.  Default is ``sim0.yml``
            (fiducial parameters from Andrews et al. 2017).
        path_out (str): Path to output dir.  Default is current dir.
    """
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')

    if config_file is None:
        config_file = join(path_flexce, 'data', 'config', 'sim0.yml')

    if path_out is None:
        path_out = os.getcwd()

    gal = ChemEvol(params=config_file)
    gal.save(path=path_out)


if __name__ == '__main__':
    main()
