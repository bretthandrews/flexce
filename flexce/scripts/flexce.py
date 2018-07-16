# @Author: Brett Andrews <andrews>
# @Date:   2018-05-30 17:05:89
# @Last modified by:   andrews
# @Last modified time: 2018-07-16 11:07:12

"""
FILE
    flexce

USAGE
    flexce config_files --path_out path/to/output

    Example::

        $ flexce sim0.yml --path_out output

DESCRIPTION
    Run flexCE from the command line using config files.
"""

import os
from os.path import join

import click

from flexce.chemevol import ChemEvol


@click.command()
@click.option('--path_out', default=None, type=click.Path())
@click.argument('config_files', default=None, nargs=-1, required=False, type=click.Path())
def main(config_files, path_out):
    """
    Args:
        config_files (str): Config file name.  Default is ``sim0.yml``
            (fiducial parameters from Andrews et al. 2017).
        path_out (str): Path to output dir.  Default is current dir.
    """
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')

    if config_files is None:
        config_files = [join(path_flexce, 'data', 'config', 'sim0.yml')]

    if path_out is None:
        path_out = os.getcwd()

    for cfile in config_files:
        gal = ChemEvol(params=cfile)
        gal.save(path=path_out)


if __name__ == '__main__':
    main()
