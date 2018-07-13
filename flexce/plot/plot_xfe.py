# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 16:07:07
# @Last modified by:   andrews
# @Last modified time: 2018-07-13 09:07:81

"""
FILE
    plot_xfe.py

USAGE
    plot_xfe [config_file]

DESCRIPTION
    Plot [X/Fe]--[Fe/H], MDF, and [X/Fe] distribution function.
"""

import os
from os.path import join

import click
import matplotlib.pyplot as plt
import seaborn as sns

from flexce.fileio import txt, yml
import flexce.plot.utils as putils


@click.command()
@click.option('--path_sim', default=None, type=click.Path())
@click.option('--path_out', default=None, type=click.Path())
@click.argument('config_file', default=None, required=False, type=click.Path())
def main(config_file, path_sim, path_out):
    """
    Args:
        config_file (str): Config file name.  Default is
            ``ofe_sim0.yml``, which plots [O/Fe]--[Fe/H] for fiducial
            parameters from Andrews et al. (2017).
        path_sim (str): Path to simulations.  Default is current dir.
        path_out (str): Path to output dir.  Default is current dir.
    """
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
    path_flexce_root = os.path.abspath(join(path_flexce, '..'))

    if config_file is None:
        config_file = join(path_flexce_root, 'examples', 'ofe_sim0.yml')

    if path_sim is None:
        path_sim = os.getcwd()

    if path_out is None:
        path_out = os.getcwd()

    params = yml.read(config_file)

    abunds = [txt.read_abundances(path_sim, sim_id=sim_id) for sim_id in params['sim_ids']]

    colors = putils.get_colors(params)

    # Make plot
    for ii, (ab, color) in enumerate(zip(abunds, colors)):
        marg_kws = {
            'norm_hist': True,
            'hist_kws': {'weights': ab.Survivors.values}
        }

        if ii == 0:
            fig = sns.jointplot('[Fe/H]', params['ab'], data=ab, stat_func=None,
                                color=color, marginal_kws=marg_kws)

        else:
            fig = putils.joint_overplot('[Fe/H]', params['ab'], data=ab,
                                        fig=fig, color=color, marg_kws=marg_kws)

    # Make legend
    p = putils.get_path_collections(fig)
    leg_args = putils.get_leg_args(params)
    leg = fig.ax_joint.legend(p, params['labels'], **leg_args)

    fout = join(path_out, os.path.splitext(os.path.basename(config_file))[0] + '.pdf')
    plt.savefig(fout)
    print(f'Wrote: {fout}')


if __name__ == '__main__':
    main()
