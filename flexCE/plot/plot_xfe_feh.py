"""
FILE
    plot_xfe_feh.py

USAGE
    python plot_xfe_feh.py [config_file]

DESCRIPTION
    Plot [X/Fe]--[Fe/H], MDF, and [X/Fe] distribution function.
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

import matplotlib.pyplot as plt
import seaborn as sns

import utils

# ----- Set paths -----
try:
    path_plot = os.path.abspath(os.path.dirname(__file__))
except NameError as e:
    path_plot = os.getcwd()

path_flexce_top = os.path.abspath(join(path_plot, '../..'))
path_output = join(path_flexce_top, 'output')
path_fileio = join(path_flexce_top, 'flexCE', 'fileio')
path_plots_top = join(path_flexce_top, 'plots')
path_config = join(path_plots_top, 'config')
path_plots = join(path_plots_top, 'plots')
# ---------------------

sys.path.append(path_fileio)
import txt_io
import cfg_io

# Read config file
fin = utils.get_filename(sys.argv, path_config)
cfg = cfg_io.read_plot_config(join(path_config, fin))
colors = utils.get_colors(cfg)
abund = cfg['General']['abundance']
labels = cfg['General']['labels']

# Read in simulation results
sims = []
for sim_id in cfg['General']['sim_ids']:
    sims.append(txt_io.load_dataframe(path_output, sim_id))

# Make plot
for i, (sim, color) in enumerate(zip(sims, colors)):
    marg_kws = dict(norm_hist=True,
                    hist_kws=dict(weights=sim.Survivors.values))
    if i == 0:
        fig = sns.jointplot('[Fe/H]', abund, data=sim, stat_func=None,
                            color=color, marginal_kws=marg_kws)
    else:
        fig = utils.joint_overplot('[Fe/H]', abund, df=sim, fig=fig,
                                   color=color, marg_kws=marg_kws)

# Make legend
p = utils.get_path_collections(fig)
leg_args = utils.get_leg_args(cfg)
leg = fig.ax_joint.legend(p, labels, **leg_args)

# Save plot
fout = ''.join((fin.strip('.cfg'), '.pdf'))
plt.savefig(join(path_plots, fout))
