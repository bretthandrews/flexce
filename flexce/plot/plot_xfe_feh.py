"""
FILE
    plot_xfe_feh.py

USAGE
    python plot_xfe_feh.py [config_file]

DESCRIPTION
    Plot [X/Fe]--[Fe/H], MDF, and [X/Fe] distribution function.
"""


raise ValueError('Decouple package and output directory structures.')  # .............................................


from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

import matplotlib.pyplot as plt
import seaborn as sns

# ----- Set paths -----
try:
    path_plot = os.path.abspath(os.path.dirname(__file__))
except NameError as e:
    path_plot = os.getcwd()

path_flexce_top = os.path.abspath(join(path_plot, '../..'))
path_flexce = join(path_flexce_top, 'flexce')
path_io = join(path_flexce_top, 'flexce', 'io')
path_plots = join(path_flexce_top, 'plots')
# ---------------------

sys.path.insert(0, path_flexce)
import utils
import plot.utils as putils
from flexce.io import cfg, txt

default_config_path = join(path_plots, 'config')
default_output_path = join(path_flexce_top, 'output')
fin, path_config = utils.set_path(sys.argv[1], default_config_path)

try:
    stem = path_config.split('config/')[1]
except IndexError:
    stem = ''

path_output = join(default_output_path, stem)

path_plot_out = utils.substitute_dir_in_path(path_config, 'config', 'plots')
if not os.path.isdir(path_plot_out):
    os.makedirs(path_plot_out)

# Read config file
cfg_in = cfg.read_plot_config(join(path_config, fin))
colors = putils.get_colors(cfg_in)
abund = cfg_in['General']['abundance']
labels = cfg_in['General']['labels']

# Read in simulation results
sims = []
for sim_id in cfg_in['General']['sim_ids']:
    sims.append(txt.load_dataframe(path_output, sim_id))

# Make plot
for i, (sim, color) in enumerate(zip(sims, colors)):
    marg_kws = dict(norm_hist=True,
                    hist_kws=dict(weights=sim.Survivors.values))
    if i == 0:
        fig = sns.jointplot('[Fe/H]', abund, data=sim, stat_func=None,
                            color=color, marginal_kws=marg_kws)
    else:
        fig = putils.joint_overplot('[Fe/H]', abund, df=sim, fig=fig,
                                    color=color, marg_kws=marg_kws)

# Make legend
p = putils.get_path_collections(fig)
leg_args = putils.get_leg_args(cfg_in)
leg = fig.ax_joint.legend(p, labels, **leg_args)

# Save plot
fout = ''.join((fin.strip('.cfg'), '.pdf'))
plt.savefig(join(path_plot_out, fout))
