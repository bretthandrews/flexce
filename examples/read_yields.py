"""Read in fiducial yields for inspection in interactive session."""

import os
from os.path import join

from fileio.cfg_io import read_sim_cfg
from utils import define_mass_bins
from utils import load_yields

path_flexce = join(os.getenv('HOME'), 'flexCE', 'flexCE')
path_config = join(path_flexce, 'config')
path_data = join(path_flexce, 'data')
file_in = join(path_config, 'sim0.cfg')

(sim_id, yld_args, initialize_args, mass_bins_args, snia_dtd_args,
 inflows_args, outflows_args, warmgasres_args, sf_args) = read_sim_cfg(file_in)

mass_bins = define_mass_bins(**mass_bins_args)
ylds = load_yields(path_data, yld_args, mass_bins)
