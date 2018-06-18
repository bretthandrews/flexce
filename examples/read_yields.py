"""Read in fiducial yields for inspection in interactive session."""

import os
from os.path import join

from flexce.io.yml import read_yml
from flexce.utils import set_mass_bins, load_yields

from flexce.yields import Yields

# Edit for docs
path_flexce = join(os.path.expanduser('~'), 'projects', 'flexce')
path_config = join(path_flexce, 'examples')
path_data = join(path_flexce, 'flexce', 'data')
file_in = join(path_config, 'sim0.yml')

params = read_yml(file_in)

mass_bins = set_mass_bins(**params['mass_bins'])
ylds = load_yields(path_data, mass_bins, params['yields'])


ylds = Yields(path_data, mass_bins, params['yields'])
