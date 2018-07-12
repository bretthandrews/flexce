"""Read in simulation output for inspection in interactive session."""

import os
from os.path import join
from flexce.fileio import pck, txt


path = join(os.expanduer('~'), 'flexce_output')

# Read in pickled sim instance
gal = pck.read_sim(path, sim_id=0)

# Read in timesteps, surviving stars, and abundances
ab = txt.read_abundances(path, sim_id=0)
