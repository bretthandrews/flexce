"""Read in simulation output for inspection in interactive session."""

import os
from os.path import join
from flexce.fileio import pickle_io
from fileio.fileio import txt_io
from flexce.chemevol import ChemEvol


path_output = join(os.expanduer('~'), 'flexce_output')

# Read in pickled box and ab objects individually...
box0 = pickle_io.box_read(path_output, sim_id=0)
ab0 = pickle_io.ab_read(path_output, sim_id=0)
# ...or simultaneously
box1, ab1 = pickle_io.sim_read(path_output, sim_id=1)

# Read in timesteps, surviving stars, and abundances of select elements
time0, survivors0, feh0, abunds0 = txt_io.txt_read(path_output, sim_id=0)
# abunds0 is a recarray with element abbreviation (e.g., 'C') as the column
# names (instead of '[C/Fe]').
