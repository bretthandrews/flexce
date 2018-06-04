"""Read in simulation output for inspection in interactive session."""

import os
from os.path import join
from flexce.io import pck
from flexce.io import txt
from flexce.chemevol import ChemEvol


path_output = join(os.expanduer('~'), 'flexce_output')

# Read in pickled box and ab objects individually...
box0 = pck.box_read(path_output, sim_id=0)
ab0 = pck.ab_read(path_output, sim_id=0)
# ...or simultaneously
box1, ab1 = pck.sim_read(path_output, sim_id=1)

# Read in timesteps, surviving stars, and abundances of select elements
time0, survivors0, feh0, abunds0 = txt.txt_read(path_output, sim_id=0)
# abunds0 is a recarray with element abbreviation (e.g., 'C') as the column
# names (instead of '[C/Fe]').
