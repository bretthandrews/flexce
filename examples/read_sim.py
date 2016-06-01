import os
from os.path import join
import sys

path_flexce_root = join(os.path.expanduser('~'), 'flexCE')
path_flexce = join(path_flexce_root, 'flexCE')
sys.path.insert(0, path_flexce)
from chemevol import ChemEvol
from fileio import pickle_io
from fileio import txt_io

# Location of your simulation output directory
path_output = join(path_flexce_root, 'output', '')

# Read in pickled box and ab objects individually...
box0 = pickle_io.box_read(path_output, sim_id=0)
ab0 = pickle_io.ab_read(path_output, sim_id=0)
# ...or simultaneously
box1, ab1 = pickle_io.sim_read(path_output, sim_id=1)

# Read in timesteps, surviving stars, and abundances of select elements
time0, survivors0, feh0, abunds0 = txt_io.txt_read(path_output, sim_id=0)
# abunds0 is a recarray with element abbreviation (e.g., 'C') as the column
# names (instead of '[C/Fe]').
