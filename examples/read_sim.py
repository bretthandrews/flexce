import os
from os.path import join
import sys

# use setup.py file to write flexCE path here?
#
# put flexCE in PYTHONPATH
# path_examples = join(os.path.abspath(os.path.dirname(__file__)), '')
# path_flexce = join('/'.join(path_examples.split('/')[:-2]), '')
# sys.path.append(path_flexce)

path_flexce = '/Users/andrews/flexCE/flexCE/'

from fileio import pickle_io
from fileio import txt_io

# Location of your simulation output directory
path_out = join(path_flexce, 'output', '')

# Read in pickled box and ab objects individually...
box0 = pickle_io.box_read(path_out, sim_id=0)
ab0 = pickle_io.ab_read(path_out, sim_id=0)
# ...or simultaneously
box1, ab1 = pickle_io.sim_read(path_out, sim_id=1)

# Read in timesteps, surviving stars, and abundances of select elements
time0, survivors0, feh0, abunds0 = txt_io.txt_read(path_out, sim_id)
# abunds0 is a recarray with element abbreviation (e.g., 'C') as the column
# names (instead of '[C/Fe]').