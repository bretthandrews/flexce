"""Create pickled dictionary of neutron-capture isotopic abundances.

Read in the isotopic abundances of neutron-capture isotopes from
Sneden et al. (2008) and create a pickled dictionary.
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

import numpy as np

#---- Set Paths -----
path_calc_yields = join(os.path.abspath(os.path.dirname(__file__)), '')
path_flexce = join('/'.join(path_calc_yields.split('/')[:-2]), '')
path_fileio = join(path_flexce, 'fileio')
path_data = join(path_flexce, 'data')
path_yields = join(path_data, 'yields')
path_yldgen = join(path_yields, 'general')
sys.path.append(path_fileio)
# -------------------

from pickle_io import pickle_write

# ----- Load Data -----
s08 = {}
fin = open(join(path_yldgen, 'sneden08.txt'), 'r')
for line in fin:
    if line.split()[0] == '#':
        pass
    else:
        cols = line.strip().split()
        if len(cols) == 5:
            name, at_num, at_mass, ns, nr = cols
            s08[name] = {'element': name, 'Z': int(at_num),
                         'Isotope': [int(at_mass)], 'N[s]': [float(ns)],
                         'N[r]': [float(nr)]}
        else:
            at_mass, ns, nr = cols
            s08[name]['Isotope'].append(int(at_mass))
            s08[name]['N[s]'].append(float(ns))
            s08[name]['N[r]'].append(float(nr))

fin.close()
# ---------------------


for e in s08.keys():
    for k in ['N[r]', 'N[s]', 'Isotope']:
        s08[e][k] = np.array(s08[e][k])


# Calculate the fraction of each element produced in the r- and s-processes
# ('fraction[r]' and 'fraction[s]') and the fraction of each isotope that
# comes from the r- or s-processes relative to the total amount of the element
# produced in the r- or s-processes('isotopic_fraction[r]' and
# 'isotopic_fraction[s]')
for e in s08.keys():
    ncap_tot = np.sum(s08[e]['N[r]'].sum() + s08[e]['N[s]'].sum())
    s08[e]['fraction[r]'] = s08[e]['N[r]'].sum() / ncap_tot
    s08[e]['fraction[s]'] = s08[e]['N[s]'].sum() / ncap_tot
    if s08[e]['N[r]'].sum() > 0.:
        s08[e]['isotopic_fraction[r]'] = s08[e]['N[r]'] / s08[e]['N[r]'].sum()
    else:
        s08[e]['isotopic_fraction[r]'] = np.zeros(len(s08[e]['N[r]']))
    if s08[e]['N[s]'].sum() > 0.:
        s08[e]['isotopic_fraction[s]'] = s08[e]['N[s]'] / s08[e]['N[s]'].sum()
    else:
        s08[e]['isotopic_fraction[s]'] = np.zeros(len(s08[e]['N[s]']))


pickle_write(s08, join(path_yldgen, 'sneden08.pck'))
