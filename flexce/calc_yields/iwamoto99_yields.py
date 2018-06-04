"""Convert Iwamoto et al. (1999) W70 SNIa yields to pickled arrays."""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

import numpy as np
import pandas as pd
import string

# ---- Set Paths -----
path_calc_yields = join(os.path.abspath(os.path.dirname(__file__)), '')
path_flexce = join('/'.join(path_calc_yields.split('/')[:-2]), '')
path_io = join(path_flexce, 'io')
path_data = join(path_flexce, 'data')
path_yields = join(path_data, 'yields')
path_yldgen = join(path_yields, 'general')
path_i99 = join(path_yields, 'iwamoto99')
sys.path.append(path_io)
# -------------------

from pickle_io import pck_write


# ---- Iwamoto et al. (1999) W70 SNIa yields ----

"""models:
w7: single-degenerate model from Nomoto et al. (1984)
w70: zero-metallicity version of w7
wdd1:
wdd2:
wdd3:
cdd1:
cdd2:
"""

# Read in yields
snia_in = join(path_i99, 'iwamoto99.txt')
cols = ['sym_in', 'w7', 'w70', 'wdd1', 'wdd2', 'wdd3', 'cdd1', 'cdd2']
i99 = pd.read_csv(snia_in, delim_whitespace=True, skiprows=2,
                  usecols=[0, 2, 3, 4, 5, 6, 7, 8], names=cols)
snia_sym = np.array([''.join((item[2:], item[:2])) for item in i99.sym_in])

# Read in isotopes
species_in = pd.read_csv(join(path_yldgen, 'species.txt'),
                         delim_whitespace=True, skiprows=1, usecols=[1],
                         names=['name'])
species = np.array(species_in['name'])
n_species = len(species)
el_name = np.array([item.strip(string.digits) for item in species])
iso_mass = np.array([int(item.strip(string.ascii_letters))
                     for item in species])
elements = []
for item in el_name:
    if item not in elements:
        elements.append(item)

elements = np.array(elements)


# Solar abundances are prevalence "by mass"
solar_isotopes = pd.read_csv(join(path_yldgen, 'Solar_isotopes.txt'),
                             delim_whitespace=True, skiprows=1,
                             usecols=[0, 1], names=['name', 'ab'])
solar_iso = np.array(solar_isotopes['name'])
solar_ab = np.array(solar_isotopes['ab'])
solar_el_name = np.array([item.strip(string.digits) for item in solar_iso])


# indices within "species" array of the elements for which CL04 give a solar
# abundance
ind_iso = []
for i in range(len(solar_iso)):
    ind_iso.append(np.where(solar_iso[i] == species)[0][0])

ind_iso = np.array(ind_iso)


# map elemental yields onto the dominant isotope
snia_yields = {}
for k in i99.keys():
    if k is not 'sym_in':
        snia_yields[k] = np.zeros(n_species)
        for j in range(n_species):
            if species[j] in snia_sym:
                snia_yields[k][j] = i99[k].ix[np.where(snia_sym == species[j])]


# write to file
for k in snia_yields.keys():
    pck_write(snia_yields[k], join(path_i99, k + '_yields.pck'))
