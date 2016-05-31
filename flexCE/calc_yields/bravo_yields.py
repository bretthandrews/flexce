"""Convert Eduardo Bravo (in prep.) SNIa yields to pickled arrays."""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import pandas as pd

# ---- Set Paths -----
path_calc_yields = join(os.path.abspath(os.path.dirname(__file__)), '')
path_flexce = join('/'.join(path_calc_yields.split('/')[:-2]), '')
path_fileio = join(path_flexce, 'fileio')
path_data = join(path_flexce, 'data')
path_yields = join(path_data, 'yields')
path_yldgen = join(path_yields, 'general')
path_bravo = join(path_yields, 'bravo')
sys.path.append(path_fileio)
# -------------------

from pickle_io import pickle_write


# ---- Eduardo Bravo (in prep.) SNIa yields ----

# models: DDTa, DDTc, DDTe, DDTf

# Read in yields
ddt_in = join(path_bravo, 'DDT_Yields_Z.dat')
data = {}
bravo_el = None
with open(ddt_in, 'r') as infile:
    for line in infile:
        if 'DDT' in line:
            model_name, model_metallicity = line.strip().split()
            if model_name not in data.keys():
                data[model_name] = {}
                data[model_name][model_metallicity] = {}
        else:
            parsed = line.strip().split('   ')
            tmp_sym = []
            tmp_yld = []
            for it in parsed:
                tmp_sym_ind, tmp_yld_ind = it.split(': ')
                tmp_sym.append(tmp_sym_ind)
                tmp_yld.append(float(tmp_yld_ind))
            if bravo_el is None:
                bravo_el = tmp_sym
            data[model_name][model_metallicity] = tmp_yld

bravo_el = np.array(bravo_el)
bravo_met = [0.00025, 0.0025, 0.01, 0.025, 0.075]
n_bravo_met = len(bravo_met)

# Read in isotopes
species_in = pd.read_csv(join(path_yldgen, 'species.txt'),
                         delim_whitespace=True, skiprows=1, usecols=[1],
                         names=['name'])
species = np.array(species_in['name'])
n_species = len(species)

# Solar abundances are prevalence "by mass"
solar_isotopes = pd.read_csv(join(path_yldgen, 'Solar_isotopes.txt'),
                             delim_whitespace=True, skiprows=1,
                             usecols=[0, 1], names=['name', 'ab'])

# Map elemental yields onto dominant isotope
snia_sym = []
for it in bravo_el:
    ind_tmp = []
    for ii, siso in enumerate(solar_isotopes.name):
        # remove numbers from isotope name
        elname = ''.join([i for i in siso if not i.isdigit()])
        if it == elname:
            ind_tmp.append(ii)
    ind_dominant = np.argmax(solar_isotopes.ab[ind_tmp])
    snia_sym.append(solar_isotopes.name[ind_dominant])

snia_sym = np.array(snia_sym)

DDTa = pd.DataFrame(data['DDTa'], index=snia_sym)
DDTc = pd.DataFrame(data['DDTc'], index=snia_sym)
DDTe = pd.DataFrame(data['DDTe'], index=snia_sym)
DDTf = pd.DataFrame(data['DDTf'], index=snia_sym)

# Bravo yields projected onto the isotope list of the CCSN yields
DDTa_iso = np.zeros((n_bravo_met, n_species))
DDTc_iso = np.zeros((n_bravo_met, n_species))
DDTe_iso = np.zeros((n_bravo_met, n_species))
DDTf_iso = np.zeros((n_bravo_met, n_species))

# Project yields onto isotope list of the CCSN yields
snia_yields = {}  # metallicity independent yields
models = [DDTa, DDTc, DDTe, DDTf]
models_iso = [DDTa_iso, DDTc_iso, DDTe_iso, DDTf_iso]
model_names = ['DDTa', 'DDTc', 'DDTe', 'DDTf']
for mod, imod, mname in zip(models, models_iso, model_names):
    for ii, bmet in enumerate(bravo_met):
        met = str(bmet)
        mzname = '_z'.join((mname, met))
        snia_yields[mzname] = np.zeros(n_species)
        for jj in range(n_species):
            if species[jj] in snia_sym:
                snia_yields[mzname][jj] = mod[met].ix[snia_sym == species[jj]]
                imod[ii][jj] = mod[met].ix[snia_sym == species[jj]]


# Metallicity dependent yields

# use same metallicity grid as Limongi & Chieffi SN yields
n_metal_bin = 1001
z_grid_lc = np.array([1e-6, 1e-4, 1e-3, 6e-3, 2e-2])
logz_grid_lc = np.log10(z_grid_lc)
z_grid_bravo = np.array([0.00025, 0.0025, 0.01, 0.025, 0.075])

# evenly sample metallicity (in log Z) between grid points
logz_final = np.zeros(n_metal_bin)
dind = (n_metal_bin - 1) / (len(z_grid_lc) - 1)
for i in range(len(z_grid_lc) - 1):
    dlogz = (logz_grid_lc[i+1] - logz_grid_lc[i]) / \
            ((n_metal_bin - 1) / (len(z_grid_lc) - 1))
    logz_final[i*dind:i*dind+dind+1] = np.arange(logz_grid_lc[i],
                                                 logz_grid_lc[i+1]+1e-9, dlogz)

# metallicity of final grid
z_final = 10.**logz_final


# Bravo yields projected onto the isotope list and metallicity grid of the CCSN
# yields
# Only goes up to solar metallicity!
# TODO Should I interpolate the SNIa yields to super-solar metallicities?
DDTa_final = np.zeros((n_metal_bin, n_species))
DDTc_final = np.zeros((n_metal_bin, n_species))
DDTe_final = np.zeros((n_metal_bin, n_species))
DDTf_final = np.zeros((n_metal_bin, n_species))
models_final = [DDTa_final, DDTc_final, DDTe_final, DDTf_final]
for imod, fmod in zip(models_iso, models_final):
    for jj in range(n_species):
        itmp = InterpolatedUnivariateSpline(z_grid_bravo, imod[:, jj], k=1)
        fmod[:, jj] = itmp(z_final)


# write Bravo yields to file for each model--metallicity combination
for k in snia_yields.keys():
    pickle_write(snia_yields[k], join(path_bravo, k + '_yields.pck'))


# write Bravo yields to file for each model interpolated onto the CCSN
# metallicity grid
for fmod, mname in zip(models_final, model_names):
    pickle_write(fmod, join(path_bravo, mname + '_yields.pck'))
