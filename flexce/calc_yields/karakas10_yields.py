"""Generate finely spaced grid of AGB yields using the Karakas (2010) yields.

Karakas & Lattanzio (2010): M = 1.0, 1.25, 1.5, 1.75, 1.9, 2.0, 2.1, 2.25,
2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 6.5 (Zsolar only); Z = 0.0001, 0.004, 0.008,
0.02
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys
import copy

import numpy as np
from scipy import interpolate
import pandas as pd


#---- Set Paths -----
path_calc_yields = join(os.path.abspath(os.path.dirname(__file__)), '')
path_flexce = join('/'.join(path_calc_yields.split('/')[:-2]), '')
path_fileio = join(path_flexce, 'fileio')
path_data = join(path_flexce, 'data')
path_yields = join(path_data, 'yields')
path_yldgen = join(path_yields, 'general')
path_k10 = join(path_yields, 'karakas10', 'iso_yields')
sys.path.append(path_fileio)
#-------------------

from pickle_io import pickle_write

#----- Read in Computed Yields -----

# K10 yields: 61 metallicity/mass combinations; 4 metallicities: Z = 1e-4,
# 4e-3, 8e-3, 2e-2; 17 masses (but not all masses at each metallicity): 1.0,
# 1.25, 1.5, 1.75, 1.9, 2.0, 2.1, 2.25, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
# 6.5 Msun; 19 elements (74 species)

# K10 yields do not include synthetic thermal pulses (see K10 Section 4),
# whereas KL07 recommend using the lambda_last = lambda_final yields (see KL07
# Section 4)that account for extra thermal pulses that likely occur near the
# tip of the AGB where the code has difficulty converging (but there is still a
# large envelope remaining).

# K10 computed yields for the 6 Msun, 1e-4 Zsun model assuming two different
# mass loss prescriptions (VW93 and Reimer's).  I chose to use the Reimer's
# mass loss because that prescription is used for the other intermediate mass
# models (3.0--5.5 Msun) at that metallicity.  See K10 Table 1 for mass loss
# prescriptions for each model.

# The K10 "yield" for a given species refers to actual mass of that species
# returned to the ISM minus the mass of that species that would have been
# returned if the wind had the initial composition of the star (i.e., the net
# yield).  When I implement mass return in the simulation, I need to calculate
# how much total mass was lost, multiply this by the initial composition of the
# star, and then add the net yield (K10 call this the "net yield" in their data
# table), which can be negative for species that are destroyed.

data = pd.read_csv(join(path_k10, 'tablea2.dat'), delim_whitespace=True,
                   usecols=[3, 4], names=['spec', 'at_mass'])
species = np.array(data.spec)[:77]
atomic_mass = np.array(data.at_mass)[:77]

cols = ['mgrid', 'zgrid', 'mrem', 'yield', 'mlost', 'minit']
kwargs = dict(delim_whitespace=True, usecols=(0, 1, 2, 5, 6, 7), names=cols)
z22tab = np.array(pd.read_csv(join(path_k10, 'tablea2.dat'), **kwargs)).T
z83tab = np.array(pd.read_csv(join(path_k10, 'tablea3.dat'), **kwargs)).T
z43tab = np.array(pd.read_csv(join(path_k10, 'tablea4.dat'), **kwargs)).T
z14tab = np.array(pd.read_csv(join(path_k10, 'tablea5.dat'), **kwargs)).T

tab = [z14tab, z43tab, z83tab, z22tab]
zm_grid = []  # metallicities and masses of the models
m_grid = []  # masses of the models
z_grid = []  # metallicities of the models
m_rem = []  # remnant (final) mass
for i in range(4):
    for j in range(tab[i].shape[1]):
        ztmp_str = '{:.0e}'.format(tab[i][1, j])[0] + \
                   '{:.0e}'.format(tab[i][1, j])[-1]
        zmtmp = ''.join(('z', ztmp_str, 'm', str(tab[i][0, j])))
        if zmtmp in zm_grid:
            pass
        else:
            zm_grid.append(zmtmp)
            m_grid.append(tab[i][0, j])
            z_grid.append(tab[i][1, j])
            m_rem.append(tab[i][2, j])

zm_grid = np.array(zm_grid)
m_grid = np.array(m_grid)
z_grid = np.array(z_grid)
m_rem = np.array(m_rem)
m_ejected = m_grid - m_rem  # ejected mass

# net yields of each species
yields_spec = np.concatenate([z14tab[3], z43tab[3], z83tab[3], z22tab[3]])
yields_spec.resize(62, 77)
# delete model with Z = 1e-4, M = 6 and VW93 mass loss
yields_spec = np.delete(yields_spec, np.s_[15], axis=0)

# mass lost of each species
m_lost_spec = np.concatenate([z14tab[4], z43tab[4], z83tab[4], z22tab[4]])
m_lost_spec.resize(62, 77)
m_lost_spec = np.delete(m_lost_spec, np.s_[15], axis=0)

# mass lost of each species if wind had the initial composition of the star
m_init_wind_spec = np.concatenate([z14tab[5], z43tab[5], z83tab[5], z22tab[5]])
m_init_wind_spec.resize(62, 77)
m_init_wind_spec = np.delete(m_init_wind_spec, np.s_[15], axis=0)

k10_m = np.array([1.0, 1.25, 1.5, 1.75, 1.9, 2.0, 2.1, 2.25, 2.5, 3.0, 3.5,
                  4.0, 4.5, 5.0, 5.5, 6.0, 6.5])


# Included 'al-6' but not 'al*6' [always = 0], which are al26, in the total
# Al mass ejected.  Renamed 'al-6' as 'al26'.  Does not keep track of 'g'
# (elements from ni64 to Bi).  s34 includes the abundances of species from s34
# to Mn, and the total S mass ejected includes these other species, although
# they almost always are sub-dominant to s32 by about a factor of 20.)

k10_iso = np.delete(copy.deepcopy(species[2:]), np.s_[39])
k10_iso[38] = 'al26'
n_sym = len(k10_iso)

# net yields for 74 isotopes
yields = np.delete(copy.deepcopy(yields_spec), np.s_[0, 1, 41], axis=1)

# mass lost for 74 isotopes
m_lost = np.delete(copy.deepcopy(m_lost_spec), np.s_[0, 1, 41], axis=1)


# -----------------------------


# ---- Extrapolated Yields -----

# --- Major grid points (filling out computed grid)

# -- (M > 6 Msun) ---

# linearly extrapolate yields up to 8 Msun

extrap_yld_z14 = np.zeros((4, n_sym))
extrap_yld_z43 = np.zeros((4, n_sym))
extrap_yld_z83 = np.zeros((4, n_sym))
extrap_yld_z22 = np.zeros((3, n_sym))

m_extrap = np.arange(6.5, 8.1, 0.5)
zm_grid_extrap = np.array(['z14m6.5', 'z14m7.0', 'z14m7.5', 'z14m8.0',
                           'z43m6.5', 'z43m7.0', 'z43m7.5', 'z43m8.0',
                           'z83m6.5', 'z83m7.0', 'z83m7.5', 'z83m8.0',
                           'z22m7.0', 'z22m7.5', 'z22m8.0'])
z_grid_extrap = np.concatenate((np.ones(4) * 1e-4, np.ones(4) * 4e-3,
                                np.ones(4) * 8e-3, np.ones(3) * 2e-2))
m_grid_extrap = np.concatenate((m_extrap, m_extrap, m_extrap, m_extrap[1:]))
mej_extrap = np.zeros(len(m_grid_extrap))
ind_extrap = [np.arange(4), np.arange(4, 8), np.arange(8, 12),
              np.arange(12, 15)]
ind_grid = [[13, 14], [28, 29], [43, 44], [59, 60]]
for i in range(4):
    itmp = interpolate.InterpolatedUnivariateSpline(m_grid[ind_grid[i]],
                                                    m_ejected[ind_grid[i]],
                                                    k=1)
    mej_extrap[ind_extrap[i]] = itmp(m_grid_extrap[ind_extrap[i]])

for k in range(n_sym):
    itmp = interpolate.InterpolatedUnivariateSpline(m_grid[13:15],
                                                    yields[13:15, k], k=1)
    extrap_yld_z14[:, k] = itmp(m_extrap)

for k in range(n_sym):
    itmp = interpolate.InterpolatedUnivariateSpline(m_grid[28:30],
                                                    yields[28:30, k], k=1)
    extrap_yld_z43[:, k] = itmp(m_extrap)

for k in range(n_sym):
    itmp = interpolate.InterpolatedUnivariateSpline(m_grid[43:45],
                                                    yields[43:45, k], k=1)
    extrap_yld_z83[:, k] = itmp(m_extrap)

for k in range(n_sym):
    itmp = interpolate.InterpolatedUnivariateSpline(m_grid[59:],
                                                    yields[59:, k], k=1)
    extrap_yld_z22[:, k] = itmp(m_extrap[1:])

# -----------------------------


# ---- chemical evolution model mass bins
# IMF
alpha = 2.35
Gamma = 1. - alpha
alpha2 = 2. - alpha
m_min = 0.1
m_max = 100.
a = alpha2 / (m_max**alpha2 - m_min**alpha2)
m_cutoff = 8.

# Bins of Stars
'''Bin lower bounds (bins, bins_low, bins_high). Bin width (dbin_low,
dbin_high).  Number of bins (n_bins, n_bins_low, n_bins_high).  Average mass
per bin (m_ave_high, m_ave_low, m_ave), fraction of total mass (f_mtot).
Fraction of total mass in a stellar generation going into each mass bin (f_m,
f_m_low, f_m_high).  '''

dbin_low = 0.1
bins_low = np.arange(m_min, m_cutoff, dbin_low)
n_bins_low = len(bins_low)
m_ave_low = (Gamma / alpha2) * \
            ((bins_low + dbin_low)**alpha2 - bins_low**alpha2) / \
            ((bins_low + dbin_low)**Gamma - bins_low**Gamma)

dbin_high = 1.
bins_high = np.arange(m_cutoff, m_max, dbin_high)
n_bins_high = len(bins_high)
m_ave_high = (Gamma / alpha2) * \
             ((bins_high + dbin_high)**alpha2 - bins_high**alpha2) / \
             ((bins_high + dbin_high)**Gamma - bins_high**Gamma)

m_ave = np.append(m_ave_low, m_ave_high)
# -----



# ---- Minor grid points (mass bins spaced in ~0.1 Msun, but at the original 4
# ---- metallicity values [Z = 1e-4, 4e-3, 8e-3, 2e-2])

# Interpolate across mass to generate yields at each mass bin of m_ave_low for
# 4 metallicity values: Z = 1e-4, 4e-3, 8e-3, 2e-2 (solar)

# computed and extrapolated yields in coarse grid
k10_grid_yld_z14 = np.concatenate((yields[:15], extrap_yld_z14), axis=0)
k10_grid_yld_z43 = np.concatenate((yields[15:30], extrap_yld_z43), axis=0)
k10_grid_yld_z83 = np.concatenate((yields[30:45], extrap_yld_z83), axis=0)
k10_grid_yld_z22 = np.concatenate((yields[45:], extrap_yld_z22), axis=0)

k10_grid_yld_list = [k10_grid_yld_z14, k10_grid_yld_z43, k10_grid_yld_z83,
                     k10_grid_yld_z22]

k10_interp_mass = np.zeros((4, n_bins_low, n_sym))
k10_interp_mej = np.zeros((4, n_bins_low))
k10_interp_rem = np.zeros((4, n_bins_low))

ind1 = [np.arange(15), np.arange(15, 30), np.arange(30, 45), np.arange(45, 61)]
ind2 = [np.arange(4), np.arange(4), np.arange(4), np.arange(1, 4)]
for i in range(4):
    m_tmp = np.zeros((n_bins_low, n_sym))
    m_gr_tmp = np.concatenate((m_grid[ind1[i]], m_grid_extrap[ind2[i]]))
    for k in range(n_sym):
        itmp = interpolate.InterpolatedUnivariateSpline(
            m_gr_tmp, k10_grid_yld_list[i][:, k], k=1)
        m_tmp[:, k] = itmp(m_ave_low)
    k10_interp_mass[i] = m_tmp
    # mass ejected
    m_ej_tmp = np.concatenate((m_ejected[ind1[i]], mej_extrap[ind2[i]]))
    itmp = interpolate.InterpolatedUnivariateSpline(m_gr_tmp, m_ej_tmp, k=1)
    k10_interp_mej[i] = itmp(m_ave_low)
    k10_interp_mej[i][np.where(k10_interp_mej[i] < 0.)] = 0.
    # remnant mass
    k10_interp_rem[i] = m_ave_low - k10_interp_mej[i]


# ---- Interpolate across metallicity to generate yields at each mass bin of
# ---- m_ave_low for N = n_metal_bin metallicity values

# use same metallicity grid as Limongi & Chieffi SN yields
n_metal_bin = 1001
z_grid2 = np.array([1e-6, 1e-4, 1e-3, 6e-3, 2e-2])
logz_grid2 = np.log10(z_grid2)
z_grid3 = np.array([1e-4, 4e-3, 8e-3, 2e-2])

# evenly sample metallicity (in log Z) between grid points
logz_final = np.zeros(n_metal_bin)
dind = (n_metal_bin - 1) / (len(z_grid2) - 1)
for i in range(len(z_grid2) - 1):
    dlogz = (logz_grid2[i+1] - logz_grid2[i]) / \
            ((n_metal_bin - 1) / (len(z_grid2) - 1))
    logz_final[i*dind:i*dind+dind+1] = np.arange(logz_grid2[i],
                                                 logz_grid2[i+1]+1e-9, dlogz)

# metallicity of final grid
z_final = 10.**logz_final


# output interpolated yields
k10_final = np.zeros((n_metal_bin, n_bins_low, n_sym))
k10_final_mej = np.zeros((n_metal_bin, n_bins_low))
k10_final_rem = np.zeros((n_metal_bin, n_bins_low))
# at each mass, interpolate each element for each metallicity
for i in range(n_bins_low):
    for j in range(n_sym):
        itmp = interpolate.InterpolatedUnivariateSpline(
            z_grid3, k10_interp_mass[:, i, j], k=1)
        k10_final[:, i, j] = itmp(z_final)
    # mass ejected
    itmp = interpolate.InterpolatedUnivariateSpline(z_grid3,
                                                    k10_interp_mej[:, i], k=1)
    k10_final_mej[:, i] = itmp(z_final)
    k10_final_mej[:, i][np.where(k10_final_mej[:, i] < 0.)] = 0.
    # remnant mass
    k10_final_rem[:, i] = (np.ones(n_metal_bin) * m_ave_low[i] -
                           k10_final_mej[:, i])


# pickle the interpolated yields array and the metallicity grid used
pickle_write(k10_final, join(path_k10, 'interp_yields.pck'))
pickle_write(k10_final_mej, join(path_k10, 'interp_meject.pck'))
pickle_write(k10_final_rem, join(path_k10, 'interp_mremnant.pck'))
# -----------------------------
