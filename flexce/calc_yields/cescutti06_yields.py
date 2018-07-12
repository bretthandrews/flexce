# @Author: Brett Andrews <andrews>
# @Date:   2018-06-21 11:06:54
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 11:07:15

"""
FILE
    cescutti06.py

DESCRIPTION
    Generate a finely spaced grid of Ba & Eu yields from SNII.

    Generate a finely spaced grid of Ba & Eu yields from SNII from the
    'empirical' r-process yields (determined by matching Ba and Eu
    abundances from a chemical evolution model that used these yields
    to observations) of Cescutti et al. (2006).
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join

import numpy as np
from scipy import interpolate
import pandas as pd

# ---- Set Paths -----
path_calc_yields = join(os.path.abspath(os.path.dirname(__file__)), '')
path_flexce = join('/'.join(path_calc_yields.split('/')[:-2]), '')
path_data = join(path_flexce, 'data')
path_yields = join(path_data, 'yields')
path_yldgen = join(path_yields, 'general')
path_c06 = join(path_yields, 'cescutti06')
# -------------------

from flexce.fileio import pck


# ---- Load Data ----
data_in = pd.read_csv(join(path_c06, 'cescutti06_rprocess.txt'),
                      delim_whitespace=True, skiprows=5,
                      names=['M', 'XBa', 'XEu'])
m_c06 = np.array(data_in['M'])
x_c06 = np.array(data_in[['XBa', 'XEu']]).T

c06_orig = x_c06 * m_c06  # convert X_Ba & X_Eu to a mass of Ba & Eu
# -------------------


# -----
# chemical evolution model mass bins
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
dbin_high = 1.
bins_high = np.arange(m_cutoff, m_max, dbin_high)
n_bins_high = len(bins_high)
m_ave_high = (Gamma / alpha2) * \
             ((bins_high + dbin_high)**alpha2 - bins_high**alpha2) / \
             ((bins_high + dbin_high)**Gamma - bins_high**Gamma)
# ----


# ---- Interpolated Yields -----

# ---- Minor grid points (interpolate across mass to generate yields at each
# ---- mass bin of m_ave_low at the original metallicity values)

# yields = 0 for M < 12 Msun and M > 30 Msun

c06_interp_metal = np.zeros((n_bins_high, len(c06_orig)))
for i in range(len(c06_orig)):
    itmp = interpolate.InterpolatedUnivariateSpline(m_c06, c06_orig[i], k=1)
    m_tmp = itmp(m_ave_high)
    m_tmp[np.where(m_ave_high < 12.)] = 0.
    m_tmp[np.where(m_ave_high > 30.)] = 0.
    m_tmp[np.where(m_tmp < 0.)] = 0.
    c06_interp_metal[:, i] = m_tmp


# project into 3D arrays (metallicities, masses, isotopes)
z_final = pck.read(join(path_yldgen, 'interp_metallicity.pck'))
n_metal_bin = len(z_final)
c06_final = np.ones((n_metal_bin, n_bins_high, len(c06_orig)))
c06_final[:] = c06_interp_metal

# -----------------------------

# pickle the interpolated yields array and the metallicity grid used
pck.write(c06_final, join(path_c06, 'cescutti06_yields.pck'))
