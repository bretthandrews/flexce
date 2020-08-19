"""Generate finely spaced grid of SN II isotopic yields.

Use a combination of the Chieffi & Limongi (2004) & Limongi & Chieffi (2006).

Chieffi & Limongi (2004): M = 13--35 Msun; Z = 0--solar
Limongi & Chieffi (2006): M = 11--120; Z = solar

Mass cut = 0.1 Msun Ni56
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

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
path_cl04 = join(path_yields, 'chieffi04', 'iso_yields')
path_lc06 = join(path_yields, 'limongi06', 'iso_yields')
path_lc06_el = join(path_yields, 'limongi06', 'el_yields')
sys.path.append(path_fileio)
# -------------------

from pickle_io import pickle_write

# ---- Read in Computed Yields -----
# CL04
cl04_files = []
cl04_nik = []
cl04_mass = []
cl04_metal = []
input_Z_cl04 = ['00', '16', '14', '13', '63', '22']
input_M_cl04 = ['13', '15', '20', '25', '30', '35']
input_q_cl04 = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                '20']
for i in range(len(input_Z_cl04)):
    for j in range(len(input_M_cl04)):
        cl04_mass.append(int(input_M_cl04[j]))
        cl04_metal.append(float(''.join((input_Z_cl04[i][0], 'e-',
                                         input_Z_cl04[i][1]))))
        cl04_nik.append(''.join(['z', input_Z_cl04[i], 'm', input_M_cl04[j],
                                 '.nik']))
        ftmp = []
        for k in range(len(input_q_cl04)):
            ftmp.append(''.join(['z', input_Z_cl04[i], 'm', input_M_cl04[j],
                                 '.q', input_q_cl04[k]]))
        cl04_files.append(ftmp)

cl04_files = np.array(cl04_files)
cl04_mass = np.array(cl04_mass)
cl04_metal = np.array(cl04_metal)


# LC06
lc06_files = []
lc06_mass = []
lc06_metal = []
input_M_lc06 = ['011', '012', '013', '014', '015', '016', '017', '020', '025',
                '030', '035', '040', '060', '080', '120']
for i in range(len(input_M_lc06)):
    lc06_files.append(''.join(['m', input_M_lc06[i], 'isoexp.difdeca']))
    lc06_mass.append(int(input_M_lc06[i]))
    lc06_metal.append(2e-2)

lc06_mass = np.array(lc06_mass)
lc06_metal = np.array(lc06_metal)

# LC06 elemental yields files (need the ni56 mass produced by each mass cut)
lc06_el_files = []
for i in range(len(input_M_lc06)):
    lc06_el_files.append(''.join(['z', '22', 'm', input_M_lc06[i],
                                  'elecum.ngms']))


# solar abundance of metals---needed to subtract the initial metal abundances
# of the stellar models (also assume Y = 0.285)---in relative amounts (not
# Msun), that is, sum(solar_ab) = 1.
solar_isotopes = pd.read_csv(join(path_yldgen, 'Solar_isotopes.txt'),
                             delim_whitespace=True, skiprows=1,
                             usecols=[0, 1], names=['name', 'ab'])
solar_iso = np.array(solar_isotopes['name'])
solar_ab = np.array(solar_isotopes['ab'])



# -- Calculate Net Yields ---

# cl04_raw: yields (1) 36 metallicity/mass combinations, (2) over all mass cuts
# for each metallicity/mass combination (number of mass cuts vary amongst
# metallicity/mass combinations), (3) 293 isotopes

# cl04_mcut: mass cuts (1) 36 metallicity/mass combinations, (2) over all mass
# cuts for each metallicity/mass combination (number of mass cuts vary amongst
# metallicity/mass combinations)

# cl04_ni56: mass of Ni56 produced (1) 36 metallicity/mass combinations, (2)
# over all mass cuts for each metallicity/mass combination (number of mass cuts
# vary amongst metallicity/mass combinations)

# cl04_abs_05ni: interpolated absolute yields for mass cut that produces 0.05
# Msun Ni56 (1) 36 metallicity/mass combinations, (2) 293 elements

# cl04_init_comp: initial composition of the stellar models (only 116 isotopes
# given, so all others assumed to be 0)

# cl04_net_05ni: interpolated net yields for mass cut that produces 0.05 Msun
# Ni56 (1) 36 metallicity/mass combinations, (2) 293 elements

# (Similar naming structure for lc06 and lc12 arrays)

mass_ni56 = 0.1 # mass of ni56 produced

# isotope names
species = np.array([])
for j in range(cl04_files.shape[1]):
    itmp = pd.read_csv(join(path_cl04, cl04_files[0, j]),
                       delim_whitespace=True)
    species = np.concatenate((species, np.array(itmp.columns)[1:]))

species = species[:-3]
n_species = len(species)

# indices within "species" array of the elements for which CL04 give a solar
# abundance
ind_iso = []
for i in range(len(solar_iso)):
    ind_iso.append(np.where(solar_iso[i] == species)[0][0])

ind_iso = np.array(ind_iso)


# CL04 absolute yields
cl04_raw = []
cl04_mcut = []
cl04_ni56 = []
cl04_abs_010ni = []
for i in range(len(cl04_files)):
    # read in mass cuts and ni56 yields
    nitmp = pd.read_csv(join(path_cl04, cl04_nik[i]), delim_whitespace=True,
                        names=['mcut', 'ni56'])
    cl04_mcut.append(np.array(nitmp.mcut))
    cl04_ni56.append(np.array(nitmp.ni56))
    # read in yields (several input files need to be stitched together)
    ytmp = np.ones((len(nitmp), n_species))
    cnt = 0
    for j in range(cl04_files.shape[1]):
        # Skip first two rows because the first row of data corresponds to the
        # surface of the star, which is not present in the nickel yield files
        # (e.g., z00m13.nik)
        qtmp0 = pd.read_csv(join(path_cl04, cl04_files[i, j]),
                            delim_whitespace=True)
        qtmp = np.array(qtmp0.iloc[1:])  # skip first row (surface of star)
        if j == cl04_files.shape[1] - 1:
            qtmp = qtmp[:, :-3]
        ytmp[:, cnt:cnt+qtmp.shape[1]-1] = qtmp[:, 1:]
        cnt += qtmp.shape[1] - 1
    cl04_raw.append(ytmp)
    # interpolated yields for the mass cut that produces the desired ni56 mass
    if i == 0:
        mass_ni56 = 9.13e-2
    yld_tmp = []
    for j in range(n_species):
        yld_tmp.append(interpolate.interp1d(nitmp.ni56,
                                            ytmp[:, j])(mass_ni56))
    cl04_abs_010ni.append(np.array(yld_tmp))
    if i == 0:
        mass_ni56 = 0.1

cl04_abs_010ni = np.array(cl04_abs_010ni)

# CL04 initial composition
cl04_init_comp = np.zeros(cl04_abs_010ni.shape)
for i in range(6):
    indt = np.arange(6*i, 6*i+6)
    if i == 0:
        cl04_init_comp[indt, 0] = (1. - 0.23) * np.sum(cl04_abs_010ni[indt],
                                                       axis=1)  # H
        cl04_init_comp[indt, 4] = (0.23 *
                                   np.sum(cl04_abs_010ni[indt], axis=1))  # He
    else:
        ztmp = cl04_metal[indt][0]
        cl04_init_comp[indt, 0] = ((1. - 0.285 - ztmp) *
                                   np.sum(cl04_abs_010ni[indt], axis=1))  # H
        cl04_init_comp[indt, 4] = (0.285 *
                                   np.sum(cl04_abs_010ni[indt], axis=1))  # He
        for j in range(len(ind_iso)):
            cl04_init_comp[indt, ind_iso[j]] = ztmp * solar_ab[j] * \
                                   np.sum(cl04_abs_010ni[indt], axis=1)  # C->Mo

# CL04 net yields = absolute yields - initial composition of stellar model
cl04_net_010ni = cl04_abs_010ni - cl04_init_comp


# CL04 mass ejected
cl04_mej = np.sum(cl04_abs_010ni, axis=1)

# CL04 remnant mass
cl04_rem = cl04_mass - cl04_mej


# The LC06 isotope files (e.g., m011isoexp.difdeca) give the mass in each shell
# and the mass fraction of each isotope.  To calculate the cumulative mass
# ejected for a given mass cut, multiply the mass of mass of each shell (i.e.,
# the difference in mass between a row and the one above it) by the mass
# fraction in that shell (columns 2-->end) and sum all of the overlying shells
# (the surface point is the first row below the column labels).

# I'm not sure what three of the last four columns (except the very last) of
# the LC06 isotopic yields files mean ('P', 'A', 'N', 'Ni56'), and the 'Ni56'
# column is set to 1e-30, so it provides no information.


# LC06 absolute yields
lc06_raw = []
lc06_mcut = []
lc06_ni56 = []
lc06_abs_010ni = []
for i in range(len(lc06_files)):
    tmp = np.array(pd.read_csv(join(path_lc06, lc06_files[i]),
                               delim_whitespace=True))
    # see note above about calculating cumulative mass ejected
    cum_yld_tmp = np.cumsum(np.multiply(np.diff(tmp[:, 0]),
                                        tmp[1:, 2:].T) * -1, axis=1).T
    cum_yld_tmp[np.where(cum_yld_tmp == 0.)] = 1e-30
    lc06_mcut.append(tmp[1:, 0])
    lc06_raw.append(cum_yld_tmp[:, :-4])  # see not above about last 4 columns
    # read in ni56 mass produced from elemental yield file (mass coordinate
    # stars at center and goes outward---opposite to isotopic yield files)
    tmp2 = np.array(pd.read_csv(join(path_lc06_el, lc06_el_files[i]),
                                delim_whitespace=True, header=None))
    lc06_ni56.append(tmp2[:, -1][::-1])
    # interpolate between the yields from different mass cuts such that 0.05
    # Msun of Ni56 produced
    yld_tmp = []
    for j in range(n_species):
        yld_tmp.append(interpolate.interp1d(tmp2[:, -1][::-1],
                                        cum_yld_tmp[:, :-4][:, j])(mass_ni56))
    lc06_abs_010ni.append(np.array(yld_tmp))

lc06_abs_010ni = np.array(lc06_abs_010ni)

# LC06 initial composition
lc06_init_comp = np.ones(lc06_abs_010ni.shape) * 1e-30
for i in range(len(lc06_abs_010ni)):
    ztmp = lc06_metal[0]
    lc06_init_comp[i, 0] = (1. - 0.285 - ztmp) * np.sum(lc06_abs_010ni[i])  # H
    lc06_init_comp[i, 4] = 0.285 * np.sum(lc06_abs_010ni[i])  # He
    for j in range(len(ind_iso)):
        lc06_init_comp[i, ind_iso[j]] = ztmp * solar_ab[j] * \
                                        np.sum(lc06_abs_010ni[i])  # C->Mo

# LC06 net yields = absolute yields - initial composition of stellar model
lc06_net_010ni = lc06_abs_010ni - lc06_init_comp
lc06_net_010ni[np.where(lc06_net_010ni == 0.)] = 1e-30


# LC06 mass ejected
lc06_mej = np.sum(lc06_abs_010ni, axis=1)

# LC06 remnant mass
lc06_rem = lc06_mass - lc06_mej
# -----------------------------


# ----
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
# ----



# ---- Interpolated Yields -----

# ---- Minor grid points (mass bins spaced in ~1 Msun, but at the original 5
#---- metallicity values [Z = 1e-6, 1e-4, 1e-3, 6e-3, 2e-2])

# Interpolate across mass to generate yields at each mass bin of m_ave_high for
# 5 metallicity values: Z = 1e-6, 1e-4, 1e-3, 6e-3, & 2e-2 (solar)

# Linearly extrapolate the mass ejected and the remnant mass at fixed
# metallicity

# Do NOT linearly extrapolate net yields.  I tried extrapolate in mass as far
# as half the difference in mass between the last two grid points, but this
# destroyed more K39 than existed (in ~solar metallicity 8.5, 9.5, & 10.5 Msun
# stars).

lc_m_subsolar = np.array([13, 15, 20, 25, 30, 35])
lc_m_solar = np.array([11, 12, 13, 14, 15, 16, 17, 20, 25, 30, 35, 40, 60, 80,
                       120])

# dm_grid_low: mass difference between two lowest grid points

# m_floor: mass of lowest grid point minus the mass difference between two
# lowest grid points (i.e., dm_grid_low); the mass below which I use the yields
# extrapolated to this mass (i.e., the mass of the yield floor)

# yld_floor: yield of a given isotope at the mass of the yield floor (m_floor)


# indices of grid points at fixed metallicity
ind_cl04_net = [[6, 7, 8, 9, 10, 11],
                [12, 13, 14, 15, 16, 17],
                [18, 19, 20, 21, 22, 23],
                [24, 25, 26, 27, 28, 29]]

# output yields
lc_interp_mass = np.zeros((5, n_bins_high, n_species))

# loop over 4 intermediate metallicity values
for i in range(4):
    m_tmp = np.zeros((n_bins_high, n_species))
    for k in range(n_species):
        itmp = interpolate.InterpolatedUnivariateSpline(
            lc_m_subsolar, cl04_net_010ni[ind_cl04_net[i], k], k=1)
        m_tmp[:, k] = itmp(m_ave_high)
        m_tmp[np.where(m_ave_high < lc_m_subsolar[0]), k] = \
                                  itmp(lc_m_subsolar[0])
        m_tmp[np.where(m_ave_high > lc_m_subsolar[-1]), k] = \
                                  itmp(lc_m_subsolar[-1])
#       dm_grid_low = lc_m_subsolar[1] - lc_m_subsolar[0]
#       dm_grid_high = lc_m_subsolar[-1] - lc_m_subsolar[-2]
#        m_floor = lc_m_subsolar[0] - (dm_grid_low / 2.)
#        yld_floor = itmp(m_floor)
#        m_ceil = lc_m_subsolar[-1] + (dm_grid_high / 2.)
#        yld_ceil = itmp(m_ceil)
#        m_tmp[np.where(m_ave_high < m_floor), k] = yld_floor
#        m_tmp[np.where(m_ave_high > m_ceil), k] = yld_ceil
    lc_interp_mass[i] = m_tmp


# solar metallicity
m_tmp = np.zeros((n_bins_high, n_species))
for k in range(n_species):
    itmp = interpolate.InterpolatedUnivariateSpline(lc_m_solar,
                                                    lc06_net_010ni[:, k], k=1)
    m_tmp[:, k] = itmp(m_ave_high)
    m_tmp[np.where(m_ave_high < lc_m_solar[0]), k] = itmp(lc_m_solar[0])

#    dm_grid_low = lc_m_solar[1] - lc_m_solar[0]
#    m_floor = lc_m_solar[0] - (dm_grid_low / 2.)
#    yld_floor = itmp(m_floor)
#    m_tmp[np.where(m_ave_high < m_floor), k] = yld_floor

lc_interp_mass[4] = m_tmp


# mass ejected
lc_interp_mej = np.zeros((5, n_bins_high))
for i in range(4):
    itmp = interpolate.InterpolatedUnivariateSpline(
        lc_m_subsolar, cl04_mej[ind_cl04_net[i]], k=1)
    lc_interp_mej[i] = itmp(m_ave_high)

itmp = interpolate.InterpolatedUnivariateSpline(lc_m_solar, lc06_mej, k=1)
lc_interp_mej[4] = itmp(m_ave_high)


# remnant mass
lc_interp_rem = np.zeros((5, n_bins_high))
for i in range(5):
    lc_interp_rem[i] = m_ave_high - lc_interp_mej[i]


# -----------------------------


# ---- Interpolate across metallicity to generate yields at each mass bin of
#---- m_ave_high for N = n_metal_bin metallicity values

n_metal_bin = 1001
z_grid = np.array([1e-6, 1e-4, 1e-3, 6e-3, 2e-2])
logz_grid = np.log10(z_grid)

# evenly sample metallicity (in log Z) between grid points
logz_final = np.zeros(n_metal_bin)
dind = (n_metal_bin - 1) / (len(z_grid) - 1)
for i in range(len(z_grid) - 1):
    dlogz = (logz_grid[i+1] - logz_grid[i]) / \
            ((n_metal_bin - 1) / (len(z_grid) - 1))
    logz_final[i*dind:i*dind+dind+1] = np.arange(logz_grid[i],
                                                 logz_grid[i+1]+1e-9, dlogz)

# metallicity of final grid
z_final = 10.**logz_final

# output interpolated yields
lc_final = np.zeros((n_metal_bin, n_bins_high, n_species))
# at each mass, interpolate each element for each metallicity
for i in range(n_bins_high):
    for j in range(n_species):
        itmp = interpolate.InterpolatedUnivariateSpline(
            z_grid, lc_interp_mass[:, i, j], k=1)
        lc_final[:, i, j] = itmp(z_final)


# mass ejected
lc_final_mej = np.zeros((n_metal_bin, n_bins_high))
for i in range(n_bins_high):
    itmp = interpolate.InterpolatedUnivariateSpline(z_grid,
                                                    lc_interp_mej[:, i], k=1)
    lc_final_mej[:, i] = itmp(z_final)

# remnant mass
lc_final_rem = np.zeros((n_metal_bin, n_bins_high))
for i in range(n_metal_bin):
    lc_final_rem[i] = m_ave_high - lc_final_mej[i]


# pickle the interpolated yields array and the metallicity grid used
pickle_write(z_final, join(path_yldgen, 'interp_metallicity.pck'))
pickle_write(lc_final, join(path_lc06, 'interp_yields.pck'))
pickle_write(lc_final_mej, join(path_lc06, 'interp_meject.pck'))
pickle_write(lc_final_rem, join(path_lc06, 'interp_mremnant.pck'))


##############################
