# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 12:06:85
# @Last modified by:   andrews
# @Last modified time: 2018-06-06 12:06:61

"""
FILE
    lifetimes.py

DESCRIPTION
    Utility functions for populating a stellar lifetimes.
"""

import numpy as np

import flexce.imf
from flexce.imf import integrate_multi_power_law


def lifetime_int(mass):
    """Compute the lifetime of intermediate mass star (M=0.6-6.6 Msun).

    Args:
        mass (array): Stellar mass.

    Returns:
        array: Stellar lifetimes.
    """
    return 10.**((1.338 - np.sqrt(1.790 - 0.2232 * (7.764 - np.log10(mass)))) /
                 0.1116 - 9.) * 1000.


def lifetime_high(mass):
    """Compute the lifetime of a high mass star (M > 6.6 Msun).

    Args:
        mass (array): Stellar mass.

    Returns:
        array: Stellar lifetimes.
    """
    return (1.2 * mass**(-1.85) + 0.003) * 1000.


def invert_lifetime(time):
    """Compute stellar masses given lifetimes (valid for <50 Gyr).

    Args:
        time (array): Lifetime in Myr.

    Returns:
        array: Stellar masses.
    """
    ind_int = (time >= 40.) & (time < 50000)
    ind_high = (time > 1.) & (time < 40.)

    mass = np.zeros(len(time))
    mass[ind_int] = invert_lifetime_int(time[ind_int])
    mass[ind_high] = invert_lifetime_high(time[ind_high])

    return mass


def invert_lifetime_int(time):
    """Compute stellar masses given lifetimes (for 40 Myr-50 Gyr).

    Args:
        time (array): Lifetime in Myr.

    Returns:
        array: Stellar masses.
    """
    return 10.**(7.764 - (1.790 - (0.3336 - 0.1116 * np.log10(time / 1000.))**2.) / 0.2232)


def invert_lifetime_high(time):
    """Compute stellar masses given lifetimes (valid for <40 Myr).

    Args:
        time (array): Lifetime in Myr.

    Returns:
        array: Stellar masses.
    """
    return (((time / 1000.) - 0.003) / 1.2)**(-1. / 1.85)


def set_lifetimes(mass_ave):
    """Stellar lifetimes adopted from Padovani & Matteucci (1993).

    See Romano et al. (2005) for motivation.
    """
    tau = 160000. * np.ones(len(mass_ave))  # [Myr]

    ind_mint = (mass_ave > 0.6) & (mass_ave <= 6.6)
    ind_mhigh = (mass_ave > 6.6)

    tau[ind_mint] = lifetime_int(mass_ave[ind_mint])
    tau[ind_mhigh] = lifetime_high(mass_ave[ind_mhigh])

    return tau


def frac_evolve(time, bins, alpha, breaks):
    """Compute fraction of stars born in a given timestep will evolve.

    Figure out which mass bins will have at least some stars evolving in
    a given timestep (ind_ev) and what fraction of the stars in that mass
    bin will evolve (frac_ev).

    Args:
        time (array): Time.
        mass_bins (array): Mass bins.
        alpha (array): IMF power law slopes.
        mass_breaks (arry): Break masses in multi-slope IMF.
    Returns:
        array, array: Index and fraction of mass bin of stars that will
           evolve in a given time step.
    """
    # lowest mass star that would evolve in a timestep
    m_ev = invert_lifetime(time)
    m_ev[0] = bins[-1]

    # integrate the IMF in each mass bin
    n_bins = len(bins) - 1
    alpha2 = 2 - alpha
    norm_factor = flexce.imf.normalize_imf(alpha, breaks)

    mass_int_tmp = np.zeros(n_bins)
    for jj in range(n_bins):
        mbin = bins[jj:jj + 2]
        mass_int_tmp[jj] = integrate_multi_power_law(mbin, alpha2, breaks, norm_factor)

    # figure out which mass bins will have at least some stars evolving in
    # a given timestep (ind_ev) and what fraction of the stars in that mass
    # bin will evolve (frac_ev)
    ind_ev = []
    frac_ev = []

    for ii in range(len(time) - 1):
        indtmp = []
        fractmp = []

        for jj in range(n_bins):
            # mass bin that spans the top end of the mass range that will
            # evolve in this timestep
            if (m_ev[ii] < bins[jj + 1]) and (m_ev[ii] >= bins[jj]):
                mlow_tmp = np.maximum(bins[jj], m_ev[ii + 1])
                mbin_tmp = np.array([mlow_tmp, m_ev[ii]])

            # mass bins fully contained within the mass range that will
            # evolve in this timestep
            elif (bins[jj] < m_ev[ii]) and (bins[jj] > m_ev[ii + 1]):
                mbin_tmp = bins[jj:jj + 2]

            # mass bin that spans bottom top end of the mass range that
            # will evolve in this timestep
            elif (m_ev[ii + 1] < bins[jj + 1]) and (m_ev[ii + 1] > bins[jj]):
                mbin_tmp = np.array([m_ev[ii + 1], bins[jj + 1]])

            else:
                continue

            indtmp.append(jj)
            mass_int2_tmp = integrate_multi_power_law(mbin_tmp, alpha2, breaks, norm_factor)
            fractmp.append(mass_int2_tmp / mass_int_tmp[jj])

        indtmp = np.array(indtmp)
        ind_ev.append(indtmp)

        fractmp = np.array(fractmp)[:, 0]
        frac_ev.append(fractmp)

    return ind_ev, frac_ev
