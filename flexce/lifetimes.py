# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 12:06:85
# @Last modified by:   andrews
# @Last modified time: 2018-06-05 12:06:49

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


