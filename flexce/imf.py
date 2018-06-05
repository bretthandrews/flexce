# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:56
# @Last modified by:   andrews
# @Last modified time: 2018-06-05 11:06:42

"""
FILE
    imf.py

DESCRIPTION
    Utility functions for populating a stellar initial mass function.
"""

import numpy as np


def integrate_power_law(exponent, bins=None):
    """Integrate a power law distribution.

    Args:
        exponent (float): power law exponent.
        bins (array): stellar mass bins. Default is ``None``.
    """
    if exponent == 0.:
        integral = bins[1:] - bins[:-1]
    elif exponent == -1.:
        integral = np.log(bins[1:]) - np.log(bins[:-1])
    else:
        integral = (1. / exponent) * (bins[1:]**(exponent) -
                                      bins[:-1]**(exponent))
    return integral


def integrate_multi_power_law(bins, exponents, breaks, mass_bins,
                              norm_factor):
    """Integrate over each section of multi-slope power law distribution.

    Args:
        bins (array): stellar mass bins.
        exponents (array): slope of power law.
        breaks (array): stellar masses of breaks in multi-slope power law.
        norm_factor (array): normalization factor of integrals.
    """
    if check_multi_slope_compatibility(exponents, breaks, mass_bins):
        integral = np.zeros(len(bins) - 1)
        for i in range(len(exponents)):
            if i == 0:
                if len(breaks) > 0:
                    ind = np.where(np.around(bins, decimals=5) <= breaks[0])[0]
                else:
                    ind = np.arange(len(bins), dtype=int)
            elif i != len(exponents) - 1:
                ind = np.where((bins >= breaks[i-1]) &
                               (bins <= breaks[i]))[0]
            else:
                ind = np.where(bins >= breaks[-1])[0]
            ind_int = ind[:-1]
            integral[ind_int] = (integrate_power_law(exponents[i],
                                                     bins[ind]) *
                                 norm_factor[i])
        return integral


def check_multi_slope_compatibility(exponents, breaks, mass_bins):
    """Check if the parameters of multi-slope power law are allowed.

    Args:
        exponents (array): Power law exponents.
        breaks (array): Stellar masses of breaks in multi-slope power law.
        mass_bins (array): Stellar mass bins.

    Returns:
        bool
    """
    for b in breaks:
        if np.around(b, decimals=5) not in \
               np.around(mass_bins, decimals=5):
            raise ValueError('Breaks in power law IMF must be located ' +
                             'at the edge of a mass bin.')
    if (len(exponents) > 1) and (len(exponents) - len(breaks) != 1):
        raise ValueError('Number of power law IMF slopes must be exactly ' +
                         'ONE greater than the number of breaks in the ' +
                         'power law.')
    return True
