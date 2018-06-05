# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:56
# @Last modified by:   andrews
# @Last modified time: 2018-06-05 12:06:56

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
        integral = (1. / exponent) * (bins[1:]**(exponent) - bins[:-1]**(exponent))

    return integral


def integrate_multi_power_law(bins, exponents, breaks, mass_bins, norm_factor):
    """Integrate over multi-slope power law distribution.

    Args:
        bins (array): stellar mass bins.
        exponents (array): slope of power law.
        breaks (array): stellar masses of breaks in multi-slope power
            law.
        norm_factor (array): normalization factor of integrals.

    Returns:
        array
    """
    for bb in breaks:
        assert np.around(bb, decimals=5) in np.around(mass_bins, decimals=5), \
            'Breaks in power law IMF must be located at the edge of a mass bin.'

    assert (len(exponents) > 1) and (len(exponents) - len(breaks) == 1), \
        ('Number of power law IMF slopes must be exactly ONE greater than the number '
         'of breaks in the power law.')

    integral = np.zeros(len(bins) - 1)

    for ii in range(len(exponents)):
        if ii == 0:
            if len(breaks) > 0:
                ind = np.where(np.around(bins, decimals=5) <= breaks[0])[0]
            else:
                ind = np.arange(len(bins), dtype=int)

        elif ii != len(exponents) - 1:
            ind = np.where((bins >= breaks[ii - 1]) & (bins <= breaks[ii]))[0]

        else:
            ind = np.where(bins >= breaks[-1])[0]

        ind_int = ind[:-1]
        integral[ind_int] = (integrate_power_law(exponents[ii], bins[ind]) *
                             norm_factor[ii])
    return integral
