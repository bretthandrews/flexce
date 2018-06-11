# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:56
# @Last modified by:   andrews
# @Last modified time: 2018-06-11 14:06:27

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


def integrate_multi_power_law(bins, exponents, breaks, norm_factor):
    """Integrate over multi-slope power law distribution.

    Args:
        bins (array): Stellar mass bins to integrate over.
        exponents (array): Power law slopes.
        breaks (array): Stellar masses of breaks in multi-slope power
            law.
        norm_factor (array): Normalization factor of integrals.

    Returns:
        array: IMF sampled over the bin sizes
    """
    for bb in breaks:
        if (bb >= np.min(bins)) and (bb <= np.max(bins)):
            assert np.around(bb, decimals=5) in np.around(bins, decimals=5), \
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


def normalize_imf(alpha, breaks):
    """Normalize stellar initial mass function.

    Args:
        alpha (array): IMF slopes.
        breaks (array): Power law break masses.

    Returns:
        array: Normalization factor for each power law segment.
    """
    norm_factor = np.ones(len(alpha))
    norm_factor[1:] = breaks**(-alpha[:1]) / breaks**(-alpha[1:])
    return norm_factor


def set_imf(imf, alpha, mass_breaks):
    """Choose IMF or input user-defined power-law IMF.

    Available IMFs:
        'salpeter': Salpeter (1955).
        'kroupa': Kroupa (2001).
        'bell': Bell & de Jong (2001) Fig 4 and Bell et al. (2003).
        'power_law': User-defined power law slopes and breaks.

    Args:
        imf (str): Stellar initial mass function to use.
        imf_alpha (array): Power law slopes of user-defined stellar initial
            mass function. Must set ``imf`` to 'power_law'.
        imf_mass_breaks (array): Mass breaks between different power law
           slopes of user-defined stellar initial mass function. Must set
           ``imf`` to 'power_law'.
    """
    # TODO convert to Warning
    if alpha is not None:
        assert imf == 'power_law', 'Setting ``imf_alpha`` only sets IMF slope for a power law IMF'

    # TODO convert to Warning
    if mass_breaks is not None:
        assert imf == 'power_law', 'Setting ``imf_mass_breaks`` only sets IMF mass breaks for a power law IMF'

    if imf == 'power_law':
        alpha = np.atleast_1d(np.array(alpha))

        if mass_breaks is None:
            mass_breaks = []

        mass_breaks = np.atleast_1d(np.array(mass_breaks))

    elif imf == 'salpeter':
        alpha = np.array([2.35])
        mass_breaks = np.array([])

    elif imf == 'kroupa':
        alpha = np.array([1.3, 2.3])
        mass_breaks = np.array([0.5])

    elif imf == 'bell':
        alpha = np.array([1., 2.35])
        mass_breaks = np.array([0.6])

    else:
        raise ValueError('Valid IMFs: "kroupa", "salpeter", "bell", or "power_law".')

    assert len(alpha) - len(mass_breaks) == 1, \
        ('Number of power law IMF slopes must be exactly ONE greater than the '
         'number of breaks in the power law.')

    params = {
        'alpha': alpha,
        'mass_breaks': mass_breaks,
    }

    return params
