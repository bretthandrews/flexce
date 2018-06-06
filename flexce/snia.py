# @Author: Brett Andrews <andrews>
# @Date:   2018-06-06 12:06:40
# @Last modified by:   andrews
# @Last modified time: 2018-06-06 12:06:46

"""
FILE
    snia.py

DESCRIPTION
    Functions for computing the SNIa delay time distribution.
"""

import traceback

import numpy as np

import flexce.utils


def snia_dtd(func='exponential', kwargs=None):
    """Set SNIa delay time distribution.

    Args:
        func (str): functional form of DTD. Default is 'exponential'.
        kwargs (dict): keyword arguments to pass to individual DTD
            functions. Default is ``None``.

    Returns:
        dict: SNIa params
    """
    kwargs = flexce.utils.none_to_empty_dict(kwargs)
    snia_param = {'func': func, 'k': kwargs}
    try:
        if func == 'exponential':
            dtd_exp(**kwargs)
        elif func == 'power_law':
            dtd_powerlaw(**kwargs)
        elif func == 'prompt_delayed':
            dtd_prompt_delayed(**kwargs)
        elif func == 'single_degenerate':
            snia_dtd_single_degenerate(**kwargs)
    except TypeError:
        print(traceback.print_exc())
        print(
            '\nValid keywords:\n'
            'exponential: timescale, min_snia_time, snia_fraction\n'
            'power_law: min_snia_time, nia_per_mstar, slope\n'
            'prompt_delayed: A, B, min_snia_time\n'
            'single_degenerate: no keywords\n'
        )

    return snia_param

def dtd_exp(min_snia_time=150., timescale=1500., snia_fraction=0.078):
    """Implement exponential SNIa delay time distribution.

    If we adopt the SNIa prescription of Schoenrich & Binney (2009a)
    and a Salpeter IMF, 7.8% of the white dwarf mass formed form stars with
    initial mass between 3.2-8.0 Msun in a stellar population explodes as a
    SNIa (once we adjust to a mass range between 3.2-8.0 Msun instead of
    7.5% of the white dwarf mass that forms from stars of initial mass
    between 3.2-8.5 Msun).  For a Kroupa (2001) IMF, 5.5% of the white
    dwarf mass will explode as SNIa.

    Args:
        min_snia_time (float): Minimum delay time for SNIa in Myr.
            Default is 150.
        timescale (float): Exponential decay timescale of delay time
            distribution in Myr. Default is 1500.
        snia_fraction (float): Fraction of white dwarf mass formed from
           stars with initial mass M=3.2-8.0 Msun that will explode in
           SNIa (see extended description). Default is 0.078.
    """
    self.snia_fraction = snia_fraction
    self.min_snia_time = min_snia_time
    self.snia_timescale = timescale
    self.dMwd = self.dt / self.snia_timescale


def dtd_powerlaw(self, min_snia_time=40., nia_per_mstar=2.2e-3, slope=-1.):
    """Implement power-law SNIa delay time distribution.

    Args:
        min_snia_time (float): Minimum delay time for SNIa in Myr.
            Defaults to 150.
        nia_per_mstar (float): number of SNIa per stellar mass formed
            that explode within 10 Gyr. Defaults to 2.2e-3.
        slope (float): power law slope. Defaults to -1.
    """
    self.min_snia_time = min_snia_time
    ind_min = np.where(self.t >= min_snia_time)
    ind10000 = np.where(self.t <= 10000.)
    ria = np.zeros(len(self.t))
    ria[ind_min] = self.t[ind_min]**slope
    norm = nia_per_mstar / ria[ind10000].sum()
    self.ria = ria * norm


def dtd_prompt_delayed(self, A=4.4e-8, B=2.6e3, min_snia_time=40.):
    """Implement prompt plus delayed SNIa delay time distribution.

    Args:
        A (float): coefficient connected to stellar mass of galaxy
            (see extended description). Defaults to 4.4e-8.
        B (float): Defaults to 2.6e3.
        min_snia_time (float): Minimum delay time for SNIa in Myr.
            Defaults to 150.

    Scannapieco & Bildstein (2005) prompt + delayed components to SNIa
    rate Equation 1\:

    N_Ia / (100 yr)^-1 = A [Mstar / 10^10 Msun] +
    B [SFR / (10^10 Msun Gyr^-1)]

    A = 4.4e-2 (errors: +1.6e-2 -1.4e-2)

    B = 2.6 (errors: +/-1.1)

    In units useful for flexCE\:
    N_Ia per timestep = {4.4e-8 [Mstar / Msun] +
    2.6e3 [SFR / (Msun yr^-1)]} * (len of timestep in Myr)
    see also Mannucci et al. (2005)
    """
    self.prob_delay = A
    self.prob_prompt = B
    self.min_snia_time = min_snia_time
    return
