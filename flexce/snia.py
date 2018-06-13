# @Author: Brett Andrews <andrews>
# @Date:   2018-06-06 12:06:40
# @Last modified by:   andrews
# @Last modified time: 2018-06-11 17:06:97

"""
FILE
    snia.py

DESCRIPTION
    Functions for computing the SNIa delay time distribution.
"""

import traceback

import numpy as np

from flexce.lifetimes import invert_lifetime


def set_snia_dtd(func='exponential', **kwargs):
    """Set SNIa delay time distribution.

    Args:
        func (str): functional form of DTD. Default is 'exponential'.
        kwargs (dict): keyword arguments to pass to individual DTD
            functions. Default is ``None``.

    Returns:
        dict: SNIa params
    """
    kwargs = kwargs if kwargs is not None else {}

    try:
        # TODO use map
        if func == 'exponential':
            params_snia = exponential(**kwargs)
        elif func == 'power_law':
            params_snia = power_law(**kwargs)
        elif func == 'prompt_delayed':
            params_snia = prompt_delayed(**kwargs)
        elif func == 'single_degenerate':
            params_snia = single_degenerate(**kwargs)

        params_snia['func'] = func

    except TypeError:
        print(traceback.print_exc())
        print(
            '\nValid keywords:\n'
            'exponential: dtime, mass, min_time, timescale, fraction\n'
            'power_law: time, min_time, nia_per_mstar, slope\n'
            'prompt_delayed: prob_prompt, prob_delay, min_time\n'
            'single_degenerate: time, dtime, alpha, breaks, num_int, mass_int, A,'
            'gam, eps=, normalize, nia_per_mstar\n'
        )

    return params_snia


def exponential(dt, mass, min_time=150., timescale=1500., fraction=0.078):
    """Implement exponential SNIa delay time distribution.

    If we adopt the SNIa prescription of Schoenrich & Binney (2009a)
    and a Salpeter IMF, 7.8% of the white dwarf mass formed form stars with
    initial mass between 3.2-8.0 Msun in a stellar population explodes as a
    SNIa (once we adjust to a mass range between 3.2-8.0 Msun instead of
    7.5% of the white dwarf mass that forms from stars of initial mass
    between 3.2-8.5 Msun).  For a Kroupa (2001) IMF, 5.5% of the white
    dwarf mass will explode as SNIa.

    Args:
        dt (float): Length of time step in Myr.
        mass (float): Mass of an individual SNIa.
        min_time (float): Minimum delay time for SNIa in Myr. Default
            is 150.
        timescale (float): Exponential decay timescale of delay time
            distribution in Myr. Default is 1500.
        snia_fraction (float): Fraction of white dwarf mass formed from
           stars with initial mass M=3.2-8.0 Msun that will explode in
           SNIa (see extended description). Default is 0.078.
    """
    dMwd = dt / timescale

    params_snia = {
        'min_time': min_time,
        'timescale': timescale,
        'fraction': fraction,
        'dMwd': dMwd,
        'mass': mass,
    }

    return params_snia


def power_law(time, min_time=40., nia_per_mstar=2.2e-3, slope=-1.):
    """Implement power-law SNIa delay time distribution.

    Args:
        min_snia_time (float): Minimum delay time for SNIa in Myr.
            Defaults to 150.
        nia_per_mstar (float): Number of SNIa per stellar mass formed
            that explode within 10 Gyr. Default is 2.2e-3.
        slope (float): Power law slope. Default is -1.
    """
    ind_min = (time >= min_time)
    ind10000 = (time <= 10000.)
    ria = np.zeros(len(time))
    ria[ind_min] = time[ind_min]**slope
    norm = nia_per_mstar / ria[ind10000].sum()
    ria = ria * norm

    params_snia = {
        'min_time': min_time,
        'ria': ria,
    }

    return params_snia


def prompt_delayed(prob_prompt=2.6e3, prob_delay=4.4e-8, min_time=40.):
    """Implement prompt plus delayed SNIa delay time distribution.

    Scannapieco & Bildstein (2005) prompt (B) + delayed (A) components
    to SNIa rate Equation 1\:
        N_Ia / (100 yr)^-1 = A [Mstar / 10^10 Msun] +
                             B [SFR / (10^10 Msun Gyr^-1)]

        A = 4.4e-2 (errors: +1.6e-2 -1.4e-2)
        B = 2.6 (errors: +/-1.1)

    In units useful for flexCE\:
        N_Ia per timestep = {4.4e-8 [Mstar / Msun] +
                             2.6e3 [SFR / (Msun yr^-1)]} *
                             (len of timestep in Myr)

    See also Mannucci et al. (2005).

    Args:
        prob_prompt (float): Default is 2.6e3.
        prob_delay (float): Coefficient connected to stellar mass of
            galaxy (see extended description). Defaults to 4.4e-8.
        min_snia_time (float): Minimum delay time for SNIa in Myr.
            Defaults to 150.


    """
    params_snia = {
        'min_time': min_time,
        'prob_prompt': prob_prompt,
        'prob_delay': prob_delay,
    }

    return params_snia


def single_degenerate(
    time,
    dt,
    alpha,
    breaks,
    num_int,
    mass_int,
    A=5e-4,
    gam=2.,
    eps=1.,
    normalize=False,
    nia_per_mstar=1.54e-3,
):
    """SNIa DTD for the single degenerate scenario.

    Solve for the SNIa rate (ria) according to Greggio (2005). The
    minimum primary mass is either
        (1) the mass of the secondary,
        (2) the mass required to form a carbon-oxygen white dwarf
            (2 Msun), or
        (3) the mass needed such that the WD mass plus the envelope of
        the secondary (accreted at an efficiency [eps]) equals the
        Chandrasekhar limit (1.4 Msun).

    Args:
        time (array): Time steps.
        dt (float): Length of time step.
        alpha (array): Power law slopes of the IMF.
        breaks (array): Mass of breaks in multi-slope IMF.
        num_int (array): Number of stars per mass bin per stellar mass
            formed.
        mass_int (array): Mass of stars per mass bin per stellar mass
            formed.
        A (float): Constant.
        gam (float): Power law index.
        eps (float): Accretion efficiency.
        normalize (bool): If ``True``, normalize the SNIa rate to match
            ``nia_per_mstar``, the number of SNIa within 10 Gyr per
            stellar mass formed. Default is ``False``.
        nia_per_mstar (float): Number of SNIa within 10 Gyr per stellar
            mass formed to normalize to if ``normalize`` is ``True``.
    """
    t2 = np.arange(29., time[-1] + 1., 1.)  # time in 1 Myr intervals
    m2 = invert_lifetime(t2)

    # calculate the envelope mass of the secondary
    m2ca = 0.3 * np.ones(len(t2))
    m2cb = 0.3 + 0.1 * (m2 - 2.)
    m2cc = 0.5 + 0.15 * (m2 - 4.)
    m2c = np.max((m2ca, m2cb, m2cc), axis=0)  # secondary core mass
    m2e = m2 - m2c  # secondary envelope mass

    mwdn = 1.4 - (eps * m2e)  # minimum WD mass
    m1na = 2. * np.ones(len(t2))
    m1nb = 2. + 10. * (mwdn - 0.6)

    # min prim. mass set by min CO WD mass
    m1n = np.max((m1na, m1nb), axis=0)  # min prim. mass
    m1low1 = invert_lifetime(t2)
    m1low = np.max((m1low1, m1n), axis=0)  # min primary mass

    m1up = 8.
    k_alpha = num_int.sum() / mass_int.sum()
    nm2 = np.zeros(len(m1low))

    for ii, aa in enumerate(alpha):
        if ii == 0:
            if len(breaks) > 0:
                ind = np.where(np.around(m1low, decimals=5) <= breaks[0])[0]
            else:
                ind = np.arange(len(m1low), dtype=int)

            ind_int = ind[:-1]

        elif ii != len(alpha) - 1:
            ind_int = np.where((m1low >= breaks[ii - 1]) &
                               (m1low <= breaks[ii]))[0][:-1]

        else:
            ind_int = np.where(m1low >= breaks[-1])[0]

        nm2[ind_int] = ((m2[ind_int]**-aa) *
                        ((m2[ind_int] / m1low[ind_int])**(aa + gam) -
                         (m2[ind_int] / m1up)**(aa + gam)))

    # from Greggio (2005): t**-1.44 approximates log(dm/dt) = log(m) -
    # log(t), which works for either the Padovani & Matteucci (1993) or the
    # Greggio (2005)/Girardi et al. (2000) stellar lifetimes
    dm2dt = 10.**4.28 * t2**1.44
    fia2 = nm2 / dm2dt
    fia = fia2 / fia2.sum()
    ria1 = k_alpha * A * fia

    ind_tbin = np.where(t2 % dt == 0.)[0]

    ria = np.zeros(len(time) - 1)
    ria[0] = ria1[:ind_tbin[0]].sum()
    for i in range(1, len(time) - 1):
        ria[i] = ria1[ind_tbin[i - 1]:ind_tbin[i]].sum()

    if normalize:
        ind10000 = np.where(time <= 10000.)
        ria = ria / ria[ind10000].sum() * nia_per_mstar

    params = {
        'ria': ria
    }

    return params


def snia_ev(params, tstep, dt, mstar, mstar_tot, sfr, Mwd_Ia):
    """Calculate the expected number of SNIa of a stellar population from
    a previous timestep.  The delay time distribution can be\:

    1. exponential
    2. empirical t^-1 power law
    3. empirical two component model with a prompt [~SFR] component and
       a delayed component [~Mstar]).
    4. theoretical DTD based on the single degenerate scenario

    Mannucci et al. (2005) find that the Rate SNIa / Rate SNII =
    0.35 +/- 0.08 in young stellar populations. Maoz et al. (2011) find
    that the time-integrated Rate SNII / Rate SNIa from a stellar
    population is about 5:1.

    Args:
        params (dict): SNIa DTD parameters.
        tstep (int): Time step.
        dt (float): Length of time step.
        mstar (float): Stellar mass formed in each time step.
        mstar_tot (float): Total stellar mass left.
        sfr (float): Star formation rate in each time step.
        Mwd_Ia (float): Mass in white dwarfs that will explode as SNIa
            (for exponetial DTD).

    Returns:
        float, float: Statistical expectation for the number of SNIa
            that should explode in a given time step, and mass in white
            dwarfs that will explode as SNIa (for exponetial DTD).
    """

    if params['func'] == 'exponential':
        ind_min_t = (tstep - np.ceil(params['min_time'] / dt).astype(int))
        if ind_min_t > 0:
            Nia_stat = np.sum(Mwd_Ia[:ind_min_t + 1] * params['dMwd'] / params['mass'])
            Mwd_Ia[:ind_min_t + 1] *= 1. - params['dMwd']
        else:
            Nia_stat = 0.

    elif params['func'] == 'power_law':
        Nia_stat = np.sum(params['ria'][:tstep] * np.sum(mstar[1:tstep + 1], axis=1)[::-1])

    elif params['func'] == 'prompt_delayed':
        ind = tstep - np.ceil(params['min_time'] / dt)
        Nia_prompt = sfr[ind] * params['prob_prompt'] if ind > 0 else 0.
        Nia_delay = mstar_tot * params['prob_delay']
        Nia_stat = (Nia_prompt + Nia_delay) * dt

    elif params['func'] == 'single_degenerate':
        Nia_stat = np.sum(params['ria'][:tstep] * np.sum(mstar[1:tstep + 1], axis=1)[::-1])

    return Nia_stat, Mwd_Ia
