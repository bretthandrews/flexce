"""Run chemical evolution model."""
from __future__ import print_function, division, absolute_import

import os
from os.path import join
import copy
import traceback
import time

import numpy as np
import pandas as pd

from flexce import utils


def integrate_power_law(exponent, bins=None):
    """Integrate a power law distribution.

    Args:
        exponent (float): power law exponent.
        bins (array): stellar mass bins. Defaults to None.
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


def lifetime_int(m):
    """Compute the lifetime of an intermediate mass star (M=0.6--6.6 Msun).

    Args:
        m (array): stellar mass.

    Returns:
        array: stellar lifetimes.
    """
    return 10.**((1.338 - np.sqrt(1.790 - 0.2232 * (7.764 - np.log10(m)))) /
                 0.1116 - 9.) * 1000.


def lifetime_high(m):
    """Compute the lifetime of a high mass star (M > 6.6 Msun)

    Args:
        m (array): stellar mass.

    Returns:
        array: stellar lifetimes.
    """
    return (1.2 * m**(-1.85) + 0.003) * 1000.


def invert_lifetime(t):
    """Compute stellar masses given lifetimes (valid for <50 Gyr).

    Args:
        t (array): lifetime in Myr.

    Returns:
        array: stellar masses.
    """
    m = np.zeros(len(t))
    ind_int = np.where((t >= 40.) & (t < 50000))
    ind_high = np.where((t > 1.) & (t < 40.))
    m[ind_int] = invert_lifetime_int(t[ind_int])
    m[ind_high] = invert_lifetime_high(t[ind_high])
    return m


def invert_lifetime_int(t):
    """Compute stellar masses given lifetimes (valid for 40 Myr-50 Gyr).

    Args:
        t (array): lifetime in Myr.

    Returns:
        array: stellar masses.
    """
    return 10.**(7.764 - (1.790 - (0.3336 - 0.1116 *
                                   np.log10(t / 1000.))**2.) / 0.2232)


def invert_lifetime_high(t):
    """Compute stellar masses given lifetimes (valid for <40 Myr).

    Args:
        t (array): lifetime in Myr.

    Returns:
        array: stellar masses.
    """
    return (((t / 1000.) - 0.003) / 1.2)**(-1./1.85)


def random_poisson(x):
    '''Draw a number from a Poisson distribution.  Used for determining
    the number of actual stars to form given an expected value.
    np.random.poisson cannot handle numbers larger than 2147483647
    (~2.14e9) because it uses the C long type.  For numbers larger than
    this value, round the expected (statistical) value to the nearest
    integer.  Since 2e9 is a large number, the Poisson fluctuations would
    have been relatively small anyway.

    This function is intended to replace this line of code:
    self.Nstar[i] = np.random.poisson(self.Nstar_stat[i])
    '''
    try:
        y = np.random.poisson(x)
    except ValueError:
        y = np.zeros(len(x), dtype=np.int64)
        for i, item in enumerate(x):
            try:
                y[i] = np.random.poisson(item)
            except ValueError:
                y[i] = np.round(item)
    return y


def evolve(yld, initialize_kws, snia_dtd_kws, inflows_kws, outflows_kws,
           warmgasres_kws, sf_kws):
    """Evolve the galaxy.

    Wrapper to initialize and run a simulation.

    Args:
        yld: Yields instance
        initialize_kws (dict): args to initialize ChemEvol instance.
        mass_bins_args (dict): args to define stellar mass bins.
        snia_dtd_kws (dict): args to set SNIa delay time distribution of
            ChemEvol instance.
        inflows_kws (dict): args to set inflow rate and composition of
            ChemEvol instance.
        outflows_kws (dict): args to set outflow rate and composition of
            ChemEvol instance.
        warmgasres_kws (dict): turn on warm ISM reservoir in ChemEvol
            instance.
        sf_kws (dict): args to set star formation rate in ChemEvol instance.

    Returns:
        ChemEvol instance
    """
    gal = ChemEvol(yld.mass_bins, **initialize_kws)
    gal.snia_dtd(**snia_dtd_kws)
    gal.inflow_rx(**inflows_kws)
    gal.outflow_rx(**outflows_kws)
    gal.warmgasres_rx(**warmgasres_kws)
    gal.star_formation(**sf_kws)
    gal.evolve_box(yields=yld)
    return gal


class ChemEvol:
    """Run chemical evolution model."""

    def __init__(self, mass_bins, radius=10., time_tot=12000., dt=30.,
                 imf='kroupa', imf_alpha=None, imf_mass_breaks=None,
                 sim_id=None):
        """Initialize model.

        Args:
            mass_bins (array): Stellar mass bins [Msun].
            radius (float): Radius of zone [kpc]. Only invoked if N_kslaw not
                equal to 1. Defaults to 10.
            time_tot (float): length of simulation [Myr]. Defaults to 12000.
            dt (float): length of time step [Myr]. Defaults to 30.
            imf (str): Stellar initial mass function. Defaults to 'kroupa'.
            imf_alpha (array): Power law slopes of user-defined stellar
                initial mass function. Must set imf to 'power_law'. Defaults
                to None.
            imf_mass_breaks (array): Mass breaks between different power law
               slopes of user-defined stellar initial mass function. Must set
               imf to 'power_law'. Defaults to None.
            sim_id (str): simulation ID number.
        """
        path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
        self.path_yldgen = join(path_flexce, 'data', 'yields', 'general', '')
        self.sim_id = 'box' + sim_id
        self.mass_bins = mass_bins
        self.n_bins = len(self.mass_bins) - 1
        self.n_bins_high = len(np.where(self.mass_bins >= 8)[0]) - 1
        self.n_bins_low = len(np.where(self.mass_bins < 8)[0])
        self.radius = radius  # kpc
        self.area = self.radius**2. * np.pi * 1e6  # pc^2
        self.timesteps(time_tot, dt)
        self.imf = imf
        self.select_imf(imf, imf_alpha, imf_mass_breaks)
        self.stellar_lifetimes()
        self.frac_evolve()

    def timesteps(self, time_tot, dt):
        """Set time steps.

        Args:
            time_tot (float): Length of simulation in Myr.
            dt (float): Size of time step in Myr.
        """
        self.time_tot = time_tot
        self.dt = dt
        self.t = np.arange(0., self.time_tot + 1., self.dt)
        self.n_steps = int(self.time_tot / self.dt + 1.)

    def select_imf(self, imf, imf_alpha, imf_mass_breaks):
        """Choose IMF or input user-defined power-law IMF.

        Args:
            imf (str): Stellar initial mass function to use.
            imf_alpha (array): Power law slopes of user-defined stellar
                initial mass function. Must set imf to 'power_law'.
            imf_mass_breaks (array): Mass breaks between different power law
               slopes of user-defined stellar initial mass function. Must set
               imf to 'power_law'.
        """
        if imf == 'power_law':
            self.powerlaw_imf(imf_alpha, imf_mass_breaks)
        elif imf == 'salpeter':
            self.salpeter()
        elif imf == 'kroupa':
            self.kroupa()
        elif imf == 'bell':
            self.bell()
        else:
            raise ValueError('Use valid IMF type: "kroupa", "salpeter",' +
                             '"bell", or "power_law".')
        self.mass_per_bin()

    def powerlaw_imf(self, alpha, mass_breaks):
        """Single or multiple slope power law IMF.

        Args:
            alpha (array): Power law slopes of user-defined stellar initial
                mass function.
            mass_breaks (array): Mass breaks between different power law
               slopes of user-defined stellar initial mass function.
        """
        self.alpha = np.atleast_1d(np.array(alpha))
        if mass_breaks is None:
            mass_breaks = []
        self.mass_breaks = np.atleast_1d(np.array(mass_breaks))
        self.imf_setup()

    def salpeter(self):
        """Set slope and mass breaks of Salpeter (1955) IMF."""
        self.alpha = np.array([2.35])
        self.mass_breaks = np.array([])
        self.imf_setup()

    def kroupa(self):
        """Set slope and mass breaks of Kroupa (2001) IMF."""
        self.alpha = np.array([1.3, 2.3])
        self.mass_breaks = np.array([0.5])
        self.imf_setup()

    def bell(self):
        """Set slope and mass breaks of Bell IMF.

        See Bell & de Jong (2001) (see Figure 4) and Bell et al. (2003) IMF.
        """
        self.alpha = np.array([1., 2.35])
        self.mass_breaks = np.array([0.6])
        self.imf_setup()

    def imf_setup(self):
        """Create reduced exponentials for IMF integration."""
        self.alpha1 = self.alpha - 1.
        self.alpha2 = self.alpha - 2.

    def mass_per_bin(self):
        """Calculate mass fraction that goes into each stellar mass bin."""
        # Normalize phi(m) for continuity between different IMF slopes
        try:
            norm_factor = self.normalize_imf()
        except IndexError:
            print(traceback.print_exc())
            raise ValueError('Number of power law IMF slopes must be ' +
                             'exactly ONE greater than the number of breaks ' +
                             'in the power law.')
        self.mass_int = integrate_multi_power_law(self.mass_bins,
                                                  self.alpha2 * -1.,
                                                  self.mass_breaks,
                                                  self.mass_bins,
                                                  norm_factor)
        self.num_int = integrate_multi_power_law(self.mass_bins,
                                                 self.alpha1 * -1.,
                                                 self.mass_breaks,
                                                 self.mass_bins,
                                                 norm_factor)
        self.mass_ave = self.mass_int / self.num_int
        a = 1. / np.sum(self.mass_int)
        self.mass_frac = a * self.mass_int
        # as a function of timestep
        self.mass_bins2 = invert_lifetime(self.t)
        self.mass_bins2[0] = self.mass_bins[-1]
        self.mass_int2 = integrate_multi_power_law(self.mass_bins2,
                                                   self.alpha2 * -1.,
                                                   self.mass_breaks,
                                                   self.mass_bins,
                                                   norm_factor * -1.)
        self.num_int2 = integrate_multi_power_law(self.mass_bins2,
                                                  self.alpha1 * -1.,
                                                  self.mass_breaks,
                                                  self.mass_bins,
                                                  norm_factor * -1.)
        self.mass_ave2 = self.mass_int2 / self.num_int2
        self.mass_frac2 = a * self.mass_int2

    def normalize_imf(self):
        """Normalize stellar initial mass function."""
        norm_factor = np.ones(len(self.alpha))
        if len(self.mass_breaks) > 0:
            for i in range(1, len(self.alpha)):
                norm_factor[i] = self.mass_breaks[i-1]**(-self.alpha[i-1]) / \
                                 self.mass_breaks[i-1]**(-self.alpha[i])
        return norm_factor

    def stellar_lifetimes(self):
        """Stellar lifetimes adopted from Padovani & Matteucci (1993).

        See Romano et al. (2005) for motivation.
        """
        self.tau_m = 160000. * np.ones(self.n_bins)  # [Myr]
        ind_mint = np.where((self.mass_ave > 0.6) & (self.mass_ave <= 6.6))[0]
        ind_mhigh = np.where(self.mass_ave > 6.6)[0]
        self.tau_m[ind_mint] = lifetime_int(self.mass_ave[ind_mint])
        self.tau_m[ind_mhigh] = lifetime_high(self.mass_ave[ind_mhigh])

    def frac_evolve(self):
        """Compute fraction of stars born in a given timestep will evolve.

        Figure out which mass bins will have at least some stars evolving in
        a given timestep (ind_ev) and what fraction of the stars in that mass
        bin will evolve (frac_ev).
        """
        self.ind_ev = []
        self.frac_ev = []
        # lowest mass star that would evolve in a timestep
        m_ev = invert_lifetime(self.t)
        m_ev[0] = self.mass_bins[-1]
        # integrate the IMF in each mass bin
        mass_int_tmp = np.zeros(self.n_bins)
        norm_factor = self.normalize_imf()
        for j in range(self.n_bins):
            mbin = self.mass_bins[j:j+2]
            mass_int_tmp[j] = integrate_multi_power_law(mbin,
                                                        self.alpha2 * -1,
                                                        self.mass_breaks,
                                                        self.mass_bins,
                                                        norm_factor)
        # figure out which mass bins will have at least some stars evolving in
        # a given timestep (ind_ev) and what fraction of the stars in that mass
        # bin will evolve (frac_ev)
        for i in range(self.n_steps - 1):
            indtmp = []
            fractmp = []
            for j in range(self.n_bins):
                # mass bin that spans the top end of the mass range that will
                # evolve in this timestep
                if (m_ev[i] >= self.mass_bins[j]) and \
                       (m_ev[i] < self.mass_bins[j+1]):
                    indtmp.append(j)
                    mlow_tmp = np.maximum(self.mass_bins[j], m_ev[i+1])
                    mbin_tmp = np.array([mlow_tmp, m_ev[i]])
                    mass_int2_tmp = integrate_multi_power_law(mbin_tmp,
                                                              self.alpha2 * -1,
                                                              self.mass_breaks,
                                                              self.mass_bins,
                                                              norm_factor)
                    fractmp.append(mass_int2_tmp / mass_int_tmp[j])
                # mass bins fully contained within the mass range that will
                # evolve in this timestep
                elif ((self.mass_bins[j] > m_ev[i+1]) and
                      (self.mass_bins[j] < m_ev[i])):
                    indtmp.append(j)
                    mbin_tmp = self.mass_bins[j:j+2]
                    mass_int2_tmp = integrate_multi_power_law(mbin_tmp,
                                                              self.alpha2 * -1,
                                                              self.mass_breaks,
                                                              self.mass_bins,
                                                              norm_factor)
                    fractmp.append(mass_int2_tmp / mass_int_tmp[j])
                # mass bin that spans bottom top end of the mass range that
                # will evolve in this timestep
                elif ((m_ev[i+1] > self.mass_bins[j]) and
                      (m_ev[i+1] < self.mass_bins[j+1])):
                    indtmp.append(j)
                    mbin_tmp = np.array([m_ev[i+1], self.mass_bins[j+1]])
                    mass_int2_tmp = integrate_multi_power_law(mbin_tmp,
                                                              self.alpha2 * -1,
                                                              self.mass_breaks,
                                                              self.mass_bins,
                                                              norm_factor)
                    fractmp.append(mass_int2_tmp / mass_int_tmp[j])
            indtmp = np.array(indtmp)
            self.ind_ev.append(indtmp)
            fractmp = np.array(fractmp)[:, 0]
            self.frac_ev.append(fractmp)
        self.frac_ev_tot = np.zeros(self.n_bins)
        for j in range(self.n_steps - 1):
            self.frac_ev_tot[self.ind_ev[j]] += self.frac_ev[j]

    def snia_dtd(self, func='exponential', kwargs=None):
        """Set SNIa delay time distribution.

        Args:
            func (str): functional form of DTD. Defaults to 'exponential'.
            kwargs (dict): keyword arguments to pass to individual DTD
                functions. Defaults to None.
        """
        kwargs = utils.none_to_empty_dict(kwargs)
        self.snia_param = dict(func=func, k=kwargs)
        self.snia_dtd_func = func
        try:
            if func == 'exponential':
                self.snia_dtd_exp(**kwargs)
            elif func == 'power_law':
                self.snia_dtd_powerlaw(**kwargs)
            elif func == 'prompt_delayed':
                self.snia_dtd_prompt_delayed(**kwargs)
            elif func == 'single_degenerate':
                self.snia_dtd_single_degenerate(**kwargs)
        except TypeError:
            print(traceback.print_exc())
            print('\nValid keywords:\n')
            print('exponential: timescale, min_snia_time, snia_fraction\n')
            print('power_law: min_snia_time, nia_per_mstar, slope\n')
            print('prompt_delayed: A, B,  min_snia_time\n')
            print('single_degenerate: no keywords\n')

    def snia_dtd_exp(self, min_snia_time=150., timescale=1500.,
                     snia_fraction=0.078):
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
                Defaults to 150.
            timescale (float): exponential decay timescale of delay time
                distribution in Myr. Defaults to 1500.
            snia_fraction (float): fraction of white dwarf mass formed from
               stars with initial mass M=3.2-8.0 Msun that will explode in
               SNIa (see extended description). Defaults to 0.078.
        """
        self.snia_fraction = snia_fraction
        self.min_snia_time = min_snia_time
        self.snia_timescale = timescale
        self.dMwd = self.dt / self.snia_timescale

    def snia_dtd_powerlaw(self, min_snia_time=40., nia_per_mstar=2.2e-3,
                          slope=-1.):
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

    def snia_dtd_prompt_delayed(self, A=4.4e-8, B=2.6e3, min_snia_time=40.):
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

    def snia_dtd_single_degenerate(self, A=5e-4, gam=2., eps=1.,
                                   normalize=False, nia_per_mstar=1.54e-3):
        '''SNIa DTD for the single degenerate scenario (SDS).

        Solve for the SNIa rate (ria) according to Greggio (2005).  The minimum
        primary mass is either (1) the mass of the secondary, (2) the mass
        required to form a carbon-oxygen white dwarf (2 Msun), or (3) the mass
        needed such that the WD mass plus the envelope of the secondary
        (accreted at an efficiency [eps]) equals the Chandrasekhar limit (1.4
        Msun).
        '''
        t2 = np.arange(29., self.t[-1]+1., 1.)  # time in 1 Myr intervals
        m2 = invert_lifetime(t2)
        # mass_int2_tmp = self.integrate_multi_power_law(
        #     m2, self.alpha2 * -1, self.mass_breaks, self.mass_bins,
        #     self.normalize_imf() * -1)
        # num_int2_tmp = self.integrate_multi_power_law(
        #     m2, self.alpha1 * -1, self.mass_breaks, self.mass_bins,
        #     self.normalize_imf() * -1)
        # mass_ave2 = mass_int2_tmp / num_int2_tmp
        # a = 1. / self.mass_int.sum()
        # a2 = 1. / np.sum(mass_int2_tmp)
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
        k_alpha = self.num_int.sum() / self.mass_int.sum()
        nm2 = np.zeros(len(m1low))
        for i in range(len(self.alpha)):
            if i == 0:
                if len(self.mass_breaks) > 0:
                    ind = np.where(np.around(m1low, decimals=5) <=
                                   self.mass_breaks[0])[0]
                else:
                    ind = np.arange(len(m1low), dtype=int)
                ind_int = ind[:-1]
            elif i != len(self.alpha) - 1:
                ind = np.where((m1low >= self.mass_breaks[i-1]) &
                               (m1low <= self.mass_breaks[i]))[0]
                ind_int = ind[:-1]
            else:
                ind = np.where(m1low >= self.mass_breaks[-1])[0]
                ind_int = ind
            nm2[ind_int] = ((m2[ind_int]**-self.alpha[i]) *
                            ((m2[ind_int]/m1low[ind_int])**(self.alpha[i]+gam) -
                             (m2[ind_int]/m1up)**(self.alpha[i] + gam)))
        # from Greggio (2005): t**-1.44 approximates log(dm/dt) = log(m) -
        # log(t), which works for either the Padovani & Matteucci (1993) or the
        # Greggio (2005)/Girardi et al. (2000) stellar lifetimes
        dm2dt = 10.**4.28 * t2**1.44
        fia2 = nm2 / dm2dt
        fia = fia2 / fia2.sum()
        ria1 = k_alpha * A * fia
        ind_tbin = np.where(t2 % self.dt == 0.)[0]
        self.ria = np.zeros(self.n_steps - 1)
        self.ria[0] = ria1[:ind_tbin[0]].sum()
        for i in range(1, self.n_steps - 1):
            self.ria[i] = ria1[ind_tbin[i-1]:ind_tbin[i]].sum()
        if normalize:
            ind10000 = np.where(self.t <= 10000.)
            self.ria = self.ria / self.ria[ind10000].sum() * nia_per_mstar

    def snia_ev(self, tstep, snia_mass, mstar_tot, sfr):
        '''Calculate the expected number of SNIa of a stellar population from
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

        snia_mass: mass of an individual SNIa

        min_snia_time: the minimum delay time from the birth of a stellar
        population
        '''
        if self.snia_dtd_func == 'exponential':
            ind_min_t = (tstep -
                         np.ceil(self.min_snia_time / self.dt).astype(int))
            if ind_min_t > 0:
                Nia_stat = np.sum(self.Mwd_Ia[:ind_min_t+1] * self.dMwd /
                                  snia_mass)
                self.Mwd_Ia[:ind_min_t+1] *= 1. - self.dMwd
            else:
                Nia_stat = 0.
        elif self.snia_dtd_func == 'power_law':
            Nia_stat = np.sum(self.ria[:tstep] *
                              np.sum(self.mstar[1:tstep+1], axis=1)[::-1])
        elif self.snia_dtd_func == 'prompt_delayed':
            ind = tstep - np.ceil(self.min_snia_time / self.dt)
            if ind < 0.:
                Nia_prompt = 0.
            else:
                Nia_prompt = sfr[ind] * self.prob_prompt
            Nia_stat = (Nia_prompt + (mstar_tot * self.prob_delay)) * self.dt
        elif self.snia_dtd_func == 'single_degenerate':
            Nia_stat = np.sum(self.ria[:tstep] *
                              np.sum(self.mstar[1:tstep+1], axis=1)[::-1])
        return Nia_stat

    def inflow_rx(self, func='double_exp', mgas_init=0., k=None,
                  inflow_rate=None, inflow_ab_pattern='bbns',
                  inflow_metallicity=1.):
        '''
        inflow_rate [=] Msun/Myr
        func: double_exp, exp, te-t, constant_mgas, or custom
        double_exp: M1, b1, M2, b2 (see Eq. 6 in Schoenrich & Binney 2009) with
        b1 & b2 in Myr
        exp: M1, b1 with b1 in Myr
        te-t: M1, b1 with b1 in Myr
        constant_mgas: inflow_rate will be dynamically defined in evolve_box
        custom: user-defined inflow rate

        inflow_comp: alpha enhanced
        '''
        self.inflow_param = dict(
            mgas_init=mgas_init, func=func, k=k, ab_pattern=inflow_ab_pattern,
            metallicity=inflow_metallicity)
        self.mgas_init = mgas_init
        self.inflow_func = func
        if func == 'double_exp':
            self.inflow_rate = ((k['M1']/k['b1']) * np.exp(-self.t/k['b1']) +
                                (k['M2']/k['b2']) * np.exp(-self.t/k['b2']))
        elif func == 'exp':
            self.inflow_rate = (k['M1']/k['b1']) * np.exp(-self.t/k['b1'])
        elif func == 'te-t':
            self.inflow_rate = ((k['M1']/k['b1']) * (self.t/k['b1']) *
                                np.exp(-self.t/k['b1']))
        elif func == 'constant_mgas':
            self.inflow_rate = np.zeros(self.n_steps)
        elif func == 'custom':
            self.inflow_rate = inflow_rate
        else:
            print('\nValid inflow functions: "double_exp", "exp", "te-t",' +
                  ' "constant_mgas", and "custom\n')
        self.inflow_ab_pattern = inflow_ab_pattern
        self.inflow_metallicity = inflow_metallicity

    def inflow_composition(self, yields, tstep):
        '''Compute the mass fraction of each element in the inflowing gas.
        "bbns": Big Bang Nucleosynthesis abundance pattern
        "alpha_enhanced": abundance pattern of a simulation before SNIa
        "scaled_solar": solar abundance pattern
        "recycled": abundance pattern of last timestep
        scaling factor is relative to solar (i.e., solar = 1)

        Set hydrogen mass fraction to 0.75 and helium mass fraction to 0.25 -
        the mass fraction of metals.  You need a hydrogen to helium mass
        fraction ratio of ~3 to avoid negative absolute yields of hydrogen.  (I
        had originally set things up to match the hydrogen/helium ratio of the
        ISM but that ran away to negative hydrogen masses).

        '''
        scaling_factor = self.inflow_metallicity  # relative to solar
        ind_h = np.where(yields.sym == 'H1')
        ind_he = np.where(yields.sym == 'He4')
        ind_metal = np.where(yields.sym_mass > 4.)
        if self.inflow_ab_pattern == 'bbns':
            inflow = yields.bbmf
            return inflow
        elif self.inflow_ab_pattern == 'alpha_enhanced':
            inftmp = pd.read_csv(self.path_yldgen + 'Z_0.1-Zsun_alpha_enhanced.txt',
                                 skiprows=6, header=None)
            inflow_init = np.array(inftmp).T
            scaling_factor *= 10.
        elif self.inflow_ab_pattern == 'scaled_solar':
            inflow_init = copy.deepcopy(yields.solar_mfrac)
        elif self.inflow_ab_pattern == 'recycled':
            ind = tstep - 1
            inflow_init = self.mgas_iso[ind] / self.mgas_iso[ind].sum()
            scaling_factor = (0.02 * scaling_factor /
                              inflow_init[ind_metal].sum())
        else:
            print('\nValid inflow compositions: "bbns", "alpha_enhanced",' +
                  ' "scaled_solar", and "recycled"\n')
        inflow = np.zeros(yields.n_sym)
        inflow[ind_metal] = inflow_init[ind_metal] * scaling_factor
        tmp = inflow.sum()
        # Set H & He mass fraction to 0.75 & 0.25 - Z, respectively.
        inflow[ind_h] = 0.75
        inflow[ind_he] = 0.25 - tmp
        return inflow

    def outflow_rx(self, outflow_source='ism', eta_outflow=1.,
                   variable_eta=None, feject=0.15):
        '''outflow_source = "ism" (ambient ISM is ejected in the wind; standard
        Mdot_wind = eta * SFR treatment) or "stellar_ejecta" (the yields from
        SNII, SNIa, and AGB stars makes up the wind; from Schoenrich & Binney
        2009). '''
        self.outflow_param = dict(outflow_source=outflow_source,
                                  eta_outflow=eta_outflow,
                                  variable_eta=variable_eta, feject=feject)
        self.outflow_source = outflow_source
        self.variable_eta = variable_eta
        if outflow_source == 'ism':
            if self.variable_eta is not None:
                self.eta_outflow = self.variable_eta
            else:
                self.eta_outflow = eta_outflow
            self.feject = 0.
        elif outflow_source == 'stellar_ejecta':
            self.feject = feject
            self.eta_outflow = 0.
        else:
            print('\nValid outflow sources: "ism" and "stellar_ejecta"\n')

    def outflow_calc(self, timestep, sfr, snii, agb, snia):
        if self.outflow_source == 'ism':
            if self.variable_eta is not None:
                return self.eta_outflow[timestep] * sfr
            else:
                return self.eta_outflow * sfr
        elif self.outflow_source == 'stellar_ejecta':
            return self.feject * (snii + agb + snia)
        else:
            print('\nValid outflow sources: "ism" and "stellar_ejecta"\n')

    def warmgasres_rx(self, mwarmgas_init=0., fdirect=0.01, tcool=1200.,
                      warmgas=True):
        '''Parameters that control gas flow into and out of the warm gas
        reservoir.
        Schoenrich & Binney (2009) fiducial values:
        mwarmgas_init = 5e8 Msun
        fdirect = 0.01 (feject=0.15 for R < 3.5 kpc and 0.04 for R > 3.5 kpc)
        tcool = 1200 Myr
        fwarm = 1 - fdirect - feject
        '''
        self.warmgas_on = warmgas
        if warmgas:
            filename = 'warmgas_abundance_pattern.txt'
            tmp = pd.read_csv(self.path_yldgen + filename,
                              delim_whitespace=True, skiprows=10,
                              names=['el', 'ab'])
            self.warmgas_ab_pattern = np.array(tmp['ab'])
            self.mwarmgas_init = mwarmgas_init
            self.fdirect = fdirect
            self.tcool = tcool
            self.fwarm = 1. - fdirect
            if self.outflow_source == 'stellar_ejecta':
                self.fwarm -= self.feject
        else:
            self.warmgas_ab_pattern = 0.
            self.mwarmgas_init = 0.
            self.fdirect = 1. - self.feject
            self.tcool = 0.
            self.fwarm = 0.
        self.warmgasres_param = dict(warmgas_on=warmgas,
                                     mwarmgas_init=mwarmgas_init,
                                     fdirect=fdirect, fwarm=self.fwarm,
                                     tcool=tcool)

    def star_formation(self, nu_kslaw=2.5e-10, N_kslaw=1.4):
        '''From Kennicutt (1998): Sigma_SFR = 2.5e-4 * Sigma_gas^1.4, where
        Sigma_SFR [=] Msun yr^-1 kpc^-2 and Sigma_gas [=] Msun yr^-1.  nu =
        2.5-10 if Sigma_SFR is converted into units of Msun yr^-1 pc^-2.  '''
        self.nu_kslaw = nu_kslaw
        self.N_kslaw = N_kslaw
        self.sf_param = dict(nu=nu_kslaw, N=N_kslaw)

    def sf_law(self, mgas):
        '''Calculate the SFR [Msun/yr] given the gas mass and the two free
        parameters that determine the SF law: nu_kslaw [Gyr^-1] and N_kslaw
        (~1.0--2.0).

        From Kennicutt (1998): Sigma_SFR = 2.5e-4 * Sigma_gas^1.4, where
        Sigma_SFR [=] Msun yr^-1 kpc^-2 and Sigma_gas [=] Msun yr^-1.  nu =
        2.5-10 if Sigma_SFR is converted into units of Msun yr^-1 pc^-2.
        '''
        return (self.nu_kslaw) * self.area * (mgas / self.area)**self.N_kslaw

    def initialize_arrays(self, yld, long_output):
        self.agb = np.zeros((self.n_steps, yld.n_sym))
        self.gas_cooling = np.zeros((self.n_steps, yld.n_sym))
        self.dm_sfr = np.zeros(self.n_steps)
        self.inflow = np.zeros((self.n_steps, yld.n_sym))
        self.metallicity = np.zeros(self.n_steps)
        self.mfrac = np.zeros((self.n_steps, yld.n_sym))
        self.mgas_iso = np.zeros((self.n_steps, yld.n_sym))
        self.mremnant = np.zeros(self.n_steps)
        self.mstar = np.zeros((self.n_steps, self.n_bins))
        self.mstar_left = np.zeros((self.n_steps, self.n_bins))
        self.mstar_stat = np.zeros((self.n_steps, self.n_bins))
        self.mwarmfrac = np.zeros((self.n_steps, yld.n_sym))
        self.mwarmgas_iso = np.zeros((self.n_steps, yld.n_sym))
        self.Mwd = np.zeros(self.n_steps)
        self.Mwd_Ia = np.zeros(self.n_steps)
        self.Mwd_Ia_init = np.zeros(self.n_steps)
        self.NIa = np.zeros(self.n_steps, dtype=int)
        self.Nstar = np.zeros((self.n_steps, self.n_bins), dtype=np.int64)
        self.Nstar_left = np.zeros((self.n_steps, self.n_bins))
        self.Nstar_stat = np.zeros((self.n_steps, self.n_bins))
        self.outflow = np.zeros((self.n_steps, yld.n_sym))
        self.sf = np.zeros((self.n_steps, yld.n_sym))
        self.sfr = np.zeros(self.n_steps)
        self.snia = np.zeros((self.n_steps, yld.n_sym))
        self.snii = np.zeros((self.n_steps, yld.n_sym))
        self.random_num_state_Nstar = []
        self.random_num_state_snia = []
        if long_output:
            self.snii_agb = np.zeros((self.n_steps, self.n_bins, yld.n_sym))
            self.snii_agb_net = np.zeros((self.n_steps, self.n_bins,
                                          yld.n_sym))
            self.snii_agb_rec = np.zeros((self.n_steps, self.n_bins,
                                          yld.n_sym))

    def evolve_box(self, yields, sfh=None, two_infall=False,
                   two_infall_args=None, set_state_Nstar=None,
                   set_state_snia=None, long_output=False):
        '''Start with a box of primordial gas.  Calculate the metallicity and
        mass fraction of each isotope (mfrac).  SFR is calculated assuming a
        constant gas depletion timescale.  Statistical expectation for total
        mass of stars formed is calculated (dm_sfr).  The statistical
        expectation of mass and number of stars formed per bin are calculated
        (mstar_stat, Nstar_stat).  The expected number of stars per bin is fed
        into a Poisson random number generator, which returns the
        stochastically sampled number of stars formed per bin (Nstar_high) and
        the mass per bin (mstar).

        High mass stars explode as SNe II, and low mass stars evolve as AGB
        stars.  Both sources return gas to the ISM according to
        metallicity-dependent yields (dmgas_iso).  The amount of mass locked up
        in remnants (neutron stars and white dwarfs) is calculated (mremnant).
        The mass of each isotope in the gas reservoir is calculated by taking
        the gas mass from the previous time-step, adding the gas mass of each
        isotope returned by SNe II and AGB stars, and then subtracting the gas
        mass of each isotope that went into stars (mgas_iso).

        The model is destroying more of certain isotopes (e.g., N15) than
        exist, so it forces isotopes with negative quantities to be a very
        small, positive quantity.

        sfh: manually set the star formation history of the zone

        two infall model: t_sf_off = time span ([0]=start, [1]=end) when no SF
        occurs (to mimic the gas surface density dropping below a threshold
        density for SF), sfe_thick = thick disk SFE / thin disk SFE
        '''
        self.initialize_arrays(yields, long_output)
        ind_yld = np.zeros(self.n_steps, dtype=int)
        ind8 = np.where(self.mass_bins == 8.)[0][0]
        ind_ia = np.where((self.mass_ave >= 3.2) & (self.mass_ave <= 8.))[0]
        start = time.time()
        # initial conditions
        self.mgas_iso[0] = yields.bbmf * self.mgas_init
        if self.warmgas_on:
            self.mwarmgas_iso[0] = (self.warmgas_ab_pattern *
                                    self.mwarmgas_init / 4.)
        ind_metal = np.where(yields.sym_mass > 4.)
        self.metallicity[0] = (np.sum(self.mgas_iso[0, ind_metal]) /
                               self.mgas_iso[0, 0])
        self.mfrac[0] = self.mgas_iso[0] / np.sum(self.mgas_iso[0])
        if np.sum(self.mwarmgas_iso[0]) > 0.:
            self.mwarmfrac[0] = (self.mwarmgas_iso[0] /
                                 np.sum(self.mwarmgas_iso[0]))
        snii_agb_yields = np.append(yields.agb_yields, yields.snii_yields,
                                    axis=1)
        for i in range(1, self.n_steps):
            if self.t[i] % 1000 < self.dt:
                print('time: {} Myr'.format(int(self.t[i])))
            self.metallicity[i] = (np.sum(self.mgas_iso[i-1, ind_metal]) /
                                   self.mgas_iso[i-1, 0])
            self.mfrac[i] = self.mgas_iso[i-1] / np.sum(self.mgas_iso[i-1])
            if np.sum(self.mwarmgas_iso[i-1]) > 0.:
                self.mwarmfrac[i] = (self.mwarmgas_iso[i-1] /
                                     np.sum(self.mwarmgas_iso[i-1]))
            if sfh is None:
                self.sfr[i] = self.sf_law(np.sum(self.mgas_iso[i-1]))
                # mimic gap in SF in the two infall model caused by a SF
                # threshold gas surface density
                if two_infall:
                    if (self.t[i] > two_infall_args['t_sf_off'][0]) & \
                           (self.t[i] < two_infall_args['t_sf_off'][1]):
                        self.sfr[i] = 0.
                    elif self.t[i] < 1000.:
                        self.sfr[i] = (self.sfr[i] *
                                       two_infall_args['sfe_thick'])
            else:
                self.sfr[i] = sfh[i]  # [=] Msun/yr
            self.dm_sfr[i] = self.sfr[i] * (self.dt * 1e6)
            # draw from IMF
            self.mstar_stat[i] = self.dm_sfr[i] * self.mass_frac
            self.Nstar_stat[i] = (self.dm_sfr[i] * self.mass_frac /
                                  self.mass_ave)
            self.random_num_state_Nstar.append(np.random.get_state())
            if set_state_Nstar is not None:
                np.random.set_state(set_state_Nstar[i-1])
            self.Nstar[i] = random_poisson(self.Nstar_stat[i])
            self.mstar[i] = self.Nstar[i] * self.mass_ave
            self.Nstar_left[i] = self.Nstar[i]
            self.mstar_left[i] = self.mstar[i]
            # SNII and AGB yields
            if self.metallicity[i] < yields.snii_agb_z[0]:
                ind_yld[i] = 0
            elif self.metallicity[i] > yields.snii_agb_z[-1]:
                ind_yld[i] = -1
            else:
                ind_yld[i] = np.where(self.metallicity[i] <
                                      yields.snii_agb_z)[0][0]
            # ind_yld[i] = -1 # uncomment for solar metallicity yields only
            # Evolve stars from previous timesteps
            snii_agb_tmp = np.zeros((self.n_bins, yields.n_sym))
            # mass_returned_tot = 0.
            mass_remnant_tot = 0.
            for j in range(1, i+1):
                # ind_ev is a list of indices of mass bins from a given birth
                # time-step that will evolve in the current time-step.
                ind = self.ind_ev[i-j]
                # abs_yld (171, 300) = net yield + (isotopic mass fraction at
                # birth) * (mass returned to ISM)
                abs_yld = (snii_agb_yields[ind_yld[j]] +
                           (self.mfrac[j] *
                            ((self.mass_ave -
                              yields.snii_agb_rem[ind_yld[j]]) *
                             np.ones((yields.n_sym, self.n_bins))).T))
                # number of stars to evolve
                N_ev = self.Nstar[j, ind] * self.frac_ev[i-j]
                snii_agb_tmp[ind] += (N_ev * abs_yld[ind].T).T
                if long_output:
                    self.snii_agb_net[i, ind] += (
                        N_ev * snii_agb_yields[ind_yld[j], ind].T).T
                # mass_returned (300); mass_returned_tot (300)
                # mass_returned = np.sum(N_ev * abs_yld[ind].T, axis=1)
                # mass_returned_tot += mass_returned
                mass_remnant_tot += np.sum(yields.snii_agb_rem[ind_yld[j], ind]
                                           * N_ev)
                self.Nstar_left[j, ind] -= N_ev
                self.mstar_left[j, ind] -= (self.mstar[j, ind] *
                                            self.frac_ev[i-j])
            self.snii[i] = np.sum(snii_agb_tmp[ind8:], axis=0)
            self.agb[i] = np.sum(snii_agb_tmp[:ind8], axis=0)
            if long_output:
                self.snii_agb[i] = snii_agb_tmp
            # SNIa
            # mass of WDs that will be formed from the stellar population that
            # is born in the current timestep
            if self.snia_dtd_func == 'exponential':
                self.Mwd[i] = np.sum(self.Nstar[i, ind_ia] *
                                     yields.agb_rem[ind_yld[i], ind_ia])
                self.Mwd_Ia[i] = self.Mwd[i] * self.snia_fraction
                self.Mwd_Ia_init[i] = self.Mwd[i] * self.snia_fraction
            self.random_num_state_snia.append(np.random.get_state())
            if set_state_snia is not None:
                np.random.set_state(set_state_snia[i-1][i-1])
            cnt_ia = np.random.poisson(
                self.snia_ev(i, yields.snia_yields.sum(),
                             self.mstar_left.sum(), self.sfr))
            self.NIa[i] = cnt_ia
            self.snia[i] = yields.snia_yields * self.NIa[i]
            self.mremnant[i] = (mass_remnant_tot - self.NIa[i] *
                                yields.snia_yields.sum())
            # gas flows
            self.sf[i] = np.sum(self.mstar[i]) * self.mfrac[i]
            self.outflow[i] = self.outflow_calc(i, self.sf[i], self.snii[i],
                                                self.agb[i], self.snia[i])
            if self.tcool > 0.:
                self.gas_cooling[i] = (self.mwarmgas_iso[i-1] * self.dt /
                                       self.tcool)
            if self.inflow_func == 'constant_mgas':
                self.inflow_rate[i] = (np.sum(
                    self.sf[i] + self.outflow[i] - self.gas_cooling[i] -
                    self.fdirect * (self.snii[i] + self.agb[i] + self.snia[i]))
                                       / self.dt)
            self.inflow[i] = (self.inflow_composition(yields, i) *
                              self.inflow_rate[i] * self.dt)
            self.mgas_iso[i] = (self.mgas_iso[i-1] +
                                (self.fdirect + self.feject) *
                                (self.snii[i] + self.agb[i] + self.snia[i]) +
                                self.gas_cooling[i] - self.sf[i] +
                                self.inflow[i] - self.outflow[i])
            self.mwarmgas_iso[i] = (
                self.mwarmgas_iso[i-1] - self.gas_cooling[i] + self.fwarm *
                (self.snii[i] + self.agb[i] + self.snia[i]))
            if (i < 4) & self.warmgas_on:
                self.mwarmgas_iso[i] += (self.warmgas_ab_pattern *
                                         self.mwarmgas_init / 4.)
        self.Nstar_left = self.Nstar_left.astype(int)
        self.mstar_left[np.where(np.abs(self.mstar_left) < -1e-8)] = 0.
        if long_output:
            self.snii_agb_rec = self.snii_agb - self.snii_agb_net
        print('Time elapsed:', time.time() - start)
        self.outflow_rate = np.sum(self.outflow, axis=1) / self.dt
        self.check_mass_conservation(yields)
        self.snii_snia_rate()
        # Set all negative masses equal to a small positive number.
        ind_neg = np.where(self.mgas_iso < 0.)
        self.mgas_iso[ind_neg] = 1e-30
        self.survivors = np.sum(self.Nstar_left[:, 1:], axis=1)
        self.survivors = self.survivors.round().astype(int)
        self.param = dict(sim_id=self.sim_id, inflow=self.inflow_param,
                          outflow=self.outflow_param,
                          warmgas=self.warmgasres_param, sf=self.sf_param,
                          snia=self.snia_param, yields=yields.sources)

    def check_mass_conservation(self, yields):
        '''Check mass conservation.

        inflow[0].sum() is fractionally larger than inflow_rate.sum() by 7.5e-5
        (due to sum(bbmf) = 1. + 7.5e-5)'''
        m_in = (self.inflow.sum() + self.mgas_iso[0].sum() +
                self.mwarmgas_iso[0].sum() * 4.)
        m_out = self.outflow.sum()
        m_gas = self.mgas_iso[-1].sum() + self.mwarmgas_iso[-1].sum()
        m_star = (self.mstar.sum() - self.snii.sum() - self.agb.sum() -
                  self.NIa.sum() * yields.snia_yields.sum())
        print('mass_in - mass_out:', m_in - m_out - m_gas - m_star)

    def snii_snia_rate(self):
        '''Instantaneous Rate(SNII) / Rate(SNIa) ~ 3:1 (Li et al.).
        Actually I think that it should be more like 5:1 (REF?)'''
        ind_snii = np.where(self.mass_ave > 8.)[0]
        self.NII = np.sum(self.Nstar[:, ind_snii], axis=1)
        rate_snia_snii = np.zeros(len(self.NII))
        ind = np.where((self.NII > 0) & (self.NIa > 0))[0]
        rate_snia_snii[ind] = self.NIa[ind].astype(float) / self.NII[ind]
        print('Rate of SNII to SNIa in last 100 timesteps:',
              1. / np.mean(rate_snia_snii[-100:]))
