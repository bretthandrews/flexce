# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:88
# @Last modified by:   andrews
# @Last modified time: 2018-06-07 20:06:78

"""
FILE
    chemevol.py

USAGE
    gal = ChemEvol(**params)

DESCRIPTION
    Main module for running the chemical evolution model.
"""

import os
from os.path import join
import time

import numpy as np

import flexce.imf
from flexce.imf import integrate_multi_power_law
import flexce.lifetimes
import flexce.utils


class ChemEvol:
    """Chemical evolution model.

    Args:
        params (dict): parameters to set up simulation run.
    """
    def __init__(self, params):
        # if any particular parameters are not specified,
        # then set them at the function call.
        self.params = params
        self.mass_bins = flexce.utils.set_mass_bins(params['mass_bins'])
        self.set_box(params['box'])
        self.set_yields(params['yields'], self.mass_bins)
        self.snia_dtd(params['snia_dtd'])  # TODO rename set_snia_dtd
        self.inflow_rx(params['inflows'])  # TODO rename set_inflow
        self.outflow_rx(params['outflows'])  # TODO rename set_outflow
        self.warmgasres_rx(params['warmgasres'])  # TODO rename set_warmgas_res
        self.star_formation(params['sf'])  # TODO rename set_star_formation

    def set_yields(self):
        # check for existing yields
        # if not, then calculate them.
        # path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
        # self.path_yldgen = join(path_flexce, 'data', 'yields', 'general', '')  # TODO fix

        pass

    def run(self):
        pass

    def set_box(
        self,
        radius=10.,
        time_tot=12000.,
        dt=30.,
        imf='kroupa',
        imf_alpha=None,
        imf_mass_breaks=None,
        sim_id=None
    ):
        """Initialize box.

        Args:
            radius (float): Radius of zone [kpc]. Only invoked if N_kslaw not
                equal to 1. Default is 10.
            time_tot (float): length of simulation [Myr]. Default is 12000.
            dt (float): length of time step [Myr]. Default is 30.
            imf (str): Stellar initial mass function. Default is 'kroupa'.
            imf_alpha (array): Power law slopes of user-defined stellar
                initial mass function. Must set ``imf`` to 'power_law'. Default
                is ``None``.
            imf_mass_breaks (array): Mass breaks between different power law
               slopes of user-defined stellar initial mass function. Must set
               ``imf`` to 'power_law'. Default is None.
            sim_id (str): Simulation ID number.
        """
        self.sim_id = 'box' + sim_id

        self.n_bins = len(self.mass_bins) - 1
        self.n_bins_high = len(np.where(self.mass_bins >= 8)[0]) - 1
        self.n_bins_low = len(np.where(self.mass_bins < 8)[0])

        self.radius = radius  # kpc
        self.area = self.radius**2. * np.pi * 1e6  # pc^2

        self.time_tot = time_tot
        self.dt = dt
        self.t = np.arange(0., self.time_tot + 1., self.dt)
        self.n_steps = int(self.time_tot / self.dt + 1.)

        self.imf = imf
        self.select_imf(imf, imf_alpha, imf_mass_breaks)

        self.tau_m = flexce.lifetimes.stellar_lifetimes(self.mass_ave)
        self.ind_ev, self.frac_ev = flexce.lifetimes.frac_evolve()

    def set_imf(self, imf, imf_alpha, imf_mass_breaks):
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
        if imf_alpha is not None:
            assert imf == 'power_law', 'Setting ``imf_alpha`` only sets IMF slope for a power law IMF'

        # TODO convert to Warning
        if imf_mass_breaks is not None:
            assert imf == 'power_law', 'Setting ``imf_mass_breaks`` only sets IMF mass breaks for a power law IMF'

        if imf == 'power_law':
            self.alpha = np.atleast_1d(np.array(imf_alpha))

            if imf_mass_breaks is None:
                imf_mass_breaks = []

            self.mass_breaks = np.atleast_1d(np.array(imf_mass_breaks))

        elif imf == 'salpeter':
            self.alpha = np.array([2.35])
            self.mass_breaks = np.array([])

        elif imf == 'kroupa':
            self.alpha = np.array([1.3, 2.3])
            self.mass_breaks = np.array([0.5])
            flexce.imf.kroupa()

        elif imf == 'bell':
            self.alpha = np.array([1., 2.35])
            self.mass_breaks = np.array([0.6])

        else:
            raise ValueError('Valid IMFs: "kroupa", "salpeter", "bell", or "power_law".')

        assert len(self.alpha) - len(self.mass_breaks) == 1, \
            ('Number of power law IMF slopes must be exactly ONE greater than the '
             'number of breaks in the power law.')

        alpha1 = (self.alpha - 1) * -1
        alpha2 = (self.alpha - 2) * -1

        norm = flexce.imf.normalize_imf()

        self.mass_int = integrate_multi_power_law(self.mass_bins, alpha2, self.mass_breaks, norm)
        self.num_int = integrate_multi_power_law(self.mass_bins, alpha1, self.mass_breaks, norm)
        self.mass_ave = self.mass_int / self.num_int

        aa = 1. / np.sum(self.mass_int)
        self.mass_frac = aa * self.mass_int

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
        self.metallicity[0] = (np.sum(self.mgas_iso[0, ind_metal]) / self.mgas_iso[0, 0])
        self.mfrac[0] = self.mgas_iso[0] / np.sum(self.mgas_iso[0])

        if np.sum(self.mwarmgas_iso[0]) > 0.:
            self.mwarmfrac[0] = (self.mwarmgas_iso[0] / np.sum(self.mwarmgas_iso[0]))

        snii_agb_yields = np.append(yields.agb_yields, yields.snii_yields, axis=1)

        for i in range(1, self.n_steps):
            if self.t[i] % 1000 < self.dt:
                print('time: {} Myr'.format(int(self.t[i])))

            self.metallicity[i] = (np.sum(self.mgas_iso[i - 1, ind_metal]) /
                                   self.mgas_iso[i - 1, 0])

            self.mfrac[i] = self.mgas_iso[i - 1] / np.sum(self.mgas_iso[i - 1])

            if np.sum(self.mwarmgas_iso[i - 1]) > 0.:
                self.mwarmfrac[i] = (self.mwarmgas_iso[i - 1] /
                                     np.sum(self.mwarmgas_iso[i - 1]))

            if sfh is None:
                self.sfr[i] = self.sf_law(np.sum(self.mgas_iso[i - 1]))

                # mimic gap in SF in the two infall model caused by a SF
                # threshold gas surface density
                if two_infall:
                    if ((self.t[i] > two_infall_args['t_sf_off'][0]) &
                            (self.t[i] < two_infall_args['t_sf_off'][1])):
                        self.sfr[i] = 0.
                    elif self.t[i] < 1000.:
                        self.sfr[i] = (self.sfr[i] * two_infall_args['sfe_thick'])
            else:
                self.sfr[i] = sfh[i]  # [=] Msun/yr
            self.dm_sfr[i] = self.sfr[i] * (self.dt * 1e6)
            # draw from IMF
            self.mstar_stat[i] = self.dm_sfr[i] * self.mass_frac
            self.Nstar_stat[i] = (self.dm_sfr[i] * self.mass_frac /
                                  self.mass_ave)
            self.random_num_state_Nstar.append(np.random.get_state())
            if set_state_Nstar is not None:
                np.random.set_state(set_state_Nstar[i - 1])
            self.Nstar[i] = flexce.utils.robust_random_poisson(self.Nstar_stat[i])
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
            for j in range(1, i + 1):
                # ind_ev is a list of indices of mass bins from a given birth
                # time-step that will evolve in the current time-step.
                ind = self.ind_ev[i - j]
                # abs_yld (171, 300) = net yield + (isotopic mass fraction at
                # birth) * (mass returned to ISM)
                abs_yld = (snii_agb_yields[ind_yld[j]] +
                           (self.mfrac[j] *
                            ((self.mass_ave -
                              yields.snii_agb_rem[ind_yld[j]]) *
                             np.ones((yields.n_sym, self.n_bins))).T))
                # number of stars to evolve
                N_ev = self.Nstar[j, ind] * self.frac_ev[i - j]
                snii_agb_tmp[ind] += (N_ev * abs_yld[ind].T).T
                if long_output:
                    self.snii_agb_net[i, ind] += (
                        N_ev * snii_agb_yields[ind_yld[j], ind].T).T
                # mass_returned (300); mass_returned_tot (300)
                # mass_returned = np.sum(N_ev * abs_yld[ind].T, axis=1)
                # mass_returned_tot += mass_returned
                mass_remnant_tot += np.sum(yields.snii_agb_rem[ind_yld[j], ind] * N_ev)
                self.Nstar_left[j, ind] -= N_ev
                self.mstar_left[j, ind] -= (self.mstar[j, ind] *
                                            self.frac_ev[i - j])
            self.snii[i] = np.sum(snii_agb_tmp[ind8:], axis=0)
            self.agb[i] = np.sum(snii_agb_tmp[:ind8], axis=0)
            if long_output:
                self.snii_agb[i] = snii_agb_tmp

            # SNIa

            if self.params['snia']['func'] == 'exponential':
                # mass of WDs that will be formed from the stellar population
                # that is born in the current timestep
                self.Mwd[i] = np.sum(self.Nstar[i, ind_ia] *
                                     yields.agb_rem[ind_yld[i], ind_ia])
                self.Mwd_Ia[i] = self.Mwd[i] * self.snia_fraction
                self.Mwd_Ia_init[i] = self.Mwd[i] * self.snia_fraction

            self.random_num_state_snia.append(np.random.get_state())
            if set_state_snia is not None:
                np.random.set_state(set_state_snia[i - 1][i - 1])

            self.NIa[i] = np.random.poisson(self.snia_ev(
                params=self.params['snia'],
                time=i,
                dt=self.dt,
                mstar=self.mstar,
                mstar_tot=self.mstar_left.sum(),
                sfr=self.sfr,
                Mwd_Ia=self.Mwd_Ia,
            ))
            self.snia[i] = yields.snia_yields * self.NIa[i]
            self.mremnant[i] = (mass_remnant_tot - self.NIa[i] * yields.snia_yields.sum())

            # gas flows
            self.sf[i] = np.sum(self.mstar[i]) * self.mfrac[i]
            self.outflow[i] = self.outflow_calc(i, self.sf[i], self.snii[i],
                                                self.agb[i], self.snia[i])

            if self.tcool > 0.:
                self.gas_cooling[i] = (self.mwarmgas_iso[i - 1] * self.dt /
                                       self.tcool)

            if self.inflow_func == 'constant_mgas':
                self.inflow_rate[i] = (np.sum(
                    self.sf[i] + self.outflow[i] - self.gas_cooling[i] -
                    self.fdirect * (self.snii[i] + self.agb[i] + self.snia[i])) / self.dt)

            self.inflow[i] = (self.inflow_composition(yields, i) *
                              self.inflow_rate[i] * self.dt)

            self.mgas_iso[i] = (self.mgas_iso[i - 1] +
                                (self.fdirect + self.feject) *
                                (self.snii[i] + self.agb[i] + self.snia[i]) +
                                self.gas_cooling[i] - self.sf[i] +
                                self.inflow[i] - self.outflow[i])

            self.mwarmgas_iso[i] = (
                self.mwarmgas_iso[i - 1] - self.gas_cooling[i] + self.fwarm *
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

        # TODO set self.param before evolving the box
        # self.param = dict(sim_id=self.sim_id, inflow=self.inflow_param,
        #                   outflow=self.outflow_param,
        #                   warmgas=self.warmgasres_param, sf=self.sf_param,
        #                   snia=self.snia_param, yields=yields.sources)
