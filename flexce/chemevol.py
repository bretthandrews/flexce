# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:88
# @Last modified by:   andrews
# @Last modified time: 2018-06-13 10:06:79

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
import flexce.inflows
import flexce.lifetimes
import flexce.outflows
import flexce.star_formation
import flexce.snia
import flexce.utils
import flexce.warm_gas_reservoir as warmgas


class ChemEvol:
    """Chemical evolution model.

    Args:
        params (dict): Parameters of simulation. Default is ``None``.
        ylds: Yields instance. Default is ``None``.
    """
    def __init__(self, params=None, ylds=None):

        params = params if params is not None else {}
        props = ['mass_bins', 'box', 'yields', 'snia_dtd', 'inflows', 'outflows',
                 'warmgasres', 'sf']
        for prop in props:
            if prop not in params.keys():
                params[prop] = {}

        self.params = params
        self.mass_bins = flexce.utils.set_mass_bins(**params['mass_bins'])
        self.set_box(**params['box'])

        if ylds is None:
            path_data = join(os.path.dirname(__file__), 'data')
            # TODO try to load existing yields. If they don't exist, then calculate them.
            ylds = flexce.utils.load_yields(path_data, self.mass_bins, params['yields'])

        self.params['snia_dtd'] = flexce.snia.set_snia_dtd(
            time=self.time,
            dtime=self.dtime,
            mass=ylds.snia_yields.sum(),
            alpha=self.params['imf']['alpha'],
            breaks=self.params['imf']['mass_breaks'],
            mass_int=self.mass_int,
            num_int=self.num_int,
            **params['snia_dtd'])

        self.params['inflows'], self.inflow_rate = flexce.inflows.set_inflows(
            time=self.time,
            **params['inflows'],
        )

        self.params['outflows'], self.eta, self.feject = flexce.outflows.set_outflows(
            **params['outflows'],
        )

        self.params['warmgas'], self.warmgas_ab_pattern = warmgas.set_warm_gas_reservoir(
            feject=self.feject,
            outflow_source=self.params['outflows']['source'],
            **params['warmgasres'],
        )

        self.params['sf'] = flexce.star_formation.set_sflaw(**params['sf'])

        self.evolve_box(
            ylds=ylds,
            sfh=None,                # TODO set as params['sf']['sfh']
            two_infall=False,        # TODO set as params['inflows']['two_infall']
            two_infall_kwargs=None,  # TODO set as params['inflows']['two_infall_kwargs']
            set_state_Nstar=None,    # TODO set as params['box']['set_state_Nstar']
            set_state_snia=None,     # TODO set as params['box']['set_state_snia']
            long_output=False,       # TODO set as params['box']['long_output']
        )

    def set_box(
        self,
        radius=10.,
        time_tot=12000.,
        dtime=30.,
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
        self.params['box']['sim_id'] = sim_id

        self.n_bins = len(self.mass_bins) - 1
        self.n_bins_high = len(np.where(self.mass_bins >= 8)[0]) - 1
        self.n_bins_low = len(np.where(self.mass_bins < 8)[0])

        self.params['box']['radius'] = radius  # kpc
        self.params['box']['area'] = radius**2. * np.pi * 1e6  # pc^2

        self.params['box']['time_tot'] = time_tot
        self.dtime = dtime
        self.time = np.arange(0., time_tot + 1., self.dtime)
        self.n_steps = int(time_tot / self.dtime + 1.)

        # IMF
        self.params['imf'] = flexce.imf.set_imf(imf, imf_alpha, imf_mass_breaks)

        alpha1 = 1 - self.params['imf']['alpha']
        alpha2 = 2 - self.params['imf']['alpha']

        norm = flexce.imf.normalize_imf(self.params['imf']['alpha'],
                                        self.params['imf']['mass_breaks'])

        self.mass_int = integrate_multi_power_law(self.mass_bins, alpha2,
                                                  self.params['imf']['mass_breaks'], norm)
        self.num_int = integrate_multi_power_law(self.mass_bins, alpha1,
                                                 self.params['imf']['mass_breaks'], norm)
        self.mass_ave = self.mass_int / self.num_int
        self.mass_frac = self.mass_int / np.sum(self.mass_int)

        self.tau_m = flexce.lifetimes.set_lifetimes(self.mass_ave)
        self.ind_ev, self.frac_ev = flexce.lifetimes.frac_evolve(
            self.time,
            self.mass_bins,
            self.params['imf']['alpha'],
            self.params['imf']['mass_breaks'],
        )

    def check_mass_conservation(self, yields):
        """Checks for mass conservation.

        ``inflow[0].sum()`` is fractionally larger than
        ``inflow_rate.sum()`` by 7.5e-5 (due to
        ``sum(bbmf)`` = 1. + 7.5e-5)

        Args:
            yields: ``Yields`` instance.
        """
        m_in = (self.inflow.sum() + self.mgas_iso[0].sum() + self.mwarmgas_iso[0].sum() * 4.)
        m_out = self.outflow.sum()
        m_gas = self.mgas_iso[-1].sum() + self.mwarmgas_iso[-1].sum()
        m_star = (self.mstar.sum() - self.snii.sum() - self.agb.sum() -
                  self.NIa.sum() * yields.snia_yields.sum())
        dm = m_in - m_out - m_gas - m_star
        assert dm < 1e-4, f'Mass not conserved!\nmass_in - mass_out: {dm}'

    def snii_snia_rate(self):
        """Prints ratio of SNII to SNIa rates in last 100 time steps.

        Instantaneous Rate(SNII) / Rate(SNIa) ~ 3:1 (Li et al.).
        Actually I think that it should be more like 5:1 (REF?)
        """
        ind_snii = (self.mass_ave > 8.)
        self.NII = np.sum(self.Nstar[:, ind_snii], axis=1)

        rate_snia_snii = np.zeros(len(self.NII))

        ind = ((self.NII > 0) & (self.NIa > 0))
        rate_snia_snii[ind] = self.NIa[ind].astype(float) / self.NII[ind]

        print('Rate of SNII to SNIa in last 100 timesteps:',
              1. / np.mean(rate_snia_snii[-100:]))

    def initialize_arrays(self, n_sym, long_output):
        """Initialize arrays for simulation.

        Args:
            n_sym (int): Number of isotopes in ``Yields`` instance.
            long_output (bool): If ``True``, record SNII and AGB
                yields.
        """
        n_steps = len(self.time)
        self.agb = np.zeros((n_steps, n_sym))
        self.gas_cooling = np.zeros((n_steps, n_sym))
        self.dm_sfr = np.zeros(n_steps)
        self.inflow = np.zeros((n_steps, n_sym))
        self.metallicity = np.zeros(n_steps)
        self.mfrac = np.zeros((n_steps, n_sym))
        self.mgas_iso = np.zeros((n_steps, n_sym))
        self.mremnant = np.zeros(n_steps)
        self.mstar = np.zeros((n_steps, self.n_bins))
        self.mstar_left = np.zeros((n_steps, self.n_bins))
        self.mstar_stat = np.zeros((n_steps, self.n_bins))
        self.mwarmfrac = np.zeros((n_steps, n_sym))
        self.mwarmgas_iso = np.zeros((n_steps, n_sym))
        self.Mwd = np.zeros(n_steps)
        self.Mwd_Ia = np.zeros(n_steps)
        self.Mwd_Ia_init = np.zeros(n_steps)
        self.NIa = np.zeros(n_steps, dtype=int)
        self.Nstar = np.zeros((n_steps, self.n_bins), dtype=np.int64)
        self.Nstar_left = np.zeros((n_steps, self.n_bins))
        self.Nstar_stat = np.zeros((n_steps, self.n_bins))
        self.outflow = np.zeros((n_steps, n_sym))
        self.sf = np.zeros((n_steps, n_sym))
        self.sfr = np.zeros(n_steps)
        self.snia = np.zeros((n_steps, n_sym))
        self.snii = np.zeros((n_steps, n_sym))
        self.random_num_state_Nstar = []
        self.random_num_state_snia = []
        if long_output:
            self.snii_agb = np.zeros((n_steps, self.n_bins, yld.n_sym))
            self.snii_agb_net = np.zeros((n_steps, self.n_bins, yld.n_sym))
            self.snii_agb_rec = np.zeros((n_steps, self.n_bins, yld.n_sym))

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
            if self.time[i] % 1000 < self.dtime:
                print('time: {} Myr'.format(int(self.time[i])))

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
                    if ((self.time[i] > two_infall_args['t_sf_off'][0]) &
                            (self.time[i] < two_infall_args['t_sf_off'][1])):
                        self.sfr[i] = 0.
                    elif self.time[i] < 1000.:
                        self.sfr[i] = (self.sfr[i] * two_infall_args['sfe_thick'])
            else:
                self.sfr[i] = sfh[i]  # [=] Msun/yr
            self.dm_sfr[i] = self.sfr[i] * (self.dtime * 1e6)
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
                dt=self.dtime,
                mstar=self.mstar,
                mstar_tot=self.mstar_left.sum(),
                sfr=self.sfr,
                Mwd_Ia=self.Mwd_Ia,
            ))
            self.snia[i] = yields.snia_yields * self.NIa[i]
            self.mremnant[i] = (mass_remnant_tot - self.NIa[i] * yields.snia_yields.sum())

            # gas flows
            self.sf[i] = np.sum(self.mstar[i]) * self.mfrac[i]
            self.outflow[i] = self.outflow_calc(
                params=self.params['outflows'],
                eta=self.eta_outflow,
                feject=self.feject,
                timestep=i,
                sfr=self.sf[i],
                stellar_ejecta=self.snii[i] + self.agb[i] + self.snia[i],
            )

            if self.tcool > 0.:
                self.gas_cooling[i] = (self.mwarmgas_iso[i - 1] * self.dtime / self.tcool)

            if self.params['inflows']['func'] == 'constant_mgas':
                self.inflow_rate[i] = (
                    np.sum(
                        self.sf[i] + self.outflow[i] - self.gas_cooling[i] -
                        self.fdirect * (self.snii[i] + self.agb[i] + self.snia[i])
                    ) / self.dtime
                )

            self.inflow[i] = (
                self.inflow_composition(self.params['inflows'], yields, self.mgas_iso[i - 1]) *
                self.inflow_rate[i] * self.dtime
            )

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

        self.outflow_rate = np.sum(self.outflow, axis=1) / self.dtime
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
