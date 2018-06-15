# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:88
# @Last modified by:   andrews
# @Last modified time: 2018-06-15 12:06:83

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
        state (dict): Random number state for drawing from the IMF and
            exploding SNIa. Default is ``None``.
    """
    def __init__(self, params=None, ylds=None, state=None):

        params = params if params is not None else {}
        props = ['mass_bins', 'box', 'yields', 'snia_dtd', 'inflows', 'outflows',
                 'warmgasres', 'sf']
        for prop in props:
            if prop not in params.keys():
                params[prop] = {}

        self.params = params

        if state is None:
            state = {'Nstar': [], 'snia': []}

        self.state = state

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
            **params['snia_dtd'],
        )

        self.params['inflows'], self.inflow_rate = flexce.inflows.set_inflows(
            time=self.time,
            **params['inflows'],
        )

        self.params['outflows'] = flexce.outflows.set_outflows(**params['outflows'])

        self.params['warmgas'], self.warmgas_ab_pattern = warmgas.set_warm_gas_reservoir(
            feject=self.params['outflows']['feject'],
            outflow_source=self.params['outflows']['source'],
            **params['warmgas'],
        )

        self.params['sf'] = flexce.star_formation.set_sflaw(**params['sf'])

        self.evolve_box(ylds=ylds)

    def set_box(
        self,
        radius=10.,
        time_tot=12000.,
        dtime=30.,
        imf='kroupa',
        imf_alpha=None,
        imf_mass_breaks=None,
        sim_id=None,
        long_output=False,
    ):
        """Initialize box.

        Args:
            radius (float): Radius of zone [kpc]. Only invoked if
                N_kslaw is not equal to 1. Default is 10.
            time_tot (float): length of simulation [Myr]. Default is
                12000.
            dt (float): length of time step [Myr]. Default is 30.
            imf (str): Stellar initial mass function. Default is
                'kroupa'.
            imf_alpha (array): Power law slopes of user-defined stellar
                initial mass function. Must set ``imf`` to 'power_law'.
                Default is ``None``.
            imf_mass_breaks (array): Mass breaks between different
                power law slopes of user-defined stellar initial mass
                function. Must set ``imf`` to 'power_law'. Default is
                ``None``.
            sim_id (str): Simulation ID number. Default is ``None``.
            long_output (bool): If ``True``, record total, net, and
                recycled SNII and AGB yields from each previous time
                step. Default is ``False``.
        """
        self.params['box']['sim_id'] = sim_id
        self.params['box']['long_output'] = long_output

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

    def initialize_arrays(self, n_sym):
        """Initialize arrays for simulation.

        Args:
            n_sym (int): Number of isotopes in ``Yields`` instance.

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
        if self.params['box']['long_output']:
            self.snii_agb = np.zeros((n_steps, self.n_bins, n_sym))
            self.snii_agb_net = np.zeros((n_steps, self.n_bins, n_sym))
            self.snii_agb_rec = np.zeros((n_steps, self.n_bins, n_sym))

    def evolve_box(self, ylds):
        """Run the chemical evolution model.

        Start with a box of primordial gas. Calculate the metallicity
        and mass fraction of each isotope (``mfrac``). SFR is
        calculated assuming a constant gas depletion timescale.
        Statistical expectation for total mass of stars formed is
        calculated (``dm_sfr``). The statistical expectation of mass
        and number of stars formed per bin are calculated
        (``mstar_stat`` and ``Nstar_stat``, respectively).  The
        expected number of stars per bin is fed into a Poisson random
        number generator, which returns the stochastically sampled
        number of stars formed per bin (``Nstar_high``) and the mass
        per bin (``mstar``).

        High mass stars explode as SNe II, and low mass stars evolve as
        AGB stars. Both sources return gas to the ISM according to
        metallicity-dependent yields (``dmgas_iso``).  The amount of
        mass locked up in remnants (neutron stars and white dwarfs) is
        calculated (``mremnant``). The mass of each isotope in the gas
        reservoir is calculated by taking the gas mass from the
        previous time-step, adding the gas mass of each isotope
        returned by SNe II and AGB stars, and then subtracting the gas
        mass of each isotope that went into stars (``mgas_iso``).

        The model is destroying more of certain isotopes (e.g., N15)
        than exist, so it forces isotopes with negative quantities to
        be a very small, positive quantity.

        Args:
            ylds: ``Yields`` instance.

        """
        self.initialize_arrays(ylds.n_sym)

        ind_yld = np.zeros(self.n_steps, dtype=int)
        ind8 = np.where(self.mass_bins == 8.)[0][0]
        ind_ia = np.where((self.mass_ave >= 3.2) & (self.mass_ave <= 8.))[0]

        preset_state_Nstar = bool(self.state['Nstar'])
        preset_state_snia = bool(self.state['snia'])

        start = time.time()

        # initial conditions
        self.mgas_iso[0] = ylds.bbmf * self.params['inflows']['mgas_init']

        if self.params['warmgas']['warmgas_on']:
            self.mwarmgas_iso[0] = (self.warmgas_ab_pattern *
                                    self.params['warmgas']['mwarmgas_init'] / 4.)

        ind_metal = (ylds.sym_mass > 4.)
        self.metallicity[0] = np.sum(self.mgas_iso[0, ind_metal]) / self.mgas_iso[0, 0]
        self.mfrac[0] = self.mgas_iso[0] / np.sum(self.mgas_iso[0])

        if np.sum(self.mwarmgas_iso[0]) > 0.:
            self.mwarmfrac[0] = self.mwarmgas_iso[0] / np.sum(self.mwarmgas_iso[0])

        snii_agb_yields = np.append(ylds.agb_yields, ylds.snii_yields, axis=1)

        for ii in range(1, self.n_steps):
            if self.time[ii] % 1000 < self.dtime:
                print('time: {} Myr'.format(int(self.time[ii])))

            self.metallicity[ii] = (np.sum(self.mgas_iso[ii - 1, ind_metal]) /
                                    self.mgas_iso[ii - 1, 0])

            self.mfrac[ii] = self.mgas_iso[ii - 1] / np.sum(self.mgas_iso[ii - 1])

            if np.sum(self.mwarmgas_iso[ii - 1]) > 0.:
                self.mwarmfrac[ii] = self.mwarmgas_iso[ii - 1] / np.sum(self.mwarmgas_iso[ii - 1])

            if self.params['sf']['sfh'] is None:
                self.sfr[ii] = flexce.star_formation.sf_law(
                    mgas=np.sum(self.mgas_iso[ii - 1]),
                    params=self.params
                )

                # mimic gap in SF in the two infall model caused by a SF
                # threshold gas surface density
                if self.params['inflows']['func'] == 'two_infall':
                    end_thick, start_thin = self.params['inflows']['coeff']['t_sf_off']

                    if (self.time[ii] > end_thick) and (self.time[ii] < start_thin):
                        self.sfr[ii] = 0.

                    elif self.time[ii] <= end_thick:
                        self.sfr[ii] *= self.params['inflows']['coeff']['sfe_thick']

            else:
                self.sfr[ii] = self.params['sf']['sfh'][ii]  # [=] Msun/yr

            self.dm_sfr[ii] = self.sfr[ii] * (self.dtime * 1e6)

            # draw from IMF
            self.mstar_stat[ii] = self.dm_sfr[ii] * self.mass_frac
            self.Nstar_stat[ii] = self.dm_sfr[ii] * self.mass_frac / self.mass_ave

            if preset_state_Nstar:
                np.random.set_state(self.state['Nstar'][ii - 1])
            else:
                self.state['Nstar'].append(np.random.get_state())

            self.Nstar[ii] = flexce.utils.robust_random_poisson(self.Nstar_stat[ii])
            self.mstar[ii] = self.Nstar[ii] * self.mass_ave
            self.Nstar_left[ii] = self.Nstar[ii]
            self.mstar_left[ii] = self.mstar[ii]

            # SNII and AGB yields
            if self.metallicity[ii] < ylds.snii_agb_z[0]:
                ind_yld[ii] = 0
            elif self.metallicity[ii] > ylds.snii_agb_z[-1]:
                ind_yld[ii] = -1
            else:
                ind_yld[ii] = np.where(self.metallicity[ii] < ylds.snii_agb_z)[0][0]

            # TODO add solar metallicity yields as a config option
            # ind_yld[i] = -1 # uncomment for solar metallicity yields only

            # Evolve stars from previous timesteps
            snii_agb_tmp = np.zeros((self.n_bins, ylds.n_sym))
            # mass_returned_tot = 0.
            mass_remnant_tot = 0.
            for jj in range(1, ii + 1):
                # ind_ev is a list of indices of mass bins from a given birth
                # time-step that will evolve in the current time-step.
                ind = self.ind_ev[ii - jj]

                # abs_yld = net yield + (isotopic mass fraction at birth) * (mass returned to ISM)
                mass_return_to_ism = self.mass_ave[ind] - ylds.snii_agb_rem[ind_yld[jj], ind]
                iso_return_to_ism = self.mfrac[jj] * np.tile(mass_return_to_ism, (ylds.n_sym, 1)).T
                abs_yld = snii_agb_yields[ind_yld[jj], ind] + iso_return_to_ism

                # number of stars to evolve
                N_ev = self.Nstar[jj, ind] * self.frac_ev[ii - jj]

                snii_agb_tmp[ind] += (N_ev * abs_yld.T).T

                if self.params['box']['long_output']:
                    self.snii_agb_net[ii, ind] += (
                        N_ev * snii_agb_yields[ind_yld[jj], ind].T).T

                # mass_returned (300); mass_returned_tot (300)
                # mass_returned = np.sum(N_ev * abs_yld[ind].T, axis=1)
                # mass_returned_tot += mass_returned

                mass_remnant_tot += np.sum(ylds.snii_agb_rem[ind_yld[jj], ind] * N_ev)
                self.Nstar_left[jj, ind] -= N_ev
                self.mstar_left[jj, ind] -= (self.mstar[jj, ind] * self.frac_ev[ii - jj])

            self.snii[ii] = np.sum(snii_agb_tmp[ind8:], axis=0)
            self.agb[ii] = np.sum(snii_agb_tmp[:ind8], axis=0)
            if self.params['box']['long_output']:
                self.snii_agb[ii] = snii_agb_tmp

            # SNIa
            if self.params['snia_dtd']['func'] == 'exponential':
                # mass of WDs that will be formed from the stellar population
                # that is born in the current timestep
                self.Mwd[ii] = np.sum(self.Nstar[ii, ind_ia] *
                                      ylds.agb_rem[ind_yld[ii], ind_ia])
                self.Mwd_Ia[ii] = self.Mwd[ii] * self.params['snia_dtd']['fraction']
                self.Mwd_Ia_init[ii] = self.Mwd[ii] * self.params['snia_dtd']['fraction']

            if preset_state_snia:
                np.random.set_state(self.state['snia'][ii - 1])
            else:
                self.state['snia'].append(np.random.get_state())

            NIa_stat, self.Mwd_Ia = flexce.snia.snia_ev(
                params=self.params['snia_dtd'],
                tstep=ii,
                dtime=self.dtime,
                mstar=self.mstar,
                mstar_tot=self.mstar_left.sum(),
                sfr=self.sfr,
                Mwd_Ia=self.Mwd_Ia,
            )
            self.NIa[ii] = np.random.poisson(NIa_stat)

            self.snia[ii] = ylds.snia_yields * self.NIa[ii]
            self.mremnant[ii] = (mass_remnant_tot - self.NIa[ii] * ylds.snia_yields.sum())

            # gas flows
            self.sf[ii] = np.sum(self.mstar[ii]) * self.mfrac[ii]
            self.outflow[ii] = flexce.outflows.outflow_calc(
                params=self.params['outflows'],
                timestep=ii,
                sfr=self.sf[ii],
                stellar_ejecta=self.snii[ii] + self.agb[ii] + self.snia[ii],
            )

            if self.params['warmgas']['tcool'] > 0.:
                self.gas_cooling[ii] = (self.mwarmgas_iso[ii - 1] * self.dtime /
                                       self.params['warmgas']['tcool'])

            if self.params['inflows']['func'] == 'constant_mgas':
                self.inflow_rate[ii] = (
                    np.sum(
                        self.sf[ii] + self.outflow[ii] - self.gas_cooling[ii] -
                        (self.params['warmgas']['fdirect'] * (self.snii[ii] + self.agb[ii] + self.snia[ii]))
                    ) / self.dtime
                )

            self.inflow[ii] = (
                flexce.inflows.inflow_composition(
                    params=self.params['inflows'],
                    yields=ylds,
                    mgas_iso_last=self.mgas_iso[ii - 1],
                ) * self.inflow_rate[ii] * self.dtime
            )

            # If outflow source is stellar_ejecta (hence feject is non-zero),
            # then we need to add the stellar ejecta lost in the outflow
            # (feject * yields) to counteract the outflow term, which is
            # non-zero for bookkeeping but in this line would correspond to
            # ejecting cold ISM.
            self.mgas_iso[ii] = (
                self.mgas_iso[ii - 1] +
                (self.params['warmgas']['fdirect'] * yields_all_sources) +
                (self.params['outflows']['feject'] * yields_all_sources) +
                self.gas_cooling[ii] -
                self.sf[ii] +
                self.inflow[ii] -
                self.outflow[ii]
            )

            self.mwarmgas_iso[ii] = (
                self.mwarmgas_iso[ii - 1] -
                self.gas_cooling[ii] +
                self.params['warmgas']['fwarm'] * yields_all_sources
            )

            if (ii < 4) and self.params['warmgas']['warmgas_on']:
                self.mwarmgas_iso[ii] += (self.warmgas_ab_pattern *
                                          self.params['warmgas']['mwarmgas_init'] / 4.)

        self.Nstar_left = self.Nstar_left.astype(int)
        self.mstar_left[np.where(np.abs(self.mstar_left) < -1e-8)] = 0.

        if self.params['box']['long_output']:
            self.snii_agb_rec = self.snii_agb - self.snii_agb_net

        print('Time elapsed:', time.time() - start)

        self.outflow_rate = np.sum(self.outflow, axis=1) / self.dtime
        self.check_mass_conservation(ylds)
        self.snii_snia_rate()

        # Set all negative masses equal to a small positive number.
        ind_neg = np.where(self.mgas_iso < 0.)
        self.mgas_iso[ind_neg] = 1e-30

        self.survivors = np.sum(self.Nstar_left[:, 1:], axis=1)
        self.survivors = self.survivors.round().astype(int)
