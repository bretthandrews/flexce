# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:88
# @Last modified by:   andrews
# @Last modified time: 2018-06-06 12:06:88

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

import numpy as np

import flexce.utils
import flexce.imf
from flexce.imf import integrate_multi_power_law
import flexce.lifetimes


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
        self.ind_ev, self.frac_ev, self.frac_ev_tot = flexce.lifetimes.frac_evolve()

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
