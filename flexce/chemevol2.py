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
