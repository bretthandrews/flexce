# @Author: Brett Andrews <andrews>
# @Date:   2018-06-05 11:06:88
# @Last modified by:   andrews
# @Last modified time: 2018-06-05 12:06:84

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

class ChemEvol:
    """Chemical evolution model.

    Args:
        params (dict): parameters to set up simulation run.
    """
    def __init__(self, params):
        # if any particular parameters are not specified,
        # then set them at the function call.
        self.params = params


    def set_yields(self):
        # check for existing yields
        # if not, then calculate them.
        pass


    def run(self):
        self.mass_bins = flexce.utils.set_mass_bins(self.params['mass_bins'])
        self.set_box(self.params['box'])
        self.set_yields(self.params['yields'], self.mass_bins)
        self.snia_dtd(self.params['snia_dtd'])  # TODO rename set_snia_dtd
        self.inflow_rx(self.params['inflows'])  # TODO rename set_inflow
        self.outflow_rx(self.params['outflows'])  # TODO rename set_outflow
        self.warmgasres_rx(self.params['warmgasres'])  # TODO rename set_warmgas_res
        self.star_formation(self.params['sf'])  # TODO rename set_star_formation


    def set_box(self, mass_bins, radius=10., time_tot=12000., dt=30.,
                 imf='kroupa', imf_alpha=None, imf_mass_breaks=None,
                 sim_id=None):
        """Initialize box.

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
