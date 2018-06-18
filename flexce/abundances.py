"""Compute abundances."""
from os.path import join
import sys
import re
import numpy as np
import pandas as pd


def calc_abundances(path, sym, mgas, survivors, time, parameters):
    """Calculate abundances of box.

    Wrapper for Abundances class.

    Args:
        path (str): data directory.
        sym (array): Isotope abbreviations.
        mgas (array): Mass of each isotope in gas-phase at each timestep.
        survivors (array): Number of stars from each timestep that survive to
            the end of the simulation.
        time (array): time in Myr.
        parameters (dict): parameters of the simulation.

    Returns:
        Abundances instance
    """
    abund = Abundances(path, sym, mgas, survivors, time, parameters)
    abund.load_solar_abund()
    abund.calc_abundances()
    apogee_el = np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
                          'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni'])
    abund.select_elements(apogee_el)
    return abund


class Abundances:
    """Compute abundances of model.
    """

    def __init__(self, path_parent, sym_iso, mgas_iso, weight, timesteps=None,
                 sim_params=None):
        """Initialize Abundances instance.

        Args:
            path_parent (str): data directory.
            sym_iso (array): isotope abbreviations.
            mgas_iso (array): Gas mass of each isotope from ChemEvol instance.
            weight (array): Number of surviving stars from each generation.
            timesteps (array): Timesteps from ChemEvol instance. Defaults to
                None.
            sim_params (dict): Parameters of ChemEvol instance.
        """
        self.path_parent = path_parent
        self.path_yields = join(path_parent, 'yields')
        self.path_yldgen = join(self.path_yields, 'general')
        self.isotope = sym_iso
        self.setup()
        self.split_element_mass()
        self.mgas_iso = mgas_iso
        self.n_steps = len(self.mgas_iso)
        self.survivors = weight
        self.t = timesteps
        self.param = sim_params
        self.sim_id = sim_params['box']['sim_id']
#        self.apogee_elements()

    def setup(self):
        """Read in atomic numbers and element abbreviations."""
        el_sym = pd.read_csv(join(self.path_yldgen, 'sym_atomicnum.txt'),
                             delim_whitespace=True, usecols=[0, 1],
                             names=['num', 'el'])
        self.all_atomic_num = np.array(el_sym['num'])
        self.all_elements = np.array(el_sym['el'])

    def split_element_mass(self):
        """Convert isotope abbreviation to element and mass.

        Takes an array of isotopes (element & mass) and creates a separate
        arrays of element symbols and isotope masses with the same length as
        the isotope array. Also creates a dictionary with the indices of each
        element in the isotope array."""
        self.n_isotope = len(self.isotope)
        if sys.version_info[0] < 3:
            self.sym = np.array(['' for i in range(self.n_isotope)], dtype='|S2')
        else:
            self.sym = np.array(['' for i in range(self.n_isotope)], dtype='<U2')
        self.isotope_mass = np.zeros(self.n_isotope, dtype=int)
        self.elements = []
        for i in range(self.n_isotope):
            match = re.match(r"([a-z]+)([0-9]+)", self.isotope[i], re.I)
            if match:
                self.sym[i], self.isotope_mass[i] = match.groups()
            if self.sym[i] not in self.elements:
                self.elements.append(self.sym[i])
        self.elements = np.array(self.elements)
        self.n_elements = len(self.elements)
        self.ind_element = {}
        for item in self.elements:
            self.ind_element[item] = np.where(self.sym == item)[0]

    def load_solar_abund(self, source='lodders'):
        """Read in solar abundances.

        Args:
            source (str): Reference for solar abundances.  Defaults to
                'lodders'.
        """
        if source == 'lodders':
            fin = join(self.path_yldgen, 'lodders03_solar_photosphere.txt')
            solar_ab = pd.read_csv(fin, delim_whitespace=True, skiprows=8,
                                   usecols=[0, 1], names=['el', 'ab'])
            self.solar_element = np.array(solar_ab['el'])
            self.solar_ab = np.array(solar_ab['ab'])
            self.solar_h = np.zeros(self.n_elements)
            self.solar_fe = np.zeros(self.n_elements)
            for i in range(self.n_elements):
                ind = np.where(self.solar_element == self.elements[i])[0]
                ind_fe = np.where(self.elements == 'Fe')
                self.solar_h[i] = self.solar_ab[ind]
                self.solar_fe[i] = np.log10(
                    10.**(self.solar_ab[ind] - 12.) /
                    10.**(self.solar_ab[ind_fe] - 12.))
        # elif source == 'asplund':
        # elif source == 'aspcap':
        #  see deprecated apogee_solar_abundances function
        # else:
        #    Raise exception

    def calc_abundances(self):
        """Calculate abundances relative to hydrogen and iron."""
        self.ngas_iso = np.divide(self.mgas_iso, self.isotope_mass)
        self.niso_h = np.array([
            self.ngas_iso[j] / self.ngas_iso[j, self.ind_element['H']].sum()
            for j in range(1, self.n_steps)])
        self.niso_fe = np.array([
            self.ngas_iso[j] / self.ngas_iso[j, self.ind_element['Fe']].sum()
            for j in range(1, self.n_steps)])
        self.xh_abs = np.log10([
            np.sum(self.niso_h[:, self.ind_element[item]], axis=1)
            for item in self.elements]) + 12.
        self.xfe_abs = np.log10([
            np.sum(self.niso_fe[:, self.ind_element[item]], axis=1)
            for item in self.elements])
        self.xh_all = np.subtract(self.xh_abs.T, self.solar_h).T
        self.feh = self.xh_all[np.where(self.elements == 'Fe')][0]
        self.xfe_all = np.subtract(self.xfe_abs.T, self.solar_fe).T

    def select_elements(self, el=np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al',
                                           'Si', 'S', 'K', 'Ca', 'Ti', 'V',
                                           'Cr', 'Mn', 'Co', 'Ni'])):
        """Downselect abundances to elements of interest.

        Args:
            el (array): array of elements. Defaults to APOGEE set of elements
                (i.e., np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
                'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni']).
        """
        ind = []
        for item in el:
            ind.append(np.where(self.elements == item)[0][0])
        self.xfe = self.xfe_all[ind]
        self.elements_out = self.elements[ind]
        ind2 = []
        for item in el:
            ind2.append(np.where(self.all_elements == item)[0][0])
        self.atomic_num_out = self.all_atomic_num[ind2]
