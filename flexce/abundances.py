# @Author: Brett Andrews <andrews>
# @Date:   2018-04-16 20:04:30
# @Last modified by:   andrews
# @Last modified time: 2018-06-26 12:06:28

"""
FILE
    abundances.py

DESCRIPTION
    Compute abundances.
"""

import os
from os.path import join
import re

import numpy as np
import pandas as pd

from flexce.yields import Yields


def calc_abundances(sym, mgas, survivors, time, parameters):
    """Calculate abundances of box.

    Wrapper for Abundances class.

    Args:
        sym (array): Isotope abbreviations.
        mgas (array): Mass of each isotope in gas-phase at each timestep.
        survivors (array): Number of stars from each timestep that survive to
            the end of the simulation.
        time (array): time in Myr.
        parameters (dict): parameters of the simulation.

    Returns:
        Abundances instance
    """
    abund = Abundances(sym, mgas, survivors, time, parameters)
    abund.load_solar_abund()
    abund.calc_abundances()
    apogee_el = np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
                          'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni'])
    abund.select_elements(apogee_el)
    return abund


def split_element_mass(isotopes):
    """Convert isotope abbreviation to element and mass.

    Takes an array of isotopes (element & mass) and creates a separate
    arrays of element symbols and isotope masses with the same length as
    the isotope array. Also creates a dictionary with the indices of each
    element in the isotope array."""

    sym = np.array(['' for _ in isotopes], dtype='<U2')
    iso_mass = np.zeros(len(isotopes), dtype=int)

    elements = []

    for ii, iso in enumerate(isotopes):
        match = re.match(r"([a-z]+)([0-9]+)", iso, re.I)

        if match:
            sym[ii], iso_mass[ii] = match.groups()

        if sym[ii] not in elements:
            elements.append(sym[ii])

    ind_element = pd.DataFrame(np.array([sym == it for it in elements]).T, columns=elements)

    return iso_mass, ind_element


class Abundances(object):
    """Compute abundances of model.
    """

    def __init__(self, box, ylds=None, params=None):
        """Initialize Abundances instance.

        Args:
            box: ChemEvol instance.
            ylds: Yields instance. Default is ``None``.
            params (dict): Abundances parameters. Default is ``None``.
        """

        if ylds is None:
            ylds = Yields(params=box.params['yields'], mass_bins=box.mass_bins)

        default = {
            'solar': {
                'source': 'lodders',
            },
            'sim_id': box.params['box']['sim_id'],
        }

        params = params if params is not None else {}
        params = {k: v if k not in params.keys() else params[k] for k, v in default.items()}

        # self.setup()
        iso_mass, ind_element = split_element_mass(ylds.sym)
        self.mgas_iso = box.mgas_iso
        self.n_steps = len(self.mgas_iso)
        self.weights = box.survivors
        self.param = box.params
        self.solar = self.load_solar_abund(params['solar']['source'])
        self.calc_abundances()
        apogee_el = np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
                              'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni'])
        self.select_elements(apogee_el)

#        self.apogee_elements()

    # def setup(self):
    #     """Read in atomic numbers and element abbreviations."""
    #     path_yldgen = join(os.path.dirname(__file__), 'data', 'yields', 'general')
    #     el_sym = pd.read_csv(
    #         join(path_yldgen, 'sym_atomicnum.txt'),
    #         delim_whitespace=True,
    #         usecols=[0, 1],
    #         names=['num', 'el']
    #     )
    #     self.elements2 = el_sym
    #
    #     self.all_atomic_num = np.array(el_sym['num'])
    #     self.all_elements = np.array(el_sym['el'])

    def load_solar_abund(self, elements, source='lodders'):
        """Read in solar abundances.

        Args:
            elements (list): Elements.
            source (str): Reference for solar abundances.  Default is
                'lodders'.
        """
        if source == 'lodders':
            path_yldgen = join(os.path.dirname(__file__), 'data', 'yields', 'general')

            solar_cat = pd.read_csv(
                join(path_yldgen, 'lodders03_solar_photosphere.txt'),
                delim_whitespace=True,
                skiprows=8,
                usecols=[0, 1],
                names=['el', 'ab'],
            )

            # self.solar_element = np.array(solar_ab['el'])
            # self.solar_ab = np.array(solar_ab['ab'])
            # solar_h = np.zeros(self.n_elements)
            # self.solar_fe = np.zeros(self.n_elements)
            xh = solar_cat['ab'][solar_cat['el'] == el]
            xfe = np.log10(
                10.**(solar_cat['ab'][solar_cat['el'] == el] - 12) /
                10.**(solar_cat['ab'][solar_cat['el'] == 'Fe'] - 12)
            )
            for ii, el in enumerate(self.elements):
                self.solar_h[ii] = self.solar_ab[self.solar_element == el]
                self.solar_fe[ii] = np.log10(
                    10.**(self.solar_ab[self.solar_element == el] - 12.) /
                    10.**(self.solar_ab[self.solar_element == 'Fe'] - 12.)
                )

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
            ind2.append(np.where(self.elements2['el'] == item)[0][0])
            # ind2.append(np.where(self.all_elements == item)[0][0])

        self.atomic_num_out = self.elements2['num'][ind2]
        # self.atomic_num_out = self.all_atomic_num[ind2]
