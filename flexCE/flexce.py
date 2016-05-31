"""Main script.

Command line args:
    config file that starts with "sim".

"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

import numpy as np

from fileio.cfg_io import read_sim_cfg
from fileio.pickle_io import pickle_write
from fileio.txt_io import txt_write

from utils import define_mass_bins
from utils import load_yields

from chemevol import ChemEvol
from abundances import Abundances


def evolve(yld, initialize_kws, snia_dtd_kws, inflows_kws, outflows_kws,
           warmgasres_kws, sf_kws):
    """Evolve the galaxy.

    Args:
        yld: Yields instance
        initialize_kws (dict): args to initialize instance of ChemEvol class.
        mass_bins_args (dict): args to define stellar mass bins.
        snia_dtd_kws (dict): args to set SNIa delay time distribution of
            ChemEvol instance.
        inflows_kws (dict): args to set inflow rate and composition of
            ChemEvol instance.
        outflows_kws (dict): args to set outflow rate and composition of
            ChemEvol instance
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


def calc_abundances(path, sym, mgas, survivors, time, parameters, sim_id):
    """Calculate abundances of box.

    Args:
        path (str): data directory.
        sym (array): Isotope abbreviations.
        mgas (array): Mass of each isotope in gas-phase at each timestep.
        survivors (array): Number of stars from each timestep that survive to
            the end of the simulation.
        time (array): time in Myr.
        parameters (dict): parameters of the simulation.
        sim_id (str): simulation ID number.

    Returns:
        Abundances instance
    """
    abund = Abundances(path, sym, mgas, survivors, time, parameters, sim_id)
    abund.load_solar_abund()
    abund.calc_abundances()
    apogee_el = np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
                          'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni'])
    abund.select_elements(apogee_el)
    return abund


def output(path, sim_id, gal, abund):
    """Write simulation results to pickle and txt files.

    Args:
        path (str): output directory.
        sim_id (str): simulation ID number.
        gal: ChemEvol instance.
        abund: Abundances instance.
    """
    path_sim = join(path, ''.join(['sim', sim_id]))
    if not os.path.isdir(path_sim):
        os.mkdir(path_sim)

    pickle_write(gal, join(path_sim, ''.join(('box', sim_id, '.pck'))))
    pickle_write(abund, join(path_sim, ''.join(('ab', sim_id, '.pck'))))

    txt_write(path_out, sim_id, gal, abund)


if __name__ == '__main__':
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
    path_flexce_top = os.path.abspath(join(path_flexce, '..'))
    path_data = join(path_flexce, 'data', '')
    path_config = join(path_flexce_top, 'config', '')
    path_out = join(path_flexce_top, 'output', '')

    argv = None
    if argv is None:
        argv = sys.argv
    try:
        # TODO autosearch for config file in the standard config directory
        file_in = argv[1]
        if file_in.split('/')[-1][:3] != 'sim':
            raise IndexError
    except IndexError:
        fname = 'sim0.cfg'
        home = os.path.expanduser('~')
        path_config = join(home, 'flexCE', 'examples', '')
        file_in = join(path_config, fname)
        print('\nUsing default parameters in \n{}'.format(file_in))

    (simulation_id, yld_args, initialize_args, mass_bins_args, snia_dtd_args,
     inflows_args, outflows_args, warmgasres_args, sf_args) = \
        read_sim_cfg(file_in)
    mass_bins = define_mass_bins(**mass_bins_args)
    ylds = load_yields(path_data, yld_args, mass_bins)
    box = evolve(ylds, initialize_args, snia_dtd_args, inflows_args,
                 outflows_args, warmgasres_args, sf_args)
    ab = calc_abundances(path_data, ylds.sym, box.mgas_iso, box.survivors,
                         box.t, box.param, box.sim_id)
    output(path_out, simulation_id, box, ab)

# TODO (specify elements_out in config file)
