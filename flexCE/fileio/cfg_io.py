from __future__ import print_function, division, absolute_import

try:
    import ConfigParser as configparser
except ImportError as e:
    import configparser


def read_sim_cfg(file_in):
    """Read in config file.

    Args:
        file_in (str): name of config file (can include path) with format
            'sim<sim_id>.cfg' where sim_id is an integer.

    Returns:
        str, dict, dict, dict, dict, dict, dict, dict, dict:
        sim_id\: simulation ID number, yld_args\: args to initialize instance
        of Yield class, initialize_args\: args to initialize instance of
        ChemEvol class, mass_bins_args\: args to define stellar mass bins,
        snia_dtd_args\: args to set SNIa delay time distribution of ChemEvol
        instance, inflows_args\: args to set inflow rate and composition of
        ChemEvol instance, outflows_args\: args to set outflow rate and
        composition of ChemEvol instance, warmgasres_args\: turn on warm ISM
        reservoir in ChemEvol instance, sf_args\: args to set star formation
        rate in ChemEvol instance.

    """
    sim_id = file_in.split('/')[-1].strip('sim').strip('.cfg')

    yld_args = {}
    initialize_args = {'sim_id': sim_id}
    mass_bins_args = {}
    snia_dtd_args = {}
    inflows_args = {}
    outflows_args = {}
    warmgasres_args = {}
    sf_args = {}

    f = open(file_in, 'r')
    for line in f:
        if (line[0] != '#') and (line != '\n'):
            if ' = ' in line:
                k = line.strip().split(' = ')[0]
                v = line.strip().split(' = ')[1]
                try:
                    v = float(v)
                except ValueError:
                    pass
                try:
                    if v[0] == '[':
                        v = v.strip('[').strip(']').split(', ')
                        try:
                            v = [float(item) for item in v]
                        except ValueError:
                            if v == ['']:
                                v = []
                except (ValueError, TypeError):
                    pass
                if v == 'True':
                    v = True
                elif v == 'False':
                    v = False
            if 'yields' in k:
                if (',' in v) and (v[-1] != ','):
                    yld_args[k.split('yields_')[1]] = v.split(', ')
                elif (',' in v) and (v[-1] == ','):
                    yld_args[k.split('yields_')[1]] = [v.strip(',')]
                else:
                    yld_args[k.split('yields_')[1]] = v
            elif 'initialize' in k:
                initialize_args[k.split('initialize_')[1]] = v
            elif 'mass_bins' in k:
                mass_bins_args[k.split('mass_bins_')[1]] = v
            elif 'snia_dtd' in k:
                if k.split('snia_dtd_')[1] == 'func':
                    snia_dtd_args[k.split('snia_dtd_')[1]] = v
                else:
                    if 'kwargs' not in snia_dtd_args:
                        snia_dtd_args['kwargs'] = {}
                    snia_dtd_args['kwargs'][k.split('snia_dtd_')[1]] = v
            elif 'inflows' in k:
                if k.split('inflows_')[1] in ['M1', 'M2', 'b1', 'b2']:
                    if 'k' not in inflows_args:
                        inflows_args['k'] = {}
                    inflows_args['k'][k.split('inflows_')[1]] = v
                else:
                    inflows_args[k.split('inflows_')[1]] = v
            elif 'outflows' in k:
                outflows_args[k.split('outflows_')[1]] = v
            elif 'warmgasres' in k:
                warmgasres_args[k.split('warmgasres_')[1]] = v
            elif 'sf' in k:
                sf_args[k.split('sf_')[1]] = v

    f.close()
    return (sim_id, yld_args, initialize_args, mass_bins_args, snia_dtd_args,
            inflows_args, outflows_args, warmgasres_args, sf_args)


def _get_section(config, section):
    """Get all options from a section in a config file.

    Args:
        config: configparser instance.
        section (str): Section name.

    Returns:
        dict: Contents of a section from config file.
    """
    dict1 = {}
    options = config.options(section)
    for option in options:
        val = config.get(section, option)
        if ',' in val:
            val = [item.strip() for item in val.split(',')]
        dict1[option] = val
    return dict1


def read_plot_config(filename):
    """Read in plotting config file.

    Args:
        filename (str): Full path and name of config file.

    Returns:
        dict: Contents of config file.
    """
    config = configparser.ConfigParser()
    config.read(filename)
    out = {}
    for section in config.sections():
        out[section] = _get_section(config, section)
    return out
