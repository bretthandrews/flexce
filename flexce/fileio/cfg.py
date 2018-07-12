# @Author: Brett Andrews <andrews>
# @Date:   2018-06-21 11:06:11
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 16:07:28

"""
FILE
    cfg.py

DESCRIPTION
    Read config files.
"""

from __future__ import print_function, division, absolute_import

try:
    import ConfigParser as configparser
except ImportError as e:
    import configparser


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
