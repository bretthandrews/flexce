# @Author: Brett Andrews <andrews>
# @Date:   2018-06-18 10:06:09
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 16:07:47

"""
FILE
    yml.py

DESCRIPTION
    Read yaml files.
"""
import re

import yaml


def read(filename):
    """Read in YAML file.

    Args:
        filename: YAML file.
    """

    # Parse scientific notation
    # https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
    loader = yaml.SafeLoader
    loader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    with open(filename, 'r') as fin:
        cfg = yaml.load(fin, Loader=loader)

    return cfg
