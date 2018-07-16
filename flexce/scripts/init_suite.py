"""
FILE
    init_suite.py

USAGE
    python init_suite.py new-suite-name

DESCRIPTION
    Creates subdirectories names "new-suite-name" within
        - flexce/config
        - flexce/output
        - flexce/plots/config
        - flexce/plots/plots
    to organize the config files and output files for a new suite of
    simulations.
"""


raise ValueError('Decouple package and output directory structures.')  # .............................................


from __future__ import print_function, division, absolute_import
import os
import sys

suite = sys.argv[1]

path_flexce = os.getcwd()
path_flexce_top = os.path.abspath(os.path.join(path_flexce, '..'))

items = [['config'], ['output'], ['plots', 'config'], ['plots', 'plots']]

for it in items:
    parts = [path_flexce_top]
    for ii in it:
        parts.append(ii)
    parts.append(suite)
    path = os.path.join(*parts)
    if os.path.isdir(path):
        print('{} already exists'.format(path))
    else:
        os.mkdir(path)
        print('Created:  {} '.format(path))
