"""
FILE
    make_yield_grids.py

DESCRIPTION
    Generate yield grids.

USAGE
    python make_yield_grids.py [ww95]
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys

# ----- Set paths -----
path_flexce = join(os.path.abspath('.'), '')
path_data = join(path_flexce, 'data')
path_yields = join(path_data, 'yields')
path_calc_yields = join(path_flexce, 'calc_yields')
# ---------------------


make_ww95 = False

if len(sys.argv) > 1:
    if sys.argv[1] == 'ww95':
        make_ww95 = True
    else:
        print(' '.join(['\nRun without keyword or with one of the following',
                        'accepted keywords:\nww95\n']))
        sys.exit(1)


print('\nGenerating pickled yield grids...\n')


ylds = {
    'busso01': {},
    'cescutti06': {},
    'iwamoto99': {},
    'karakas10': {},
    'limongi06': {},
    'sneden08': {},
    'ww95': {}
    }

for k in sorted(ylds.keys()):
    py_file = k + '_yields.py'
    pck_file = 'interp_yields.pck'
    path = k + '/'
    if k == 'busso01':
        pck_file = 'busso01_yields.pck'
    elif k == 'cescutti06':
        pck_file = 'cescutti06_yields.pck'
    elif k == 'iwamoto99':
        pck_file = 'w70_yields.pck'
    elif k == 'karakas10':
        path = path + 'iso_yields/'
    if k == 'limongi06':
        py_file = 'limongi_chieffi_yields.py'
        path = path + 'iso_yields/'
    elif k == 'sneden08':
        py_file = 'sneden08.py'
        pck_file = 'sneden08.pck'
        path = 'general/'
    elif k == 'ww95':
        path = path + 'half_fe/'
    ylds[k]['script'] = py_file
    ylds[k]['pck'] = pck_file
    ylds[k]['path'] = path


keys = ['limongi06', 'ww95']
for k in ylds.keys():
    if k not in ['limongi06', 'ww95']:
        keys.append(k)

for k in keys:
    if os.path.isfile(join(path_yields, ylds[k]['path'], ylds[k]['pck'])):
        print('yield grids already exist: ', k)
    else:
        if k == 'ww95':
            if not make_ww95:
                continue
        print('\ngenerating yields: ', k)
        os.system('python ' + join(path_calc_yields, ylds[k]['script']))

print()
