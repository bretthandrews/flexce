# @Author: Brett Andrews <andrews>
# @Date:   2018-04-16 20:04:47
# @Last modified by:   andrews
# @Last modified time: 2018-06-15 16:06:86

"""
FILE
    make_yield_grids.py

USAGE
    make_yield_grids [--ww95]

DESCRIPTION
    Generate yield grids.
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join

import click


@click.command()
@click.option('--ww95', is_flag=True, default=False,
              help='Create Woosley & Weaver (1995) yield grid.')
def main(ww95):
    path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
    path_data = join(path_flexce, 'data')
    path_yields = join(path_data, 'yields')
    path_calc_yields = join(path_flexce, 'calc_yields')

    make_ww95 = ww95

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

    for kk in sorted(ylds.keys()):
        py_file = kk + '_yields.py'
        pck_file = 'interp_yields.pck'

        path = kk + '/'
        if kk == 'busso01':
            pck_file = 'busso01_yields.pck'

        elif kk == 'cescutti06':
            pck_file = 'cescutti06_yields.pck'

        elif kk == 'iwamoto99':
            pck_file = 'w70_yields.pck'

        elif kk == 'karakas10':
            path = join(path, 'iso_yields')

        if kk == 'limongi06':
            py_file = 'limongi_chieffi_yields.py'
            path = join(path, 'iso_yields')

        elif kk == 'sneden08':
            py_file = 'sneden08.py'
            pck_file = 'sneden08.pck'
            path = 'general'

        elif kk == 'ww95':
            path = join(path, 'half_fe')

        ylds[kk]['script'] = py_file
        ylds[kk]['pck'] = pck_file
        ylds[kk]['path'] = path

    keys = ['limongi06', 'ww95']
    for kk in ylds.keys():
        if kk not in ['limongi06', 'ww95']:
            keys.append(kk)

    for kk in keys:
        if os.path.isfile(join(path_yields, ylds[kk]['path'], ylds[kk]['pck'])):
            print('yield grids already exist: ', kk)

        else:
            if kk == 'ww95':
                if not make_ww95:
                    continue
            print('\ngenerating yields: ', kk)
            os.system('python ' + join(path_calc_yields, ylds[kk]['script']))


if __name__ == '__main__':
    main()
