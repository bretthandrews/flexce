# @Author: Brett Andrews <andrews>
# @Date:   2018-06-07 21:06:85
# @Last modified by:   andrews
# @Last modified time: 2018-06-07 21:06:79

"""
FILE
    warm_gas_reservoir.py

DESCRIPTION
    Functions for controlling gas flow into and out of a warm gas
    reservoir.
"""

import pandas as pd


def warmgasres_rx(self, mwarmgas_init=0., fdirect=0.01, tcool=1200.,
                  warmgas=True):
    '''Parameters that control gas flow into and out of the warm gas
    reservoir.
    Schoenrich & Binney (2009) fiducial values:
    mwarmgas_init = 5e8 Msun
    fdirect = 0.01 (feject=0.15 for R < 3.5 kpc and 0.04 for R > 3.5 kpc)
    tcool = 1200 Myr
    fwarm = 1 - fdirect - feject
    '''
    self.warmgas_on = warmgas
    if warmgas:
        filename = 'warmgas_abundance_pattern.txt'
        tmp = pd.read_csv(self.path_yldgen + filename,
                          delim_whitespace=True, skiprows=10,
                          names=['el', 'ab'])
        self.warmgas_ab_pattern = np.array(tmp['ab'])
        self.mwarmgas_init = mwarmgas_init
        self.fdirect = fdirect
        self.tcool = tcool
        self.fwarm = 1. - fdirect
        if self.outflow_source == 'stellar_ejecta':
            self.fwarm -= self.feject
    else:
        self.warmgas_ab_pattern = 0.
        self.mwarmgas_init = 0.
        self.fdirect = 1. - self.feject
        self.tcool = 0.
        self.fwarm = 0.
    self.warmgasres_param = dict(warmgas_on=warmgas,
                                 mwarmgas_init=mwarmgas_init,
                                 fdirect=fdirect, fwarm=self.fwarm,
                                 tcool=tcool)
