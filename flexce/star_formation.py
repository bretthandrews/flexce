# @Author: Brett Andrews <andrews>
# @Date:   2018-06-11 13:06:00
# @Last modified by:   andrews
# @Last modified time: 2018-06-11 13:06:30

"""
FILE
    star_formation.py

DESCRIPTION
    Functions for setting the star formation rate.
"""


def star_formation(self, nu_kslaw=2.5e-10, N_kslaw=1.4):
    '''From Kennicutt (1998): Sigma_SFR = 2.5e-4 * Sigma_gas^1.4, where
    Sigma_SFR [=] Msun yr^-1 kpc^-2 and Sigma_gas [=] Msun yr^-1.  nu =
    2.5-10 if Sigma_SFR is converted into units of Msun yr^-1 pc^-2.  '''
    self.nu_kslaw = nu_kslaw
    self.N_kslaw = N_kslaw
    self.sf_param = dict(nu=nu_kslaw, N=N_kslaw)


def sf_law(self, mgas):
    '''Calculate the SFR [Msun/yr] given the gas mass and the two free
    parameters that determine the SF law: nu_kslaw [Gyr^-1] and N_kslaw
    (~1.0--2.0).

    From Kennicutt (1998): Sigma_SFR = 2.5e-4 * Sigma_gas^1.4, where
    Sigma_SFR [=] Msun yr^-1 kpc^-2 and Sigma_gas [=] Msun yr^-1.  nu =
    2.5-10 if Sigma_SFR is converted into units of Msun yr^-1 pc^-2.
    '''
    return (self.nu_kslaw) * self.area * (mgas / self.area)**self.N_kslaw
