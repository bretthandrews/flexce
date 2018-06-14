# @Author: Brett Andrews <andrews>
# @Date:   2018-06-11 13:06:00
# @Last modified by:   andrews
# @Last modified time: 2018-06-14 14:06:35

"""
FILE
    star_formation.py

DESCRIPTION
    Functions for setting the star formation rate.
"""


def set_sflaw(nu_kslaw=1e-9, N_kslaw=1., sfh=None):
    """Set parameters of the star formation law.

    From Kennicutt (1998):
        Sigma_SFR = nu * Sigma_gas^N, where
            Sigma_SFR [=] Msun yr^-1 kpc^-2,
            Sigma_gas [=] Msun yr^-1,
            N = 1.4, and
            nu = 2.5e-4 (= 2.5e-10 if Sigma_SFR is converted into units
                of Msun yr^-1 pc^-2).
    Args:
        nu_kslaw (float): Normalization constant of Kennicutt-Schmidt
            Law. Default is 1e-9 (for Sigma_SFR in units of
            Msun yr^-1 pc^-2).
        N_kslaw (float): Power law slope of Kennicutt-Schmidt Law.
            Default is 1.
        sfh (array): User-defined star formation history. Default is
            ``None``.

    Returns:
        dict: Star formation law parameters.
    """
    params = {'nu': nu_kslaw, 'N': N_kslaw, 'sfh': sfh}

    return params


def sf_law(mgas, params):
    """Calculate star formation rate.

    Calculate the SFR [Msun/yr] given the gas mass and the two free
    parameters that determine the SF law: ``nu_kslaw`` [Gyr^-1] and
    ``N_kslaw`` (~1.0--2.0).

    Args:
        mgas (float): Gas mass.

    """
    return (params['sf']['nu'] * params['box']['area'] *
            (mgas / params['box']['area'])**params['sf']['N'])
