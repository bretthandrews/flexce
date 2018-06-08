# @Author: Brett Andrews <andrews>
# @Date:   2018-06-07 21:06:39
# @Last modified by:   andrews
# @Last modified time: 2018-06-07 21:06:90

"""
FILE
    inflows.py

DESCRIPTION
    Functions for setting the inflow rate and composition.
"""


def set_outflows(
    source='ism',
    eta=1.,
    variable_eta=None,
    feject=0.15
):
    """Set outflow parameters.

    Available outflow sources:
        ism: ambient ISM is ejected in the wind (Mdot_out = eta * SFR)
        stellar_ejecta: the yields from SNII, SNIa, and AGB stars make
            up the wind as used in Schoenrich & Binney (2009).

    Args:
        source (str): Source of outflowing gas. Default is 'ism'.
        eta (float): Outflow mass-loading parameter. Only used if
            ``source`` is 'ism'. Default is 1.
        variable_eta (array): Time variable outflow mass-loading
            parameter. Only used if ``source`` is 'ism'. Default is
            ``None``.
        feject (float): Fraction of new stellar ejecta that leaves the
            galaxy in the outflow (if ``source`` is 'stellar_ejecta').
            Default is 0.15.

    """
    params = {
        'source': source,
        'eta': eta,
        'variable_eta': variable_eta,
        'feject': feject,
    }

    if source == 'ism':
        feject = 0.

        if variable_eta is not None:
            eta = variable_eta

    elif source == 'stellar_ejecta':
        eta = 0.

    else:
        raise ValueError('Valid outflow sources: "ism" and "stellar_ejecta".')

    return params, eta, feject


def outflow_calc(params, eta, feject, timestep, sfr, stellar_ejecta):

    if params['source'] == 'ism':
        if params['variable_eta'] is not None:
            return eta[timestep] * sfr

        else:
            return eta * sfr

    elif params['source'] == 'stellar_ejecta':
        return feject * (stellar_ejecta)

    else:
        raise ValueError('Valid outflow sources: "ism" and "stellar_ejecta".')

    return eta
