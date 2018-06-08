# @Author: Brett Andrews <andrews>
# @Date:   2018-06-07 21:06:39
# @Last modified by:   andrews
# @Last modified time: 2018-06-07 21:06:09

"""
FILE
    inflows.py

DESCRIPTION
    Functions for setting the inflow rate and composition.
"""


def set_outflows(
    outflow_source='ism',
    eta_outflow=1.,
    variable_eta=None,
    feject=0.15
):
    """outflow_source = "ism" (ambient ISM is ejected in the wind; standard
    Mdot_wind = eta * SFR treatment) or "stellar_ejecta" (the yields from
    SNII, SNIa, and AGB stars makes up the wind; from Schoenrich & Binney
    2009).
    """
    params = {
        'outflow_source': outflow_source,
        'eta_outflow': eta_outflow,
        'variable_eta': variable_eta,
        'feject': feject,
    }

    if outflow_source == 'ism':
        feject = 0.

        if variable_eta is not None:
            eta_outflow = variable_eta

    elif outflow_source == 'stellar_ejecta':
        eta_outflow = 0.

    else:
        raise ValueError('Valid outflow sources: "ism" and "stellar_ejecta".')

    return params, eta_outflow, feject


def outflow_calc(self, timestep, sfr, snii, agb, snia):
    if self.outflow_source == 'ism':
        if self.variable_eta is not None:
            return self.eta_outflow[timestep] * sfr
        else:
            return self.eta_outflow * sfr
    elif self.outflow_source == 'stellar_ejecta':
        return self.feject * (snii + agb + snia)
    else:
        print('\nValid outflow sources: "ism" and "stellar_ejecta"\n')
