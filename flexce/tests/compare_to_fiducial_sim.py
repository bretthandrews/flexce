import pickle
import sys

from flexce.chemevol import ChemEvol

sys.path.append('/Users/andrews/projects/pcaa_chemevol/')
import chemevol_main

# load fiducial simulation
with open('/Users/andrews/projects/pcaa_chemevol/sims/paper_masscut01/runs/box0.pck', 'rb') as fin:
    box0 = pickle.load(fin)

# run flexCE with fiducial parameters
gal = ChemEvol()

# Parameters
# Inflows
for k, v in gal.params['inflows'].items():
    if k == 'coeff':
        k = 'k'
    assert box0.param['inflow'][k] == v, k

# Outflows
for k, v in gal.params['outflows'].items():
    if k == 'eta':
        k = 'eta_outflow'
    elif k == 'source':
        k = 'outflow_source'
    if k != 'feject':
        assert box0.param['outflow'][k] == v, k

# Warm gas reservoir
for k, v in gal.params['warmgas'].items():
    if k not in ['fdirect', 'tcool']:
        assert box0.param['warmgas'][k] == v, k

# Star formation
for k, v in gal.params['sf'].items():
    assert box0.param['sf'][k] == v, k

# SNIa
for k, v in gal.params['snia_dtd'].items():
    if k == 'func':
        assert box0.param['snia'][k] == v, k
    elif k == 'min_time':
        k = 'min_snia_time'
    elif k == 'fraction':
        k = 'snia_fraction'
    assert box0.param['snia']['k'][k] == v, k
