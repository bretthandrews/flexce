flexCE's Change Log
===================

[2.0.0] - Future Date
---------------------

Added
^^^^^
- Ability to import ``flexce``.
- Requirements file.
- Ability to evolve box without passing in arguments (defaults to fiducial simulation from Andrews et al. 2017).

- Ability to pass in config file name to ``ChemEvol`` instance.
- Two infall scenario can be specified in config.
- Solar metallicity only can be specified in config.

- Save options "slim" (disk space efficient), "state" (random number state), and "yields" (full yields used).
- Ability to load and use the random number state from a previous simulation.

- Test to check that simulation output is the same as fiducial simulation from Andrews et al. (2017).
- Significant amount to docstrings.

Removed
^^^^^^^
- Python 2.7 compatibility.
- ``config``, ``plots``, and ``output`` directories.
- ``flexce.utils.load_yields()`` because initializing a ``Yields`` instance directly is simple.
- ``flexce.io.cfg.read_sim_cfg()`` because config file are now YAML files.

Changed
^^^^^^^
- Package name to flexce (all lowercase).
- ``run_flexce`` to be an entry_point script.

- Config file type to YAML (.yml).

- ``ChemEvol.params`` carries all config options.

- Refactored IMF functions into ``imf.py``.
- Refactored inflow functions into ``inflows.py``.
- Refactored stellar lifetimes functions into ``lifetimes.py``.
- Refactored outflow functions into ``outflows.py``.
- Refactored SNIa DTD functions into ``snia.py``.
- Refactored star formation functions into ``star_formation.py``.
- Refactored warm gas reservoir functions into ``warm_gas_reservoir.py``.

- Converted existing tests to ``pytest``.

- Decreased run time by a factor of 6 by reducing the size of arrays created when calculating absolute yields from previous time steps.

- Simplifies call of Abundances.

- ``make_yield_grids.py`` to a console script.

- Increased readbility of code (PEP-8 compliance).

Fixed
^^^^^
- Reading in random state for SNIa draws.
- Issue #4 (over-consumption of gas) by limiting the SFR + outflow rate to be less than the gas mass / time step length. If the Poisson draw from the IMF over-consumes gas, then deterministically populate the IMF by rounding the expectation value of the number of stars in each mass bin down to the nearest integer. This should not make much of a difference in the vast majority of simulations because the numbers will be much larger than 1.


[1.0.1] - 2018/06/04
--------------------

Added
^^^^^
- Dependencies to README.
- Example config file with all options.
- Script to initialize a suite of simulations.


Changed
^^^^^^^
- Utility scripts ``read_sim.py`` and ``read_yields.py``.
- Plotting script to generate [X/Fe]-[Fe/H] plots.
