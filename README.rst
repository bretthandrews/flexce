======
flexCE
======
A flexible Chemical Evolution model
-----------------------------------

Choose a location for flexCE to live. We'll assume that it's in your home
directory ("~/").



Installation
^^^^^^^^^^^^
::
    
    cd ~/flexCE
    python setup.py install


**Generate yields**::

    cd ~/flexCE/flexCE/
    python make_yield_grids.py


Run the code
^^^^^^^^^^^^
::

    cd ~/flexCE/flexCE/
    python flexce.py ../config/sim0.cfg



**Change the parameters of the code**::

1. create a new `simN.cfg` file (where **N** is a number)
2. edit `simN.cfg`
3. re-run the code::

    python flexce.py ../config/simN.cfg


**Plot [O/Fe]-[Fe/H]**::

    cd flexCE/flexCE/plot
    python plot_xfe_feh.py ofe_feh_sim0.cfg


**Go to output plot directory**::

    cd flexCE/plots/plots



License
^^^^^^^
Copyright 2016 Brett Andrews.

flexCE is free software made available under the MIT License. For details see
the LICENSE file.
