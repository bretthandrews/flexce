======
flexCE
======
A flexible Chemical Evolution model in python
---------------------------------------------

|ascl:1612.006|

Documentation
^^^^^^^^^^^^^

For the latest documenation please see `here <http://bretthandrews.github.io/flexce>`_.


Python Versions and Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- `Python 2.7 or 3.3+ <https://www.python.org/>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <http://scipy.org/>`_
- `matplotlib <http://matplotlib.org/>`_
- `pandas <http://pandas.pydata.org/>`_
- `seaborn <http://web.stanford.edu/~mwaskom/software/seaborn/index.html>`_

Installation
^^^^^^^^^^^^

Choose a location for flexCE to live. We'll assume that it's in your home
directory ("~/").

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

1. create a new :code:`simN.cfg` file (where **N** is a number)
2. edit :code:`simN.cfg`
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


Attribution
^^^^^^^^^^^
Please cite `Andrews et al. (2017)
<https://ui.adsabs.harvard.edu/abs/2017ApJ...835..224A>`_ if you find
this code helpful in your research. You may also want to check out
`Weinberg, Andrews, & Freudenburg (2017)
<https://ui.adsabs.harvard.edu/abs/2017ApJ...837..183W>`_ for a companion
analytic model.

The BibTeX entry for Andrews et al. (2017)::

    @ARTICLE{2017ApJ...835..224A,
           author = {{Andrews}, B.~H. and {Weinberg}, D.~H. and {Sch{\"o}nrich}, R. and {Johnson}, J.~A.},
            title = "{Inflow, Outflow, Yields, and Stellar Population Mixing in Chemical Evolution Models}",
          journal = {\apj},
    archivePrefix = "arXiv",
           eprint = {1604.08613},
         keywords = {galaxies: ISM, Galaxy: evolution, Galaxy: formation, Galaxy: general, Galaxy: stellar content, stars: abundances},
             year = 2017,
            month = feb,
           volume = 835,
              eid = {224},
            pages = {224},
              doi = {10.3847/1538-4357/835/2/224},
           adsurl = {http://adsabs.harvard.edu/abs/2017ApJ...835..224A},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }



.. |ascl:1612.006| image:: https://img.shields.io/badge/ascl-1612.006-blue.svg?colorB=262255
   :target: http://ascl.net/1612.006
