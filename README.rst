======
flexCE
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

    cd ~/flexce
     python setup.py install


**Generate yields**::

    cd ~/flexce/flexce/
     python make_yield_grids.py


Run the code
^^^^^^^^^^^^
::

    cd ~/flexce/flexce/
    cd ~/flexce/flexce/
     python flexce.py ../config/sim0.cfg



**Change the parameters of the code**::

1. create a new :code:`simN.cfg` file (where **N** is a number)
2. edit :code:`simN.cfg`
3. re-run the code::

    python flexce.py ../config/simN.cfg


**Plot [O/Fe]-[Fe/H]**::

    cd flexce/flexce/plot
    cd flexce/flexce/plot
     python plot_xfe_feh.py ofe_feh_sim0.cfg


**Go to output plot directory**::

    cd flexce/plots/plots
    cd flexce/plots/plots



License
^^^^^^^
Copyright 2016 Brett Andrews.

flexCE is free software made available under the MIT License.  For details see
the LICENSE file.


Attribution
^^^^^^^^^^^
Please cite `Andrews et al. (2016)
<http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1604.08613>`_ if you find
this code helpful in your research. You may also want to check out
`Weinberg, Andrews, & Freudenburg (2016)
<http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1604.07435>`_ for a companion
analytic model.

The BibTeX entry for Andrews et al. (2016)::

    @ARTICLE{2016arXiv160408613A,
           author = {{Andrews}, B.~H. and {Weinberg}, D.~H. and {Sch{\"o}nrich}, R. and {Johnson}, J.~A.},
            title = "{Inflow, Outflow, Yields, and Stellar Population Mixing in Chemical Evolution Models}",
          journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
           eprint = {1604.08613},
         keywords = {Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
             year = 2016,
            month = apr,
           adsurl = {http://adsabs.harvard.edu/abs/2016arXiv160408613A},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

.. |ascl:1612.006| image:: https://img.shields.io/badge/ascl-1612.006-blue.svg?colorB=262255
   :target: http://ascl.net/1612.006
