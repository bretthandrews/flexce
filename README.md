# flexCE
A flexible Chemical Evolution model for all.

Choose a location for flexCE to live. We'll assume that it's in your home
directory ("~/")



### Installation
cd ~/flexCE
python setup.py install


#### Generate yields:
cd ~/flexCE/flexCE/
python make_yield_grids.py


### Run the code
cd ~/flexCE/flexCE/
python flexce.py ../config/sim0.cfg


#### Change the parameters of the code
1. create a new sim<N>.cfg file (where <N> is a number)
2. edit sim\<N\>.cfg
3. re-run the code:
python flexce.py ../config/sim<N>.cfg


# Plot [O/Fe]-[Fe/H]
cd flexCE/flexCE/plot
python plot_xfe_feh.py ofe_feh_sim0.cfg

# Go to output plot directory
cd flexCE/plots/plots
