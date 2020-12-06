# Black hole dynamics for globular clusters using NBODY6

This repository contains tools to compute the temporal evolution of black hole populations using the output of NBODY6 globular cluster simulations.
We use a suite of simulations on Swinburne's OzStar supercomputer, run by Ruggero de Vita for the following [publications](https://ui.adsabs.harvard.edu/search/filter_author_facet_hier_fq_author=AND&filter_author_facet_hier_fq_author=author_facet_hier%3A%220%2Fde%20Vita%2C%20R%22&filter_author_facet_hier_fq_author=author_facet_hier%3A%220%2FMacLeod%2C%20M%22&filter_property_fq_property=AND&filter_property_fq_property=property%3A%22refereed%22&fq=%7B!type%3Daqp%20v%3D%24fq_author%7D&fq=%7B!type%3Daqp%20v%3D%24fq_property%7D&fq_author=(author_facet_hier%3A%220%2Fde%20Vita%2C%20R%22%20AND%20author_facet_hier%3A%220%2FMacLeod%2C%20M%22)&fq_property=(property%3A%22refereed%22)&p_=0&q=%20%20author%3A%22de%20Vita%22%20%20author%3A%22Trenti%22%20%20author%3A%22MacLeod%22&sort=date%20desc%2C%20bibcode%20desc).

## simulations

number of stars range from 50,000-200,000 stars, with varying metallicity, physics prescriptions (e.g. neutron star kicks) etc.

## data

Simulations are organised in the convention explained here.
For a single cluster simulation, there are n_output.dat (e.g. 0_output.dat ... 4_output.dat)
files which give a summary of all noteworthy events per timestep.
There is also a single snapdata.hdf5 file which contains e.g position and velocity vectors of every single star. This file is ~100GBs so these scripts must be run on OzStar.

## analysis-scripts

contains:

-  Analysis for largest cluster (200k stars, canonical)
  - [x] ipython notebook
  - [ ] .py script
    - [x] quick local density computations with kdTree, using nearest neighbour structures.
  - Results
    - [x] .csv files for important events e.g. black hole ejections, exchanges
    -  Movies for time evolution of:
    - [ ] local density around each black hole
    - [x] mass density as a function of (log) radius
    - [x] number density as a function of (log) radius
    - [ ] mass distribution (histograms)
    - [ ] root-mean-square distance to centre of mass, as a function of mass
