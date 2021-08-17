# Black hole dynamics in globular clusters using NBODY6

This repository contains tools to compute the temporal evolution of black hole populations using the output of NBODY6 globular cluster simulations.
We use a suite of simulations on Swinburne's OzStar supercomputer, run by Ruggero de Vita for the following [publications](https://ui.adsabs.harvard.edu/search/filter_property_fq_property=AND&filter_property_fq_property=property%3A%22refereed%22&fq=%7B!type%3Daqp%20v%3D%24fq_property%7D&fq_property=(property%3A%22refereed%22)&p_=0&q=%20%20author%3A%22de%20Vita%22%20%20author%3A%22Trenti%22%20%20author%3A%22MacLeod%22&sort=date%20desc%2C%20bibcode%20desc).

## Project overview

Recently, [Antonini and Gieles 2020a](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2936A/abstract) created an efficient model for the evolution of the black hole population in globular clusters.
They determine parameters in a coupled system of equations for the time evolution of e.g. the total black hole mass.
Results for this "background" are consistent with NBODY6 across many orders of magnitude in e.g. density, star number etc. It can therefore be used in place of running NBODY simulations for a given cluster model, if we are not necessarily interested in individually resolving encounters.

A "foreground" model for the black hole population is run alongside the background model, which ultimately determines the merger rates in the cluster. This depends only on the bulk properties determined by the background model, and several simplifying assumptions about the black hole population.

A [Mathematica notebook](https://github.com/mcneilllucy/mcneilllucy.github.io/blob/master/assets/mathematica-notebooks/GC-background-foreground-evolution-AG2020.pdf) for this model is available for [download](https://github.com/mcneilllucy/mcneilllucy.github.io/blob/master/assets/mathematica-notebooks/GC-background-foreground-evolution-AG2020.nb).

This foreground model has been used to make inferences about merger rates in [Antonini and Gieles 2020b](https://ui.adsabs.harvard.edu/abs/2020PhRvD.102l3016A/abstract).

We will make use of the background model, but we are analysing NBODY6 simulations to create a different formulation of the black hole foreground, in the form of a simple monte carlo simulation that can be run in seconds. For example, it includes [Ginat and Perets 2020](https://ui.adsabs.harvard.edu/abs/2021PhRvX..11c1020G/citations) to determine outcomes of three body encounters.

The purpose is not to rival or make redundant cluster monte carlo or NBODY. Instead, we would like to do inference on cluster simulations generally and understand how different prescriptions for e.g. supernova kicks, inital mass functions impact black hole merger rates in dense environments.

## simulations (on OzStar)

number of stars range from 50,000-200,000 stars, with varying metallicity, physics prescriptions (e.g. neutron star kicks) etc.

## data (on OzStar)

Simulations are organised in the convention explained here.
For a single cluster simulation, there are n_output.dat (e.g. 0_output.dat ... 4_output.dat)
files which give a summary of all noteworthy events per timestep.
There is also a single snapdata.hdf5 file which contains e.g position and velocity vectors of every single star. This file is ~100GBs so these scripts must be run on OzStar.

## analysis-scripts (here, soon to be posted in package form)

contains:

-  Analysis for largest cluster (200k stars, canonical)
  - [x] ipython notebook
  - [x] .py script
    - [x] quick local density computations with kdTree, using nearest neighbour structures.
  - Results
    - [x] .csv files for important events e.g. black hole ejections, exchanges
    -  Movies for time evolution of:
    - [x] local density around each black hole
    - [x] mass density as a function of (log) radius
    - [x] number density as a function of (log) radius
    - [ ] mass distribution (histograms)
    - [ ] root-mean-square distance to centre of mass, as a function of mass

### local density around each black hole

We make use of kdTree to optimize local density calculations around each black hole.
Out of the 200,000 stars in the example simulation, at most ~100 of them are black holes.
We probe the mass dependence on local density for all the black holes over time, shown in the following .gif:

![local-dens-200k-can](/analysis-scripts/Example-with-largest-200k-cluster/results/local-densities-radius-log0.041no-self.gif)
