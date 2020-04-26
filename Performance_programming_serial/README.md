The aim of this piece of work  was to perform performance optimisation for a relatively
simple application code carrying out a molecular dynamics (MD).

# Contents of scripts

The programme reads an initial state from a file (i.e. input.dat) and then performs
5 blocks of 100 timesteps worth of simulation writing an output file after each block.
The code automatically reports timing information for each block of 100 timesteps
as well as the total time for running the entire simulation.

The MD.c file contains the computations carried out in order to simulate the collisions.
The util.c file contains the declaration of some helper functions that are used
in the MD.c file. The coord.h defines some pointers to some dynamically allocated
arrays contained in the control.c files. The main variables and arrays used in the
code and their associated physical concept are listed in Table 1. The control.c file
controls the update of the MD simulation and is responsible for reading the input
file and writing out the output to disk.

# Overview of code optimisations

Code optimisations included compiler flags optimisations, memory structures optimisation, 
data arrays alignment and hand optimisations.

# Conclusions
During this piece of work the performance of a simple application code used for carrying
out a MD simulation was improved such that the total application runtime was
decreased to 35.34 seconds, which is equivalent to a 99.32% decrease relative to the
original one. The optimised code performs approximately 25 times less calls to the
function required for calculating forces, carries out approximately 1000 times less
loop entries and also exhibits reduced data cache miss rates while still generating the
correct simulation output.