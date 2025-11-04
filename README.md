This repository contains the Fortran source codes and input files used to simulate mixed ionic–electronic conducting (MIEC) electrodes voxelized into three-dimensional resistor networks. The numerical model resolves the coupled transport of electrons, oxygen ions, and gaseous species, incorporating the oxygen reduction reaction (ORR) through an interfacial resistance formalism. These simulations were developed in the context of the manuscript “Methodology for the simulation of voxelized MIEC electrodes using resistor networks: coupling between electronic, ionic, and gaseous transport.” The provided scripts allow the reproduction of conductance–pO₂ curves and the exploration of the impact of microstructural parameters on the overall electrochemical performance.

Instructions for running the simulations:

The codes were developed and tested using the Force Fortran compiler on Windows, which is freely available for download from its official website. Any standard Fortran compiler (e.g., gfortran) can also be used to compile and execute the programs.

To run the simulations, create a working directory containing the required input data files:
ait1.dat, alado.dat, av.dat, avG.dat, and avIO.dat.
These files define the voxelized electrode structure and the corresponding physical parameters for the resistor-network simulation. The main program reads these files automatically during initialization.

Once the files are in place, compile the Fortran source code and execute the resulting binary from the same directory. Output files containing the potential fields, current distributions, and total conductance are generated upon completion.

Description of input and output files:

The simulation requires several .dat files that define both the input parameters and the voxelized structure of the electrode:

alado.dat – Contains the dimensions of the 3D matrix corresponding to the synthetic tomography used to generate the voxelized electrode.

ait1.dat – Stores all relevant simulation parameters (e.g., conductivities, boundary conditions, iteration settings) along with the computed global results such as total conductance.

av.dat – Represents a 2D cross-section of the electronic potential field extracted from the 3D tomography.

avIO.dat – Represents a 2D cross-section of the ionic potential field.

avG.dat – Represents a 2D cross-section of the gaseous phase concentration field.

ac.dat – Corresponds to a 2D structural slice of the voxelized microstructure, useful for visualization or debugging.

These files are generated or updated within the subroutine guardsimple, which handles the reading and writing of the simulation data during each computational cycle.

