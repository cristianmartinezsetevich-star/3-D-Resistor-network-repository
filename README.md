
This repository contains the Fortran source codes and input files used to simulate mixed ionic–electronic conducting (MIEC) electrodes voxelized into three-dimensional resistor networks. The numerical model resolves the coupled transport of electrons, oxygen ions, and gaseous species, incorporating the oxygen reduction reaction (ORR) through an interfacial resistance . These simulations were developed in the context of the manuscript “Methodology for the simulation of voxelized MIEC electrodes using resistor networks: coupling between electronic, ionic, and gaseous transport.” The provided scripts allow the reproduction of conductance–pO₂ curves and the exploration of the impact of microstructural parameters on the overall electrochemical performance. Simplified 3D Resistor Network Simulation Code. This repository contains a simplified version of the 3D resistor network model used to simulate Mixed Ionic–Electronic Conductor (MIEC) electrodes. It is intended as a practical example that preserves the essential physical coupling mechanisms of the full model while remaining easy to read and modify.

----------------------------------------------------------------------------------------------------------------------------------
Overview
-------

The code implements a voxel-based representation of the electrode microstructure, where each phase (gas, ionic conductor, and electronic conductor) is described by an independent resistive network in the three spatial directions.

It includes routines for:

Generating or importing a random or tomographic porous structure

Building resistor matrices for each transport network

Solving the potential and current distributions through iterative relaxation

Computing and saving total and partial currents (electronic, ionic, and gas)

Exporting data and potential maps for post-processing or visualization

Purpose

This version has been simplified for demonstration purposes. However, it retains the core logic of the coupled transport processes and can be used as a foundation for further development or adaptation.

Author

Cristian Martínez-Setevich
CITEDEF / CONICET — Argentina
2025

------------------------------------------------------------------------------------------------------------------------------------

Instructions for running the simulations:

The codes were developed and tested using the Force Fortran compiler on Windows, which is freely available for download from its official website. Any standard Fortran compiler (e.g., gfortran) can also be used to compile and execute the programs.

To run the simulations, create a working directory containing the required input data files:
ait1.dat, alado.dat, av.dat, avG.dat, and avIO.dat.
These files define the voxelized electrode structure and the corresponding physical parameters for the resistor-network simulation. The main program reads these files automatically during initialization.

Once the files are in place, compile the Fortran source code and execute the resulting binary from the same directory. Output files containing the potential fields, current distributions, and total conductance are generated upon completion.

------------------------------------------------------------------------------------------------------------------------------------

Description of input and output files:

The simulation requires several .dat files that define both the input parameters and the voxelized structure of the electrode:

alado.dat – Contains the dimensions of the 3D matrix corresponding to the synthetic tomography used to generate the voxelized electrode.

ait1.dat – Stores all relevant simulation parameters (e.g., conductivities, boundary conditions, iteration settings) along with the computed global results such as total conductance.

av.dat – Represents a 2D cross-section of the electronic potential field extracted from the 3D tomography.

avIO.dat – Represents a 2D cross-section of the ionic potential field.

avG.dat – Represents a 2D cross-section of the gaseous phase concentration field.

ac.dat – Corresponds to a 2D structural slice of the voxelized microstructure, useful for visualization or debugging.

These files are generated or updated within the subroutine guardsimple, which handles the reading and writing of the simulation data during each computational cycle.

------------------------------------------------------------------------------------------------------------------------------------

Code structure and main variables:

The main program, tresredes, initializes the data structures required to simulate a voxelized MIEC (Mixed Ionic–Electronic Conductor) electrode using three coupled resistor networks representing electronic, ionic, and gaseous transport.

The arrays declared at the beginning of the program define the 3D computational domain and the local transport properties in each direction (x, y, z):

r* arrays (e.g., rh, rv, rz, rih, riv, riz, rho, rvo, rzo) – Represent the resistive elements of the three interpenetrating networks (electronic, ionic, and gaseous) along the three spatial directions.

rcv, rco – Contain the local charge-transfer resistances coupling the electronic and ionic networks, corresponding to the oxygen reduction reaction (ORR).

i* arrays (ih, iv, iz, iih, iiv, iiz, ioh, iov, ioz) – Store the current (or flux) components along each direction for electrons (i), oxygen ions (ii), and gaseous species (iO).

v, vIO, vOO – Represent the electrochemical potentials of the electronic, ionic, and gaseous phases at each voxel, respectively.

c – 3D matrix defining the synthetic tomography of the electrode (voxel phase distribution).

per – Auxiliary matrix used for percolation analysis and connectivity verification of the conducting phases.

bo – Used to identify boundary voxels and apply boundary conditions during the iterative solution.

im – Indicates voxels corresponding to infiltrated regions or composite electrodes.

it, itt, ittt – Variables that store the total resulting current from the simulation, integrated over the entire electrode.

The header section also opens the main output file ait1.dat, which records the key simulation parameters and the resulting global quantities such as total conductance, porosity, and interface resistance.

The program allocates fixed-size 3D arrays of dimensions (200 × 200 × 200) for all field variables (resistances, potentials, and currents). This configuration provides a balance between spatial resolution and computational cost, and can be adjusted in the source code if a finer or coarser discretization is required.

To ensure numerical stability and accuracy in the solution of the coupled resistor networks, different levels of floating-point precision are employed:

Variables related to potential and current fields (v, vIO, vOO, ih, iv, iz, etc.) are defined as real*10, providing extended precision for iterative convergence and current balance.

Variables defining resistances, diffusion coefficients, and structural properties (rh, rv, rz, sigmae, sigmaO, etc.) are defined as real (single precision), as their variations are typically less sensitive to round-off errors.

The program uses a structured memory layout where each physical quantity is stored in a separate 3D array, allowing efficient indexing of neighboring voxels during the iterative computation of currents and potentials in the three transport networks. This organization simplifies the implementation of coupling terms and facilitates post-processing or visualization of any individual field.

-----------------------------------------------------------------------------------------------------------------------------------

Initialization of structural and simulation parameters:

The section labeled “Initial values” defines the structural parameters of the synthetic electrode and the iteration settings used in the relaxation method that solves the coupled resistor networks.

rmf, rmi, rmv – Define the range and current value of the pore diameter in voxels. The comment indicates that odd numbers yield better symmetry in the discretized structure.

rof, roi, rov – Control the porosity level or the number of pores in the generated microstructure.

ref, rei, rev – Define the diameter of infiltrated particles, used to represent composite or impregnated cathodes.

porf, pori, pov – Set the percentage of infiltration relative to the maximum possible impregnation.

b and a – Define the spatial dimensions of the working matrix; in this configuration, b = 40 determines the size of the cubic voxel domain, while a = 200 is used for indexing and boundary handling.

kk – Specifies the number of relaxation iterations performed in the solver loop.

Before starting the computation, the 3D matrices are initialized with uniform values:

c = 3 assigns an initial material label to all voxels.

bo = 0 and im = 0 reset auxiliary matrices related to porosity and impregnation.

The following scalar variables compute derived geometric or physical quantities:

ro = rov / rof → normalized porosity (gas fraction).

mr = rmv → current pore diameter.

rmelect = rev → size of the infiltrated particle.

porcen = pov / porf → fractional impregnation relative to the maximum.

These parameters define the initial state of the 3D synthetic microstructure prior to assigning resistances and solving the potential and current fields.

-----------------------------------------------------------------------------------------------------------------------------------

Physical parameters and unit definitions

This section defines the physical parameters governing electronic, ionic, and gaseous transport in the BSCF (Ba₀.₅Sr₀.₅Co₀.₈Fe₀.₂O₃₋δ) cathode, as well as the corresponding units and temperature dependencies used in the resistor-network model.

Temperature and oxygen partial pressure

The example simulation temperature is fixed at 750 °C (temp = 750).

The oxygen partial pressure is given by podos = 0.21, corresponding to air at 1 atm (can be varied by powers of 10 through diii).

The inverse temperature utemp = 1000 / (temp + 273) is used in Arrhenius-type relations.

Voxel geometry

lpix = 0.3 µm → voxel side length.

apix = lpix² → voxel cross-sectional area.

vpix = lpix × apix → voxel volume.

fara = 4 × 96485 C/mol, accounting for the four-electron transfer in the oxygen reduction reaction (ORR).

Electrical conductivities (BSCF)

Electronic conductivity:

𝜎𝑒=0.8500×10E(2(𝑇+273)/1275)×𝑝𝑂2E(0.25)

Ionic conductivity:

𝜎𝑂2−=0.0035×10E(2.5(𝑇+273)/1275)×𝑝𝑂2E(0.25)


These empirical relations reproduce the temperature and 𝑝𝑂2 dependence reported for BSCF at 750 °C.

Interfacial resistances (Ω·cm²)

resCT – Charge-transfer resistance associated with the ORR and activation energy (Baumann, 2006).

resinter – Interfacial resistance between the MIEC electrode and the electrolyte (Baumann, 2006).

rabsor – Resistance associated with surface adsorption and dissociative oxygen incorporation; values typically range from 0.01–0.1 Ω·cm².

Gas-phase parameters

Ccat – Gas concentration at the cathode surface, computed from the ideal gas law.

Di – Molecular diffusion coefficient in air, scaled with 𝑇E(3/2).

Dn – Knudsen diffusivity based on the voxel size and molecular mass of oxygen.

Def = 1 / (1/Di + 1/Dn) – Effective diffusivity combining molecular and Knudsen mechanisms.

Dimensionless coupling and transport coefficients

ctgm, ctmt, and ct – Represent the combined inverse resistances for charge transfer and adsorption processes, converted to consistent units using the Faraday constant and voxel area.

me and ls – Effective electronic and ionic conductance coefficients (1/Ω), scaled by voxel geometry and interfacial coupling terms.

ctml – Coupling term between the MIEC electrode and the electrolyte interface.

Diu – Dimensionless ionic diffusivity incorporating gas concentration (Ccat) and voxel geometry.

ga – Effective gaseous parameter, set equal to Diu.

finalunit – Total interfacial conductance considering both charge-transfer (resCT) and adsorption (rabsor) resistances in series.

The block concludes by printing the initialized values to the console, providing a summary of the main parameters (σ_e, σ_O, resCT, resinter, T, pO₂, ga, me, ls, etc.) before the iterative solver begins.

-----------------------------------------------------------------------------------------------------------------------------------

After computing the main physical quantities (conductivities, interfacial resistances, concentrations, and diffusivities), the code prints the initial parameters for verification. This diagnostic block ensures consistency of units and magnitudes before building the resistor networks.

rmv, rov, b, lpix*b, kk: mean resistance, porosity, voxel domain size, electrode length (µm), and simulation index.

sigmae, sigmaO: effective electronic and ionic conductivities of BSCF.

resCT, resinter: charge-transfer and MIEC–electrolyte interfacial resistances.

temp, podos: operating temperature and oxygen partial pressure.

ga, me, ls, ctml, ct, ctgm, ctmt: coupling and scaling coefficients for the electronic, ionic, and gas-phase sub-networks.

The write(*,*) commands generate a console output summarizing these values, allowing quick validation of input and scaling parameters before the simulation starts.

----------------------------------------------------------------------------------------------------------------------------------

Subroutine Calls: Network Generation

After initializing all physical parameters, the code calls a series of subroutines to construct the 3D resistor networks representing the mixed-conducting (MIEC) and gas phases:

porosidad(a,b,c,mr,ro) – Generates the porous structure based on the porosity (ro) and pore size (mr).

front(a,b,c,bo,im,frontera,impeso,rmelect) – Defines the boundary conditions and infiltration region for the electrode.

Resis(...) – Builds the resistor network for the electronic phase.

ResisIO(...) – Builds the resistor network for the ionic phase.

ResisOO(...) – Constructs the gas-phase diffusion network.

ResisCT(...)** and **ResisCO(...)` – Define the charge-transfer coupling resistances between electronic/ionic and gas/electrode interfaces.

Each subroutine fills the corresponding resistance matrices (e.g., rh, rv, rz for electrons, rih, riv, riz for ions), which are later used to solve the coupled transport and reaction processes in the MIEC electrode.

---------------------------------------------------------------------------------------------------------------------------------

Relaxation Method and Current Calculation

After constructing the resistance networks, the program computes the electrochemical potentials and current distributions using an iterative relaxation method (kir, kirB).
The corrientes subroutine calculates local current densities across all three sub-networks (electronic, ionic, and gaseous).
The total cell current and individual phase contributions are then obtained through itotal.

A convergence criterion is applied (if(abs(100*eit/it) >= 5) go to 100) to ensure current balance between electronic and ionic paths.
Finally, the results and key simulation parameters are saved with guardit, and optional routines can store 2D slices or current maps (guardsimple, guardcorrien).

The program concludes by closing all loops and recording the final timestamp.

----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE porosidad

This routine generates the 3D porous structure used as the synthetic tomography input for the simulation.
It randomly places spherical pores of radius mr/2 within a cubic domain until the target porosity ro is reached.

Periodic boundary conditions are applied on all faces to avoid edge artifacts during relaxation.
The electrolyte region (value = 2) is converted to the MIEC phase (value = 3) where appropriate.
The resulting 3D matrix c(i,j,k) defines the spatial distribution of solid and pore phases in the voxelized electrode.

----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE front

This subroutine identifies the interfacial (boundary) voxels between the MIEC and the pore phases in the 3D structure generated previously.
It scans all neighboring voxels with periodic boundary conditions to detect phase transitions.

For each boundary voxel, it:

Marks it in the auxiliary matrix bo,

Increments the total boundary count (frontera),

Optionally generates a spherical impregnation region around it in the matrix im, simulating infiltration of a secondary phase of radius rmelect.

Finally, it computes the maximum impregnation fraction (impeso = impre / (a*a*b)), which represents the relative volume of infiltrated regions.

----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE Resis, ResisIO, ResisOO, ResisCT, and ResisCO

These subroutines generate the 3D resistor networks corresponding to each transport path within the voxelized MIEC electrode structure.
Each network represents a specific physical process, assigning a local conductivity (or resistance) to every voxel connection based on the material phases stored in matrix c:

Resis → Electronic network (rh, rv, rz) for the metallic phase (σₑ = me).

ResisIO → Ionic network (rih, riv, riz) for transport through ionic and mixed regions (ls, ctml).

ResisOO → Gas-phase network (rho, rvo, rzo) for diffusive transport through pores (ga).

ResisCT → Charge-transfer network (rcv) coupling gas and metal interfaces via ctgm.

ResisCO → Charge-transfer network (rco) for metal–gas coupling governed by ct.

Each routine loops over all neighboring voxels (x, y, z directions) and fills the corresponding resistance matrices using phase-dependent coupling rules defined in con(u,v).

----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE kir, corrientes, corrientORR, and itotal

These routines perform the potential relaxation and current evaluation across the three coupled resistor networks (electronic, ionic, and gaseous).

kir → Solves the coupled potential fields (v, vIO, vOO) through iterative relaxation over kk cycles.
It includes dynamic adjustment of the charge-transfer resistances (rcv, rco) at mixed (MIEC) sites to ensure consistent current balance between gas–ion–electron pathways.
Periodic boundary conditions are applied laterally, while fixed potentials define the top and bottom electrode boundaries.

corrientes → Computes local current densities along the three spatial directions (x, y, z) for each transport network:

iv, ih, iz → electronic

iiv, iih, iiz → ionic

iov, ioh, ioz → gaseous

corrientORR → Evaluates the oxygen reduction reaction (ORR) current at each plane of the domain by integrating interfacial fluxes through rcv and rco.
The results (electronic and gaseous contributions) are written to the output file ORR.dat.

itotal → Integrates all current contributions to obtain the total current density and consistency check between electronic and gaseous fluxes at the current collector and electrolyte interface.

---------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE guardit

This routine saves the global simulation parameters and total current results to the file ait1.dat.
It stores a single line containing the most relevant output variables for post-processing or parametric analysis, such as:

Geometry and porosity parameters: a, b, lpix, podos, rmv, ro

Material and interface parameters: ga, me, ls, ctml, ctgm, ctmt, ct, sigmae, sigmao, resCT, resinter, rabsor, Ccat, Di

Temperature and total current results: temp, it, eit, itt, ittt

This file (ait1.dat) is appended after each simulation loop, building a compact summary dataset for later plotting or statistical evaluation.

SUBROUTINE guardcorrien

Writes the 3D current maps for each direction of the resistor network to separate files:

acorrien1.dat → current along the x-direction (iv)

acorrien2.dat → current along the y-direction (ih)

acorrien3.dat → current along the z-direction (iz)

Each file contains the current values for all planes (b × a × a), allowing visualization of current distribution within the MIEC electrode.

SUBROUTINE guardsimple

Exports 2D slices of the simulation domain and potential fields for visualization.
At a given cross-section (corte = 10), it writes:

ac.dat → material phase map (c)

av.dat → electronic potential (v)

avIO.dat → ionic potential (vIO)

avG.dat → gaseous potential (vOO)

The subroutine also handles value normalization (e.g., replacing boundary potentials of 1 or 0 with 0.5 for display purposes).
Additionally, the file alado.dat stores the grid dimensions (a, b) for reference.

----------------------------------------------------------------------------------------------------------------------------------
