# Ferret Brain Morphogenesis

The repository contains numerical simulation codes and data for our paper:

G. P. T. Choi, C. Liu, S. Yin, G. Séjourné, R. S. Smith, C. A. Walsh, L. Mahadevan,
"[Biophysical basis for brain folding and misfolding patterns in ferrets and humans.](https://doi.org/10.1101/2025.03.05.641682)"
Preprint, bioRxiv 2025.03.05.641682. 

===============================================================

Numerical simulations:
* All codes are in the folder `code_simulation`
* `Brains_simulation.cpp`: Simulating ferret brain growth
* `Brains_simulation_modified_thickness.cpp`: Simulating ferret brain growth with modified cortical layer thickness
* `Brains_simulation_modified_growth_rate.cpp`: Simulating ferret brain growth with modified growth rate
* `Brains_simulation_modified_local.cpp`: Simulating ferret brain growth with modified thickness/growth rate at a localized region

Steps for running the numerical simulations:
* Prepare a tetrahedral mesh data file (with `.mesh` format) and change the file name in the cpp code
* Run the following command for Brains_simulation (similar for the other codes):
  
`cd your_code_directory`

`g++ -Ofast -fopenmp eig3.cpp Brains_simulation.cpp -o Brains_simulation`

`./Brains_simulation`

* The outputs will be the simulation results at different time points, e.g. `B0.off`, `B500.off` 
  
=============================================================

Morphometric analyses:
* All codes are in the folder `code_morphometrics`
* Landmark-aligned spherical parameterization
* Shape index quantification
* Spherical harmonic analysis

Steps for running the morphometric analyses:
* Use `read_off.m` to read the simulated .off files
* Label the sulcal landmark positions on the meshes
* Run `demo_morphometrics.m` with the required vertex, face, and landmark information of two meshes to be compared
* Obtain the parameterization results, shape indices and the spherical harmonic representations
