MADNESS Skyrme-Hartree-Fock code
--------------------------------

Branch: mshf\_stable: This branch uses the master branch at commit 
commit ba1218958388532cba5891b3b4fa5fd9235e928e


For general information, see:

Sagert, I.; Fann, G. I.; Fattoyev, F. J.; Postnikov, S.; Horowitz, C. J.</br>
Physical Review C, Volume 93, Issue 5, id.055801, 2016


To compile the code after building MADNESS:
```
make clean
make mshf
```

When running the code, make sure that the input files ```mshf_skyrme``` and ```mshf_input``` are in the same directory as the executable. Otherwise, mshf will run a default calculation of Oxygen-16. ```mshf_skyrme``` contains the parameters of the Skyrme nuclear force while ```mshf_input``` contains general input parameters, e.g. the size of the simulation volume and the number of neutrons and protons. For more information on the code, see the manual in the ```mshf_manual``` folder. 


The ```mshf``` folder contains the following files
- ```input.h```: Manages input parameters. 
- ```mshf.cc```: Main code file. Currently only it calls ```ground.cc``` 
- ```ground.cc```: Calculates the ground state of nuclear matter, e.g. ground states of nuclei or nuclear pasta at zero temperature. 
- ```densities.cc```: Calculates the nuclear matter density and all density terms that are required to set up the Skyrme potentials and energies
- ```potentials.cc```: Sets up Skyrme potentials for Hartree-Fock iterations
- ```energies.cc```: Sets up Skyrme energies for the Hartree-Fock iterations
- ```iterate.cc```: Main iteration routine
- ```orthonorm.cc```: Orthonormalize the nucleon single particle states. Is called during each iteration 
- ```states.cc```: Provides various manipulation routines for the single particle states, including initial setup of the wavefunctions, truncation and normalization. 
- ```loadbalance.cc```: Routines for loadbalancing
- ```output.cc```: Creates output files
- ```mshf_manual```: Manual for mshf code
- ```mshf_skyrme``` and ```mshf_input```: Input files for running simulations
