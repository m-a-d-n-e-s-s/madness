MADNESS Skyrme-Hartree-Fock code
--------------------------------

Branch: mshf\_stable: This branch uses the master branch at commit ba1218958388532cba5891b3b4fa5fd9235e928e


For general information, see:

Sagert, I.; Fann, G. I.; Fattoyev, F. J.; Postnikov, S.; Horowitz, C. J.</br>
Physical Review C, Volume 93, Issue 5, id.055801, 2016

Please also cite the paper if you download and use the code. 


To compile the code after building MADNESS:
```
make clean
make mshf
```

When running the code, make sure that the input files ```skyrme.input``` and ```mshf.input``` are in the same directory as the executable. Otherwise, mshf will run a default calculation of Oxygen-16. ```skyrme.input``` contains the parameters of the Skyrme nuclear force while ```mshf.input``` contains general input parameters, e.g. the size of the simulation volume and the number of neutrons and protons. For more information on the code, see the manual. 


The ```mshf``` folder contains the following files: 

- ```mshf.cc```: Main code file. Currently only it calls ```ground.cc``` 
- ```ground.cc, ground.h```: For ground state of nuclear matter, e.g. ground states of nuclei or nuclear pasta at zero temperature. 
- ```densities.cc, densities.h```: Nuclear matter density and all density terms to set up the Skyrme potentials and energies
- ```potentials.cc, potentials.h```: Skyrme potentials for Hartree-Fock iterations
- ```energies.cc, energies.h```: Skyrme energies for the Hartree-Fock iterations
- ```iterate.cc, iterate.h```: Main iteration routine
- ```orthonorm.cc, orthonorm.h```: Orthonormalization of the nucleon single particle states. Called in each iteration step
- ```states.cc, states.h```: Various manipulation routines for the single particle states, including initial setup of the wavefunctions, truncation and normalization. 
- ```loadbalance.cc, loadbalance.h```: Loadbalancing
- ```output.cc, output.h```: Output files
- ```input.h```: Manages input parameters. 
- ```skyrme.input``` and ```mshf.input```: Input files for running simulations
