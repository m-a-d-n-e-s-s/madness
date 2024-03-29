# QCSchema Json

Integration of Madness into QCArchive requires to ability to translate data
generated by Madness into a standard quantum chemistry variable defined in 
`QCSchema`.  Within QCEngine, the madness `harvester.py` code is responible 
for reading in madness output files and harvesting relevant output data.  
This is done using the `regex` python api which pattern matches.  

`json.hpp` is included to allow applications to write json files with
relevant Quantum Chemistry data.  The ultimate goal is to be able
to generate data requested by QCArchive according to the `QCSchema` as
explained (QCSchema)['https://molssi-qc-schema.readthedocs.io/en/latest/auto_props.html'].

This ways we generate `json` output files that can easily be read within 
QCA.  


## Calculation Information 

`calc_info.json`

A list of fields that involve basic information of the requested computation.

- `calcinfo_nbasis`	The number of basis functions for the computation.	number
    - Not relevant in MADNESS
- `calcinfo_nmo`	The number of molecular orbitals for the computation.	number
    - MOLDFT nmo_alpha + nmo_beta ... derived
- `calcinfo_nalpha`	The number of alpha electrons in the computation.	number
    - nalpha
- `calcinfo_nbeta`	The number of beta electrons in the computation.	number
    - nbeta
- `calcinfo_natom`	The number of atoms in the computation.	number
    - molecule .size
- `return_energy`	The energy of the requested method, identical to return_value for energy computations.	number

## Self-Consistent Field Information

A list of fields added at the SCF level.  (HF and DFT)

`scf_info.json`

- scf_one_electron_energy	The one-electron (core Hamiltonian) energy contribution to the total SCF energy.	number
- scf_two_electron_energy	The two-electron energy contribution to the total SCF energy.	number
- nuclear_repulsion_energy	The nuclear repulsion energy contribution to the total SCF energy.	number
- scf_vv10_energy	The VV10 functional energy contribution to the total SCF energy.	number
- scf_xc_energy	The functional energy contribution to the total SCF energy.	number
- scf_dispersion_correction_energy	The dispersion correction appended to an underlying functional when a DFT-D method is requested.	number
- scf_dipole_moment	The X, Y, and Z dipole components.	array[number]
- scf_total_energy	The total electronic energy of the SCF stage of the calculation. This is represented as the sum of the … quantities.	number
- scf_iterations

## Wavefunction Schema

A list of valid quantum chemistry wavefunction properties tracked byu the schema.  
Matrices are in column-major order. AO basis functions are ordered according to the CCA standard as implemented in libint.




