#ifndef MADNESS_VIBANAL_INCLUDED
#define MADNESS_VIBANAL_INCLUDED

#include <madness/tensor/tensor.h>
#include <chem/molecule.h>

madness::Tensor<double> compute_frequencies(const madness::Molecule& molecule,
                                            const madness::Tensor<double>& hessian, madness::Tensor<double>& normalmodes,
                                            const bool project_tr=true, const bool print_hessian=false);


madness::Tensor<double> compute_reduced_mass(const madness::Molecule& molecule, const madness::Tensor<double>& normalmodes);


#endif //  MADNESS_VIBANAL_INCLUDED



 
