


/// returns the vibrational frequencies

/// @param[in]  hessian the hessian matrix (not mass-weighted)
/// @param[out] normalmodes the normal modes
/// @param[in]  project_tr whether to project out translation and rotation
/// @param[in]  print_hessian   whether to print the hessian matrix
/// @return the frequencies in atomic units
Tensor<double> compute_frequencies(const Molecule& molecule,
                                   const Tensor<double>& hessian, Tensor<double>& normalmodes,
                                   const bool project_tr=true, const bool print_hessian=false) {
    
    // compute mass-weighing matrices
    Tensor<double> M=molecule.massweights();
    Tensor<double> Minv(3*molecule.natom(),3*molecule.natom());
    for (int i=0; i<3*molecule.natom(); ++i) Minv(i,i)=1.0/M(i,i);
    
    // mass-weight the hessian
    Tensor<double> mwhessian=inner(M,inner(hessian,M));
    
    // remove translation and rotation
    if (project_tr) MolecularOptimizer::remove_external_dof(mwhessian,molecule);
    
    if (print_hessian) {
        if (project_tr) {
            print("mass-weighted hessian with translation and rotation projected out");
        } else {
            print("mass-weighted unprojected hessian");
        }
        Tensor<double> mmhessian=inner(Minv,inner(mwhessian,Minv));
        print(mwhessian);
        print("mass-weighted unprojected hessian; mass-weighing undone");
        print(mmhessian);
    }
    
    Tensor<double> freq;
    syev(mwhessian,normalmodes,freq);
    for (long i=0; i<freq.size(); ++i) {
        if (freq(i)>0.0) freq(i)=sqrt(freq(i)); // real frequencies
        else freq(i)=-sqrt(-freq(i));           // imaginary frequencies
    }
    return freq;
}


Tensor<double> compute_reduced_mass(const Molecule& molecule, const Tensor<double>& normalmodes) {
    
    Tensor<double> M=molecule.massweights();
    Tensor<double> D=MolecularOptimizer::projector_external_dof(molecule);
    Tensor<double> L=copy(normalmodes);
    Tensor<double> DL=inner(D,L);
    Tensor<double> MDL=inner(M,DL);
    Tensor<double> mu(3*molecule.natom());
    
    for (int i=0; i<3*molecule.natom(); ++i) {
        double mu1=0.0;
        for (int j=0; j<3*molecule.natom(); ++j) mu1+=MDL(j,i)*MDL(j,i);
        if (mu1>1.e-14) mu(i)=1.0/(mu1*constants::atomic_mass_in_au);
    }
    return mu;
}
