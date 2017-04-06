#include <chem/vibanal.h>
#include <madness/tensor/tensor_lapack.h>

using namespace madness;

/// returns the vibrational frequencies

/// @param[in]  hessian the hessian matrix (not mass-weighted)
/// @param[out] normalmodes the normal modes
/// @param[in]  project_tr whether to project out translation and rotation
/// @param[in]  print_hessian   whether to print the hessian matrix
/// @return the frequencies in atomic units
Tensor<double> compute_frequencies(const Molecule& molecule,
                                   const Tensor<double>& hessian, Tensor<double>& normalmodes,
                                   const bool project_tr, const bool print_hessian) {
    
    // compute mass-weighing matrices
    Tensor<double> M=molecule.massweights();
    Tensor<double> Minv(3*molecule.natom(),3*molecule.natom());
    for (int i=0; i<3*molecule.natom(); ++i) Minv(i,i)=1.0/M(i,i);
    
    // mass-weight the hessian
    Tensor<double> mwhessian=inner(M,inner(hessian,M));
    
    // remove translation and rotation
    if (project_tr) remove_external_dof(mwhessian,molecule);
    
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
    Tensor<double> D=projector_external_dof(molecule);
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


/// compute the projector to remove transl. and rot. degrees of freedom

/// taken from http://www.gaussian.com/g_whitepap/vib.htm
/// I don't really understand the concept behind the projectors, but it
/// seems to work, and it is not written down explicitly anywhere.
/// NOTE THE ERROR IN THE FORMULAS ON THE WEBPAGE !
Tensor<double> projector_external_dof(const Molecule& mol) {
    
    // compute the translation vectors
    Tensor<double> transx(3*mol.natom());
    Tensor<double> transy(3*mol.natom());
    Tensor<double> transz(3*mol.natom());
    for (int i=0; i<mol.natom(); ++i) {
        transx[3*i  ]=sqrt(mol.get_atom(i).get_mass_in_au());
        transy[3*i+1]=sqrt(mol.get_atom(i).get_mass_in_au());
        transz[3*i+2]=sqrt(mol.get_atom(i).get_mass_in_au());
    }
    
    // compute the rotation vectors
    
    // move the molecule to its center of mass and compute
    // the moment of inertia tensor
    Tensor<double> com=mol.center_of_mass();
    Molecule mol2=mol;
    mol2.translate(-1.0*com);
    Tensor<double> I=mol2.moment_of_inertia();
    I.scale(constants::atomic_mass_in_au);
    
    // diagonalize the moment of inertia
    Tensor<double> v,e;
    syev(I, v, e);  // v being the "X" tensor on the web site
    v=transpose(v);
    
    //        Tensor<double> B(e.size());
    //        for (long i=0; i<e.size(); ++i) B(i)=1.0/(2.0*e(i));
    //        print("rotational constants in cm-1");
    //        print(constants::au2invcm*B);
    
    // rotation vectors
    Tensor<double> rotx(3*mol.natom());
    Tensor<double> roty(3*mol.natom());
    Tensor<double> rotz(3*mol.natom());
    
    for (int iatom=0; iatom<mol.natom(); ++iatom) {
        
        // coordinates wrt the center of mass ("R" on the web site)
        Tensor<double> coord(3);
        coord(0l)=mol.get_atom(iatom).x-com(0l);
        coord(1l)=mol.get_atom(iatom).y-com(1l);
        coord(2l)=mol.get_atom(iatom).z-com(2l);
        
        // note the wrong formula on the Gaussian website:
        // multiply with sqrt(mass), do not divide!
        coord.scale(sqrt(mol.get_atom(iatom).get_mass_in_au()));
        
        // p is the dot product of R and X on the web site
        Tensor<double> p=inner(coord,v);
        
        // Eq. (5)
        rotx(3*iatom + 0)=p(1)*v(0,2)-p(2)*v(0,1);
        rotx(3*iatom + 1)=p(1)*v(1,2)-p(2)*v(1,1);
        rotx(3*iatom + 2)=p(1)*v(2,2)-p(2)*v(2,1);
        
        roty(3*iatom + 0)=p(2)*v(0,0)-p(0l)*v(0,2);
        roty(3*iatom + 1)=p(2)*v(1,0)-p(0l)*v(1,2);
        roty(3*iatom + 2)=p(2)*v(2,0)-p(0l)*v(2,2);
        
        rotz(3*iatom + 0)=p(0l)*v(0,1)-p(1)*v(0,0);
        rotz(3*iatom + 1)=p(0l)*v(1,1)-p(1)*v(1,0);
        rotz(3*iatom + 2)=p(0l)*v(2,1)-p(1)*v(2,0);
        
    }
    
    // move the translational and rotational vectors to a common tensor
    Tensor<double> ext_dof(6,3*mol.natom());
    ext_dof(0l,_)=transx;
    ext_dof(1l,_)=transy;
    ext_dof(2l,_)=transz;
    ext_dof(3l,_)=rotx;
    ext_dof(4l,_)=roty;
    ext_dof(5l,_)=rotz;
    
    // normalize
    for (int i=0; i<6; ++i) {
        double norm=ext_dof(i,_).normf();
        if (norm>1.e-14) ext_dof(i,_).scale(1.0/norm);
        else ext_dof(i,_)=0.0;
    }
    
    // compute overlap to orthonormalize the projectors
    Tensor<double> ovlp=inner(ext_dof,ext_dof,1,1);
    syev(ovlp,v,e);
    ext_dof=inner(v,ext_dof,0,0);
    
    // normalize or remove the dof if necessary (e.g. linear molecules)
    for (int i=0; i<6; ++i) {
        if (e(i)<1.e-14) {
            ext_dof(i,_).scale(0.0);      // take out this degree of freedom
        } else {
            ext_dof(i,_).scale(1.0/sqrt(e(i)));   // normalize
        }
    }
    
    // construct projector on the complement of the rotations
    Tensor<double> projector(3*mol.natom(),3*mol.natom());
    for (int i=0; i<3*mol.natom(); ++i) projector(i,i)=1.0;
    
    // compute the outer products of the projectors
    // 1- \sum_i | t_i >< t_i |
    projector-=inner(ext_dof,ext_dof,0,0);
    
    // construct random tensor for orthogonalization
    Tensor<double> D(3*mol.natom(),3*mol.natom());
    D.fillrandom();
    for (int i=0; i<6; ++i) D(i,_)=ext_dof(i,_);
    
    // this works only for non-linear molecules -- projector seems simpler
    //        ovlp=inner(D,D,1,1);
    //        cholesky(ovlp); // destroys ovlp
    //        Tensor<double> L=transpose(ovlp);
    //        Tensor<double> Linv=inverse(L);
    //        D=inner(Linv,D,1,0);
    //        ovlp=inner(D,D,1,1);
    //
    //        for (int i=0; i<6; ++i) D(i,_)=0.0;
    //        D=copy(D(Slice(6,8),_));
    
    //        return transpose(D);
    
    return projector;
    
}

/// remove translational degrees of freedom from the hessian
void remove_external_dof(Tensor<double>& hessian,
                                const Molecule& mol) {
    
    // compute the translation of the center of mass
    Tensor<double> projector_ext=projector_external_dof(mol);
    
    // this is P^T * H * P
    hessian=inner(projector_ext,inner(hessian,projector_ext),0,0);
}
