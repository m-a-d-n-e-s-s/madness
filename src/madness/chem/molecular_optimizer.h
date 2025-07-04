/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

*/


/// \file chem/molecular_optimizer.h
/// \brief optimize the geometrical structure of a molecule
#ifndef MADNESS_CHEM_MOLECULAR_OPTIMIZER_H__INCLUDED
#define MADNESS_CHEM_MOLECULAR_OPTIMIZER_H__INCLUDED

#include <madness/tensor/solvers.h>
#include<madness/chem/molecule.h>
#include"madness/mra/QCCalculationParametersBase.h"

namespace madness {


struct MolecularOptimizationTargetInterface : public OptimizationTargetInterface {

    /// return the molecule of the target
    virtual Molecule& molecule() {
        MADNESS_EXCEPTION("you need to return a molecule",1);
        return *(new Molecule());   // this is a memory leak silencing warnings
    }
};


class MolecularOptimizationParameters : public QCCalculationParametersBase {
public:

	/// ctor reading out the input file
    MolecularOptimizationParameters() {
        initialize < std::string > ("update", "bfgs", "Quasi-Newton update method",{"bfgs","sr1"});
        initialize < int > ("maxiter", 20);
        initialize < double > ("tol", 1.e-4,"geometry convergence threshold for the gradient norm");
        initialize < double > ("value_precision", 1.e-12);
        initialize < double > ("gradient_precision", 1.e-12);
        initialize < bool > ("printtest", false);
        initialize < std::string > ("cg_method", "polak_ribiere","conjugate gradients update",{"polak_ribiere","fletcher-reeves"});
        initialize < std::vector<std::string> >
        ("remove_dof", {"tx", "ty", "tz", "rx", "ry", "rz"}, "degree of freedom projected out: translation/rotation");
    }

    MolecularOptimizationParameters(World& world, const commandlineparser& parser)
            : MolecularOptimizationParameters() {
		// read input file
        read_input_and_commandline_options(world,parser,"geoopt");
        set_derived_values();
	}
    std::string get_tag() const override {
        return std::string("geoopt");
    }

	void set_derived_values() {
	}

	std::string update() const {return get<std::string>("update");}
	std::string cg_method() const {return get<std::string>("cg_method");}
	int maxiter() const {return get<int>("maxiter");}
	double tol() const {return get<double>("tol");}
	double value_precision() const {return get<double>("value_precision");}
	double gradient_precision() const {return get<double>("gradient_precision");}
	bool printtest() const {return get<bool>("printtest");}
	std::vector<std::string> remove_dof() const {return get<std::vector<std::string> >("remove_dof");}

};




/// Molecular optimizer derived from the QuasiNewton optimizer

/// Essentially the QuasiNewton optimizer, but with the additional feature
/// of projecting out rotational and translational degrees of freedom
class MolecularOptimizer : public OptimizerInterface {
private:
    /// How to update the hessian: BFGS or SR1
    std::shared_ptr<MolecularOptimizationTargetInterface> target;
    Tensor<double> h;
    double f=1.e10;
    double gnorm=1.e10;

public:
    MolecularOptimizationParameters parameters;

    /// same ctor as the QuasiNewton optimizer
    MolecularOptimizer(World& world, const commandlineparser& parser,
    		const std::shared_ptr<MolecularOptimizationTargetInterface>& tar)
			: target(tar), parameters(world, parser) {
    }

    /// optimize the underlying molecule

    /// @param[in]  x   the coordinates to compute energy and gradient
    bool optimize(Tensor<double>& x) {
        bool converge;
        converge=optimize_quasi_newton(x);
//        converge=optimize_conjugate_gradients(x);
        return converge;
    }

    bool converged() const {
        return gradient_norm();
    }

    bool converged(const Tensor<double>& displacement) const {
        return (gradient_norm()<parameters.tol() and (displacement.normf()/displacement.size())<parameters.tol());
    }

    double value() const {return 0.0;}

    double gradient_norm() const {return gnorm;}

    /// set an (initial) hessian
    void set_hessian(const Tensor<double>& hess) {
        h=copy(hess);
    }


private:

    bool optimize_quasi_newton(Tensor<double>& x) {

        if(parameters.printtest())  target->test_gradient(x, parameters.value_precision());

        bool h_is_identity = (h.size() == 0);
        if (h_is_identity) {
            int n = x.dim(0);
            h = Tensor<double>(n,n);

            // mass-weight the initial hessian
            for (size_t i=0; i<target->molecule().natom(); ++i) {
                h(3*i  ,3*i  )=1.0/(target->molecule().get_atom(i).mass);
                h(3*i+1,3*i+1)=1.0/(target->molecule().get_atom(i).mass);
                h(3*i+2,3*i+2)=1.0/(target->molecule().get_atom(i).mass);
            }
            h*=10.0;
            madness::print("using the identity as initial Hessian");
        } else {
            Tensor<double> normalmodes;
            Tensor<double> freq=MolecularOptimizer::compute_frequencies(
                    target->molecule(),h,normalmodes,parameters.remove_dof());
            madness::print("\ngopt: projected vibrational frequencies (cm-1)\n");
            printf("frequency in cm-1   ");
            for (int i=0; i<freq.size(); ++i) {
                printf("%10.1f",constants::au2invcm*freq(i));
            }

        }

        remove_external_dof(h,target->molecule(),parameters.remove_dof());

        // the previous gradient
        Tensor<double> gp;
        // the displacement
        Tensor<double> dx;

        for (int iter=0; iter<parameters.maxiter(); ++iter) {
            Tensor<double> gradient;

            target->value_and_gradient(x, f, gradient);
//            print("gopt: new energy",f);
//            const double rawgnorm = gradient.normf()/sqrt(gradient.size());
//            print("gopt: raw gradient norm ",rawgnorm);

            // remove external degrees of freedom (translation and rotation)
            Tensor<double> project_ext=projector_external_dof(target->molecule(),parameters.remove_dof());
            gradient=inner(gradient,project_ext);
            gnorm = gradient.normf()/sqrt(gradient.size());
//            print("gopt: projected gradient norm ",gnorm);
//            const double gradratio=rawgnorm/gnorm;


//            if (iter == 1 && h_is_identity) {
//                // Default initial Hessian is scaled identity but
//                // prefer to reuse any existing approximation.
//                h.scale(gradient.trace(gp)/gp.trace(dx));
//            }

            if (iter > 0) {
                if (parameters.update() == "bfgs") QuasiNewton::hessian_update_bfgs(dx, gradient-gp,h);
                else QuasiNewton::hessian_update_sr1(dx, gradient-gp,h);
            }


            remove_external_dof(h,target->molecule(),parameters.remove_dof());
            Tensor<double> v, e;
            syev(h, v, e);
            Tensor<double> normalmodes;
            Tensor<double> freq=MolecularOptimizer::compute_frequencies(
                    target->molecule(),h,normalmodes,parameters.remove_dof());
            madness::print("\ngopt: projected vibrational frequencies (cm-1)\n");
            printf("frequency in cm-1   ");
            for (int i=0; i<freq.size(); ++i) {
                printf("%10.3f",constants::au2invcm*freq(i));
            }
            printf("\n");

            // this will invert the hessian, multiply with the gradient and
            // return the displacements
            dx = new_search_direction2(gradient,h);

            // line search only for large displacements > 0.01 bohr = 2pm
            double step=1.0;
            double maxdx=dx.absmax();
            if (h_is_identity and (maxdx>0.01)) step = QuasiNewton::line_search(1.0, f,
                    dx.trace(gradient), x, dx,target, parameters.value_precision());

            dx.scale(step);
            x += dx;
            gp = gradient;

            double disp2=dx.normf()/dx.size();
            printf(" QuasiNewton iteration %2d energy: %16.8f gradient %.2e displacement %.2e \n", iter,f,gnorm, disp2);
            if (converged(dx)) break;
        }

        if (parameters.printtest()) {
            print("final hessian");
            print(h);
        }
        return converged(dx);
    }

    bool optimize_conjugate_gradients(Tensor<double>& x) {

//        Tensor<double> project_ext=projector_external_dof(target->molecule());


        // initial energy and gradient gradient
        double energy=0.0;
        Tensor<double> gradient;

        // first step is steepest descent
        Tensor<double> displacement(x.size());
        Tensor<double> oldgradient;
        Tensor<double> old_displacement(x.size());
        old_displacement.fill(0.0);

        for (int iter=1; iter<parameters.maxiter(); ++iter) {

            // displace coordinates
            if (iter>1) x+=displacement;

            // compute energy and gradient
            target->value_and_gradient(x, energy, gradient);
//            print("gopt: new energy",energy);
            gnorm = gradient.normf()/sqrt(gradient.size());
//            print("gopt: raw gradient norm ",gnorm);

            // remove external degrees of freedom (translation and rotation)
            Tensor<double> project_ext=projector_external_dof(target->molecule(),parameters.remove_dof());
            gradient=inner(gradient,project_ext);
            gnorm = gradient.normf()/sqrt(gradient.size());
//            print("gopt: projected gradient norm ",gnorm);

            // compute new displacement
            if (iter==1) {
                displacement=-gradient;
            } else {
                double beta=0.0;
                if (parameters.cg_method()=="fletcher_reeves")
                    beta=gradient.normf()/oldgradient.normf();
                if (parameters.cg_method()=="polak_ribiere")
                    beta=gradient.normf()/(gradient-oldgradient).normf();
                displacement=-gradient + beta * old_displacement;
            }

            // save displacement for the next step
            old_displacement=displacement;

            if (converged(displacement)) break;
        }

        return converged(displacement);
    }

    /// effectively invert the hessian and multiply with the gradient
    Tensor<double> new_search_direction2(const Tensor<double>& g,
            const Tensor<double>& hessian) const {
        Tensor<double> dx, s;
        double tol = parameters.gradient_precision();
        double trust = 1.0; // This applied in spectral basis

        // diagonalize the hessian:
        // VT H V = lambda
        // H^-1   = V lambda^-1 VT
        Tensor<double> v, e;
        syev(hessian, v, e);

        // Transform gradient into spectral basis
        // H^-1 g = V lambda^-1 VT g
        Tensor<double> gv = inner(g,v); // this is VT g == gT V == gv

        // Take step applying restriction
        int nneg=0, nsmall=0, nrestrict=0;
        for (int i=0; i<e.size(); ++i) {
            if (e[i] < -tol) {
                if (parameters.printtest()) printf(
                        "   forcing negative eigenvalue to be positive %d %.1e\n", i, e[i]);
                nneg++;
                //e[i] = -2.0*e[i]; // Enforce positive search direction
                e[i] = -0.1*e[i]; // Enforce positive search direction
            }
            else if (e[i] < tol) {
                if (parameters.printtest()) printf(
                        "   forcing small eigenvalue to be zero %d %.1e\n", i, e[i]);
                nsmall++;
                e[i] = tol;
                gv[i]=0.0;   // effectively removing this direction
            }

            // this is the step -lambda^-1 gv
            gv[i] = -gv[i] / e[i];
            if (std::abs(gv[i]) > trust) { // Step restriction
                double gvnew = trust*std::abs(gv(i))/gv[i];
                if (parameters.printtest()) printf(
                        "   restricting step in spectral direction %d %.1e --> %.1e\n",
                        i, gv[i], gvnew);
                nrestrict++;
                gv[i] = gvnew;
            }
        }
        if (nneg || nsmall || nrestrict)
            printf("   nneg=%d nsmall=%d nrestrict=%d\n", nneg, nsmall, nrestrict);

        // Transform back from spectral basis to give the displacements
        // disp = -V lambda^-1 VT g = V lambda^-1 gv
        return inner(v,gv);
    }

public:

    /// compute the projector to remove transl. and rot. degrees of freedom

    /// taken from http://www.gaussian.com/g_whitepap/vib.htm
    /// I don't really understand the concept behind the projectors, but it
    /// seems to work, and it is not written down explicitly anywhere.
    /// NOTE THE ERROR IN THE FORMULAS ON THE WEBPAGE !
    /// @param[in]	do_remove_dof	which dof to remove: x,y,z,Rx,Ry,Rz (transl/rot)
    static Tensor<double> projector_external_dof(const Molecule& mol,
    		const std::vector<std::string>& remove_dof) {

        // compute the translation vectors
        Tensor<double> transx(3*mol.natom());
        Tensor<double> transy(3*mol.natom());
        Tensor<double> transz(3*mol.natom());
        for (size_t i=0; i<mol.natom(); ++i) {
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

        for (size_t iatom=0; iatom<mol.natom(); ++iatom) {

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
        auto tmp=remove_dof;
        for (auto& t : tmp) t=commandlineparser::tolower(t);
        bool remove_Tx=std::find(tmp.begin(), tmp.end(), "tx")!=tmp.end();
        bool remove_Ty=std::find(tmp.begin(), tmp.end(), "ty")!=tmp.end();
        bool remove_Tz=std::find(tmp.begin(), tmp.end(), "tz")!=tmp.end();
        bool remove_Rx=std::find(tmp.begin(), tmp.end(), "rx")!=tmp.end();
        bool remove_Ry=std::find(tmp.begin(), tmp.end(), "ry")!=tmp.end();
        bool remove_Rz=std::find(tmp.begin(), tmp.end(), "rz")!=tmp.end();
        Tensor<double> ext_dof(6,3*mol.natom());
        if (remove_Tx) ext_dof(0l,_)=transx;
        if (remove_Ty) ext_dof(1l,_)=transy;
        if (remove_Tz) ext_dof(2l,_)=transz;

        if (remove_Rx) ext_dof(3l,_)=rotx;
        if (remove_Ry) ext_dof(4l,_)=roty;
        if (remove_Rz) ext_dof(5l,_)=rotz;
        print("removing dof ",remove_Tx, remove_Ty, remove_Tz, remove_Rx, remove_Ry, remove_Rz);

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
        for (size_t i=0; i<3*mol.natom(); ++i) projector(i,i)=1.0;

        // compute the outer products of the projectors
        // 1- \sum_i | t_i >< t_i |
        projector-=inner(ext_dof,ext_dof,0,0);

        return projector;

    }

    /// remove translational degrees of freedom from the hessian

    /// @param[in]	do_remove_dof	which dof to remove: x,y,z,Rx,Ry,Rz (transl/rot)
    static void remove_external_dof(Tensor<double>& hessian, const Molecule& mol,
    		const std::vector<std::string>& remove_dof) {

        // compute the translation of the center of mass
        Tensor<double> projector_ext=projector_external_dof(mol,remove_dof);

        // this is P^T * H * P
        hessian=inner(projector_ext,inner(hessian,projector_ext),0,0);
    }


    /// returns the vibrational frequencies

    /// @param[in]  hessian the hessian matrix (not mass-weighted)
    /// @param[out] normalmodes the normal modes
    /// @param[in]  project_tr whether to project out translation and rotation
    /// @param[in]  print_hessian   whether to print the hessian matrix
    /// @return the frequencies in atomic units
    static Tensor<double> compute_frequencies(const Molecule& molecule,
            const Tensor<double>& hessian, Tensor<double>& normalmodes,
            const std::vector<std::string>& remove_dof={}, const bool print_hessian=false) {

        // compute mass-weighing matrices
        Tensor<double> M=molecule.massweights();
        Tensor<double> Minv(3*molecule.natom(),3*molecule.natom());
        for (size_t i=0; i<3*molecule.natom(); ++i) Minv(i,i)=1.0/M(i,i);

        // mass-weight the hessian
        Tensor<double> mwhessian=inner(M,inner(hessian,M));

        // remove translation and rotation
        if (remove_dof.size()>0) MolecularOptimizer::remove_external_dof(mwhessian,molecule,remove_dof);

        if (print_hessian) {
            if (remove_dof.size()>0) {
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


    static Tensor<double> compute_reduced_mass(const Molecule& molecule,
            const Tensor<double>& normalmodes) {

        Tensor<double> M=molecule.massweights();
        Tensor<double> D=MolecularOptimizer::projector_external_dof(molecule,{"Tx","Ty","Tz","Rx","Ry","Rz"});
        Tensor<double> L=copy(normalmodes);
        Tensor<double> DL=inner(D,L);
        Tensor<double> MDL=inner(M,DL);
        Tensor<double> mu(3*molecule.natom());

        for (size_t i=0; i<3*molecule.natom(); ++i) {
            double mu1=0.0;
            for (size_t j=0; j<3*molecule.natom(); ++j) mu1+=MDL(j,i)*MDL(j,i);
            if (mu1>1.e-14) mu(i)=1.0/(mu1*constants::atomic_mass_in_au);
        }
        return mu;
    }

};

}

#endif //MADNESS_CHEM_MOLECULAR_OPTIMIZER_H__INCLUDED
