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
#include <chem/molecule.h>

namespace madness {


struct MolecularOptimizationTargetInterface : public OptimizationTargetInterface {

    /// return the molecule of the target
    virtual Molecule& molecule() {
        MADNESS_EXCEPTION("you need to return a molecule",1);
        return *(new Molecule());   // this is a memory leak silencing warnings
    }
};



/// Molecular optimizer derived from the QuasiNewton optimizer

/// Essentially the QuasiNewton optimizer, but with the additional feature
/// of projecting out rotational and translational degrees of freedom
class MolecularOptimizer : public OptimizerInterface {

public:
    /// same ctor as the QuasiNewton optimizer
    MolecularOptimizer(const std::shared_ptr<MolecularOptimizationTargetInterface>& tar,
            int maxiter = 20, double tol = 1e-6, double value_precision = 1e-12,
            double gradient_precision = 1e-12)
        : update("BFGS")
        , target(tar)
        , maxiter(maxiter)
        , tol(tol)
        , value_precision(value_precision)
        , gradient_precision(gradient_precision)
        , f(tol*1e16)
        , gnorm(tol*1e16)
        , printtest(false)
        , cg_method("polak_ribiere") {
    }

    /// optimize the underlying molecule

    /// @param[in]  x   the coordinates to compute energy and gradient
    bool optimize(Tensor<double>& x) {
        bool converge;
//        converge=optimize_quasi_newton(x);
        converge=optimize_conjugate_gradients(x);
        return converge;
    }

    bool converged() const {return gradient_norm()<tol;}

    double value() const {return 0.0;}

    double gradient_norm() const {return gnorm;}

private:

    /// How to update the hessian: BFGS or SR1
    std::string update;
    std::shared_ptr<MolecularOptimizationTargetInterface> target;
    const int maxiter;
    const double tol;       // the gradient convergence threshold
    const double value_precision;  // Numerical precision of value
    const double gradient_precision; // Numerical precision of each element of residual
    double f;
    double gnorm;
    Tensor<double> h;
    bool printtest;

    /// conjugate_gradients method
    std::string cg_method;

    bool optimize_quasi_newton(Tensor<double>& x) {


        if(printtest)  target->test_gradient(x, value_precision);

        bool h_is_identity = (h.size() == 0);
        if (h_is_identity) {
            int n = x.dim(0);
            h = Tensor<double>(n,n);
            for (int i=0; i<n; ++i) h(i,i) = 1.0;

            // mass-weight the initial hessian
            for (int i=0; i<target->molecule().natom(); ++i) {
                h(i  ,i  )/=(target->molecule().get_atom(i).mass);
                h(i+1,i+1)/=(target->molecule().get_atom(i).mass);
                h(i+2,i+2)/=(target->molecule().get_atom(i).mass);
            }
        }

        // the previous gradient
        Tensor<double> gp;
        // the displacement
        Tensor<double> dx;

        for (int iter=0; iter<maxiter; ++iter) {
            Tensor<double> g;

            target->value_and_gradient(x, f, g);
            print("new energy, corresponding coords and gradient",f);
            print(x);
            print(g);
            gnorm = g.normf();
            printf(" QuasiNewton iteration %2d value %.12e gradient %.2e\n",iter,f,gnorm);
            if (converged()) break;

            if (iter == 1 && h_is_identity) {
                // Default initial Hessian is scaled identity but
                // prefer to reuse any existing approximation.
                h.scale(g.trace(gp)/gp.trace(dx));
            }

            if (iter > 0) {
                if (update == "BFGS") QuasiNewton::hessian_update_bfgs(dx, g-gp,h);
                else QuasiNewton::hessian_update_sr1(dx, g-gp,h);
            }

            Tensor<double> v, e;
//            syev(h, v, e);
//            print("hessian eigenvalues",e);
            remove_translation(h,target->molecule());
            syev(h, v, e);
            print("hessian eigenvalues",e);

            // this will invert the hessian, multiply with the gradient and
            // return the displacements
            dx = new_search_direction2(g,h);

//            double step = line_search(1.0, f, dx.trace(g), x, dx);
            double step=0.5;

            dx.scale(step);
            x += dx;
            gp = g;
        }

        if (printtest) {
            print("final hessian");
            print(h);
        }
        return converged();
    }

    bool optimize_conjugate_gradients(Tensor<double>& x) {


        // initial energy and gradient gradient
        double energy=0.0;
        Tensor<double> gradient;
        target->value_and_gradient(x, energy, gradient);

        // first step is steepest descent
        Tensor<double> displacement=-gradient;
        Tensor<double> oldgradient=gradient;
        Tensor<double> old_displacement=displacement;


        for (int iter=1; iter<maxiter; ++iter) {

            // displace coordinates
            x+=displacement;

            Tensor<double> com=center_of_mass(target->molecule());
            print("current coordinates and center of mass",com);
            print(x);
            target->value_and_gradient(x, energy, gradient);
            print("new energy and gradient",energy);
            print(gradient);
            gnorm = gradient.normf();
            print("raw gradient norm ",gnorm);
            Tensor<double> project_T=projector_translation(target->molecule());
            gradient=inner(gradient,project_T);
            gnorm = gradient.normf();
            print("projected gradient norm ",gnorm);

            // compute new displacement (Fletcher-Reeves)
            double beta=0.0;
            if (cg_method=="fletcher_reeves") beta=gradient.normf()/oldgradient.normf();
            if (cg_method=="polak_ribiere") beta=gradient.normf()/(gradient-oldgradient).normf();

            displacement=-gradient + beta * old_displacement;
            old_displacement=displacement;

            if (converged()) break;
        }

        return converged();
    }

    /// effectively invert the hessian and multiply with the gradient
    Tensor<double> new_search_direction2(const Tensor<double>& g,
            const Tensor<double>& hessian) const {
        Tensor<double> dx, s;
        double tol = gradient_precision;
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
                if (printtest) printf("   forcing negative eigenvalue to be positive %d %.1e\n", i, e[i]);
                nneg++;
                //e[i] = -2.0*e[i]; // Enforce positive search direction
                e[i] = -0.1*e[i]; // Enforce positive search direction
            }
            else if (e[i] < tol) {
                if (printtest) printf("   forcing small eigenvalue to be zero %d %.1e\n", i, e[i]);
                nsmall++;
                e[i] = tol;
                gv[i]=0.0;   // effectively removing this direction
            }

            // this is the step -lambda^-1 gv
            gv[i] = -gv[i] / e[i];
            if (std::abs(gv[i]) > trust) { // Step restriction
                double gvnew = trust*std::abs(gv(i))/gv[i];
                if (printtest) printf("   restricting step in spectral direction %d %.1e --> %.1e\n", i, gv[i], gvnew);
                nrestrict++;
                gv[i] = gvnew;
            }
        }
        if (nneg || nsmall || nrestrict) printf("   nneg=%d nsmall=%d nrestrict=%d\n", nneg, nsmall, nrestrict);

        // Transform back from spectral basis to give the displacements
        // disp = -V lambda^-1 VT g = V lambda^-1 gv
        return inner(v,gv);
    }

    /// compute the projector to remove translational degrees of freedom
    Tensor<double> projector_translation(const Molecule& mol) const {
        Tensor<double> transx(3*mol.natom());
        Tensor<double> transy(3*mol.natom());
        Tensor<double> transz(3*mol.natom());
        for (int i=0; i<3*mol.natom(); i+=3) {
            transx[i]=1.0/sqrt(mol.natom());
            transy[i+1]=1.0/sqrt(mol.natom());
            transz[i+2]=1.0/sqrt(mol.natom());
        }

        Tensor<double> identity(3*mol.natom(),3*mol.natom());
        for (int i=0; i<3*mol.natom(); ++i) identity(i,i)=1.0;

        Tensor<double> project_T=identity-outer(transx,transx)
                - outer(transy,transy) - outer(transz,transz);
        return project_T;
    }

    /// compute the projector to remove rotational degrees of freedom
    Tensor<double> projector_rotation(const Molecule& mol) const {

        Tensor<double> rotx, roty, rotz;

        // diagonalize the moment of inertia
        Tensor<double> I=mol.moment_of_inertia();
        Tensor<double> v,e;
        syev(I, v, e);

        // compute the center of mass
        Tensor<double> com=center_of_mass(mol);

        for (int iatom=0; iatom<mol.natom(); ++iatom) {

            // coordinates wrt the center of mass
            Tensor<double> coord(3);
            coord(0l)=mol.get_atom(iatom).x-com(0l);
            coord(1l)=mol.get_atom(iatom).y-com(1l);
            coord(2l)=mol.get_atom(iatom).z-com(2l);

//
//            for (int i=0; i<3; ++i) {
//                rotx=
//            }
        }


    }


    /// remove translational degrees of freedom from the hessian
    void remove_translation(Tensor<double>& hessian,
            const Molecule& mol) const {

        print("projecting out translational degrees of freedom");
        // compute the translation of the center of mass
        Tensor<double> project_T=projector_translation(mol);

        // this is P^T * H * P
        hessian=inner(project_T,inner(hessian,project_T),0,0);
    }

    /// compute the center of mass
    Tensor<double> center_of_mass(const Molecule& molecule) const {
        Tensor<double> com(3);
        double xx=0.0, yy=0.0, zz=0.0, qq=0.0;
        for (unsigned int i=0; i<molecule.natom(); ++i) {
            xx += molecule.get_atom(i).x*molecule.get_atom(i).mass;
            yy += molecule.get_atom(i).y*molecule.get_atom(i).mass;
            zz += molecule.get_atom(i).z*molecule.get_atom(i).mass;
            qq += molecule.get_atom(i).mass;
        }
        com(0l)=xx/qq;
        com(1l)=yy/qq;
        com(2l)=zz/qq;
        return com;
    }


};

}

#endif //MADNESS_CHEM_MOLECULAR_OPTIMIZER_H__INCLUDED
