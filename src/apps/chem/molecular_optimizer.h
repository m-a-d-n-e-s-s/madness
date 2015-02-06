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

#include <madness/tensor/solvers.h>
#include <chem/molecule.h>

namespace madness {

/// Molecular optimizer derived from the QuasiNewton optimizer

/// Essentially the QuasiNewton optimizer, but with the additional feature
/// of projecting out rotational and translational degrees of freedom
class MolecularOptimizer : public QuasiNewton {

public:
    /// same ctor as the QuasiNewton optimizer
    MolecularOptimizer(const std::shared_ptr<OptimizationTargetInterface>& tar,
            const Molecule& mol,
            int maxiter = 20, double tol = 1e-6, double value_precision = 1e-12,
            double gradient_precision = 1e-12)
        : QuasiNewton(tar,maxiter,tol,value_precision,gradient_precision)
        , molecule(mol), cg_method("polak_ribiere") {
    }

    /// optimize the underlying molecule

    /// @param[in]  x   the coordinates to compute energy and gradient
    bool optimize(Tensor<double>& x) {
        bool converge=optimize_conjugate_gradients(x);
        return converge;
    }

    bool optimize_quasi_newton(Tensor<double>& x) {
        if (n != x.dim(0)) {
            n = x.dim(0);
            h = Tensor<double>();
        }


       if(printtest)  target->test_gradient(x, value_precision);

        bool h_is_identity = (h.size() == 0);
        if (h_is_identity) {
            h = Tensor<double>(n,n);
            for (int i=0; i<n; ++i) h(i,i) = 1.0;

            // mass-weight the initial hessian
            for (int i=0; i<molecule.natom(); ++i) {
                h(i  ,i  )/=(molecule.get_atom(i).mass);
                h(i+1,i+1)/=(molecule.get_atom(i).mass);
                h(i+2,i+2)/=(molecule.get_atom(i).mass);
            }
        }

        // the previous gradient
        Tensor<double> gp;
        // the displacement
        Tensor<double> dx;

        for (int iter=0; iter<maxiter; ++iter) {
            Tensor<double> g;
            target->value_and_gradient(x, f, g);
            gnorm = g.normf();
            printf(" QuasiNewton iteration %2d value %.12e gradient %.2e\n",iter,f,gnorm);
            if (converged()) break;

            if (iter == 1 && h_is_identity) {
                // Default initial Hessian is scaled identity but
                // prefer to reuse any existing approximation.
                h.scale(g.trace(gp)/gp.trace(dx));
            }

            if (iter > 0) {
                if (update == "BFGS") hessian_update_bfgs(dx, g-gp,h);
                else hessian_update_sr1(dx, g-gp,h);
            }

//            Tensor<double> v, e;
//            syev(h, v, e);
//            print("hessian eigenvalues",e);
//            remove_translation(h,molecule);
//            syev(h, v, e);
//            print("hessian eigenvalues",e);
//
//            print("gradient",g);
//            // project gradients onto purified hessian
//            g=inner(v,g,0,0);
//            print("gradient (proj)",g);


            dx = new_search_direction(g);

            double step = line_search(1.0, f, dx.trace(g), x, dx);

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

            print("current coordinates and energy",energy);
            print(x);
            target->value_and_gradient(x, energy, gradient);
            print("new energy and gradient",energy);
            print(gradient);
            gnorm = gradient.normf();

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

private:

    /// need the molecule for projecting rotation and translation
    const Molecule& molecule;

    /// conjugate_gradients method
    std::string cg_method;

    /// remove translational degrees of freedom from the hessian
    void remove_translation(Tensor<double>& hessian,
            const Molecule& mol) const {

        print("projecting out translational degrees of freedom");
        // compute the translation of the center of mass
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

        print("hessian");
        print(hessian);
        print("project_T");
        print(project_T);
        // this is P^T * H * P
        hessian=inner(project_T,inner(hessian,project_T),0,0);

        print("hessian (proj)");
        print(hessian);

    }


};

}
