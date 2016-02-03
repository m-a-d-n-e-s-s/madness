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

  $Id$
*/

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file examples/nemo.cc
  \brief solve the HF equations using numerical exponential MOs

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/nemo.cc>here</a>.

*/

#include <chem/nemo.h>
#include <chem/molecular_optimizer.h>
#include <chem/cheminfo.h>
#include <chem/SCFOperators.h>

using namespace madness;

namespace madness {

const static double a=2.0;
double exp_r_inv(const coord_3d& r) {
    return exp(-a*sqrt(r.normf()));
}
double exp_r(const coord_3d& r) {
    return exp(a*sqrt(r.normf()));
}
double dx_exp_r(const coord_3d& r) {
    double rr=sqrt(r.normf());
    return a*r[0]*exp(a*rr);
}

double dy_exp_r(const coord_3d& r) {
    double rr=sqrt(r.normf());
    return a*r[1]*exp(a*rr);
}

double dz_exp_r(const coord_3d& r) {
    double rr=sqrt(r.normf());
    return a*r[2]*exp(a*rr);
}

double d2_exp_r(const coord_3d& r) {
    double rr=sqrt(r.normf());
    return a*a*exp(a*rr) + 2.0*a/(rr+0.00001)*exp(rr);
}

struct dsmooth : public FunctionFunctorInterface<double,3> {
    typedef std::shared_ptr<NuclearCorrelationFactor> ncf_ptr;
    ncf_ptr ncf;
    std::shared_ptr<SCF> calc;

    dsmooth(World& world, std::shared_ptr<SCF> calc) : calc(calc) {
        ncf=ncf_ptr(new Slater(world, calc->molecule, 2.0));
    }

    double operator()(const coord_3d& xyz) const {
        return ncf->dsmoothed_unitvec(xyz,0,calc->molecule.get_eprec())[2];
    }
};



void Nemo::do_stuff() {
    dsmooth ds(world,calc);

    plot_plane<3,dsmooth>(world,ds,"ds01");
    real_function_3d f=real_factory_3d(world).functor(ds);
    save(f,"dsmooth");


    return;
    const vecfuncT& nemo=calc->amo;
    const real_function_3d rhonemo=2.0*make_density(calc->aocc, calc->amo);
    save(rhonemo,"rhonemo");
    const real_function_3d rho = (R_square*rhonemo);
    save(rho,"rho");

    std::vector<real_function_3d> ddens(3);
    ddens[0]=this->make_ddensity(rhonemo,0);
    ddens[1]=this->make_ddensity(rhonemo,1);
    ddens[2]=this->make_ddensity(rhonemo,2);

    real_function_3d vsigaa=real_factory_3d(world);
    load(vsigaa,"vsigaa");
    std::vector<real_function_3d> ddel=mul(world,4.0*vsigaa,ddens);

    // case 1: apply nabla first, then convolve with the Poisson operator
    real_function_3d gga_pot1=real_factory_3d(world).compressed();
    for (int axis=0; axis<3; axis++) {
        Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
        functionT vxc2=D(ddel[axis]);
        gga_pot1-=vxc2;//.truncate();
    }
    real_function_3d Ggga_pot1=(*poisson)(gga_pot1);
    save(Ggga_pot1,"Ggga_pot1");

    // case 2: apply the derivative Coulomb operator, accumulate its results
    const double thresh=FunctionDefaults<3>::get_thresh();
    std::vector<real_convolution_3d_ptr> g = GradCoulombOperator(world, 1e-3, thresh);
    std::vector<real_function_3d> Gddel=apply(world,g,ddel);
    real_function_3d Ggga_pot2=real_factory_3d(world).compressed();
    for (int axis=0; axis<3; ++axis) {
        Ggga_pot2-=Gddel[axis];
    }
    save(Ggga_pot2,"Ggga_pot2");


    return;

    Laplacian<double,3> DD(world);
    real_function_3d d2nemo=DD(rhonemo);
    save(d2nemo,"d2rhonemo");

//    typedef NonlinearSolver solverT;
//    solverT solver;
//
//    // solve the inverse Poisson equation
//    for (int i=0; i<20; ++i) {
//        real_function_3d rho_rec=-1.0/(4.0*M_PI)*(*poisson)(d2nemo);
//        real_function_3d res=(rho_rec-rhonemo).truncate();
//        d2nemo = (solver.update(d2nemo, res));
//        save(d2nemo,"d2nemo"+stringify(i));
//        save(res,"res"+stringify(i));
//        double rnorm=res.norm2();
//        print("finished iteration ",i," with residual norm", rnorm);
//    }
//


    real_function_3d d2rho=this->make_laplacian_density(rhonemo);
    save(d2rho,"d2rho");

    XCOperator xcop(world,this);
    xcop.make_xc_potential();
    xcop.apply_xc_kernel(rho);
    throw;

//
//    Derivative<double,3> D0 = free_space_derivative<double,3>(world, 0);
//    Derivative<double,3> D1 = free_space_derivative<double,3>(world, 1);
//    Derivative<double,3> D2 = free_space_derivative<double,3>(world, 2);
//    real_function_3d ddens0b=D0(rho);
//    real_function_3d ddens1b=D1(rho);
//    real_function_3d ddens2b=D2(rho);
//    save(ddens0a,"ddens0a");
//    save(ddens1a,"ddens1a");
//    save(ddens2a,"ddens2a");
//    save(ddens0b,"ddens0b");
//    save(ddens1b,"ddens1b");
//    save(ddens2b,"ddens2b");
//
//    real_function_3d drhonemo0=D0(rhonemo);
//    real_function_3d drhonemo1=D1(rhonemo);
//    real_function_3d drhonemo2=D2(rhonemo);
//    save(drhonemo0,"drhonemo0");
//    save(drhonemo1,"drhonemo1");
//    save(drhonemo2,"drhonemo2");

    real_function_3d precond=real_factory_3d(world).f(exp_r);
    real_function_3d precond_inv=real_factory_3d(world).f(exp_r_inv);

    std::vector<real_function_3d> dprecond(3);
    dprecond[0]=real_factory_3d(world).f(dx_exp_r);
    dprecond[1]=real_factory_3d(world).f(dy_exp_r);
    dprecond[2]=real_factory_3d(world).f(dz_exp_r);

    real_function_3d d2_precond=real_factory_3d(world).f(d2_exp_r);

    real_function_3d rhof=rhonemo*precond;
    save(rhof,"rhof");
    real_function_3d d2rhof=DD(rhof);
    save(d2rhof,"d2rhof");

    real_function_3d tmp=dot(world,ddens,dprecond);
    real_function_3d test=precond_inv*(d2rhof - 2.0*tmp - rho*d2_precond);
    save(test,"d2rho_precond");

    throw;

}


}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  NEMO -- Hartree-Fock using numerical exponential molecular orbitals \n");
    	printf("starting at time %.1f\n", wall_time());

    }
    startup(world,argc,argv);
    std::cout.precision(6);

#ifdef MADNESS_GITREVISION
    const  char* gitrev =  MADNESS_GITREVISION;
    const std::string gitrevision(gitrev);
    if (world.rank()==0) {
    	print("    main() git revision ...",gitrevision);
    }
#endif

    if (world.rank()==0) {
        print("     main() compiled at ...",__TIME__," on ",__DATE__);
        const std::string gitrevision(info::cheminfo_git_commit());
        print("   chemlib git revision ...",gitrevision);
    }

    try {

        const std::string input="input";
        std::shared_ptr<SCF> calc(new SCF(world,input.c_str()));
        if (world.rank()==0) {
            calc->molecule.print();
            print("\n");
            calc->param.print(world);
        }

        std::shared_ptr<Nemo> nemo(new Nemo(world,calc));

        // optimize the geometry if requested
        if (calc->param.gopt) {
            print("\n\n Geometry Optimization                      ");
            print(" ----------------------------------------------------------\n");
            calc->param.gprint(world);

            Tensor<double> geomcoord = calc->molecule.get_all_coords().flat();
//            MolecularOptimizer geom(std::shared_ptr<MolecularOptimizationTargetInterface>(new Nemo(world, calc)),
            MolecularOptimizer geom(nemo,
                    calc->param.gmaxiter,
                    calc->param.gtol,  //tol
                    calc->param.gval,  //value prec
                    calc->param.gprec); // grad prec
//            geom.set_update(calc->param.algopt);
//            geom.set_test(calc->param.gtest);

            // compute initial hessian
            if (calc->param.ginitial_hessian) {
                nemo->value();
                Tensor<double> hess=nemo->hessian(calc->molecule.get_all_coords());
                geom.set_hessian(hess);
            }
            geom.optimize(geomcoord);
        } else {

            // compute the energy to get converged orbitals
//            Nemo nemo(world,calc);
            const double energy=nemo->value();
            if (world.rank()==0) {
                printf("final energy   %12.8f\n", energy);
                printf("finished at time %.1f\n", wall_time());
            }

        }

        // compute the hessian
        if (calc->param.hessian) nemo->hessian(calc->molecule.get_all_coords());


    } catch (const SafeMPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (char* s) {
        print(s);
        error("caught a string exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }


    finalize();
    return 0;
}
