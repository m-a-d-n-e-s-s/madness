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

/// \file testopdir.cc
/// \brief test different operators in different directions

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/funcplot.h>
#include <constants.h>

using namespace madness;

template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;

    Gaussian(const coordT& center, double exponent, T coefficient)
            : center(center), exponent(exponent), coefficient(coefficient) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    };
};

// The result of convolving a normalized Gaussian of expnt centered
// at the origin with a product of normalized Gaussians with
// different exponents in directions x, y, z.  The order of
// the derivative of the kernel in direction d is m[d]
class OpFExact : public FunctionFunctorInterface<double,3> {
    double expnts[3], fac;
    int m[3];

public:
    OpFExact(double expnt, const double (&expnts)[3], const int (&m)[3]) {
        this->fac = 1.0;
        for (int d=0; d<3; d++) {
            this->m[d] = m[d];
            this->expnts[d] = expnt*expnts[d]/(expnt+expnts[d]);

            if (m[d] == 0)
                fac *= this->expnts[d]/constants::pi;
            else if (m[d] == 1)
                fac *= std::pow(this->expnts[d],3.0)/constants::pi;
            else if (m[d] == 2)
                fac *= std::pow(this->expnts[d],5.0)/constants::pi;

        }
        this->fac = sqrt(fac);
    }

    double operator()(const coord_3d& r) const {
        double e = fac*exp(-(expnts[0]*r[0]*r[0] + expnts[1]*r[1]*r[1] + expnts[2]*r[2]*r[2]));
        for (int d=0; d<3; d++) {
            if (m[d] == 1)
                e *= -2.0*r[d];
            else if (m[d] == 2)
                e *=4.0*r[d]*r[d] - 2.0/expnts[d];
        }
        return e;
    }
};


void test_opdir(World& world) {
    const coord_3d origin(0.0);
    const double expnt=1.0, coeff=pow(expnt/constants::pi,1.5); // normalized gaussian with unit exponent
    const double expnts[3] = {5.0,6.0,7.0};

    if (world.rank() == 0)
        print("Test different operations in different dirs");

    real_tensor cell(3,2);
    cell(0,0)=-20; cell(0,1)=20; // Deliberately have different width and range in each dimension
    cell(1,0)=-20; cell(1,1)=30;
    cell(2,0)=-40; cell(2,1)=40;
    //FunctionDefaults<3>::set_cubic_cell(-20,20);
    FunctionDefaults<3>::set_cell(cell);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1e-6);
    FunctionDefaults<3>::set_initial_level(3);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(true);
    FunctionDefaults<3>::set_truncate_mode(0);
    FunctionDefaults<3>::set_truncate_on_project(false);

    real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new Gaussian<double,3>(origin, expnt, coeff)));
    //f.truncate(); // Comment this out to get 20x reduction in error
    f.reconstruct();

    double norm = f.trace();
    double ferr = f.err(Gaussian<double,3>(origin, expnt, coeff));
    if (world.rank() == 0) print("norm of initial function", norm, ferr);

    const real_tensor& width = FunctionDefaults<3>::get_cell_width();
    const int k = FunctionDefaults<3>::get_k();

    // These from previous computation with k=8 thresh=1e-6
    // (error is consistently reduced as compute with higher accuracy)
    const double errs[] = {5.0e-07,1.2e-06,1.4e-05,7.4e-07,1.4e-06,1.9e-05,
                           1.2e-05,3.8e-05,4.3e-05,6.8e-07,1.3e-06,6.6e-05,
                           1.0e-06,2.8e-05,6.3e-05,1.2e-05,7.9e-05,6.2e-05,
                           5.9e-06,1.9e-04,2.0e-04,6.0e-06,1.5e-04,1.3e-04,
                           1.8e-04,2.5e-04,2.1e-04};
    const char* msg[] = {"FAIL <<<<<<<<<<<<<","PASS"};
    int inderr = 0;
    for (int mx=0; mx<=2; mx++) {
        for (int my=0; my<=2; my++) {
            for (int mz=0; mz<=2; mz++) {
                const int m[3] = {mx,my,mz};
                std::vector< ConvolutionND<double,3> > ops(1);
                for (int d=0; d<3; d++) {
                    double e = expnts[d]*width[d]*width[d];              // Exponent in sim coords
                    double c = sqrt(expnts[d]/constants::pi)*width[d];   // Coeff of user-coords normalized gaussian scaled to sim coords
                    c *= pow(width[d],-m[d]);
                    ops[0].setop(d, std::shared_ptr< Convolution1D<double> >(new GaussianConvolution1D<double>(k, c, e, m[d], false, 0.0)));
                }

                real_convolution_3d op(world, ops);

                real_function_3d opf = op(f);
                //double oval=opf(origin), ovalexact=OpFExact(expnt,expnts,m)(origin);
                //if (world.rank() == 0) print("opf at origin", oval, ovalexact);
                double opfnorm = opf.trace();
                double opferr = opf.err(OpFExact(expnt,expnts,m));
                if (world.rank() == 0)
                    print("m =", m, ", norm =", opfnorm, ", err =", opferr, msg[opferr < 1.1*errs[inderr++]]);

                // This stuff useful for diagnosing problems
                // for (int i=-10; i<=10; i++) {
                //     coord_3d r(i*0.1);
                //     double num = opf(r);
                //     double anl = OpFExact(expnt,expnts,m)(r);
                //     double rat = anl/num;
                //     double err = anl-num;
                //     print("     r =", r, ", numeric =", num, ", analytic =", anl, ", anl-num =", err, ", anl/num =", rat);
                // }
            }
        }
    }

    world.gop.fence();
}

/// Factory function generating operator for convolution with grad(1/r) in 3D

/// Returns a 3-vector containing the convolution operator for the x,
/// y, and z components of grad(1/r)
// static
// inline
// std::vector< std::shared_ptr< SeparatedConvolution<double,3> > >
// GradCoulombOperator(World& world,
//                     double lo,
//                     double eps,
//                     const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
//                     int k=FunctionDefaults<3>::get_k())
// {
//     const double pi = constants::pi;
//     const Tensor<double>& width = FunctionDefaults<3>::get_width();
//     double hi = width.normf(); // Diagonal width of cell
//     const bool isperiodicsum = (bc(0,0)==BC_PERIODIC);
//     if (isperiodicsum) hi *= 100; // Extend range for periodic summation

//     // bsh_fit generates representation for 1/4Pir but we want 1/r
//     // so have to scale eps by 1/4Pi
//     Tensor<double> coeff, expnt;
//     bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);

//     if (bc(0,0) == BC_PERIODIC) {
//         truncate_periodic_expansion(coeff, expnt, width.max(), true);
//     }

//     coeff.scale(4.0*pi);
//     int rank = coeff.dim(0);

//     std::vector< std::shared_ptr< SeparatedConvolution<double,3> > > gradG(3);

//     for (int dir=0; dir<3; dir++) {
//         std::vector< ConvolutionND<double,3> > ops(rank);
//         for (int mu=0; mu<rank; mu++) {
//             // We cache the normalized operator so the factor is the value we must multiply
//             // by to recover the coeff we want.
//             Q c = std::pow(sqrt(expnt(mu)/pi),3); // Normalization coeff
//             ops[mu].setfac(coeff(mu)/c);

//             for (int d=0; d<3; d++) {
//                 if (d != dir)
//                     ops[mu].setop(d,GaussianConvolution1DCache<Q>::get(k, expnt(mu)*width[d]*width[d], 0, isperiodicsum));
//             }
//             ops[mu].setop(dir,GaussianConvolution1DCache<Q>::get(k, expnt(mu)*width[dir]*width[dir], 1, isperiodicsum));
//         }
//         ops(dir) = new SeparatedConvolution<double,3>(world, ops);
//     }

//     return gradG;
// }

void testgradG(World& world) {


}


int main(int argc, char**argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);

        test_opdir(world);

    }
    catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    finalize();

    return 0;
}

