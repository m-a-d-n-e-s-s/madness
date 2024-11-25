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

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/funcplot.h>
#include <madness/constants.h>
#include <madness/world/test_utilities.h>

using namespace madness;

bool smalltest = false;

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

class Charge : public FunctionFunctorInterface<double,3> {
public:
    double operator()(const coord_3d& r) const {
        static const double fac = std::pow(constants::pi,-1.5);
        const double rsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        return fac*exp(-rsq);
    }
};

class DPot : public FunctionFunctorInterface<double,3> {
    const int dir;
public:
    DPot(int dir) : dir(dir) {}

    double operator()(const coord_3d& r) const {
        const double rsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        const double R = sqrt(rsq);
        double dudr;

        if (R>7.0) {
            dudr = -1.0/rsq;
        }
        else if (R<0.01) {
            dudr = (-.75225277806367504925+(.45135166683820502956+(-.16119702387078751056+0.41791821003537502737e-1*rsq)*rsq)*rsq)*R;
        }
        else {
            dudr = 2*exp(-rsq)/(sqrt(constants::pi)*R) - erf(R)/rsq;
        }

        return dudr * r[dir] / R; // might need taylor expansion around r=0
    }
};


int test_opdir(World& world) {
    int success=0;
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
    //FunctionDefaults<3>::set_k(10);
    //FunctionDefaults<3>::set_thresh(1e-8);
    FunctionDefaults<3>::set_initial_level(3);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(true);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_truncate_on_project(true);

    real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new Gaussian<double,3>(origin, expnt, coeff)));
    f.truncate();

    double norm = f.trace();
    double ferr = f.err(Gaussian<double,3>(origin, expnt, coeff));
    if (world.rank() == 0) print("norm of initial function", norm, ferr);

    const real_tensor width = FunctionDefaults<3>::get_cell_width();
    const int k = FunctionDefaults<3>::get_k();

    // These from previous computation with k=8 thresh=1e-6
    // (error is consistently reduced as compute with higher accuracy)
    const double errs[] = {5.8e-06, 5.2e-06, 2.1e-05, 5.0e-06, 8.1e-06, 2.0e-05,
                           1.2e-05, 3.1e-05, 4.5e-05, 1.2e-05, 2.8e-05, 5.7e-05,
                           1.0e-05, 4.2e-05, 5.0e-05, 3.2e-05, 6.0e-05, 5.4e-05,
                           2.3e-05, 1.3e-04, 1.2e-04, 3.1e-05, 1.5e-04, 1.1e-04,
                           1.5e-04, 2.0e-04, 2.5e-04 };

    const char* msg[] = {"FAIL <<<<<<<<<<<<<","PASS"};
    int inderr = 0;

    int mhi=1;
    if (!smalltest) mhi=2;
    for (int mx=0; mx<=mhi; mx++) {
        for (int my=0; my<=mhi; my++) {
            for (int mz=0; mz<=mhi; mz++) {
                if (smalltest && mx+my+mz>1) break;
                
                const int m[3] = {mx,my,mz};
                std::vector< ConvolutionND<double,3> > ops(1);
                for (int d=0; d<3; d++) {
                    double e = expnts[d]*width[d]*width[d];              // Exponent in sim coords
                    double c = sqrt(expnts[d]/constants::pi)*width[d];   // Coeff of user-coords normalized gaussian scaled to sim coords
                    c *= pow(width[d],-m[d]);
                    ops[0].setop(d, std::shared_ptr< Convolution1D<double> >(
                              new GaussianConvolution1D<double>(k, c, e, m[d], false, 0.0)));
                }

                real_convolution_3d op(world, ops);
                op.doleaves=true;

                real_function_3d opf = op(f);
                //double oval=opf(origin), ovalexact=OpFExact(expnt,expnts,m)(origin);
                //if (world.rank() == 0) print("opf at origin", oval, ovalexact);
                double opfnorm = opf.norm2();
                double opferr = opf.err(OpFExact(expnt,expnts,m));
                const real_function_3d exact = real_factory_3d(world)
                		.functor(real_functor_3d(new OpFExact(expnt,expnts,m)));
                double exactnorm=exact.norm2();


				bool correct=opferr < 1.1 * errs[inderr++];
				if (world.rank() == 0)
					print("m =", m, ", norm =", opfnorm, ", norm exact =",
							exactnorm, ", err =", opferr, msg[correct]);
				if (not correct) success++;

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
    return success;
}


int testgradG(World& world) {
    // The potential due to a normalized gaussian with exponent 1.0 is
    //
    // u(r) = erf(r)/r
    //
    // du/dx = du/dr * dr/dx
    //
    // du/dr = 2*exp(-r^2)/(sqrt(Pi)*r)-erf(r)/r^2
    //
    // dr/dx = x/r

	int success=0;
    real_tensor cell(3,2);
    cell(0,0)=-20; cell(0,1)=20; // Deliberately have different width and range in each dimension
    cell(1,0)=-20; cell(1,1)=30;
    cell(2,0)=-40; cell(2,1)=40;

    // This will give a uniform box
    //cell(_,0)=-200;
    //cell(_,1)= 200;

    //FunctionDefaults<3>::set_cubic_cell(-20,20);

    const int k = 10;
    const double thresh = 1e-8;

    FunctionDefaults<3>::set_cell(cell);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_initial_level(3);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(true);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_truncate_on_project(false);

    std::vector<real_convolution_3d_ptr> g = GradCoulombOperator(world, 1e-3, thresh);

    real_function_3d q = real_factory_3d(world).functor(real_functor_3d(new Charge()));
    //q.truncate();
    double chargeerr = q.trace()-1.0;
    if (world.rank() == 0) print("err in Q",chargeerr);

    for (int d=0; d<3; d++) {
        real_function_3d dq = (*g[d])(q);
        double err = dq.err(DPot(d));
        if (world.rank() == 0) {
            if (err < 1.1e-05) print(d,err,"PASS");
            else {
            	print(d,err,"FAIL");
            	success++;
            }
        }
    }
    return success;
}



/// test the various operators following op_type
template<std::size_t NDIM>
int test_combined_operators(World& world) {
    FunctionDefaults<NDIM>::set_cubic_cell(-20,20);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(1.e-5);

    test_output t("testing combined operators");
//        t.set_cout_to_terminal();

    // parameters for the convolutions
    double mu=1.0;
    double lo=1.e-6;
    double thresh=FunctionDefaults<NDIM>::get_thresh();

    // we assume these three are actually working and correct
    auto slater=SeparatedConvolution<double,NDIM>(world,OperatorInfo(mu,lo,thresh,OT_SLATER));
    auto slater2=SeparatedConvolution<double,NDIM>(world,OperatorInfo(2.0*mu,lo,thresh,OT_SLATER));
    auto bsh=SeparatedConvolution<double,NDIM>(world,OperatorInfo(mu,lo,thresh,OT_BSH));
    auto bsh2=SeparatedConvolution<double,NDIM>(world,OperatorInfo(2.0*mu,lo,thresh,OT_BSH));
    auto coulomb=SeparatedConvolution<double,NDIM>(world,OperatorInfo(mu,lo,thresh,OT_G12));

    for (int itype=int(OT_ONE); itype<int(OT_SIZE); ++itype) {
        OpType type=OpType(itype);
        if (type==OT_ONE) continue;
        if (type==OT_G12) continue;
        if (type==OT_BSH) continue;
        if (type==OT_GAUSS) continue;
        if (type==OT_SLATER) continue;

        std::stringstream ss;
        ss << type;
        std::string type_str=ss.str();

        // test construction
        double mu=1.0;
        double lo=1.e-6;
        double thresh=FunctionDefaults<NDIM>::get_thresh();
        OperatorInfo info(mu,lo,thresh,type);
        SeparatedConvolution<double,NDIM> op(world,info);

        // argument is a Gaussian function exp(-r^2) in any dimension
        Function<double,NDIM> arg=FunctionFactory<double,NDIM>(world)
                .functor([](const Vector<double,NDIM>& r){return exp(-inner(r,r));});

        // test application
        Function<double,NDIM> result=op(arg);
        Function<double,NDIM> ref;
        print("result norm",result.norm2());
        double fourpi=4.0*constants::pi;

        // numerical checks
        if (type==OT_F12) { // (1 - Slater)/(2 mu)
            ref=0.5*(arg.trace() - slater(arg));
        } else if (type==OT_F212) { // ((1 - Slater)/(2 mu) )^2  = 1/(4 mu^2) (1 - 2 Slater + Slater2 )
            ref=0.25*(arg.trace() - 2*slater(arg) + slater2(arg));
        } else if (type==OT_F2G12) { // ((1 - Slater)/(2 mu))^2 1/g  = 1/(4 mu^2) (g12 - 2 bsh + bsh2 )
            ref=0.25*(coulomb(arg)- 2*fourpi*bsh(arg) + fourpi*bsh2(arg));
        } else if (type==OT_FG12) {  // (1 - Slater)/(2 mu) 1/g  = 1/(2 mu) (g12 - bsh )
            ref=0.5*(coulomb(arg)- fourpi*bsh(arg));
        }
        double error=(ref-result).norm2()/ref.norm2();
        print("refnorm",ref.norm2());
        print("diff norm  ",(ref-result).norm2());

        t.checkpoint(error,thresh,type_str);
    }

    return t.end();

}


int main(int argc, char**argv) {
    World& world=initialize(argc,argv);
    int success=0;

    try {
        startup(world,argc,argv);

        if (getenv("MAD_SMALL_TESTS")) smalltest=true;
        for (int iarg=1; iarg<argc; iarg++) if (strcmp(argv[iarg],"--small")==0) smalltest=true;
        std::cout << "small test : " << smalltest << std::endl;

        success+=test_opdir(world);
        success+=test_combined_operators<3>(world);
        if (!smalltest) success+=testgradG(world);
    }
    catch (const SafeMPI::Exception& e) {
        print(e);
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
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
//    catch (const char* s) {
//        print(s);
//        error("caught a c-string exception");
//    }
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

    return success;
}

