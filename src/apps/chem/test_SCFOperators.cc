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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness.h>
#include <apps/chem/SCFOperators.h>
#include <apps/chem/molecule.h>

using namespace madness;

/// an N-dimensional real-valued Gaussian function

/// the function looks like
/// \[
/// f(r) = x^i y^j .. z^k exp(-alpha r^2)
/// \]
template<std::size_t NDIM>
class GaussianGuess : FunctionFunctorInterface<double,NDIM> {
    typedef Vector<double,NDIM> coordT;

public:

    /// ctor

    /// @param[in]  origin  the origin of the Gauss function
    /// @param[in]  alpha   the exponent exp(-alpha r^2)
    /// @param[in]  ijk     the monomial x^i y^j z^k exp(-alpha r^2) (for NDIM)
    GaussianGuess(const coordT& origin, const double alpha,
            const std::vector<int> ijk=std::vector<int>(NDIM))
            : origin(origin), exponent(alpha), ijk(ijk) {
    }

    coordT origin;
    double exponent;        ///< exponent of the guess
    std::vector<int> ijk;   ///< cartesian exponents

    double operator()(const coordT& xyz) const {
        double arg=0.0, prefac=1.0;
        for (int i=0; i<NDIM;++i) {
            arg+=(xyz[i]-origin[i])*(xyz[i]-origin[i]);
            prefac*=pow(xyz[i],ijk[i]);
        }
        const double e=exponent*arg;
        return prefac*exp(-e);
    }
};

struct refpotfunctor {
    double arg;
    refpotfunctor(double arg) : arg(arg) {}
    double operator()(const coord_3d& xyz) const {
        double r=xyz.normf();
        return erf(arg*r)/r;
    }
};


template<typename opT, std::size_t NDIM>
int test_hermiticity(World& world, const opT& op, double thresh) {

    Vector<double,NDIM> origin(0.0);
    const std::vector<int> ijk(NDIM);   // s-symmetry

    // test hermiticity of the T operator
    std::vector<Function<double,NDIM> > amo(2);
    amo[0]=FunctionFactory<double,NDIM>(world)
            .functor2(GaussianGuess<NDIM>(origin,1.0,ijk)).truncate_on_project();
    amo[1]=FunctionFactory<double,NDIM>(world)
            .functor2(GaussianGuess<NDIM>(origin,2.0,ijk)).truncate_on_project();

    Tensor<double> tmat;
    tmat=op(amo,amo);
    Tensor<double> tmp=tmat-transpose(tmat);
    double err=tmp.normf()/tmp.size();
    print("hermiticity error",err);
    if (err>thresh) return 1;
    return 0;
}

template<typename opT, std::size_t NDIM>
int test_asymmetric(World& world, const opT& op, double thresh) {

    Vector<double,NDIM> origin(0.0);
    const std::vector<int> ijk(NDIM);   // s-symmetry
    std::vector<int> ijk_p(NDIM);       // p-symmetry
    ijk_p[NDIM-1]=1;
    std::vector<int> ijk_d(NDIM);       // d-symmetry (not quite actually..)
    ijk_d[NDIM-1]=2;

    std::vector<Function<double,NDIM> > amo(2);
    amo[0]=FunctionFactory<double,NDIM>(world)
            .functor2(GaussianGuess<NDIM>(origin,1.0,ijk)).truncate_on_project();
    amo[1]=FunctionFactory<double,NDIM>(world)
            .functor2(GaussianGuess<NDIM>(origin,2.0,ijk)).truncate_on_project();

    // test asymmetric T operator
    std::vector<Function<double,NDIM> > bmo(5);
    bmo[0]=copy(amo[0]);
    bmo[1]=copy(amo[1]);
    bmo[2]=FunctionFactory<double,NDIM>(world).functor2(GaussianGuess<NDIM>(origin,1.0,ijk_p))
                .truncate_on_project();
    bmo[3]=FunctionFactory<double,NDIM>(world).functor2(GaussianGuess<NDIM>(origin,1.0,ijk_d))
                .truncate_on_project();
    bmo[4]=FunctionFactory<double,NDIM>(world).functor2(GaussianGuess<NDIM>(origin,0.5,ijk))
                .truncate_on_project();

    Tensor<double> tmat=op(amo,amo);
    Tensor<double> tmat1=op(amo,bmo);
    Tensor<double> tmp=tmat1(_,Slice(0,1))-tmat;
    double err=tmp.normf()/tmp.size();
    print("a/symmetric algorithm error",err);
    if (err>thresh) return 1;
    return 0;

}


template<std::size_t NDIM>
int test_kinetic(World& world) {
    if (world.rank()==0) print("entering test_kinetic with dimension ",NDIM);

    FunctionDefaults<NDIM>::set_cubic_cell(-10, 10);
    double thresh=FunctionDefaults<NDIM>::get_thresh();

    Kinetic<double,NDIM> T(world);

    Vector<double,NDIM> origin(0.0);
    const std::vector<int> ijk(NDIM);   // s-symmetry
    std::vector<int> ijk_p(NDIM);       // p-symmetry
    ijk_p[NDIM-1]=1;
    std::vector<int> ijk_d(NDIM);       // d-symmetry (not quite actually..)
    ijk_d[NDIM-1]=2;

    // compare T operator with integration by parts
    double expo=2.0;
    Function<double,NDIM> f=FunctionFactory<double,NDIM>(world)
            .functor2(GaussianGuess<NDIM>(origin,expo,ijk)).truncate_on_project();
    double ip_amo=0.0;
    for (std::size_t i=0; i<NDIM; ++i) {
        std::vector<int> dijk(NDIM);
        dijk[i]=1;
        Function<double,NDIM> df=FunctionFactory<double,NDIM>(world).truncate_on_project()
                .functor2(GaussianGuess<NDIM>(origin,expo,dijk));
        df.scale(2.0*expo);
        ip_amo+=0.5*inner(df,df);
    }
    double T_amo=T(f,f);
    double err=fabs(ip_amo-T_amo);
    print("<Damo|Damo>",ip_amo,T_amo,err);
    if (err>thresh) return 1;

    // test hermiticity of the T operator
    int success=test_hermiticity<Kinetic<double,NDIM>,NDIM>(world, T, thresh);
    if (success>0) return 1;

    success=test_asymmetric<Kinetic<double,NDIM>,NDIM>(world, T, thresh);
    if (err>thresh) return 1;

    return 0;
}

int test_coulomb(World& world) {
    if (world.rank()==0) print("\nentering test_coulomb");

    FunctionDefaults<3>::set_cubic_cell(-10, 10);
    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("thresh",thresh);

    // A Gaussian charge density yields an erf potential
    // \rho(r) = Q/(sigma (2\pi)^3/2 ) exp(-r^2/(2 sigma^2)
    // V(r)    = 1/(4\pi)Q/r erf(r/(sqrt(2)*sigma))

    const double sigma=2.0;
    const double Q=3.0;
    const double pi=constants::pi;
    Vector<double,3> origin(0.0);

    // compute a trial density
    real_function_3d density=real_factory_3d(world).truncate_on_project()
            .functor2(GaussianGuess<3>(origin,0.5/(sigma*sigma))).thresh(thresh*0.1);
    double prefac=Q/pow(2.0*pi,1.5)/(sigma*sigma*sigma);
    density.scale(prefac);

    // compute the reference potential
    real_function_3d refpot=real_factory_3d(world).truncate_on_project()
            .functor2(refpotfunctor(1.0/(sqrt(2)*sigma))).thresh(thresh*0.1);
    refpot.scale(Q);
    double refpotnorm=refpot.norm2();
    print("refpotnorm",refpotnorm);

    // compute the potential from the trial density
    Coulomb J(world);
    J.potential()=J.compute_potential(density,1e-5);
    double Jpotnorm=J.potential().norm2();
    print("Jpotnorm  ",Jpotnorm);

    // compare potentials
    real_function_3d diffdensity=J.potential()-refpot;
    double err=diffdensity.norm2()/Jpotnorm;
    print("relative error in the densities: ",err);
    if (err>thresh) return 1;       // some safety

    // test matrix element
    real_function_3d g=real_factory_3d(world).truncate_on_project()
                    .functor2(GaussianGuess<3>(origin,3.0));
    real_function_3d g2=copy(g).square();
    double refelement=inner(g2,refpot);
    double element=J(g,g);
    err=fabs(element-refelement);
    print("element, refelement, error",element,refelement,err);
    if (err>thresh) return 1;


    // test hermiticity of the T operator
    int success=test_hermiticity<Coulomb,3>(world, J, thresh);
    if (success>0) return 1;

    success=test_asymmetric<Coulomb,3>(world, J, thresh);
    if (err>thresh) return 1;

    return 0;
}

int test_exchange(World& world) {
    return 0;
}


int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);

    int result=0;
    result+=test_kinetic<1>(world);
//    result+=test_kinetic<2>(world);
//    result+=test_kinetic<3>(world);
//    result+=test_kinetic<4>(world);

    result+=test_coulomb(world);
    result+=test_exchange(world);

    if (world.rank()==0) {
        if (result==0) print("\ntests passed\n");
        else print("\ntests failed\n");
    }
    madness::finalize();
    return result;
}
