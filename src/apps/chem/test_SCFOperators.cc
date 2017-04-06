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
#include <apps/chem/SCF.h>
#include <apps/chem/nemo.h>
#include <apps/chem/correlationfactor.h>

using namespace madness;

bool similar(double val1, double val2, double thresh=1.e-6) {
    return std::fabs(val1-val2)<thresh;
}

/// an N-dimensional real-valued Gaussian function

/// the function looks like
/// \[
/// f(r) = x^i y^j .. z^k exp(-alpha r^2)
/// \]
template<std::size_t NDIM>
class GaussianGuess : public FunctionFunctorInterface<double,NDIM> {
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
        for (std::size_t i=0; i<NDIM;++i) {
            arg+=(xyz[i]-origin[i])*(xyz[i]-origin[i]);
            prefac*=pow(xyz[i],ijk[i]);
        }
        const double e=exponent*arg;
        return prefac*exp(-e);
    }
};

/// reference potential for the Poisson solver of a Gaussian

/// Given a unnormalized Gaussian return the result of the Poisson equation:
/// \[
///     V(r) = \int dr' 1/|r-r'| exp(-alpha r'^2)
/// \]
/// Formulas taken from wikipedia:
/// http://en.wikipedia.org/wiki/Poisson%27s_equation#Potential_of_a_Gaussian_charge_density
/// \[
///     \nabla^2 \phi = \rho
///     \rho(r) = Q * exp(-alpha r^2)
///     \phi(r) = Q * 1/r erf(r/\sqrt{2} \sigma)
/// \]
/// where the Gaussian function is NOT normalized, and alpha = 1/(2 sigma^2)
struct refpotfunctor {
    double exponent;    // exponent of the corresponding Gaussian
    double sigma;       // parameter in the formulas
    double arg;         // argument of the error function
    double prefac;      // prefactor to account for unnormalized Gaussian

    /// ctor takes the exponent of the corresponding Gaussian exp(-alpha r^2)
    refpotfunctor(double alpha) : exponent(alpha) {
        // convert exponent into parameters of the above equations: sigma
        sigma=1.0/sqrt(2.0*exponent);
        arg=1./sqrt(2.)/sigma;
        prefac=pow(sigma*sqrt(2.*constants::pi),3.);
    }

    double operator()(const coord_3d& xyz) const {
        double r=xyz.normf();
        return prefac * erf(r*arg)/r;
    }
};


/// will write a test input and remove it from disk upon destruction
struct write_test_input {

    double eprec=FunctionDefaults<3>::get_thresh()*0.1;

    std::string filename_;
    write_test_input(std::string mol="lih") : filename_("test_input") {
        std::ofstream of(filename_);
        of << "dft\n";
        of << "xc hf\n";
        of << "no_orient\n";
        of << "k 8\n";
        of << "protocol 1.e-5 \n";
        of << "nuclear_corrfac  slater 2.0\n";
        of << "end\n";

        if (mol=="lih") {
            of << "geometry\n";
            of << "eprec " << eprec << std::endl;
            of << "Li 0.0    0.0 0.0\n";
            of << "H  1.4375 0.0 0.0\n";
            of << "end\n";
        } else if (mol=="hf") {
            double eprec=1.e-5;
            of << "geometry\n";
            of << "eprec " << eprec << std::endl;
            of << "F  0.1    0.0 0.2\n";
            of << "H  1.4375 0.0 0.0\n";
            of << "end\n";
        }
        of.close();
    }

    ~write_test_input() {
        std::remove(filename_.c_str());
    }

    std::string filename() const {return filename_;}
};


bool check_err(double err, double thresh, std::string msg) {
    if (fabs(err)>thresh) {
        print("\nfailing test:",msg,"\n");
        return true;
    }
    return false;
}

/// test the hermiticity of the operator op and its translational invariance

/// @param[in]  world   the world
/// @param[in]  op      the operator to be tested, must implement op(vecfuncT,vecfuncT)
/// @param[in]  thresh  the accuracy threshold
/// @return     0 if test passes, 1 if test fails
template<typename opT, std::size_t NDIM>
int test_hermiticity(World& world, const opT& op, double thresh) {

    print("hermiticity error/translational invariance");
    print("abs err       op(ij).normf()");
    for (int i=0; i<10; ++i) {
        Vector<double,NDIM> origin(double(i)*0.2);
        const std::vector<int> ijk(NDIM);   // s-symmetry

        // test hermiticity of the T operator
        std::vector<Function<double,NDIM> > amo(2);
        amo[0]=FunctionFactory<double,NDIM>(world)
                .functor(GaussianGuess<NDIM>(origin,1.0,ijk)).truncate_on_project();
        amo[1]=FunctionFactory<double,NDIM>(world)
                .functor(GaussianGuess<NDIM>(origin,2.0,ijk)).truncate_on_project();

        Tensor<double> tmat;
        tmat=op(amo,amo);
        Tensor<double> tmp=tmat-transpose(tmat);
        double err=tmp.normf()/tmp.size();
        print(err,tmat.normf());
        if (check_err(err,thresh*2.0,"hermiticity error")) return 1;  // tolerance 2
    }
    return 0;
}

/// test if the operator return correct matrix elements if bra!=ket

/// @param[in]  world   the world
/// @param[in]  op      the operator to be tested, must implement op(vecfuncT,vecfuncT)
/// @param[in]  thresh  the accuracy threshold
/// @return     0 if test passes, 1 if test fails
template<typename opT, std::size_t NDIM>
int test_asymmetric(World& world, const opT& op, double thresh) {

    Vector<double,NDIM> origin(1.2);
    const std::vector<int> ijk(NDIM);   // s-symmetry
    std::vector<int> ijk_p(NDIM);       // p-symmetry
    ijk_p[NDIM-1]=1;
    std::vector<int> ijk_d(NDIM);       // d-symmetry (not quite actually..)
    ijk_d[NDIM-1]=2;

    std::vector<Function<double,NDIM> > amo(2);
    amo[0]=FunctionFactory<double,NDIM>(world)
            .functor(GaussianGuess<NDIM>(origin,1.0,ijk)).truncate_on_project();
    amo[1]=FunctionFactory<double,NDIM>(world)
            .functor(GaussianGuess<NDIM>(origin,2.0,ijk)).truncate_on_project();

    // test asymmetric T operator
    std::vector<Function<double,NDIM> > bmo(5);
    bmo[0]=copy(amo[0]);
    bmo[1]=copy(amo[1]);
    bmo[2]=FunctionFactory<double,NDIM>(world).functor(GaussianGuess<NDIM>(origin,1.0,ijk_p))
                .truncate_on_project();
    bmo[3]=FunctionFactory<double,NDIM>(world).functor(GaussianGuess<NDIM>(origin,1.0,ijk_d))
                .truncate_on_project();
    bmo[4]=FunctionFactory<double,NDIM>(world).functor(GaussianGuess<NDIM>(origin,0.5,ijk))
                .truncate_on_project();

    Tensor<double> tmat=op(amo,amo);
    Tensor<double> tmat1=op(amo,bmo);
    Tensor<double> tmp=tmat1(_,Slice(0,1))-tmat;
    double err=tmp.normf()/tmp.size();
    print("a/symmetric algorithm error",err);
    if (check_err(err,thresh*2.0,"symmetry error")) return 1;
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
            .functor(GaussianGuess<NDIM>(origin,expo,ijk)).truncate_on_project();
    double ip_amo=0.0;
    for (std::size_t i=0; i<NDIM; ++i) {
        std::vector<int> dijk(NDIM);
        dijk[i]=1;
        Function<double,NDIM> df=FunctionFactory<double,NDIM>(world).truncate_on_project()
                .functor(GaussianGuess<NDIM>(origin,expo,dijk));
        df.scale(2.0*expo);
        ip_amo+=0.5*inner(df,df);
    }
    double T_amo=T(f,f);
    double err=fabs(ip_amo-T_amo);
    print("<Damo|Damo>",ip_amo,T_amo,err);
    if (check_err(err,thresh*2.0,"kinetic error")) return 1;

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

    const double alpha=1.5;
    Vector<double,3> origin(0.0);

    // compute a trial density and the reference potential
    real_function_3d density=real_factory_3d(world).truncate_on_project()
                .functor(GaussianGuess<3>(origin,alpha)).thresh(thresh*0.1);
    real_function_3d refpot=real_factory_3d(world).truncate_on_project()
            .functor(refpotfunctor(alpha)).thresh(thresh*0.1);
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
    if (check_err(err,thresh,"Coulomb density error")) return 1;

    // test matrix element
    real_function_3d g=real_factory_3d(world).truncate_on_project()
                    .functor(GaussianGuess<3>(origin,3.0));
    real_function_3d g2=copy(g).square();
    double refelement=inner(g2,refpot);
    double element=J(g,g);
    err=fabs(element-refelement);
    print("element, refelement, error",element,refelement,err);
    if (check_err(err,thresh,"Coulomb matrix element error"))  return 1;


    // test hermiticity of the T operator
    int success=test_hermiticity<Coulomb,3>(world, J, thresh);
    if (success>0) return 1;

    success=test_asymmetric<Coulomb,3>(world, J, thresh);
    if (err>thresh) return 1;

    return 0;
}

/// anchor test for the exchange operator -- partially hardwired
int exchange_anchor_test(World& world, Exchange& K, const double thresh) {

    const int nmo=2;
    Tensor<double> alpha(nmo);
    alpha(0l)=1.0;
    alpha(1l)=2.0;
    Vector<double,3> origin(0.0);

    // construct two dummy orbitals
    vecfuncT amo(nmo);
    for (int i=0; i<nmo; ++i) {
        amo[i]=real_factory_3d(world).truncate_on_project()
                .functor(GaussianGuess<3>(origin,alpha(i))).thresh(thresh*0.1);
    }
    Tensor<double> aocc(nmo);
    aocc.fill(1.0);

    // compute the reference result functions
    vecfuncT Kamo=zero_functions_compressed<double,3>(world,nmo,true);
    for (int i=0; i<nmo; ++i) {
        for (int k=0; k<nmo; ++k) {
            const double sum_alpha=alpha(i)+alpha(k);
            real_function_3d refpot=real_factory_3d(world).truncate_on_project()
                    .functor(refpotfunctor(sum_alpha)).thresh(thresh*0.1);
            Kamo[i]+=amo[k]*refpot;
        }
    }

    // compute the result functions with the exchange operator
    vecfuncT Kamo1=K(amo);

    vecfuncT diff=sub(world,Kamo,Kamo1);
    std::vector<double> norms=norm2s(world,diff);
    print("diffnorm in K",norms);
    int ierr=0;
    for (double& n : norms) {
        if (n>thresh*5.0) ierr++;           // tolerance 5.0 in the density
    }

    Tensor<double> Kmat=matrix_inner(world,amo,Kamo);
    Tensor<double> Kmat1=matrix_inner(world,amo,Kamo1);
//    print("kmat, kmat1");
//    print(Kmat);
//    print(Kmat1);
    double err=(Kmat-Kmat1).normf()/Kmat.size();
    print("diff in Kmat compared to reference potential",err);
    if (check_err(err,thresh,"Exchange potential error")) return 1;

    // hard-wire!
    Tensor<double> hardwire(2,2);
    hardwire(0,0)=5.96038972;
    hardwire(0,1)=3.70974661;
    hardwire(1,0)=3.70974661;
    hardwire(1,1)=2.36014231;
    err=(hardwire-Kmat1).normf()/Kmat.size();
    print("diff in Kmat compared to hardwired result ",err);
    if (check_err(err,thresh,"Exchange matrix element error")) return 1;
    return 0;
}

int test_exchange(World& world) {

    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("\nentering test_exchange",thresh);
    FunctionDefaults<3>::set_cubic_cell(-10, 10);

    // construct exchange operator
    Exchange K(world);

    const int nmo=2;
    Tensor<double> alpha(nmo);
    alpha(0l)=1.0;
    alpha(1l)=2.0;
    Vector<double,3> origin(0.0);

    // construct two dummy orbitals
    vecfuncT amo(nmo);
    for (int i=0; i<nmo; ++i) {
        amo[i]=real_factory_3d(world).truncate_on_project()
                .functor(GaussianGuess<3>(origin,alpha(i))).thresh(thresh*0.1);
    }
    Tensor<double> aocc(nmo);
    aocc.fill(1.0);

    K.set_parameters(amo,amo,aocc);

    // compare the exchange operator to precomputed reference values
    int success=exchange_anchor_test(world, K, thresh);
    if (success>0) return 1;

    // test hermiticity of the K operator
    success=test_hermiticity<Exchange,3>(world, K, thresh);
    if (success>0) return 1;

    // test bra/ket sets being different
    success=test_asymmetric<Exchange,3>(world, K, thresh);
    if (success>0) return 1;

    return 0;
}

int test_XCOperator(World& world) {

    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("\nentering test_XCOperator",thresh);
    FunctionDefaults<3>::set_cubic_cell(-10, 10);

    const double alpha=1.5;
    const double beta=10.0;
    Vector<double,3> origin1(0.0);
    Vector<double,3> origin2={1.0,3.0,-1.5};

    // test refine_to_common_level
    real_function_3d empty;
    real_function_3d f1=real_factory_3d(world).functor(GaussianGuess<3>(origin1,alpha));
    real_function_3d f2=real_factory_3d(world).functor(GaussianGuess<3>(origin2,beta));
    vector_real_function_3d v(3);
    v[0]=empty;
    v[1]=f1;
    v[2]=f2;
    refine_to_common_level(world,v);
    MADNESS_ASSERT(v[1].tree_size()==v[2].tree_size()); // should be identical
    MADNESS_ASSERT(v[0].tree_size()==0);    // no change here

    real_function_3d arho=copy(f1);
    arho.square();

    // hardwire-check several functionals for their energy, potential, and kernels
    // for Slater exchange their ratios must be 3/4 and 9/4
    std::vector<double> refvalues{
        -1.435302e+00, -1.884266e+00, -5.937228e-01,    // lda energy, pot, kernel
        -1.295368e+00, -1.727158e+00, -5.757192e-01,    // LDA_X energy, pot, kernel
        -1.515767e+00, -1.954087e+00, -5.796946e-01,    // pbe energy, pot, kernel
        -1.541359e+00, -1.972863e+00, -5.764567e-01     // bp energy, pot, kernel
    };

    std::vector<std::string> xcfuncs{"lda","LDA_X","pbe","bp"};
#ifndef MADNESS_HAS_LIBXC
    xcfuncs=std::vector<std::string>(1,"lda");
#endif
    int i=0;
    for (std::string xcfunc : xcfuncs) {
        /// custom ctor with information about the XC functional
        XCOperator xc(world,xcfunc,false,arho,arho);
        print("xc functional ",xcfunc);

        double a0=xc.compute_xc_energy();
        print("energy ",a0);
        MADNESS_ASSERT(similar(a0,refvalues[i++]));

        // compare xc potential to hardwired results
        double a1=2.0*inner(f1,xc(f1)); // factor 2 for RHF
        real_function_3d lda_pot=xc.make_xc_potential();
        double a11=2.0*inner(arho,lda_pot);
        print("potential ",a1);
        print("potential ",a11);
        print("ratio ",a0,a1*3.0/4.0);
        if (xcfunc=="LDA_X") MADNESS_ASSERT(std::fabs(a0-a1*3.0/4.0)<1.e-6);
        MADNESS_ASSERT(similar(a1,refvalues[i++]));

        // compare xc kernel to hardwired results
        double a2=inner(2.0*f1*f1,xc.apply_xc_kernel(2.0*arho)); // factors 2 for RHF
        print("kernel ",a2);
        print("ratio ",a0,a2*9.0/4.0);
        if (xcfunc=="LDA_X") MADNESS_ASSERT(std::fabs(a0-a2*9.0/4.0)<1.e-6);
        MADNESS_ASSERT(similar(a2,refvalues[i++]));

        // do spin-polarized
        for (int ispin=0; ispin<2   ; ++ispin) {
            XCOperator xc1(world,xcfunc,true,arho,arho);
            xc1.set_ispin(ispin);

            double a0a=xc1.compute_xc_energy();
            print("energy ", a0a);
            MADNESS_ASSERT(similar(a0,a0a));

            double a1a=2.0*inner(f1,xc(f1)); // factor 2 for RHF
            real_function_3d lda_pot=xc.make_xc_potential();
            double a11a=2.0*inner(arho,lda_pot);
            print("potential ",a1a);
            print("potential ",a11a);
            print("ratio ",a0a,a1a*3.0/4.0);
            if (xcfunc=="LDA_X") MADNESS_ASSERT(std::fabs(a0a-a1a*3.0/4.0)<1.e-6);
            MADNESS_ASSERT(similar(a1,a1a));

        }
        print("\n");

    }
    return 0;
}

int nuclear_anchor_test(World& world) {
    double thresh=FunctionDefaults<3>::get_thresh();
    write_test_input test_input;
    SCF calc(world,test_input.filename().c_str());
    calc.make_nuclear_potential(world);

    // test ncf=none
    std::shared_ptr<NuclearCorrelationFactor> ncf_none=
            std::shared_ptr<NuclearCorrelationFactor>(
                new PseudoNuclearCorrelationFactor(world,
                calc.molecule,calc.potentialmanager,1.0));
    ncf_none->initialize();

    Nuclear Vnuc(world,ncf_none);

    std::vector<int> ijk(3);
    Vector<double,3> origin{0,0.1,1.0};

    real_function_3d gaussian=FunctionFactory<double,3>(world)
                    .functor(GaussianGuess<3>(origin,2.0,ijk)).truncate_on_project();
    real_function_3d gaussian2=copy(gaussian).square();

    double V=Vnuc(gaussian,gaussian);
    double Vref=inner(calc.potentialmanager->vnuclear(),gaussian2);
    double err=fabs(V-Vref);
    print("V,Vref",V,Vref,err);
    if (check_err(err,thresh,"Nuclear matrix element error 1")) return 1;

    // test ncf=slater
    std::shared_ptr<NuclearCorrelationFactor> ncf=
    create_nuclear_correlation_factor(world, calc);
    ncf->initialize();

    Nuclear Vnuc1(world,ncf);
    Kinetic<double,3> T(world);
    real_function_3d R2gaussian=(gaussian*ncf->square());
    real_function_3d Rgaussian=(gaussian*ncf->function());
    double V0=inner(copy(Rgaussian).square(),calc.potentialmanager->vnuclear());
    double V1=Vnuc1(R2gaussian,gaussian);
    double T1=T(R2gaussian,gaussian);
    double V2=Vnuc(Rgaussian,Rgaussian);
    double T2=T(Rgaussian,Rgaussian);

    // < nemo | R (T + V) R | nemo > = < nemo | R^2 (T + U) | nemo >
    err=fabs(T2+V2 - T1-V1);
    print("T1,V1,T2,V2,err",T1,V1,T2,V2,err);
    print("V0",V0);
    if (check_err(err,thresh,"Nuclear matrix element error ncf")) return 1;

    return 0;
}


int test_nuclear(World& world) {

    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("\nentering test_nuclear",thresh);
    FunctionDefaults<3>::set_cubic_cell(-10, 10);

    int ierr=0;
    ierr+=nuclear_anchor_test(world);
    return ierr;
}

int dnuclear_anchor_test(World& world) {
    double thresh=FunctionDefaults<3>::get_thresh();
    write_test_input test_input("hf");
    SCF calc(world,test_input.filename().c_str());
    calc.make_nuclear_potential(world);

    // derivative of atom wrt axis
    int iatom=0;

    std::vector<int> ijk(3);
    Vector<double,3> origin{0.7,0.2,0.0};

    const real_function_3d gaussian=FunctionFactory<double,3>(world)
                    .functor(GaussianGuess<3>(origin,2.0,ijk)).truncate_on_project();
    const real_function_3d gaussian2=copy(gaussian).square();

    // test ncf=none
    std::shared_ptr<NuclearCorrelationFactor> ncf_none=
            std::shared_ptr<NuclearCorrelationFactor>(
                new PseudoNuclearCorrelationFactor(world,
                calc.molecule,calc.potentialmanager,1.0));
    ncf_none->initialize();

    for (int iaxis=0; iaxis<3; ++iaxis) {
        // compute matrix element and reference matrix element
        DNuclear DVnuc(world,ncf_none,iatom,iaxis);
        double V=DVnuc(gaussian,gaussian);
        MolecularDerivativeFunctor mdf(calc.molecule, iatom, iaxis);
        double Vref=inner(gaussian2,mdf);

        double err=fabs(V-Vref);
        print("V,Vref",V,Vref,err);
//        if (check_err(err,thresh,"DNuclear matrix element error 1")) return 1;
    }

    // reference values for ncf=slater, HF molecule
    const double u2ref=-23.758103991373474;
    const double u3ref=-0.0005781745459595776;
    Vector<double,9> u1xref{            // 3 * u1axis + derivativeaxis  // f
        0.002332606981637648,  0.002860180008999843, -0.0028602426401998825,
        0.002860180008999843, -0.00529473396843896,  -0.0009534186090795483,
       -0.002860242640199882, -0.0009534186090795483,-0.005294723600284016};
    Vector<double,6> u2xref{
        -13.737236204145585, -4.5791011175397855, 4.5790826623248995,   // f
         0.8528620408203597, -0.2312846148777162, 0.0                   // h
    };
    Vector<double,6> u3xref{
        -0.001241050130475841,-0.0003208126486983482,0.0003624837425800077, // hf: f
         0.000767304127536402,-0.0001325268685442912,-0.00004166157789034493 // hf: h
    };

//    // reference values for ncf=slater, LiH molecule
//    const double u2ref=-8.8837766341272;
//    const double u3ref=-0.008762753661850198;
//    Vector<double,9> u1xref{            // 3 * u1axis + derivativeaxis  // li
//        0.04211909838428108,  0.037540023223296595, 0.0,
//        0.037540023223296595,-0.07854524189957111,  0.0,
//        0.0,                  0.0,                  -0.08927218015277515};
//    Vector<double,6> u2xref{
//        -4.661029872588109, -1.3317228199701656,0.0,  // li
//        0.8528620408203597,-0.2312846148777162,0.0
//    };
//    Vector<double,6> u3xref{
//        -0.012460779625385277,-0.003938479758437209,0.0, // lih: li
//        0.008163434617542605,-0.002572827376788745,0.0 // lih: h
//    };


    // test ncf=slater

    // test U2 and U3
    std::shared_ptr<NuclearCorrelationFactor> ncf=
    create_nuclear_correlation_factor(world, calc);
    ncf->initialize();
    NuclearCorrelationFactor::U2_functor u2f(ncf.get());
    const double u2=inner(gaussian,u2f);
    double err1=fabs(u2-u2ref);
    print("u2,u2ref",u2,u2ref,err1);
    if (check_err(err1,thresh,"DNuclear matrix element error U2")) return 1;

    NuclearCorrelationFactor::U3_functor u3f(ncf.get());
    const double u3=-inner(gaussian,u3f);
    err1=fabs(u3-u3ref);
    print("u3,u3ref",u3,u3ref,err1);
    if (check_err(err1,thresh,"DNuclear matrix element error U3")) return 1;
    print("");

    // test U1X
    for (int iaxis=0; iaxis<3; ++iaxis) {
        for (int u1axis=0; u1axis<3; ++u1axis) {
            NuclearCorrelationFactor::U1X_functor u1xf(ncf.get(),iatom,u1axis,iaxis);
            const double u1x=-inner(gaussian,u1xf);     // note the sign
            double err2=fabs(u1x-u1xref[3*u1axis + iaxis]);
            print("u1x,u1xref, u1axis=",u1axis,"  ",u1x,u1xref[3*u1axis + iaxis],err2);
            if (check_err(err2,thresh,"DNuclear matrix element error u1x")) return 1;
        }
    }
    print("");

    // test U2X
    for (int iatom=0; iatom<2; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            NuclearCorrelationFactor::U2X_functor u2xf(ncf.get(),iatom,iaxis);
            const double u2x=inner(gaussian,u2xf);
            double err2=fabs(u2x-u2xref[3*iatom+iaxis]);
            print("u2x,u2xref, daxis=",iaxis,"  ",u2x,u2xref[3*iatom+iaxis],err2);
            if (check_err(err2,thresh*3.0,"DNuclear matrix element error u2x")) return 1;
        }
    }
    print("");

    // test U3X
    for (int iatom=0; iatom<2; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            NuclearCorrelationFactor::U3X_functor u3xf(ncf.get(),iatom,iaxis);
            const double u3x=inner(gaussian,u3xf);     // note the sign
            double err2=fabs(u3x-u3xref[3*iatom+iaxis]);
            print("u3x,u3xref, daxis=",iaxis,"  ",u3x,u3xref[3*iatom+iaxis],err2);
        if (check_err(err2,thresh,"DNuclear matrix element error u3x")) return 1;
        }
    }

//    int iaxis=0;
//    DNuclear DVnuc(world,ncf,iatom,iaxis);
//    double V=DVnuc(gaussian,gaussian);
//    double Vref=0.0;       // computed by mathematica
//    double err=fabs(V-Vref);
//    print("V,Vref",V,Vref,err);
//
//    if (check_err(err,thresh,"Nuclear matrix element error ncf")) return 1;

    return 0;
}

int test_dnuclear(World& world) {
    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("\nentering test_dnuclear",thresh);
    FunctionDefaults<3>::set_cubic_cell(-10, 10);

    int ierr=0;
    ierr+=dnuclear_anchor_test(world);
    return ierr;
}

int test_SCF(World& world) {
    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("\nentering test_SCF",thresh);

    write_test_input test_input;
    SCF calc(world,test_input.filename().c_str());
    MolecularEnergy ME(world, calc);
    double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
    print("energy(LiH)",energy);
    // hard-wire test
    print("energy, hard-wire, diff",energy,7.703833,energy+7.703833e+00);
    if (check_err(energy+7.703833e+00,thresh,"SCF error")) return 1;
    return 0;
}

int test_nemo(World& world) {
    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("\nentering test_nemo",thresh);

    write_test_input test_input;
    std::shared_ptr<SCF> calc_ptr(new SCF(world,test_input.filename().c_str()));
    Nemo nemo(world,calc_ptr);
    double energy=nemo.value(calc_ptr->molecule.get_all_coords().flat()); // ugh!
    print("energy(LiH)",energy);
    // hard-wire test
    if (check_err(energy+7.703832e+00,thresh,"nemo error")) return 1;
    return 0;
}


int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);

    int result=0;

    result+=test_kinetic<1>(world);
    result+=test_kinetic<2>(world);
    result+=test_kinetic<3>(world);
//    result+=test_kinetic<4>(world);

    result+=test_coulomb(world);
    result+=test_exchange(world);
    result+=test_XCOperator(world);
    result+=test_nuclear(world);
    result+=test_dnuclear(world);
    result+=test_SCF(world);
    result+=test_nemo(world);

    if (world.rank()==0) {
        if (result==0) print("\ntests passed\n");
        else print("\ntests failed\n");
    }
    madness::finalize();
    return result;
}
