/*
 * test_exchangeoperator.cc
 *
 *  Created on: Sep 4, 2020
 *      Author: fbischoff
 */

#include <madness.h>
#include <chem/SCFOperators.h>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/correlationfactor.h>
#include<typeinfo>

using namespace madness;
namespace madness {

/// anchor test for the exchange operator -- partially hardwired
template<typename T>
int exchange_anchor_test(World& world, Exchange<T,3>& K, const double thresh) {

    const int nmo=2;
    Tensor<double> alpha(nmo);
    alpha(0l)=1.0;
    alpha(1l)=2.0;
    Vector<double,3> origin(0.0);

    // construct two dummy orbitals
    std::vector<Function<T,3> > amo(nmo);
    for (int i=0; i<nmo; ++i) {
        amo[i]=FunctionFactory<T,3>(world).truncate_on_project()
                .functor(GaussianGuess<T,3>(origin,alpha(i))).thresh(thresh*0.1);
    }
    Tensor<double> aocc(nmo);
    aocc.fill(1.0);

    // compute the reference result functions
    std::vector<Function<T,3> > Kamo=zero_functions_compressed<T,3>(world,nmo,true);
    for (int i=0; i<nmo; ++i) {
        for (int k=0; k<nmo; ++k) {
            const double sum_alpha=alpha(i)+alpha(k);
            real_function_3d refpot=real_factory_3d(world).truncate_on_project()
                    .functor(refpotfunctor(sum_alpha)).thresh(thresh*0.1);
            Kamo[i]+=amo[k]*refpot;
        }
    }

    // compute the result functions with the exchange operator
    std::vector<Function<T,3> > Kamo1=K(amo);

    std::vector<Function<T,3> > diff=sub(world,Kamo,Kamo1);
    std::vector<double> norms=norm2s(world,diff);
    print("diffnorm in K",norms);
    int ierr=0;
    for (double& n : norms) {
        if (n>thresh*5.0) ierr++;           // tolerance 5.0 in the density
    }

    Tensor<T> Kmat=matrix_inner(world,amo,Kamo);
    Tensor<T> Kmat1=matrix_inner(world,amo,Kamo1);
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





template<typename T>
int test_exchange(World& world) {

    FunctionDefaults<3>::set_thresh(1.e-5);
    double thresh=FunctionDefaults<3>::get_thresh();
    if (world.rank()==0) print("\nentering test_exchange",thresh,typeid(T).name());
    FunctionDefaults<3>::set_cubic_cell(-10, 10);

    // construct exchange operator
    Exchange<T,3> K(world);

    const int nmo=2;
    Tensor<double> alpha(nmo);
    alpha(0l)=1.0;
    alpha(1l)=2.0;
    Vector<double,3> origin(0.0);

    // construct two dummy orbitals
    std::vector<Function<T,3> > amo(nmo);
    for (int i=0; i<nmo; ++i) {
        amo[i]=FunctionFactory<T,3>(world).truncate_on_project()
                .functor(GaussianGuess<T,3>(origin,alpha(i))).thresh(thresh*0.1);
    }
    Tensor<double> aocc(nmo);
    aocc.fill(1.0);

    K.set_parameters(conj(world,amo),amo,aocc);

    // compare the exchange operator to precomputed reference values
    int success=0;
    if (typeid(T)==typeid(double)) success+=exchange_anchor_test(world, K, thresh);
    if (success>0) return 1;

    if (!smalltest) {
    	// test hermiticity of the K operator
    	success=test_hermiticity<T,Exchange<T,3> ,3>(world, K, thresh);
    	if (success>0) return 1;

    	// test bra/ket sets being different
    	success=test_asymmetric<T,Exchange<T,3> ,3>(world, K, thresh);
    	if (success>0) return 1;
	}

    return 0;
}


int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);
    srand (time(NULL));
    bool smalltest=false;
    if (getenv("MAD_SMALL_TESTS")) smalltest=true;
    for (int iarg=1; iarg<argc; iarg++) if (strcmp(argv[iarg],"--small")==0) smalltest=true;
    std::cout << "small test : " << smalltest << std::endl;

    FunctionDefaults<3>::set_k(8); // needed for XC test to work

    int result=0;
    result+=test_exchange<double>(world);
#ifndef HAVE_GENTENSOR
    result+=test_exchange<double_complex>(world);
#endif

    if (world.rank()==0) {
        if (result==0) print("\ntests passed\n");
        else print("\ntests failed\n");
    }
    madness::finalize();
    return result;
}
} /* namespace madness */
