/*
 * testSepRep.cc
 *
 *  Created on: Jan 17, 2011
 *      Author: fbischoff
 */

/*
 * test the SeparatedRepresentation (SepRep) for representing matrices
 */

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "mra/SepRep.h"
#include "mra/mra.h"

using namespace madness;

int testLowRankTensor(const long& k, const long& dim, const double& eps) {




	return 0;
}

int testFullTensor(const long& k, const long& dim, const double& eps) {


	Tensor<double> t0=Tensor<double>(3,3,3,3);
	t0.fillrandom();

	// default ctor
	FullTensor<double> t1;

	// ctor with rhs=Tensor
	FullTensor<double> t2(t0);
	print("ctor with rhs=Tensor      ", (t0-t2.fullTensor()).normf());

	// ctor with rhs=FullTensor
	FullTensor<double> t3(t1);
	print("ctor with rhs=FullTensor/1", (t0-t3.fullTensor()).normf());
	FullTensor<double> t4(t2);
	print("ctor with rhs=FullTensor/2", (t0-t4.fullTensor()).normf());

	// assignment with rhs=Tensor
	t3=t0;
	print("assignment with rhs=Tensor", (t0-t3.fullTensor()).normf());

	// assignment with rhs=FullTensor
	t3=t1;
	print("assignment with rhs=FullTensor", (t0-t1.fullTensor()).normf());

	// virtual ctor
	SepRepTensor<double>* sr1=t3.clone();
	print("virtual ctor with rhs=FullTensor", (sr1->fullTensor()-t3.fullTensor()).normf());

	return 0;
}


int main(int argc, char**argv) {

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    srand(time(NULL));

    // the parameters
    const long k=5;
    const unsigned int dim=3;
    double eps=1.e-3;

#if 0
    // construct random SepRep and try to represent it
    SepRep<double> sr1(TT_3D,k,dim);
    sr1.fillWithRandom(2);
    Tensor<double> d1=sr1.reconstructTensor();

    SepRep<double> sr2(TT_3D,k,dim);
    sr2.fillWithRandom(3);
    Tensor<double> d2=sr2.reconstructTensor();

    SepRep<double> sr3(TT_3D,k,dim);
    sr3=sr1;
    sr3+=sr2;
    Tensor<double> d3=sr3.reconstructTensor();

    sr3.reduceRank(eps);
    Tensor<double> d4=sr3.reconstructTensor();

    print("residual for addition      ", (d1+d2-d3).normf());
    print("residual for rank reduction", (d3-d4).normf());

    // construct SepRep from tensor
    Tensor<double> t(3,3,3,3,3,3);
    t.fillrandom();

    SepRep<double> sr(t,eps,TT_3D);
    print("3d SR rank",sr.rank());
    Tensor<double> t2=sr.reconstructTensor();
    print("residual for rank reduction", (t-t2).normf());

    SepRep<double> sr4(t,eps,TT_2D);
    print("2d SR rank",sr4.rank());
    Tensor<double> t4=sr4.reconstructTensor();
    print("residual for rank reduction", (t-t4).normf());
#endif

#if 0
    // do some benchmarking
    Tensor<double> t5(5,5,5,5,5,5);
    t5.fillindex();
    t5.fillrandom();
    t5=0.0;
    Tensor<double> t6;
    if(world.rank() == 0) print("starting at time", wall_time());
    for (unsigned int i=0; i<1000; i++) {
        SepRep<double> sr5(t5,eps,TT_2D);
        SepRep<double> tmp(sr5);
        tmp+=sr5;
        tmp.reduceRank(eps);
        print("tmp.rank()",tmp.rank());
        t6=t5-sr5.reconstructTensor();
        print("error norm",t6.normf());
        t6=sr5.reconstructTensor();
    }
    if(world.rank() == 0) print("ending at time  ", wall_time());
    print(t6);
#endif

    print("hello world");

    testFullTensor(k,dim,eps);

    world.gop.fence();
    finalize();

    return 0;
}
