
/// \file testSepRep.cc
/// \brief test the SeparatedRepresentation (SepRep) for representing matrices

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "mra/seprep.h"
#include "mra/mra.h"

using namespace madness;

std::string is_small(const double& val, const double& eps) {
	if (val<eps) return "ok   ";
	return "fail ";
}

int testLowRankTensor(const long& k, const long& dim, const double& eps) {




	return 0;
}

int testFullTensor_ctor(const long& k, const long& dim, const double& eps) {


	Tensor<double> t0=Tensor<double>(3,3,3,3);
	t0.fillrandom();

	// default ctor
	FullTensor<double> t1;

	// ctor with rhs=Tensor
	FullTensor<double> t2(t0);
	print("ctor with rhs=Tensor      ", (t0-t2.fullTensor()).normf());

	// ctor with rhs=FullTensor
	FullTensor<double> t3(t1);
	print("ctor with rhs=FullTensor/1", (t1.fullTensor()-t3.fullTensor()).normf());
	FullTensor<double> t4(t2);
	print("ctor with rhs=FullTensor/2", (t0-t4.fullTensor()).normf());

	// assignment with rhs=Tensor
	t3=t0;
	print("assignment with rhs=Tensor", (t0-t3.fullTensor()).normf());

	// assignment with rhs=FullTensor
	t3=t2;
	print("assignment with rhs=FullTensor", (t3.fullTensor()-t2.fullTensor()).normf());

	// virtual ctor
	SepRepTensor<double>* sr1=t3.clone();
	print("virtual ctor with rhs=FullTensor", (sr1->fullTensor()-t3.fullTensor()).normf());

	print("all done");
	return 0;
}

int testFullTensor_assignment(const long& k, const long& dim, const double& eps) {


	Tensor<double> t0=Tensor<double>(4,4,4);
	t0.fillrandom();
	std::vector<Slice> s(3,Slice(0,2,1));

	// ctor with rhs=Tensor
	FullTensor<double> t2(copy(t0));
	print("ctor with rhs=Tensor      ", (t0-t2.fullTensor()).normf());

	// assignment to a number
	t2=1.4;
	t0=1.4;
	print("assignment to a number    ", (t0-t2.fullTensor()).normf());


	// assignment with SliceTensor on both sides
	t0(s)=3.1;
	t2(s)=t0(s);
	print("sliced assignment/1       ", (t0-t2.fullTensor()).normf());

	// assignment with SliceTensor on one side
	Tensor<double> t1=t0(s);
	t2(s)=t1;
	Tensor<double> t11=t2(s);

	// test Slice
	print("sliced assignment/2       ", (t1-t11).normf());

	// test rest
	t0=1.1;
	t2=copy(t0);
	t1=copy(t0);
	t2(s)=2.8*t0(s);
	t1(s)=2.8*t0(s);
	print("sliced assignment/3       ", (t1-t2.fullTensor()).normf());

	print("all done");
	return 0;
}

int testLowRankTensor_ctor(const long& k, const long& dim, const double& eps) {


	Tensor<double> t0=Tensor<double>(7,7,7);
	Tensor<double> t1=Tensor<double>(7,7,7,7);
	t0.fillindex();
	std::vector<Slice> s(4,Slice(0,5));
	t1(s).fillrandom();
	// t1 is now sort of low-rank

	// default ctor
	LowRankTensor<double> lrt1();

	// ctor with a regular tensor
	LowRankTensor<double> lrt2(t0,eps,TT_3D);
	LowRankTensor<double> lrt3(t1,eps,TT_2D);
	print("LRT, 3-way ctor   ", lrt2.rank(), (lrt2.reconstructTensor()-t0).normf());
	print("LRT, 2-way ctor   ", lrt3.rank(), (lrt3.reconstructTensor()-t1).normf());

	// copy ctor
	LowRankTensor<double> lrt4(lrt2);
	print("copy ctor         ", (lrt2.reconstructTensor()-lrt4.reconstructTensor()).normf());



	print("all done");
	return 0;
}


int testLowRankTensor_assignment(const long& k, const long& dim, const double& eps) {

	std::vector<long> d(6,k);
	Tensor<double> t0=Tensor<double>(d);
	t0.fillindex();
	std::vector<Slice> s(6,Slice(0,k/2));

	// regular assignment
	LowRankTensor<double> lrt0(t0,eps,TT_3D);
	LowRankTensor<double> lrt0a;
	lrt0a=lrt0;
	print("LRT assignment      ", (lrt0.reconstructTensor()-lrt0a.reconstructTensor()).normf());


	// sliced assignment
	LowRankTensor<double> lrt1(t0,eps,TT_3D);
	LowRankTensor<double> lrt2=lrt1(s);
	LowRankTensor<double> lrt2a=lrt1;
	lrt2a=lrt1(s);

	Tensor<double> t1=lrt1.reconstructTensor();
	Tensor<double> t2=lrt2.reconstructTensor();
	Tensor<double> t2a=lrt2a.reconstructTensor();
	const Tensor<double> t3=t1(s);
	print("sliced assignment/1 ", (t3-t2).normf(),t2.normf());
	print("sliced assignment/1a", (t3-t2a).normf(),t2.normf());

	// sliced assignment
	LowRankTensor<double> lrt3(t0,eps,TT_2D);
	LowRankTensor<double> lrt4=lrt3(s);
	t1=lrt3.reconstructTensor();
	t2=lrt4.reconstructTensor();
	print("sliced assignment/2 ", (t1(s)-t2).normf(),t2.normf());

	// sliced assignment
	t0.reshape(k*k,k*k,k*k);
	LowRankTensor<double> lrt5(t0,eps,TT_3D);
	LowRankTensor<double> lrt6=lrt3(s);
	t1=lrt5.reconstructTensor();
	t2=lrt6.reconstructTensor();
	print("sliced assignment/3 ", (t1(s)-t2).normf()>1.e-3*eps,t2.normf());


	print("all done");
	return 0;
}

int testLowRankTensor_algebra(const long& k, const long& dim, const double& eps) {

	std::vector<Slice> s(3,Slice(0,k/2));

	Tensor<double> t0=Tensor<double>(k,k,k);
	Tensor<double> t1=Tensor<double>(k,k,k);
	t0.fillrandom();
	t1.fillindex();

	LowRankTensor<double> lrt1(t0,eps,TT_3D);
	LowRankTensor<double> lrt2(t1,eps,TT_3D);

	t0+=t1;
	lrt1+=lrt2;
	print(is_small((t0-lrt1.reconstructTensor()).normf(),eps),"inplace addition LRT   ",
			(t0-lrt1.reconstructTensor()).normf());

	FullTensor<double> f0(t0);
	FullTensor<double> f1(t1);
	t0+=t1;
	f0+=f1;
	print(is_small((t0-f0.fullTensor()).normf(),1.e-14), "inplace addition full   ", (t0-f0.fullTensor()).normf());

	// lhs slice testing
	LowRankTensor<double> lrt3(copy(t0),eps,TT_3D);
	LowRankTensor<double> lrt4=lrt3(s);
	lrt3.scale(-1.0);

	lrt3(s)+=lrt4;

	Tensor<double> p4=lrt3.reconstructTensor();
	Tensor<double> p4slice=p4(s);

	print(is_small(p4slice.normf(),1.e-13), "inplace addition slice LRT", p4slice.normf());

	print("all done\n");
	return 0;
}

int main(int argc, char**argv) {

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    srand(time(NULL));

    // the parameters
    const long k=5;
    const unsigned int dim=6;
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

    testFullTensor_ctor(k,dim,eps);
    testFullTensor_assignment(k,dim,eps);

    testLowRankTensor_ctor(k,dim,eps);
    testLowRankTensor_assignment(k,dim,eps);

    testLowRankTensor_algebra(k,dim,1.e-5);

    world.gop.fence();
    finalize();

    return 0;
}
