
/// \file testSepRep.cc
/// \brief test the SeparatedRepresentation (SepRep) for representing matrices

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "mra/seprep.h"
#include "mra/mra.h"

using namespace madness;

bool is_small(const double& val, const double& eps) {
	return (val<eps);
}

std::string ok(const bool b) {if (b) return "ok   "; return "fail ";};

bool is_large(const double& val, const double& eps) {
	return (val>eps);
}

int testGenTensor_ctor(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering ctor");
	Tensor<double> t0=Tensor<double>(3,3,3,3,3,3);
	t0.fillrandom();
//	t0.fillindex();
	double norm=0.0;
	int nerror=0;

	// default ctor
	GenTensor<double> t1;


	// ctor with rhs=Tensor, deep
	GenTensor<double> t2(t0,eps,tt);
	norm=(t0-t2.full_tensor_copy()).normf();
	print(ok(is_small(norm, eps)), "ctor with rhs=Tensor/1 ",t2.what_am_i(),t2.rank(),norm);
	if (!is_small(norm,eps)) nerror++;

	// test deepness
	t0.scale(2.0);
	norm=(t0-t2.full_tensor_copy()).normf();
	print(ok(is_large(norm, eps)), "ctor with rhs=Tensor/2 ",t2.what_am_i(),t2.rank(),norm);
	if (!is_large(norm,eps)) nerror++;

	{
		// ctor with rhs=GenTensor, shallow
		GenTensor<double> t3(t0,eps,tt);
		GenTensor<double> t4(t3);
		norm=(t3.full_tensor_copy()-t4.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"ctor with rhs=GenTensor/1",t4.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

		// test deepness
		t3.scale(2.0);
		norm=(t3.full_tensor_copy()-t4.full_tensor_copy()).normf();
		print(ok(is_small(norm, eps)), "ctor with rhs=GenTensor/2 ",t4.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;
	}

	{
		// deep ctor using copy()
		GenTensor<double> t3(t0,eps,tt);
		GenTensor<double> t4(copy(t3));
		norm=(t3.full_tensor_copy()-t4.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"ctor with rhs=GenTensor, using copy()",t4.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

		// test deepness
		t3.scale(3.0);
		norm=(t3.full_tensor_copy()-t4.full_tensor_copy()).normf();
		print(ok(is_large(norm, eps)), "ctor with rhs=GenTensor, using copy()",t4.what_am_i(),norm);
		if (!is_large(norm,eps)) nerror++;
	}


	print("all done\n");
	return nerror;
}


int testGenTensor_assignment(const long& k, const long& dim, const double& eps, const TensorType& tt) {


	print("entering assignment");
	Tensor<double> t0=Tensor<double>(3,3,3,3,3,3);
	Tensor<double> t1=Tensor<double>(3,3,3,3,3,3);
	std::vector<Slice> s(dim,Slice(0,1));
	t0.fillrandom();
	t1.fillindex();

//	t0.fillindex();
	double norm=0.0;
	int nerror=0;

	// default ctor
	const GenTensor<double> g0(t0,eps,tt);
	{
		GenTensor<double> g1(copy(g0));
		g1.scale(2.0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_large(norm,eps)),"pre-assignment check",g1.what_am_i(),norm);
		if (!is_large(norm,eps)) nerror++;

	}


	// regular assignment: g1=g0
	{

		GenTensor<double> g1(t1,eps,tt);
		g1=g0;
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"assignment with rhs=GenTensor/1",g1.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

		// test deepness
		g1.scale(2.0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"assignment with rhs=GenTensor/2",g1.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

	}

	// regular assignment w/ copy: g1=copy(g0)
	{

		GenTensor<double> g1(t1,eps,tt);
		g1=copy(g0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"copy assignment with rhs=GenTensor/1",g1.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

		// test deepness
		g1.scale(2.0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_large(norm,eps)),"copy assignment with rhs=GenTensor/2",g1.what_am_i(),norm);
		if (!is_large(norm,eps)) nerror++;

	}

	// regular assignment: g1=number
	{

	}


	// sliced assignment: g1=g0(s)
	{
		GenTensor<double> g1(t1,eps,tt);
		g1=g0;
		g1.scale(34.0);
		g1=g0(s);
		norm=(g0.full_tensor_copy()(s)-g1.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"sliced assignment with rhs=GenTensor/1",g1.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;


		// test deepness
		g1.scale(2.0);
		norm=(g0.full_tensor_copy()(s)-g1.full_tensor_copy()).normf();
		print(ok(is_large(norm,eps)),"sliced assignment with rhs=GenTensor/2",g1.what_am_i(),norm,
				"SHOULD FAIL FOR FULLRANK");
//		if (!is_large(norm,eps)) nerror++;

	}

	// sliced assignment: g1(s)=g0
	{

	}

	// sliced assignment: g1(s)=g0(s)
	{

	}

	// sliced assignment: g1(s)=number
	{

	}
	print("all done\n");
	return nerror;

}

int testGenTensor_algebra(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering algebra");
	Tensor<double> t0=Tensor<double>(k,k,k,k,k,k);
	Tensor<double> t1=Tensor<double>(k,k,k,k,k,k);

	std::vector<Slice> s(dim,Slice(0,k/2));
	std::vector<Slice> s2(dim,Slice(k/2+1,k-1));
	t0.fillrandom();
	t1.fillindex();

	Tensor<double> t2=copy(t0(s));

	double norm=0.0;
	int nerror=0;

	// default ctor
	GenTensor<double> g0(t0,eps,tt);
	GenTensor<double> g1(t1,eps,tt);

	// test inplace add: g0+=g1
	{
		GenTensor<double> g0(t0,eps,tt);
		g0+=g1;
		t0+=t1;
		norm=(g0.full_tensor_copy()-t0).normf();
		print(ok(is_small(norm,eps)),"algebra g0+=g1      ",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

	}


	// test inplace add: g0+=g1(s)
	{
		GenTensor<double> g0(t2,eps,tt);
		g0+=g1(s);
		t2+=t1(s);
		norm=(g0.full_tensor_copy()-t2).normf();
		print(ok(is_small(norm,eps)),"algebra g0+=g1(s)   ",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

	}


	// test inplace add: g0(s)+=g1(s)
	{
		GenTensor<double> g0(t0,eps,tt);
		g0(s)+=g1(s);
		t0(s)+=t1(s);
		norm=(g0.full_tensor_copy()-t0).normf();
		print(ok(is_small(norm,eps)),"algebra g0(s)+=g1(s)",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

	}



	print("all done\n");
	return nerror;
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
	Tensor<double> t11=t2(s).fullTensor();

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

//int testLowRankTensor_ctor(const long& k, const long& dim, const double& eps) {
//
//	const double tol=1.e-13;
//	const double safety=1.e-3;
//
//	Tensor<double> t0=Tensor<double>(7,7,7);
//	Tensor<double> t1=Tensor<double>(7,7,7,7);
//	t0.fillindex();
//	std::vector<Slice> s(4,Slice(0,5));
//	t1(s).fillrandom();
//	// t1 is now sort of low-rank
//
//	// default ctor
//	LowRankTensor<double> lrt1();
//
//	// ctor with a regular tensor
//	LowRankTensor<double> lrt2(t0,eps,TT_3D);
//	LowRankTensor<double> lrt3(t1,eps,TT_2D);
//	double norm=(lrt2.reconstructTensor()-t0).normf();
//	print(is_small(norm,eps*safety), "LRT, 3-way ctor   ", norm, lrt2.rank());
//	norm=(lrt3.reconstructTensor()-t1).normf();
//	print(is_small(norm,eps*safety), "LRT, 2-way ctor   ", norm, lrt3.rank());
//
//	// copy ctor
//	LowRankTensor<double> lrt4(lrt2);
//	norm=(lrt2.reconstructTensor()-lrt4.reconstructTensor()).normf();
//	print(is_small(norm,tol), "copy ctor         ", norm);
//
//	/// virtual ctor should be shallow
//	SepRepTensor<double>* sr;
//	sr=lrt4.clone();
//	lrt4.scale(0.5);
//	norm=(lrt4.reconstructTensor()-sr->reconstructTensor()).normf();
//	print(is_small(norm,tol), "LRT, virtual ctor ",norm);
//
//
//
//
//	print("all done");
//	return 0;
//}

//int testLowRankTensor_assignment(const long& k, const long& dim, const double& eps) {
//
//	std::vector<long> d(6,k);
//	Tensor<double> t0=Tensor<double>(d);
//	t0.fillindex();
//	std::vector<Slice> s(6,Slice(0,k/2));
//
//	// regular assignment
//	LowRankTensor<double> lrt0(t0,eps,TT_3D);
//	LowRankTensor<double> lrt0a;
//	lrt0a=lrt0;
//	print("LRT assignment      ", (lrt0.reconstructTensor()-lrt0a.reconstructTensor()).normf());
//
//
//	// sliced assignment
//	LowRankTensor<double> lrt1(t0,eps,TT_3D);
//	LowRankTensor<double> lrt2=lrt1(s);
//	LowRankTensor<double> lrt2a=lrt1;
//	lrt2a=lrt1(s);
//
//	Tensor<double> t1=lrt1.reconstructTensor();
//	Tensor<double> t2=lrt2.reconstructTensor();
//	Tensor<double> t2a=lrt2a.reconstructTensor();
//	const Tensor<double> t3=t1(s);
//	print("sliced assignment/1 ", (t3-t2).normf(),t2.normf());
//	print("sliced assignment/1a", (t3-t2a).normf(),t2.normf());
//
//	// sliced assignment
//	LowRankTensor<double> lrt3(t0,eps,TT_2D);
//	LowRankTensor<double> lrt4=lrt3(s);
//	t1=lrt3.reconstructTensor();
//	t2=lrt4.reconstructTensor();
//	print("sliced assignment/2 ", (t1(s)-t2).normf(),t2.normf());
//
//	// sliced assignment
//	t0.reshape(k*k,k*k,k*k);
//	LowRankTensor<double> lrt5(t0,eps,TT_3D);
//	LowRankTensor<double> lrt6=lrt3(s);
//	t1=lrt5.reconstructTensor();
//	t2=lrt6.reconstructTensor();
//	print("sliced assignment/3 ", (t1(s)-t2).normf()>1.e-3*eps,t2.normf());
//
//
//	print("all done");
//	return 0;
//}


int main(int argc, char**argv) {

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    srand(time(NULL));

    // the parameters
    const long k=3;
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

    int error=0;
    print("hello world");

#if 0
    {
		Tensor<double> t0=Tensor<double>(3,3);
		Tensor<double> t1=Tensor<double>(2,2);
		std::vector<Slice> s(2,Slice(0,1));
		t0=2.0;
		t1=1.0;

		t1=t0(s);
		t0.scale(2.0);

		print(t1);
		print(t0);
    }

    {
		Tensor<double> t0=Tensor<double>(3,3);
		Tensor<double> t1=Tensor<double>(2,2);
		std::vector<Slice> s(2,Slice(0,1));
		t0=2.0;
		t1=1.0;

		t1(s)=t0(s);
		t0.scale(2.0);

		print(t1);
		print(t0);
    }
#endif

    error+=testGenTensor_ctor(k,dim,eps,TT_FULL);
    error+=testGenTensor_ctor(k,dim,eps,TT_3D);
    error+=testGenTensor_ctor(k,dim,eps,TT_2D);

    error+=testGenTensor_assignment(k,dim,eps,TT_FULL);
    error+=testGenTensor_assignment(k,dim,eps,TT_3D);
    error+=testGenTensor_assignment(k,dim,eps,TT_2D);

    error+=testGenTensor_algebra(k,dim,eps,TT_FULL);
    error+=testGenTensor_algebra(k,dim,eps,TT_3D);
    error+=testGenTensor_algebra(k,dim,eps,TT_2D);


    print(ok(error==0),error,"finished test suite\n");


//    testFullTensor_ctor(k,dim,eps);
//    testFullTensor_assignment(k,dim,eps);
//
//    testLowRankTensor_ctor(k,dim,eps);
//    testLowRankTensor_assignment(k,dim,eps);
//
//    testLowRankTensor_algebra(k,dim,1.e-5);

    world.gop.fence();
    finalize();

    return 0;
}
