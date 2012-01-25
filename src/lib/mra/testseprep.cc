
/// \file testSepRep.cc
/// \brief test the SeparatedRepresentation (SepRep) for representing matrices

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "tensor/gentensor.h"

using namespace madness;

bool is_small(const double& val, const double& eps) {
	return (val<eps);
}

std::string ok(const bool b) {if (b) return "ok   "; return "fail ";};

bool is_large(const double& val, const double& eps) {
	return (val>eps);
}
#if HAVE_GENTENSOR

int testGenTensor_ctor(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering ctor");

	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[2]=Tensor<double>(d).fillrandom();

	double norm=0.0;
	int nerror=0;

	// ctor with rhs=Tensor, deep
	for (int i=0; i<3; i++) {
		Tensor<double> t0=copy(t[i]);

		GenTensor<double> g0(t0,eps,tt);
		norm=(t0-g0.full_tensor_copy()).normf();
		print(ok(is_small(norm, eps)), "ctor with rhs=Tensor/1 ",g0.what_am_i(),g0.rank(),norm);
		if (!is_small(norm,eps)) nerror++;

		// test deepness
		t0+=(2.0);
		t0.scale(2.0);
		norm=(t0-g0.full_tensor_copy()).normf();
		print(ok(is_large(norm, eps)), "ctor with rhs=Tensor/2 ",g0.what_am_i(),g0.rank(),norm);
		if (!is_large(norm,eps)) nerror++;
	}

	// ctor with rhs=GenTensor, shallow
	for (int i=0; i<3; i++) {

		Tensor<double> t0=copy(t[i]);
		GenTensor<double> g0(t0,eps,tt);
		GenTensor<double> g1(g0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"ctor with rhs=GenTensor/1 ",g1.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

		// test deepness
		g0.scale(2.0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_small(norm, eps)), "ctor with rhs=GenTensor/2 ",g1.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;
	}

	// deep ctor using copy()
	for (int i=0; i<3; i++) {
		Tensor<double> t0=copy(t[i]);
		GenTensor<double> g0(t0,eps,tt);
		GenTensor<double> g1(copy(g0));
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"ctor with rhs=GenTensor, using copy()",g1.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

		// test deepness
		g0.scale(3.0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		double norm2=g0.normf();	// if t0 is empty
		if (norm2<eps) norm=1.0;
		print(ok(is_large(norm, eps)), "ctor with rhs=GenTensor, using copy()",g1.what_am_i(),norm);
		if (!is_large(norm,eps)) nerror++;
	}


	print("all done\n");
	return nerror;
}

int testGenTensor_assignment(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering assignment");

	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[2]=Tensor<double>(d).fillrandom();

	std::vector<Slice> s(dim,Slice(0,1));

	double norm=0.0;
	int nerror=0;

	// default ctor
	for (int i=0; i<3; i++) {
		Tensor<double> t0=t[i];

		GenTensor<double> g0(t0,eps,tt);
		GenTensor<double> g1(copy(g0));
		g1.scale(2.0);
		norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
		if (t0.normf()>eps) {
			print(ok(is_large(norm,eps)),"pre-assignment check",g1.what_am_i(),norm);
			if (!is_large(norm,eps)) nerror++;
		}
	}


	// regular assignment: g1=g0
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=t[i];
			Tensor<double> t1=t[j];

			GenTensor<double> g0(t0,eps,tt);
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
	}

	// regular assignment w/ copy: g1=copy(g0)
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=t[i];
			Tensor<double> t1=t[j];

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);
			g1=copy(g0);
			norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
			print(ok(is_small(norm,eps)),"copy assignment with rhs=GenTensor/1",g1.what_am_i(),norm);
			if (!is_small(norm,eps)) nerror++;

			// test deepness
			g1.scale(2.0);
			norm=(g0.full_tensor_copy()-g1.full_tensor_copy()).normf();
			// if t0 is zero
			if (t0.normf()>eps) {
				print(ok(is_large(norm,eps)),"copy assignment with rhs=GenTensor/2",g1.what_am_i(),norm);
				if (!is_large(norm,eps)) nerror++;
			}
		}
	}

	// regular assignment: g1=number
	{

	}


	// sliced assignment: g1=g0(s)
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=t[i];
			Tensor<double> t1=t[j];

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g1=g0;
			g1.scale(34.0);
			g1=g0(s);
			norm=(g0.full_tensor_copy()(s)-g1.full_tensor_copy()).normf();
			print(ok(is_small(norm,eps)),"sliced assignment g1=g0(s) /1",g1.what_am_i(),norm);
			if (!is_small(norm,eps)) nerror++;


			// test deepness
			g1.scale(2.0);
			norm=(g0.full_tensor_copy()(s)-g1.full_tensor_copy()).normf();
			if (t0.normf()>eps) {
				print(ok(is_large(norm,eps)),"sliced assignment g1=g0(s) /2",g1.what_am_i(),norm);
				if (!is_large(norm,eps)) nerror++;
			}
		}
	}

	// sliced assignment: g1(s)=g0
	{
	}




	// sliced assignment: g1(s)=g0(s)
	{
	}

	// sliced assignment: g1(s)=number
	for (int i=0; i<3; i++) {

		Tensor<double> t0=t[i];

		GenTensor<double> g0(t0,eps,tt);
		Tensor<double> t2=copy(t0);
		g0(s)=0.0;
		t2(s)=0.0;
		norm=(t2-g0.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"sliced assignment with rhs=0.0",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;
	}
	print("all done\n");
	return nerror;

}

int testGenTensor_algebra(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering algebra");

	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[2]=Tensor<double>(d).fillrandom();

	std::vector<Slice> s(dim,Slice(0,1));

//	Tensor<double> t2=copy(t0(s));

	double norm=0.0;
	int nerror=0;

	// test inplace add: g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1;
			t0+=t1;
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"algebra g0+=g1      ",g0.what_am_i(),norm);
			if (!is_small(norm,eps)) nerror++;
		}
	}


	// test inplace add: g0+=g1(s)
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i](s));
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1(s);
			t0+=t1(s);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"algebra g0+=g1(s)   ",g0.what_am_i(),norm);
			if (!is_small(norm,eps)) nerror++;
		}
	}


	// test inplace add: g0(s)+=g1(s)
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0(s)+=g1(s);
			t0(s)+=t1(s);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"algebra g0(s)+=g1(s)",g0.what_am_i(),norm,g0.rank());
			if (!is_small(norm,eps)) nerror++;
		}
	}

//	// test inplace scale: g=g0*=fac
//	{
//		GenTensor<double> g0(t0,eps,tt);
//		GenTensor<double> g2=g0.scale(2.0);
//		Tensor<double> t2=t0.scale(2.0);
//		norm=(g0.full_tensor_copy()-t0).normf();
//		print(ok(is_small(norm,eps)),"algebra scale",g0.what_am_i(),norm);
//		if (!is_small(norm,eps)) nerror++;
//
//	}



	print("all done\n");
	return nerror;
}

int testGenTensor_rankreduce(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering rank reduce");

	double norm=0.0;
	int nerror=0;


	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[1].scale(1.0/t[1].normf());
	t[2]=Tensor<double>(d).fillrandom();

	// test rank reduction g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1;
			t0+=t1;
			g0.config().orthonormalize(eps);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"rank reduction orthonormalize   ",g0.what_am_i(),norm,g0.rank());
			if (!is_small(norm,eps)) nerror++;
		}
	}


	// test rank reduction g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0.config().orthonormalize(eps*0.5/std::max(1.0,g1.normf()));
			g1.config().orthonormalize(eps*0.5/std::max(1.0,g0.normf()));
			g0.add_SVD(g1,eps);
			t0+=t1;
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"add SVD   ",g0.what_am_i(),norm,g0.rank());
			if (!is_small(norm,eps)) nerror++;
		}
	}

	// test rank reduction g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1;
			t0+=t1;
			g0.reconstruct_and_decompose(eps/1.e-3);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"rank reduction reconstruct   ",g0.what_am_i(),norm,g0.rank());
			if (!is_small(norm,eps)) nerror++;
		}
	}


#if 0
	// test rank reduction g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);
//			t1.scale(-1.0);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1;
			t0+=t1;
			g0.config().right_orthonormalize(eps);
			g0.config().right_orthonormalize(eps);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"rank reduction right_orthonormalize   ",g0.what_am_i(),norm,g0.rank(),t0.normf());
			if (!is_small(norm,eps)) nerror++;
		}
	}
#endif

#if 0
	// test rank reduction g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1;
			t0+=t1;
			g0.config().undo_structure();
			g0.config().rank_revealing_modified_gram_schmidt2(eps);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"rank reduction RR/MGS2   ",g0.what_am_i(),norm,g0.rank());
			if (!is_small(norm,eps)) nerror++;
		}
	}
#endif

	// test rank reduction g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1;
			t0+=t1;
			g0.config().divide_and_conquer_reduce(eps);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"rank reduction divide&conquer",OrthoMethod(),
					g0.what_am_i(),norm,g0.rank());
			if (!is_small(norm,eps)) nerror++;
		}
	}




	// test rank reduction g0+=g1
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			Tensor<double> t0=copy(t[i]);
			Tensor<double> t1=copy(t[j]);

			GenTensor<double> g0(t0,eps,tt);
			GenTensor<double> g1(t1,eps,tt);

			g0+=g1;
			t0+=t1;
			g0.reduceRank(eps);
			norm=(g0.full_tensor_copy()-t0).normf();
			print(ok(is_small(norm,eps)),"rank reduction reduceRank   ",OrthoMethod(),g0.what_am_i(),norm);
			if (!is_small(norm,eps)) nerror++;
		}
	}

	print("all done\n");
	return nerror;
}

//int testGenTensor_transform(const long& k, const long& dim, const double& eps, const TensorType& tt) {
//
//	print("entering transform");
//	Tensor<double> t0=Tensor<double>(k,k,k,k,k,k);
//	Tensor<double> c=Tensor<double> (k,k);
//	Tensor<double> cc[dim];
//	for (unsigned int idim=0; idim<dim; idim++) {
//		cc[idim]=Tensor<double>(k,k);
//		cc[idim].fillrandom();
//	}
//
//	t0.fillrandom();
//	c.fillindex();
//	c.scale(1.0/c.normf());
//
//	double norm=0.0;
//	int nerror=0;
//
//	// default ctor
//	GenTensor<double> g0(t0,eps,tt);
//
//	// test transform_dir
//	{
//		const long ndim=t0.ndim();
//
//		for (long idim=0; idim<ndim; idim++) {
////		for (long idim=0; idim<1; idim++) {
//			GenTensor<double> g1=transform_dir(g0,c,idim);
//			Tensor<double> t1=transform_dir(t0,c,idim);
//			norm=(g1.full_tensor_copy()-t1).normf();
//			print(ok(is_small(norm,eps)),"transform_dir",idim,g0.what_am_i(),norm);
//			if (!is_small(norm,eps)) nerror++;
//		}
//	}
//
//	// test transform with tensor
//	{
//		GenTensor<double> g1=transform(g0,c);
//		Tensor<double> t1=transform(t0,c);
//		norm=(g1.full_tensor_copy()-t1).normf();
//		print(ok(is_small(norm,eps)),"transform.scale",g0.what_am_i(),norm);
//		if (!is_small(norm,eps)) nerror++;
//
//	}
//
//
//	// test general_transform
//	{
//		GenTensor<double> g1=general_transform(g0,cc);
//		Tensor<double> t1=general_transform(t0,cc);
//		norm=(g1.full_tensor_copy()-t1).normf();
//		print(ok(is_small(norm,eps)),"general_transform",g0.what_am_i(),norm);
//		if (!is_small(norm,eps)) nerror++;
//
//	}
//
//
//
//	// test general_transform with scale
//	{
//		GenTensor<double> g1=general_transform(g0,cc).scale(1.9);
//		Tensor<double> t1=general_transform(t0,cc).scale(1.9);
//		norm=(g1.full_tensor_copy()-t1).normf();
//		print(ok(is_small(norm,eps)),"general_transform.scale",g0.what_am_i(),norm);
//		if (!is_small(norm,eps)) nerror++;
//
//	}
//
//
//	print("all done\n");
//	return nerror;
//
//}

int testGenTensor_reconstruct(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering reconstruct");
	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[2]=Tensor<double>(d).fillrandom();

	double norm=0.0;
	int nerror=0;

	// reconstruct
	for (int i=0; i<3; i++) {
		Tensor<double> t0=copy(t[i]);

		GenTensor<double> g0(t0,eps,tt);
		Tensor<double> t=g0.full_tensor_copy();
		norm=(t-t0).normf();
		print(ok(is_small(norm,eps)),"reconstruct",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;
	}

	// reconstruct
//	for (int i=0; i<3; i++) {
//		Tensor<double> t0=copy(t[i]);
//		GenTensor<double> g0(t0,eps,tt);
//
//		Tensor<double> t=g0.full_tensor_copy();
//		g0.accumulate_into(t,-1.0);
//
//		norm=(t).normf();
//		print(ok(is_small(norm,eps)),"accumulate_into tensor",g0.what_am_i(),norm);
//		if (!is_small(norm,eps)) nerror++;
//	}

#if 0
	// reconstruct
	for (int i=0; i<3; i++) {
		Tensor<double> t0=copy(t[i]);
		GenTensor<double> g0(t0,eps,tt);

		GenTensor<double> g=copy(g0);
		g.right_orthonormalize(eps);
		g0.accumulate_into(g,eps,-1.0);

		norm=(g.full_tensor_copy()).normf();
		print(ok(is_small(norm,eps)),"accumulate_into gentensor",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;
	}
#endif

	print("all done\n");
	return nerror;

}



int main(int argc, char**argv) {

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    srand(time(NULL));

    // the parameters
    long k=5;
    const unsigned int dim=6;
    double eps=1.e-3;
    print("k    ",k);
    print("eps  ",eps);

    // some test
#if 0

	int nloop=3*30;
	print("nloop",nloop);
	long k3=2*k*2*k*2*k;
	Tensor<double> s(k3,k3);
	Tensor<double> v, e;
	Tensor<double> ss(2*k,k3*2*k*2*k);
	Tensor<double> t(2*k,2*k);
	s.fillrandom();
	s=inner(s,s,1,1);
    double cpu0, cpu1;
    cpu0=wall_time();

    for (int i=0; i<nloop; i++) {
//        syev(s,v,e);
//        cholesky(s);
//        v=inner(s,s,1,1);
        v=inner(t,ss);
    }

	cpu1=wall_time();
	if(world.rank() == 0) print("chunk Brand at time ", cpu1-cpu0);
	return 0;

#endif



#if 0

    // do some benchmarking
	long nloop=1.e5;

    print("benchmarking with k ",k);
    print("dimensions          ",dim);
    print("nloop               ",nloop);

   	for (int R=4; R<16; R+=4) {
   	    for (int r=4; r<16; r+=4) {

   	    	// do some prep
   	    	print("r",r);
    		print("R",R);
    		double cpu0, cpu1;
    		GenTensor<double> g0(TT_2D,k,dim);
    		GenTensor<double> g1(TT_2D,k,dim);
    		GenTensor<double> g2(TT_2D,k,dim);
    		Tensor<double> t0(k,k,k,k,k,k);

    		g0.fillrandom(r);
    		g1.fillrandom(R);
    		g2.fillrandom(r+R);

    		// start benchmark
    		cpu0=wall_time();
    		for (int i=0; i<nloop; i++) {
    			GenTensor<double> g00=g0;
    			GenTensor<double> g11=g1;
    		}

    		cpu1=wall_time();
    		if(world.rank() == 0) print("baseline at time ", cpu1-cpu0);
//    		cpu0=wall_time();
//
//    		for (int i=0; i<nloop; i++) {
//    			GenTensor<double> g00=g0;
//    			g00.accumulate_into(t0);
//    		}
//
//    		cpu1=wall_time();
//    		if(world.rank() == 0) print("accumulate_into tensor at time ", cpu1-cpu0);
//    		cpu0=wall_time();
//
////			GenTensor<double> g00=g0;
////   			GenTensor<double> g11=g1;
//   			g1.config().right_orthonormalize();
//    		for (int i=0; i<nloop; i++) {
//    			GenTensor<double> g00=g0;
//    			GenTensor<double> g11=g1;
//    			g00.accumulate_into(g11);
//    		}
//
//    		cpu1=wall_time();
//    		if(world.rank() == 0) print("accumulate_into gentensor at time ", cpu1-cpu0);


    		cpu0=wall_time();
    		for (int i=0; i<nloop; i++) {
    			GenTensor<double> g00=g2;
//    			g00.config().rank_revealing_modified_gram_schmidt2(eps);
    		}
    		cpu1=wall_time();
    		if(world.rank() == 0) print("RR/MGS at time ", cpu1-cpu0);


    		cpu0=wall_time();
    		for (int i=0; i<nloop; i++) {
    			GenTensor<double> g00=g2;
    			g00.config().right_orthonormalize(eps);
    		}
    		cpu1=wall_time();
    		if(world.rank() == 0) print("right_orthonormalize at time ", cpu1-cpu0);

    		cpu0=wall_time();
    		for (int i=0; i<nloop; i++) {
    			GenTensor<double> g00=g2;
    			g00.config().sequential_orthogonalization(eps);
    		}
    		cpu1=wall_time();
    		if(world.rank() == 0) print("RR/LRT add at time ", cpu1-cpu0);

    		cpu0=wall_time();
    		for (int i=0; i<nloop; i++) {
    			GenTensor<double> g00=g2;
    			g00.config().orthonormalize(eps);
    		}
    		cpu1=wall_time();
    		if(world.rank() == 0) print("orthonormalize at time ", cpu1-cpu0);

   	    }
    }



#if 0

	Tensor<double> t0=Tensor<double>(k,k,k,k,k,k);
	Tensor<double> t1=Tensor<double>(k,k,k,k,k,k);

	t0.fillindex();
	t0=2.5;
	t1=2.3;
	GenTensor<double> g1(t1,eps,TT_2D);
	g1.fillrandom();
	SepRep<double> sr0(TT_2D,k,dim);
	SepRep<double> sr1(TT_2D,k,dim);
	double tim=wall_time();

    if(world.rank() == 0) print("starting at time", wall_time()-tim);
    for (unsigned int i=0; i<nloop; i++) {
//		GenTensor<double> g0(t0,eps,TT_2D);
//    	g1.fillrandom(4);
    	sr1.fillrandom(4);

    }
    if(world.rank() == 0) print("baseline at time  ", wall_time()-tim);
    double base1=wall_time()-tim;
    tim=wall_time();

    if(world.rank() == 0) print("starting at time", wall_time()-tim);
    for (unsigned int i=0; i<nloop; i++) {
//		GenTensor<double> g0(t0,eps,TT_2D);
      	sr1.fillrandom(8);
  }
    if(world.rank() == 0) print("baseline at time  ", wall_time()-tim);
    double base2=wall_time()-tim;
    tim=wall_time();

    for (unsigned int i=0; i<nloop; i++) {
    	sr1.fillrandom(4);
    	sr1.config().ortho3();
//		g0.finalize_accumulate();
    }
    if(world.rank() == 0) print("ortho3/4 at time  ", wall_time()-tim-base1);
    tim=wall_time();

    for (unsigned int i=0; i<nloop; i++) {
    	sr1.fillrandom(8);
    	sr1.config().ortho3();
//		g0.finalize_accumulate();
    }
    if(world.rank() == 0) print("ortho3/8 at time  ", wall_time()-tim-base2);
    tim=wall_time();

//    for (unsigned int i=0; i<nloop; i++) {
//    	sr1.fillrandom(16);
//    	sr1.config().ortho3();
////		g0.finalize_accumulate();
//    }
//    if(world.rank() == 0) print("ortho3/16 at time  ", wall_time()-tim-base2);
//    tim=wall_time();


    sr0.fillrandom(10);
	sr1.fillrandom(1);
    for (unsigned int i=0; i<nloop; i++) {
		sr0.update_by(sr1);
    }
    if(world.rank() == 0) print("update_by/1 at time  ", wall_time()-tim);
    tim=wall_time();

	sr1.fillrandom(4);
    for (unsigned int i=0; i<nloop; i++) {
		sr0.update_by(sr1);
    }
    if(world.rank() == 0) print("update_by/4 at time  ", wall_time()-tim);
    tim=wall_time();

    GenTensor<double> g0(t0,eps,TT_2D);
    for (int i=0; i<8; i++) g1+=g1;
    for (unsigned int i=0; i<nloop; i++) {
//    	g1.fillrandom();
		g1.accumulate_into(t0,1.0);
    }
    if(world.rank() == 0) print("accumulate_into at time  ", wall_time()-tim);
    tim=wall_time();
#endif

//    GenTensor<double> g2(t0,eps,TT_2D);
//    GenTensor<double> g3(t1,eps,TT_2D);
//    for (unsigned int i=0; i<nloop; i++) {
//    	g1.fillrandom();
//    	g2+=g3;
//    }
//    if(world.rank() == 0) print("inplace_add at time  ", wall_time()-tim);
//    tim=wall_time();

//    MADNESS_EXCEPTION("end benchmark",0);
#endif

    int error=0;
    print("hello world");
#if 1
    error+=testGenTensor_ctor(k,dim,eps,TT_FULL);
//    error+=testGenTensor_ctor(k,dim,eps,TT_3D);
    error+=testGenTensor_ctor(k,dim,eps,TT_2D);

    error+=testGenTensor_assignment(k,dim,eps,TT_FULL);
//    error+=testGenTensor_assignment(k,dim,eps,TT_3D);
    error+=testGenTensor_assignment(k,dim,eps,TT_2D);

    error+=testGenTensor_algebra(k,dim,eps,TT_FULL);
//    error+=testGenTensor_algebra(k,dim,eps,TT_3D);
    error+=testGenTensor_algebra(k,dim,eps,TT_2D);

    error+=testGenTensor_rankreduce(k,dim,eps,TT_FULL);
//    error+=testGenTensor_rankreduce(k,dim,eps,TT_3D);
    error+=testGenTensor_rankreduce(k,dim,eps,TT_2D);

//    error+=testGenTensor_transform(k,dim,eps,TT_FULL);
//    error+=testGenTensor_transform(k,dim,eps,TT_3D);
//    error+=testGenTensor_transform(k,dim,eps,TT_2D);

    error+=testGenTensor_reconstruct(k,dim,eps,TT_FULL);
//    error+=testGenTensor_reconstruct(k,dim,eps,TT_3D);
    error+=testGenTensor_reconstruct(k,dim,eps,TT_2D);

    print(ok(error==0),error,"finished test suite\n");
#endif

    world.gop.fence();
    finalize();

    return 0;
}

#else
int main(int argc, char** argv) {

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    print("no testseprep without having a GenTensor");
    world.gop.fence();
    finalize();
    return 0;
}

#endif


