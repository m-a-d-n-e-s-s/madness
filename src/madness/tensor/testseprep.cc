
/// \file testSepRep.cc
/// \brief test the SeparatedRepresentation (SepRep) for representing matrices

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <madness/tensor/gentensor.h>

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

	// ctor with slices
    for (int i=0; i<3; i++) {

        std::vector<Slice> s(dim,Slice(0,1));

        Tensor<double> t0=copy(t[i]);
        Tensor<double> t1=copy(t0(s));
        GenTensor<double> g0(t0,eps,tt);
        GenTensor<double> g1=g0(s);
        norm=(g1.full_tensor_copy()-t1).normf();
        print(ok(is_small(norm,eps)),"ctor with rhs=SliceGenTensor",g1.what_am_i(),norm);
        if (!is_small(norm,eps)) nerror++;

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

	// assign to a number
    // sliced assignment: g1(s)=number
    for (int i=0; i<3; i++) {

        Tensor<double> t0=copy(t[i]);
        GenTensor<double> g0(t0,eps,tt);
        g0=2.0;
        t0=2.0;
        norm=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(norm,eps)),"assignment with a number",g0.what_am_i(),norm);
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
	double t1norm=t[1].normf();
    double t2norm=t[2].normf();
    t[1].scale(1./t1norm)*2.1375;
    t[2].scale(1./t2norm)*137.72;

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

	// test inplace scale: g=g0*=fac
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {

            Tensor<double> t0=copy(t[i]);
            Tensor<double> t1=copy(t[j]);
            GenTensor<double> g0(t0,eps,tt);
            GenTensor<double> g2=g0.scale(2.0);
            Tensor<double> t2=t0.scale(2.0);
            norm=(g0.full_tensor_copy()-t0).normf();
            print(ok(is_small(norm,eps)),"algebra scale",g0.what_am_i(),norm);
            if (!is_small(norm,eps)) nerror++;
        }
	}

    // test elementwise multiplication emul:
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {

            Tensor<double> t0=copy(t[i]);
            Tensor<double> t1=copy(t[j]);
            GenTensor<double> g0(t0,eps,tt);
            GenTensor<double> g1(t1,eps,tt);
            g1.emul(g0);
            t1.emul(t0);
            norm=(g1.full_tensor_copy()-t1).normf();
            print(ok(is_small(norm,eps)),"algebra emul",g1.what_am_i(),norm);
            if (!is_small(norm,eps)) nerror++;
        }
    }




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
            g0.reduce_rank(eps);
            norm=(g0.full_tensor_copy()-t0).normf();
            print(ok(is_small(norm,eps)),"rank reduction reduceRank   ",g0.what_am_i(),norm);
            if (!is_small(norm,eps)) nerror++;
        }
    }

    // the other tests make only sense for TT_2D

    if (tt==TT_2D) {
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
    }

	print("all done\n");
	return nerror;
}

int testGenTensor_transform(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering transform");
	Tensor<double> t0=Tensor<double>(k,k,k,k,k,k);
	Tensor<double> c=Tensor<double> (k,k);
	Tensor<double> cc[TENSOR_MAXDIM];
	for (unsigned int idim=0; idim<dim; idim++) {
		cc[idim]=Tensor<double>(k,k);
		cc[idim].fillrandom();
	}

	t0.fillrandom();
	c.fillindex();
	c.scale(1.0/c.normf());

	double norm=0.0;
	int nerror=0;

	// default ctor
	GenTensor<double> g0(t0,eps,tt);

	// test transform_dir
	{
		const long ndim=t0.ndim();

		for (long idim=0; idim<ndim; idim++) {
//		for (long idim=0; idim<1; idim++) {
			GenTensor<double> g1=transform_dir(g0,c,idim);
			Tensor<double> t1=transform_dir(t0,c,idim);
			norm=(g1.full_tensor_copy()-t1).normf();
			print(ok(is_small(norm,eps)),"transform_dir",idim,g0.what_am_i(),norm);
			if (!is_small(norm,eps)) nerror++;
		}
	}

	// test transform with tensor
	{
		GenTensor<double> g1=transform(g0,c);
		Tensor<double> t1=transform(t0,c);
		norm=(g1.full_tensor_copy()-t1).normf();
		print(ok(is_small(norm,eps)),"transform.scale",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

	}


	// test general_transform
	{
		GenTensor<double> g1=general_transform(g0,cc);
		Tensor<double> t1=general_transform(t0,cc);
		norm=(g1.full_tensor_copy()-t1).normf();
		print(ok(is_small(norm,eps)),"general_transform",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

	}



	// test general_transform with scale
	{
		GenTensor<double> g1=general_transform(g0,cc).scale(1.9);
		Tensor<double> t1=general_transform(t0,cc).scale(1.9);
		norm=(g1.full_tensor_copy()-t1).normf();
		print(ok(is_small(norm,eps)),"general_transform.scale",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;

	}


	print("all done\n");
	return nerror;

}

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

int testGenTensor_deepcopy(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering deep copy");
	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[2]=Tensor<double>(d).fillrandom();

	std::vector<Slice> s0(dim,Slice(0,k-2));

	double norm=0.0;
	int nerror=0;

	// reconstruct
	for (int i=0; i<3; i++) {
		Tensor<double> t0=copy(t[i]);

		GenTensor<double> g0(t0,eps,tt);
		GenTensor<double> g1=copy(g0);
		GenTensor<double> g2=copy(g0(s0));
		norm=0.0;
		print(ok(is_small(norm,eps)),"deep copy",g0.what_am_i(),norm);
		if (!is_small(norm,eps)) nerror++;
	}

	print("all done\n");
	return nerror;

}

/// test the reduce algorithm for adding many tensors
int testGenTensor_reduce(const long& k, const long& dim, const double& eps, const TensorType& tt) {

	print("entering reduce");
	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);
	std::list<GenTensor<double> > addends;

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[1].scale(1.0/t[1].normf());
	t[2]=Tensor<double>(d).fillrandom();
	t[2].scale(1.0/t[2].normf());

	double norm=0.0;
	int nerror=0;

	// reconstruct
	for (int i=0; i<3; i++) {
		addends.clear();

		addends.push_back(GenTensor<double>(t[0],eps,tt));
		addends.push_back(GenTensor<double>(t[1],eps,tt));
		addends.push_back(GenTensor<double>(t[2],eps,tt));

		GenTensor<double> result=reduce(addends,eps,true);
		Tensor<double> tresult=t[0]+t[1]+t[2];

		norm=(tresult-result.full_tensor_copy()).normf();

		print(ok(is_small(norm,eps)),"reduce",tt,norm);
		if (!is_small(norm,eps)) nerror++;
	}

	print("all done\n");
	return nerror;

}

/// test conversion from one gentensor representation to another
int testGenTensor_convert(const long k, const long dim, const TensorArgs targs) {

    print("entering testGenTensor convert");
    int nerror=0;
    double eps=targs.thresh;

    std::vector<long> d(dim,k);
    Tensor<double> t0(d),t1(d),t2(d);

    // set up a tensor of full rank
    t0.fillrandom();

    // set up a tensor of low rank
    t1.fillindex();
    t1.scale(1.0/t1.normf());

    {
        GenTensor<double> g0(t0,targs);

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        double error1=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t0.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_FULL));
        error1=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_FULL; k=",k,"dim=",dim,"error=",error1,"size=",t0.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_TENSORTRAIN));
        error1=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_TENSORTRAIN; k=",k,"dim=",dim,"error=",error1,"size=",t0.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        error1=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t0.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_TENSORTRAIN));
        error1=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_TENSORTRAIN; k=",k,"dim=",dim,"error=",error1,"size=",t0.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_FULL));
        error1=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_FULL; k=",k,"dim=",dim,"error=",error1,"size=",t0.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        error1=(t0-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t0.size());

    }

    // set up a tensor of low rank
    {
        GenTensor<double> g0(t1,targs);

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        double error1=(t1-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t1.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_FULL));
        error1=(t1-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_FULL; k=",k,"dim=",dim,"error=",error1,"size=",t1.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_TENSORTRAIN));
        error1=(t1-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_TENSORTRAIN; k=",k,"dim=",dim,"error=",error1,"size=",t1.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        error1=(t1-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t1.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_TENSORTRAIN));
        error1=(t1-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_TENSORTRAIN; k=",k,"dim=",dim,"error=",error1,"size=",t1.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_FULL));
        error1=(t1-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_FULL; k=",k,"dim=",dim,"error=",error1,"size=",t1.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        error1=(t1-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t1.size());

    }

    // set up a tensor of zero rank
    {
        GenTensor<double> g0(t2,targs);

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        double error1=(t2-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t2.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_FULL));
        error1=(t2-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_FULL; k=",k,"dim=",dim,"error=",error1,"size=",t2.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_TENSORTRAIN));
        error1=(t2-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_TENSORTRAIN; k=",k,"dim=",dim,"error=",error1,"size=",t2.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        error1=(t2-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t2.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_TENSORTRAIN));
        error1=(t2-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_TENSORTRAIN; k=",k,"dim=",dim,"error=",error1,"size=",t2.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_FULL));
        error1=(t2-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_FULL; k=",k,"dim=",dim,"error=",error1,"size=",t2.size());

        g0=g0.convert(TensorArgs(targs.thresh,TT_2D));
        error1=(t2-g0.full_tensor_copy()).normf();
        print(ok(is_small(error1,eps)),"TT_2D; k=",k,"dim=",dim,"error=",error1,"size=",t2.size());

    }

    return nerror;
}

/// test the tensor train representation
int testTensorTrain(const long k, const long dim, const TensorArgs targs) {

	print("entering testTensorTrain");
	int nerror=0;
	double eps=targs.thresh;

	// set up a tensor of full rank
	{
		std::vector<long> d(dim,k);
		Tensor<double> t(d);
		t.fillrandom();

		TensorTrain<double> tt1(t,targs.thresh);
		double error1=(t-tt1.reconstruct()).normf();
		print(ok(is_small(error1,eps)),"full rank; k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
		if (!is_small(error1,eps)) nerror++;

		// test fusedim/splitdim
		TensorTrain<double> tt2=copy(tt1);
		tt2.fusedim(0);
		tt2=tt2.splitdim(0,k,k,eps);
		tt2-=tt1;
		double error2=tt2.normf();
        print(ok(is_small(error2,eps)),"split/fusedim; k=",k,"dim=",dim,"error=",error2,"size=",tt1.size());
        if (!is_small(error2,eps)) nerror++;

	}

	// set up a tensor of low rank
	{
		std::vector<long> d(dim,k);
		Tensor<double> t(d);
		t.fillindex();
		t.scale(1.0/t.normf());

		TensorTrain<double> tt1(t,targs.thresh);
		double error1=(t-tt1.reconstruct()).normf();
		print(ok(is_small(error1,eps)),"low rank;  k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
		if (!is_small(error1,eps)) nerror++;
	}

	// set up a tensor of zero rank
	{
		std::vector<long> d(dim,k);
		Tensor<double> t(d);

		TensorTrain<double> tt1(t,targs.thresh);
		double error1=(t-tt1.reconstruct()).normf();
		print(ok(is_small(error1,eps)),"low rank;  k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
		if (!is_small(error1,eps)) nerror++;
	}

	return nerror;
}

/// test arbitrarily
void test(const long& k, const long& dim, const TensorArgs& targs) {

	print("entering arbitrary test");
	// set up three tensors (zeros, rank 2, full rank)
	std::vector<long> d(dim,k);
	std::vector<Tensor<double> > t(3);

	t[0]=Tensor<double>(d);
	t[1]=Tensor<double>(d).fillindex();
	t[2]=Tensor<double>(d).fillrandom();
	for (std::vector<Tensor<double> >::iterator it=t.begin(); it!=t.end(); ++it) {
		*it *=1.0/it->normf();
	}

	GenTensor<double> r(t[2],targs);
	int maxrank=r.rank();
	GenTensor<double> rr=copy(r.get_configs(0,maxrank/2));
	Tensor<double> reftensor=rr.full_tensor_copy();

	TensorTrain<double> tt(reftensor,targs.thresh);
	GenTensor<double> gt(reftensor,targs);
	print("tt",tt.size(),tt.real_size());
	print("gt",gt.size(),gt.real_size(),gt.rank());

	Tensor<double> diff=reftensor-gt.full_tensor_copy();
	print("error",diff.normf());
	Tensor<double> diff1=reftensor-tt.reconstruct();
	print("error",diff1.normf());

}


int test_TT_truncate(const long k, const long dim, const TensorArgs targs) {
	print("entering test_TT_truncate");

	int nerror=0;
	double eps=targs.thresh;

	// set up a tensor of low rank
	{
		std::vector<long> d(dim,k);
		Tensor<double> t(d);
		t.fillindex();
		t.scale(1.0/t.normf());

		TensorTrain<double> tt1(t,targs.thresh);
		double error1=(t-tt1.reconstruct()).normf();
		print(ok(is_small(error1,eps)),"low rank;  k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
		if (!is_small(error1,eps)) nerror++;
		print("tt ranks",tt1.ranks());

		tt1+=tt1;
		t+=t;
//		TensorTrain<double> tt1(t,targs.thresh);
		error1=(t-tt1.reconstruct()).normf();
		print(ok(is_small(error1,eps)),"low rank;  k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
		if (!is_small(error1,eps)) nerror++;
		print("tt ranks",tt1.ranks());

		tt1.truncate(eps);
		error1=(t-tt1.reconstruct()).normf();
		print(ok(is_small(error1,eps)),"low rank;  k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
		if (!is_small(error1,eps)) nerror++;
		print("tt ranks",tt1.ranks());


	}
	return nerror;
}


int test_TT_operator_application(const long k, const long dim, const TensorArgs targs) {
    print("entering test_TT_operator_application");
    print("k, dim, thresh ",k,dim,targs.thresh);

    int nerror=0;
    double eps=targs.thresh;

    // set up a tensor of low rank
    {
        std::vector<long> d(dim,k);
        Tensor<double> t(d);
        t.fillindex();
        t.scale(1.0/t.normf());

        TensorTrain<double> tt1(t,targs.thresh);
        double error1=(t-tt1.reconstruct()).normf();
        print(ok(is_small(error1,eps)),"low rank;  k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
        if (!is_small(error1,eps)) nerror++;

        // set up identity operator
        TensorTrain<double> id=tt_identity<double>(dim,k);
        TensorTrain<double> tt2=apply(id,tt1,targs.thresh);
        tt2-=tt1;
        double error2=tt2.normf();
        print(ok(is_small(error1,eps)),"low rank;  k=",k,"dim=",dim,"error=",error1,"size=",tt1.size());
        if (!is_small(error2,eps)) nerror++;

    }
    return nerror;
}

int main(int argc, char**argv) {

//    initialize(argc,argv);
//    World world(MPI::COMM_WORLD);

    srand(time(nullptr));
    std::cout << std::scientific;

    // the parameters
    long k=4;
    const unsigned int dim=6;
    double eps=1.e-3;
    print("k    ",k);
    print("eps  ",eps);


    int error=0;
    print("hello world");

////    test(k,dim,TensorArgs(eps,TT_2D));
    for (int kk=4; kk<k+1; ++kk) {
    	for (int d=3; d<dim+1; ++d) {
    	    error+=testTensorTrain(kk,d,TensorArgs(eps,TT_2D));
    	    error+=test_TT_truncate(kk,d,TensorArgs(eps,TT_2D));
            error+=test_TT_operator_application(kk,d,TensorArgs(eps,TT_2D));
    	}
    }


#if 1
    error+=testGenTensor_ctor(k,dim,eps,TT_FULL);
    error+=testGenTensor_ctor(k,dim,eps,TT_2D);
    error+=testGenTensor_ctor(k,dim,eps,TT_TENSORTRAIN);

    error+=testGenTensor_assignment(k,dim,eps,TT_FULL);
    error+=testGenTensor_assignment(k,dim,eps,TT_2D);
    error+=testGenTensor_assignment(k,dim,eps,TT_TENSORTRAIN);

    error+=testGenTensor_algebra(k,dim,eps,TT_FULL);
    error+=testGenTensor_algebra(k,dim,eps,TT_2D);
    error+=testGenTensor_algebra(k,dim,eps,TT_TENSORTRAIN);

    error+=testGenTensor_convert(k,dim,TensorArgs(eps,TT_FULL));
    error+=testGenTensor_convert(k,dim,TensorArgs(eps,TT_2D));
    error+=testGenTensor_convert(k,dim,TensorArgs(eps,TT_TENSORTRAIN));

    error+=testGenTensor_rankreduce(k,dim,eps,TT_FULL);
    error+=testGenTensor_rankreduce(k,dim,eps,TT_2D);
    error+=testGenTensor_rankreduce(k,dim,eps,TT_TENSORTRAIN);

    error+=testGenTensor_transform(k,dim,eps,TT_FULL);
    error+=testGenTensor_transform(k,dim,eps,TT_2D);
    error+=testGenTensor_transform(k,dim,eps,TT_TENSORTRAIN);

    error+=testGenTensor_reconstruct(k,dim,eps,TT_FULL);
    error+=testGenTensor_reconstruct(k,dim,eps,TT_2D);
    error+=testGenTensor_reconstruct(k,dim,eps,TT_TENSORTRAIN);

    error+=testGenTensor_deepcopy(k,dim,eps,TT_FULL);
    error+=testGenTensor_deepcopy(k,dim,eps,TT_2D);
    error+=testGenTensor_deepcopy(k,dim,eps,TT_TENSORTRAIN);

    error+=testGenTensor_reduce(k,dim,eps,TT_2D);

    print(ok(error==0),error,"finished test suite\n");
#endif

//    world.gop.fence();
//    finalize();

    return 0;
}

#else
int main(int argc, char** argv) {

    print("no testseprep without having a GenTensor");
    return 0;
}

#endif


