/*
 * LRTensor_test.cc
 *
 *  Created on: Oct 2, 2019
 *      Author: fbischoff
 */

#ifndef ENABLE_GENTENSOR
#if ENABLE_GENTENSOR==1
#define ENABLE_GENTENSOR
#endif
#endif

#include <madness/tensor/gentensor.h>
#include <madness/tensor/RandomizedMatrixDecomposition.h>


using namespace madness;

template<typename T>
double compute_difference(const GenTensor<T>& lrt1, const GenTensor<T>& lrt2) {\
	return (lrt1.full_tensor_copy()-lrt2.full_tensor_copy()).normf();
}

template<typename T>
double compute_difference(const GenTensor<T>& lrt1, const Tensor<T>& t2) {\
	return (lrt1.full_tensor_copy()-t2).normf();
}

template<typename T>
double compute_difference(const Tensor<T>& t2, const GenTensor<T>& lrt1) {\
	return (lrt1.full_tensor_copy()-t2).normf();
}

std::vector<long> make_dimensions(bool even=false) {
	madness::Random(10);

	long ndim=(madness::RandomValue<long>() %5+2);	// anything from 2 to 6
    if (even and (ndim%2==1)) ndim++;
	std::vector<long> dim(ndim);
	for (long& d : dim) d=(madness::RandomValue<long>() %6+1);// anything from 1 to 6
	print("dimensions",dim);
	return dim;
}


template<typename T>
int test_constructor() {
	print("\nentering test_construct");
	std::vector<long> dim=make_dimensions(true);
	Tensor<T> tensor(dim);
	double thresh=1.e-5;
	for (int i=0; i<2; ++i) {
		print("repeating with random tensors");
		SVDTensor<T> svd(tensor,thresh);
		TensorTrain<T> tt(tensor,thresh);

		GenTensor<T> lrt;
		GenTensor<T> lrt1(tensor);
		GenTensor<T> lrt2(tensor,thresh,TT_2D);
		GenTensor<T> lrt3(tensor,thresh,TT_TENSORTRAIN);
		GenTensor<T> lrt4(tensor,thresh,TT_FULL);
		GenTensor<T> lrt5(lrt4);

		GenTensor<T> lrt6(svd);
		GenTensor<T> lrt7(tt);

		GenTensor<T> lrt8(dim,TT_FULL);
		GenTensor<T> lrt9(dim,TT_2D);
		GenTensor<T> lrt10(dim,TT_TENSORTRAIN);
		GenTensor<T> lrt11(TT_2D,dim[0],dim.size());
		tensor.fillrandom();

	}

	// not throwing is success
	return 0;
}

template<typename T>
int test_reconstruct() {
	print("\nentering test_reconstruct");

	std::vector<long> dim=make_dimensions();
	Tensor<T> tensor(dim);
	tensor.fillrandom();
	double thresh=1.e-5;
	double error=0.0;

	GenTensor<T> lrt1(tensor);
	GenTensor<T> lrt2(tensor,thresh,TT_2D);
	GenTensor<T> lrt3(tensor,thresh,TT_TENSORTRAIN);
	GenTensor<T> lrt4(tensor,thresh,TT_FULL);

	double error1=compute_difference(lrt1,tensor);
	double error2=compute_difference(lrt2,tensor);
	double error3=compute_difference(lrt3,tensor);
	double error4=compute_difference(lrt4,tensor);

	print("error1",error1);
	print("error2",error2);
	print("error3",error3);
	print("error4",error4);

	error=error1+error2+error3+error4;
	return 0;
}



template<typename T>
int test_addition(const TensorType& tt) {

	print("\nentering test_addition", tt);
	std::vector<long> dim=make_dimensions(tt==TT_2D);
	Tensor<T> tensor(dim);
	tensor.fillrandom();
	double error=0.0;

	GenTensor<T> lrt1(tensor,TensorArgs(1.e-4,tt));
	GenTensor<T> lrt2(tensor,TensorArgs(1.e-4,tt));
	Tensor<T> tensor2=copy(tensor);
	lrt2+=lrt1;
	tensor2+=tensor;
	double error1=compute_difference(lrt2,tensor2);
	print("error1",error1);
	error+=error1;

	Tensor<T> tensor3=tensor+tensor;
	GenTensor<T> lrt3=lrt1+lrt1;
	double error2=compute_difference(lrt3,tensor3);
	print("error2",error2);
	error+=error2;

	Tensor<T> tensor4=tensor-tensor2;
	GenTensor<T> lrt4=lrt1-lrt2;
	double error4=compute_difference(lrt4,tensor4);
	print("error4",error4);
	error+=error4;

	if (tt==TT_2D) {
		tensor4+=tensor2;
		lrt4.reduce_rank(1.e-4);
		lrt2.reduce_rank(1.e-4);
		lrt4.get_svdtensor().add_SVD(lrt2.get_svdtensor(),1.e-4);
		double error5=compute_difference(lrt4,tensor4);
		print("error5",error5);
		error+=error5;


	}


	return (error>1.e-8);
}


template<typename T>
int test_sliced_assignment(const TensorType& tt) {

	print("\nentering test_addition", tt);
	Tensor<T> tensor(10,10,10,10);
	tensor.fillrandom();
	double error=0.0;
	std::vector<Slice> s0={_,_,Slice(0,8),Slice(3,5)};
	std::vector<Slice> s1={_,_,Slice(1,9),Slice(2,4)};

	return (error>1.e-8);
}

template<typename T>
int test_sliced_addition(const TensorType& tt) {

	print("\nentering test_slicing_addition", tt);
	Tensor<T> tensor1(10,8,11,10);
	Tensor<T> tensor2=3.0*copy(tensor1);
	tensor1.fillrandom();
	double error=0.0;
	std::vector<Slice> s1={_,_,Slice(0,8),Slice(3,5)};
	std::vector<Slice> s2={_,_,Slice(1,9),Slice(2,4)};

	GenTensor<T> lrt1(tensor1,TensorArgs(1.e-4,tt));
	GenTensor<T> lrt2(tensor2,TensorArgs(1.e-4,tt));

	tensor1(s1)+=tensor2(s2);
	lrt1(s1)+=lrt2(s2);
	double error2=compute_difference(lrt1,tensor1);
	print("error2",error2);
	error+=error2;


	return (error>1.e-8);
}


template<typename T>
int test_reduce_rank(const TensorType& tt) {

	print("\nentering test_reduce_rank", tt);
	std::vector<long> dim=make_dimensions(tt==TT_2D);
	Tensor<T> tensor1(dim);
	Tensor<T> tensor2(dim);
	tensor1.fillrandom();
	tensor2.fillrandom();
	double error=0.0;
	double thresh=1.e-5;

	GenTensor<T> lrt1(tensor1,TensorArgs(1.e-4,tt));
	GenTensor<T> lrt2(tensor2,TensorArgs(1.e-4,tt));

	tensor1+=tensor2;
	lrt1+=lrt2;
	lrt1.reduce_rank(thresh);

	double error2=compute_difference(lrt1,tensor1);
	print("error2",error2);
	error+=error2;


	return (error>1.e-8);
}



template<typename T>
int test_emul(const TensorType& tt) {

	print("\nentering test_emul", tt);
	std::vector<long> dim=make_dimensions(tt==TT_2D);
	Tensor<T> tensor1(dim);
	Tensor<T> tensor2(dim);
	tensor1.fillrandom();
	tensor2.fillrandom();
	double error=0.0;
	// double thresh=1.e-5;

	GenTensor<T> lrt1(tensor1,TensorArgs(1.e-4,tt));
	GenTensor<T> lrt2(tensor2,TensorArgs(1.e-4,tt));

	tensor1.emul(tensor2);
	lrt1.emul(lrt2);

	double error2=compute_difference(lrt1,tensor1);
	print("error2",error2);
	error+=error2;


	return (error>1.e-8);
}

template<typename T>
int test_convert() {
	print("\nentering test_convert");

	Tensor<T> tensor(10,10,10,10);
	tensor.fillrandom();
	double thresh=1.e-5;

	GenTensor<T> lrt(tensor);
	GenTensor<T> lrt1=lrt.convert(TensorArgs(TT_2D,thresh));
	GenTensor<T> lrt2=lrt.convert(TensorArgs(TT_TENSORTRAIN,thresh));

	GenTensor<T> lrt3=lrt1.convert(TensorArgs(TT_FULL,thresh));
	GenTensor<T> lrt4=lrt2.convert(TensorArgs(TT_FULL,thresh));

	GenTensor<T> lrt5=lrt2.convert(TensorArgs(TT_2D,thresh));

	double error1=compute_difference(lrt1,tensor);
	double error2=compute_difference(lrt2,tensor);
	double error3=compute_difference(lrt3,tensor);
	double error4=compute_difference(lrt4,tensor);
	double error5=compute_difference(lrt5,tensor);
	print("error1", error1);
	print("error2", error2);
	print("error3", error3);
	print("error4", error4);
	print("error5", error5);

	return 0;

}

template<typename T>
int test_general_transform() {
	print("\nentering test_general_transform");

	Tensor<T> tensor(10,10,10,10);
	Tensor<T> c(10,10);
	tensor.fillrandom();
	c.fillrandom();
	double thresh=1.e-5;
	Tensor<double> cc[TENSOR_MAXDIM];
	for (long idim=0; idim<4; idim++) {
		cc[idim]=Tensor<double>(10,10);
		cc[idim].fillrandom();
	}

	GenTensor<T> lrt1(tensor);
	GenTensor<T> lrt2(tensor,thresh,TT_TENSORTRAIN);
	GenTensor<T> lrt3(tensor,thresh,TT_2D);
	GenTensor<T> result1=general_transform(lrt1,cc);
	GenTensor<T> result2=general_transform(lrt2,cc);
	GenTensor<T> result3=general_transform(lrt3,cc);

	double error2=(result2.full_tensor_copy()-result1.full_tensor_copy()).normf();
	double error3=(result2.full_tensor_copy()-result1.full_tensor_copy()).normf();
	print("error2",error2);
	print("error3",error3);

	return 0;
}


template<typename T>
int test_stuff() {
	double thresh=1.e-5;
	long k=5;
	Tensor<T> matrix(std::pow(k,3l),std::pow(k,3l));

	print("k, thresh",k,thresh);
	matrix.fillindex();
	matrix=matrix*(T(1.0)/matrix.normf());

	printf("%15s %12s %12s %12s\n","method","range","error ratio","time");

	RandomizedMatrixDecomposition<T> rmd=RMDFactory().oversampling(10);
	rmd.compute_range(matrix,thresh,{0,0});

	Tensor<T> Q;
	if (not rmd.exceeds_maxrank()) Q=rmd.get_range();

	double error=RandomizedMatrixDecomposition<T>::check_range(matrix,Q);
//	printf("dim(range), error ratio %d %12.8f\n",Q.dim(1),error/thresh);

	double wall0=wall_time();
	matrix.fillrandom();
	matrix=RandomizedMatrixDecomposition<T>::make_SVD_decaying_matrix(matrix,3);
	matrix=matrix*(T(1.0)/matrix.normf());
	double wall1=wall_time();
//	printf("wall  %12.8f\n",wall1-wall0);


	wall0=wall_time();
	Q=rmd.compute_range(matrix,thresh*0.1,{0,0});
	wall1=wall_time();
	error=RandomizedMatrixDecomposition<T>::check_range(matrix,Q);
	printf("%15s %12ld %12.4f %12.8f\n","compute_range",Q.dim(1),error/thresh,wall1-wall0);

	wall0=wall_time();
	SVDTensor<T> st=SVDTensor<T>::compute_svd_from_range(Q,matrix.reshape(k,k,k,k,k,k));
	wall1=wall_time();
	error=(st.reconstruct().flat()-matrix.flat()).normf();
	printf("%15s %12ld %12.4f %12.8f\n","svd from range",Q.dim(1),error/thresh,wall1-wall0);

	wall0=wall_time();
	SVDTensor<T> st1(matrix,thresh);
	wall1=wall_time();
	error=(st1.reconstruct().flat()-matrix.flat()).normf();
	printf("%15s %12ld %12.4f %12.8f\n","svdtensor",st1.rank(),error/thresh,wall1-wall0);

	wall0=wall_time();
	TensorTrain<T> tt(matrix.reshape(k,k,k,k,k,k),thresh);
	wall1=wall_time();
	auto ranks=tt.ranks();
	double tt_error=(matrix.flat()-tt.reconstruct(true)).normf();
	printf("%15s %12ld %12.4f %12.8f\n","TT",Q.dim(1),tt_error/thresh,wall1-wall0);

	Tensor<T> U,VT;
	Tensor<typename Tensor<T>::scalar_type> s;
	wall0=wall_time();
	svd(matrix,U,s,VT);
	wall1=wall_time();
	int rank=SRConf<T>::max_sigma(thresh, s.size(), s);
	printf("%15s %12d %12s %12.8f\n","svd",rank,"---",wall1-wall0);
	return 0;
}

int
main(int argc, char* argv[]) {

	madness::default_random_generator.setstate(int(cpu_time())%4149);
	int success=0;
	success+=test_stuff<double>();
#if 1
    success += test_constructor<double>();
    success += test_constructor<double_complex>();

//    success +=  test_reconstruct<double>();
//    success +=  test_reconstruct<double_complex>();

    success +=  test_convert<double>();
    success +=  test_convert<double_complex>();

    std::vector<TensorType> tt={TT_FULL,TT_2D,TT_TENSORTRAIN};
    for (auto& t : tt) {
    	success += test_addition<double>(t);
    	success += test_sliced_addition<double>(t);
    	success += test_reduce_rank<double>(t);
    	success += test_emul<double>(t);
    }

    success += test_general_transform<double>();
    success += test_general_transform<double_complex>();
#endif

    std::cout << "Test " << ((success==0) ? "passed" : "did not pass") << std::endl;

    return success;
}
