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

  $Id: srconf.h 2816 2012-03-23 14:59:52Z 3ru6ruWu $
 */



#include <madness/tensor/SVDTensor.h>
#include <madness/constants.h>
#include <madness/tensor/RandomizedMatrixDecomposition.h>

using namespace madness;

namespace madness {

static std::string reduction_alg="divide_conquer";

template<typename T>
void SVDTensor<T>::set_reduction_algorithm(const std::string alg) {
	reduction_alg=alg;
}

template<typename T>
std::string SVDTensor<T>::reduction_algorithm() {

	if (getenv("REDUCTION_ALGORITHM")){
		char* calg = getenv("REDUCTION_ALGORITHM");
		return std::string(calg);
	}
	return reduction_alg;
}

/// reduce the rank using SVD
template<typename T>
SVDTensor<T> SVDTensor<T>::compute_svd(const Tensor<T>& values,
		const double& eps, std::array<long,2> vectordim) {

	SVDTensor<T> result(values.ndim(),values.dims());

	// fast return if possible
	if (values.normf()<eps) return result;

	Tensor<T> values_eff=RandomizedMatrixDecomposition<T>::resize_to_matrix(values,vectordim);

	// output from svd
	Tensor<T> U;
	Tensor<T> VT;
	Tensor< typename Tensor<T>::scalar_type > s;

	svd(values_eff,U,s,VT);

	// find the maximal singular value that's supposed to contribute
	// singular values are ordered (largest first)
	const double thresh=eps;
	long i=SRConf<T>::max_sigma(thresh,s.dim(0),s);

	// convert SVD output to our convention
	if (i>=0) {
		// copy to have contiguous and tailored singular vectors
		result.set_vectors_and_weights(copy(s(Slice(0,i))), copy(transpose(U(_,Slice(0,i)))),
				copy(VT(Slice(0,i),_)));
	}
	return result;
}


/// reduce the rank using SVD
template<typename T>
SVDTensor<T> SVDTensor<T>::compute_randomized_svd(const Tensor<T>& tensor,
		const double& eps, std::array<long,2> vectordim) {

	// 1/8 is 0.125, so 0.2 should be good for all tensors other than the full rank ones
	long maxrank=std::max(125.0,floor(0.2*sqrt(tensor.size())));

	double wall0=wall_time();
	RandomizedMatrixDecomposition<T> rmd=RMDFactory().maxrank(maxrank);
	Tensor<T> Q=rmd.compute_range(tensor,eps*0.1,vectordim);
	double wall1=wall_time();
	double attempt=wall1-wall0;
	SVDTensor<T> result(tensor.ndim(),tensor.dims());
	if (rmd.exceeds_maxrank()) {
		// fallback option
		wall0=wall_time();
		result=compute_svd(tensor,eps,vectordim);
		wall1=wall_time();
		print("fall back into full SVD: rank, maxrank",Q.dim(1),maxrank, attempt, wall1-wall0);
	} else if (Q.size()>0) {
		result=compute_svd_from_range(Q,tensor);
		result.truncate_svd(eps);

	} else {
		MADNESS_ASSERT(Q.size()==0);
	}

	// TODO: fixme
//	result.orthonormalize(eps);
	return result;

}

template<typename T>
SVDTensor<T> SVDTensor<T>::compute_svd_from_range(const Tensor<T>& Q, const Tensor<T>& tensor) {
	typedef typename Tensor<T>::scalar_type scalar_type;

	MADNESS_ASSERT(tensor.size()%Q.dim(0) == 0);
	std::array<long,2> vectordim={Q.dim(0),tensor.size()/Q.dim(0)};
	const Tensor<T> matrix=RandomizedMatrixDecomposition<T>::resize_to_matrix(tensor,vectordim);

	Tensor<T> B=inner(conj(Q),matrix,0,0);
	Tensor<T> U,VT;
	Tensor<scalar_type> s;
	svd(B,U,s,VT);
	U=copy(conj_transpose(inner(Q,U)));
	SVDTensor<T> result(tensor.ndim(),tensor.dims());
	result.set_vectors_and_weights(s,U,VT);
	return result;
}

template<typename T>
void SVDTensor<T>::recompute_from_range(const Tensor<T>& Q) {

	// check dimensions
	// tensor = A = Q * Q(T) A = Q * Q(T) * left(T) * right
	MADNESS_ASSERT(Q.dim(0)==this->flat_vector(0).dim(1));

	const Tensor<T> U_ri=this->make_left_vector_with_weights();
	const Tensor<T> V_rj=this->flat_vector(1);

	Tensor<T> B=inner(inner(conj(Q),U_ri,0,1),V_rj,1,0);

	Tensor<T> U,VT;
	typedef typename Tensor<T>::scalar_type scalar_type;
	Tensor<scalar_type> s;
	svd(B,U,s,VT);

	U=copy(conj_transpose(inner(Q,U)));
	this->set_vectors_and_weights(s,U,VT);

}


/// reduce the rank using a divide-and-conquer approach
template<typename T>
void SVDTensor<T>::divide_and_conquer_reduce(const double& thresh) {

	if (this->has_no_data() or rank()==0) return;
	if (rank()==1) {
		this->normalize();
		return;
	}

	// divide the SRConf into two
	const long chunksize=8;
	if (rank()>chunksize) {
		SVDTensor<T> chunk1=this->get_configs(0,rank()/2);
		SVDTensor<T> chunk2=this->get_configs(rank()/2+1,rank()-1);
		chunk1.divide_and_conquer_reduce(thresh*0.5);
		chunk2.divide_and_conquer_reduce(thresh*0.5);

		// collect the two SRConfs
		*this=chunk1;
		this->add_SVD(chunk2,thresh);

	} else {

		// and reduce the rank
		this->orthonormalize(thresh);
	}
    MADNESS_ASSERT(this->has_structure());
}

template<typename T>
void SVDTensor<T>::orthonormalize(const double& thresh) {

	if (this->has_no_data() or rank()==0) return;
	this->normalize();
	if (rank()==1) return;

	Tensor< typename Tensor<T>::scalar_type > weights=this->weights_(Slice(0,rank()-1));

    Tensor<T> v0=this->flat_vector(0);
    Tensor<T> v1=this->flat_vector(1);
	ortho3(v0,v1,weights,thresh);
    std::swap(this->vector_[0],v0);
    std::swap(this->vector_[1],v1);
    std::swap(this->weights_,weights);
	this->make_structure();
}


template<typename T>
void SVDTensor<T>::orthonormalize_random(const double& eps) {

	if (this->has_no_data() or rank()==0) return;
	this->normalize();
	if (rank()==1) return;

	long maxrank=std::min(this->kVec(0),this->kVec(1));

	RandomizedMatrixDecomposition<T> rmd=RMDFactory().maxrank(maxrank);
	Tensor<T> scr=this->make_left_vector_with_weights();
	Tensor<T> Q=rmd.compute_range(scr,this->flat_vector(1),eps*0.1);

	recompute_from_range(Q);
	truncate_svd(eps);


}

template<typename T>
void SVDTensor<T>::truncate_svd(const double& thresh) {

	// compute the numerical rank
	long maxrank=SRConf<T>::max_sigma(thresh, rank(), this->weights_);
	*this=this->get_configs(0,maxrank);
}




// explicit instantiation
template class SVDTensor<float>;
template class SVDTensor<double>;
template class SVDTensor<float_complex>;
template class SVDTensor<double_complex>;

} // namespace madness

