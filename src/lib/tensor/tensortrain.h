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

#ifndef TENSORTRAIN_H_
#define TENSORTRAIN_H_

#include "tensor/tensor.h"
#include "tensor/srconf.h"
#include <linalg/clapack.h>
#include <linalg/tensor_lapack.h>


/**
  \file tensortrain.h
  \brief Defines and implements the tensor train decomposition as described in
         I.V. Oseledets, Siam J. Sci. Comput. 33, 2295 (2011).
  \ingroup tensor
  \addtogroup tensor

*/


namespace madness {

	/**
	 * A tensor train is a multi-modal representation of a tensor t
	 * \code
	 *  t(i,j,k,l) = \sum G^1_{a1,i,a2} G^2_{a2,j,a3} G^3_{a3,k,a4} G^4_{a4,l,a5}
	 * \endcode
	 * The "core" tensors G are connected via a linear index network, where the
	 * first index a1 and the last index a5 are boundary indices and are set to 1.
	 *
	 */
	template<typename T>
	class TensorTrain {

		/// holding the core tensors of a tensor train
		std::vector<Tensor<T> > core;

	public:

		/// ctor for a TensorTrain, with the tolerance eps

		/// the tensor train will represent the input tensor with
		/// accuracy || t - this ||_2 < eps
		/// @param[in]	t	full representation of a tensor
		/// @param[in]	eps	the accuracy threshold
		TensorTrain(const Tensor<T>& t, double eps)
			: core(std::vector<Tensor<T> >(t.ndim())) {

			const long k=t.dim(0);
			eps=eps/sqrt(t.ndim()-1);	// error is relative

			Tensor<T> u,vt;
			Tensor< typename Tensor<T>::scalar_type > s;

			Tensor<T> c=t;

			// tentative ranks
			std::vector<long> r(t.ndim()+1,k);
			r[0]=1;

			for (long d=1; d<t.ndim(); ++d) {

				const long d1=r[d-1]*k;
				c=c.reshape(d1,c.size()/d1);

				svd(c,u,s,vt);

				r[d]=SRConf<T>::max_sigma(eps,s.dim(0),s)+1;


				u=copy(u(_,Slice(0,r[d]-1)));
				vt=vt(Slice(0,r[d]-1),_);

				for (int i=0; i<vt.dim(0); ++i) {
					for (int j=0; j<vt.dim(1); ++j) {
						vt(i,j)*=s(i);
					}
				}
				core[d-1]=u.reshape(r[d-1],k,r[d]);
				c=copy(vt);
			}
			core[t.ndim()-1]=c;
			core[0]=core[0].fusedim(0);
		}

		/// reconstruct this to a full representation

		/// @return	this in full rank representation
		Tensor<T> reconstruct() const {
			Tensor<T> result=core.front();
			typename std::vector<Tensor<T> >::const_iterator it;

			for (it=++core.begin(); it!=core.end(); ++it) {
				result=inner(result,*it);
			}
			print("result.ndim", result.ndim());
			return result;
		}

		/// construct a two-mode representation (aka unnormalized SVD)

		/// @param[out] U({i},rank) the left singular vectors
		/// @param[out] VT(rank,{i}) the right singular vectors
		/// @param[out]	s(rank) vector holding 1's
		void two_mode_representation(Tensor<T>& U, Tensor<T>& VT,
				Tensor< typename Tensor<T>::scalar_type >& s) {

			/// number of dimensions needs to be even
			MADNESS_ASSERT(ndim()%2==0);

			typename std::vector<Tensor<T> >::const_iterator it1, it2;

			U=core.front();
			VT=core.back();

			for (it1=++core.begin(), it2=--(--core.end()); it1<it2; ++it1, --it2) {
				U=inner(U,*it1);
				VT=inner(*it2,VT);
			}
			s=Tensor< typename Tensor<T>::scalar_type >(VT.dim(0));
			s=1.0;
		}

		/// return the number of dimensions
		long ndim() const {return core.size();}

	};


}

#endif /* TENSORTRAIN_H_ */
