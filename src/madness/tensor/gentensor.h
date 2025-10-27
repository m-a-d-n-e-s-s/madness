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

  $Id$
*/

#ifndef GENTENSOR_H_
#define GENTENSOR_H_


/*!
 *
 * \file GenTensor.h
 * \brief Provides a tensor with taking advantage of possibly low rank
 *
 ****************************************************************
 * MAIN DIFFERENCES (Tensor t; GenTensor g)
 *  t=t1(s) is shallow
 *  g=g1(s) is deep
 ****************************************************************
 *
 * a GenTensor is a generalized form of a Tensor
 * for now only little functionality shall be implemented; feel free to extend
 * a consequence is that we (usually) can't directly access individual
 * matrix elements
 * note that a GenTensor might have zero rank, but is still a valid
 * tensor, and should therefore return a FullTensor with zeros in it.
 *
 * \par Slicing in GenTensors:
 *
 * A GenTensor differs from a Tensor in that we can't directly access
 * individual matrix elements, and thus can't directly assign or manipulate slices
 * as lvalues. For rvalues we simply provide a slices of the constituent vectors
 * in SRConf, which is a valid GenTensor by itself
 * \code
 *  lhs = rhs(s)
 * \endcode
 *
 * Manipulations of slices of a GenTensor are heavily restricted, but should
 * cover the most important cases:
 * - assignment to zero; done by inplace subtraction of the slice
 *    \code
 *    GenTensor lrt1(t);
 *    lrt1(s) = 0.0;
 *    \endcode
 * - in-place addition
 *    \code
 *    GenTensor lrt1(3,3,3,3), lrt2(t);
 *    lrt1(s) += lrt2;
 *    \endcode
 *
 * Note that *all* of these operation increase the rank of lhs
 *
 * \par Addition in GenTensors
 *
 * Addition in a GenTensor is a fairly complicated issue, and there are several
 * algorithms you can use, each with certain pros and cons
 *
 * - append()
 *   plain concatenating of the configurations will increase the rank and will
 *   introduce linear dependencies and redundancies. Also, you might have to
 *   reallocate memory and copy a lot of data around. However, no accuracy is lost
 *   and it is fairly fast. Final rank reduction might be expensive, since you
 *   might have accumulated a huge amount of terms.
 *
 * - low_rank_add_sequential()
 *   only for SVD (2-way decomposition)
 *   take the 2nd vector as a basis of a vector space, project the overlap of the
 *   old and new terms on the lhs basis, perform modified Gram-Schmidt
 *   orthogonalization on the rhs basis and increase the basis if necessary.
 *   This will require the left hand side to be right-orthonormal, and the caller is
 *   responsible for doing that. If a new GenTensor is created and only additions
 *   using this method are performed the tensor will automatically be orthonormal.
 *   Cost depends on the number of terms of the rhs and the lhs
 *
 * - addition of slices
 *   as of now we can add slices only using append()
 *
 * - addition in full rank form
 *   at times it might be sensible to accumulate your result in a full rank tensor,
 *   and after collecting all the contribution you can transform it into a low
 *   rank tensor. This will also maintain "all" digits, and cost depends not on
 *   the number of terms on the lhs.
 */

#include <madness/tensor/tensor.h>


namespace madness {


/// low rank representations of tensors (see gentensor.h)
enum TensorType {TT_FULL, TT_2D, TT_TENSORTRAIN};

static
inline
std::ostream& operator << (std::ostream& s, const TensorType& tt) {
	std::string str="confused tensor type";
	if (tt==TT_FULL) str="full rank tensor";
	if (tt==TT_2D) str="low rank tensor 2-way";
	if (tt==TT_TENSORTRAIN) str="tensor train";
//	if (tt==TT_DYNAMIC) str="dynamic";
	s << str.c_str();
	return s;
}
	/// TensorArgs holds the arguments for creating a LowRankTensor
	struct TensorArgs {
		double thresh;
		TensorType tt;
	TensorArgs() : thresh(-1.0), tt(TT_FULL) {}
		TensorArgs(const double& thresh1, const TensorType& tt1)
			: thresh(thresh1)
			, tt(tt1) {
		}
	TensorArgs(const TensorType& tt1, const double& thresh1)
	: thresh(thresh1)
	, tt(tt1) {
	}

		static std::string what_am_i(const TensorType& tt) {
			if (tt==TT_2D) return "TT_2D";
			if (tt==TT_FULL) return "TT_FULL";
            if (tt==TT_TENSORTRAIN) return "TT_TENSORTRAIN";
			return "unknown tensor type";
		}
		template <typename Archive>
		void serialize(const Archive& ar) {
		    int i=int(tt);
		    ar & thresh & i;
		    tt=TensorType(i);
		}
	};
}

// you can use low-rank tensors only when you use gentensor
#if HAVE_GENTENSOR
#include <madness/tensor/lowranktensor.h>

#else

#include <madness/tensor/SVDTensor.h>
#include <madness/tensor/tensortrain.h>

namespace madness {

	template<typename T> class SRConf;

	template <typename T>
	class GenTensor : public Tensor<T> {

	public:

		GenTensor() : Tensor<T>() {}

		GenTensor(const Tensor<T>& t1) : Tensor<T>(t1) {}
		GenTensor(const Tensor<T>& t1, const TensorArgs& targs) : Tensor<T>(t1) {}
		GenTensor(const Tensor<T>& t1, double eps, const TensorType tt) : Tensor<T>(t1) {}
		GenTensor(const TensorType tt): Tensor<T>() {}
		GenTensor(std::vector<long> v, const TensorType& tt) : Tensor<T>(v) {}
		GenTensor(std::vector<long> v, const TensorArgs& targs) : Tensor<T>(v) {}
		GenTensor(const SRConf<T>& sr1) : Tensor<T>() {MADNESS_EXCEPTION("no ctor with SRConf: use HAVE_GENTENSOR",1);}
		GenTensor(long nd, const long d[], const TensorType& tt) : Tensor<T>(nd,d){};

        /// Type conversion makes a deep copy
        template <class Q> operator GenTensor<Q>() const { // type conv => deep copy
            Tensor<Q> result = Tensor<Q>(this->_ndim,this->_dim,false);
            BINARY_OPTIMIZED_ITERATOR(Q, result, const T, (*this), *_p0 = (Q)(*_p1));
            return result;
        }

        GenTensor convert(const TensorArgs& targs) const {return copy(*this);}
        Tensor<T> reconstruct_tensor() const {return *this;}
        const Tensor<T>& full_tensor() const {return *this;}
        Tensor<T>& full_tensor() {return *this;}

        const Tensor<T>& get_tensor() const {return *this;}
        Tensor<T>& get_tensor() {return *this;}

		Tensor<T> full_tensor_copy() const {return copy(*this);}
		Tensor<T> full_tensor_copy() {return copy(*this);}

        bool is_assigned() const {return this->size()>0;};
        bool has_data() const {return this->size()>0;};
        bool has_no_data() const {return not has_data();};
		long rank() const {return -1;}
		double svd_normf() const {return this->normf();}
		size_t real_size() const {return this->size();}
		size_t nCoeff() const {return this->size();}

        void reduce_rank(const double& eps) {return;};
        void normalize() {return;}

        std::string what_am_i() const {return "GenTensor, aliased to Tensor";};
		TensorType tensor_type() const {return TT_FULL;}
        constexpr bool is_svd_tensor() const {return false;}
        constexpr bool is_tensortrain() const {return false;}
        constexpr bool is_full_tensor() const {return true;}
        bool is_of_tensortype(const TensorType& tt) const {return (tt==TT_FULL);}


        SVDTensor<T>& get_svdtensor() {MADNESS_EXCEPTION("no SVDTensor in aliased GenTensor",1);}
        SVDTensor<T>& get_tensortrain() {MADNESS_EXCEPTION("no SVDTensor in aliased GenTensor",1);}
        const SVDTensor<T>& get_svdtensor() const {MADNESS_EXCEPTION("no SVDTensor in aliased GenTensor",1);}
        const SVDTensor<T>& get_tensortrain() const {MADNESS_EXCEPTION("no SVDTensor in aliased GenTensor",1);}



		void add_SVD(const GenTensor<T>& rhs, const double& eps) {*this+=rhs;}

		SRConf<T> config() const {MADNESS_EXCEPTION("no SRConf in complex GenTensor",1);}
        SRConf<T> get_configs(const int& start, const int& end) const {MADNESS_EXCEPTION("no SRConf in complex GenTensor",1);}

		/// return the additional safety for rank reduction
		static double fac_reduce() {return -1.0;};

	};

    template <class T>
	GenTensor<T> reduce(std::list<GenTensor<T> >& addends, double eps, bool are_optimal=false) {
    	typedef typename std::list<GenTensor<T> >::iterator iterT;
    	GenTensor<T> result=copy(addends.front());
    	for (iterT it=++addends.begin(); it!=addends.end(); ++it) {
    		result+=*it;
    	}
	    return result;
    }

    /// Outer product ... result(i,j,...,p,q,...) = left(i,k,...)*right(p,q,...)

    /// \ingroup tensor
    template <class T>
    GenTensor<T> outer(const GenTensor<T>& left, const GenTensor<T>& right,
            const TensorArgs final_tensor_args) {
        return madness::outer(static_cast<Tensor<T> >(left),static_cast<Tensor<T> >(right));
    }

    /// Outer product ... result(i,j,...,p,q,...) = left(i,k,...)*right(p,q,...)

     /// \ingroup tensor
     template <class T>
     GenTensor<T> outer(const Tensor<T>& left, const Tensor<T>& right,
             const TensorArgs final_tensor_args) {
         return madness::outer(left,right);
     }

     template <typename R, typename Q>
     GenTensor<TENSOR_RESULT_TYPE(R,Q)> general_transform(
             const GenTensor<R>& t, const Tensor<Q> c[]) {
    	 const Tensor<R>& tensor=static_cast<const Tensor<R>& >(t);
         return GenTensor<TENSOR_RESULT_TYPE(R,Q)>(general_transform(tensor,c));
     }



     /// change representation to targ.tt
     template<typename T>
     void change_tensor_type(GenTensor<T>& t, const TensorArgs& targs) {
     	MADNESS_ASSERT(targs.tt==TT_FULL);
     }

	namespace archive {
	/// Serialize a tensor
	template <class Archive, typename T>
	struct ArchiveStoreImpl< Archive, GenTensor<T> > {
		static void store(const Archive& s, const GenTensor<T>& t) {
			if (t.iscontiguous()) {
				s & t.size() & t.id();
				if (t.size()) s & t.ndim() & wrap(t.dims(),TENSOR_MAXDIM) & wrap(t.ptr(),t.size());
			}
			else {
				s & copy(t);
			}
		};
	};


	/// Deserialize a tensor ... existing tensor is replaced
	template <class Archive, typename T>
	struct ArchiveLoadImpl< Archive, GenTensor<T> > {
		static void load(const Archive& s, GenTensor<T>& t) {
			long sz = 0l, id =0l;
			s & sz & id;
			if (id != t.id()) throw "type mismatch deserializing a tensor";
			if (sz) {
				long _ndim=0l, _dim[TENSOR_MAXDIM];
				s & _ndim & wrap(_dim,TENSOR_MAXDIM);
				t = Tensor<T>(_ndim, _dim, false);
				if (sz != t.size()) throw "size mismatch deserializing a tensor";
				s & wrap(t.ptr(), t.size());
			}
			else {
				t = Tensor<T>();
			}
		};
	};

	}

}   // namespace madness

#endif /* HAVE_GENTENSOR */
#endif /* GENTENSOR_H_ */
