//
// Ghaly
//
// Compresses the function, transforming into wavelet basis. 
// Possible non-blocking comm.
//
// By default fence=true meaning that this oepration completes before returning,
// othewise if fence=false it returns without fencing and the user must invoke 
// workd.gop.fence() to assure global completion before using the function
// for other purposes.
//
// Noop if already compressed or if not initialized.
//
// Since reconstruction/compression do not discard information we define them
// as const ... "logical constness" not "bitwise contness".
//
#ifndef __MADNESS_MRA_FUSET_COMPRESS_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_COMPRESS_OP__INCLUDED__

#include "PrimitiveOp.h"
#include "CompressParameters.h"
#include "../mra.h"
#include "../function_common_data.h"
#include "../../world/MADworld.h"
#include "../../tensor/tensor.h"
#include "../../tensor/gentensor.h"

namespace madness 
{
	template<typename T, std::size_t NDIM>
	class CompressOp : public PrimitiveOp<T,NDIM> 
	{
	typedef Function<T,NDIM>									KTREE;
	typedef FunctionNode<T,NDIM>								KNODE;	// identical to nodeT
	typedef Key<NDIM>											keyT;
	typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM>>	dcT;
	typedef WorldObject< FunctionImpl<T,NDIM> >					woT;	///< Base class world object type

	typedef GenTensor<T>										coeffT;	//Type of tensor used to hold coeffs
	typedef Tensor<T>											tensorT;

	///< Type of container holding the nodes
	    
    public:
									CompressOp		(string opName, KTREE* output, const KTREE* i1);
		void						compute			(const keyT& key) { }

		//!Blocking
		Future <GenTensor<T>>		computeB		(const keyT& key);
		Future <GenTensor<T>>		afterComputeB	(const keyT& key, const std::vector<Future<coeffT>> &v);
		
		bool						isDone			(const keyT& key) const;
		bool						isPre			() const { return true; }

	public:	// for CompressOp
		coeffT						filter			(const coeffT& s) const;
		std::vector<Slice>			child_patch		(const keyT& child) const;


    private:
		//!Points to operand trees
		const KTREE*						_i1;
    
		//!Points to operand nodes of the tree
		KNODE								*_t1, *_t2;

		//!Variables for CompressOp
		dcT&								_coeffs;
		const FunctionCommonData<T,NDIM>&	_cdata; 
		TensorArgs							_targs;	

		int									_k;				// Wavelet order	
		bool								_nonstandard;
		bool								_keepleaves;
		bool								_redundant;
		bool								_root;
    };

		
	// Constructor
	// World is needed for communication in the "compute" function
    template<typename T, std::size_t NDIM>
	CompressOp<T,NDIM>::CompressOp(string opName, KTREE* output, const KTREE* i1)
		: PrimitiveOp<T,NDIM>(opName, output, false)
		, _i1(i1)
		, _cdata(FunctionCommonData<T,NDIM>::get(i1->k()))
		, _targs(i1->get_impl()->get_tensor_args())
		, _nonstandard(false)
		, _keepleaves(false)
		, _redundant(false)
		, _root(false)
		, _coeffs(i1->get_impl()->get_coeffs())
		, _k(i1->get_impl()->get_k())
	{
		// output is itself.

	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));

		woT(i1->world());
	}

	// Compute  <--- compress in mra.h
	//
	// "compute" function is called by each node through traversal tree.
	// : it means that "compute" function is identical with compress_spawn 
	//   but plays "compress()" role as well.
	//
	//
	// This function should be composed of compress() and compress_spawn()
    template<typename T, std::size_t NDIM>
	Future<GenTensor<T>>
	CompressOp<T,NDIM>::computeB(const keyT& key)  
	{	
		if ( key == _i1->world().rank())
		{
			if (_i1->is_compressed() || _i1->get_impl() == NULL) 
			{
				exit(0);
			}

			_i1->get_impl()->compressed		= true;
			_i1->get_impl()->nonstandard	= this->_nonstandard;
			_i1->get_impl()->redundant		= this->_redundant;
		}
	
		// Getting Current Unit (Node) by using key.
		typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();
        if (it == _i1->get_impl()->get_coeffs().end()) 
		{
			cerr<<"This should not have happenned"<<endl;
			exit(0);
        }

		// get fetches remote data (here actually local)
		KNODE& node = _coeffs.find(key).get()->second;
		cout << "key: " << key << endl;
		cout << "node: " << node << endl;
		if (node.has_children())
		{
			// task will be created in traverseTreeB.
			// Itermediate nodes do work after children have done something.
		}
		else
		{
			// Leaf node
			Future<coeffT> result(node.coeff());
			if (!_keepleaves)
				node.clear_coeff();

			//printf ("[%s] this node (", __func__);
			//std::cout << key << ") is a leaf node\n" << std::endl;
			return result;
		}

		//return *castedReturn;
		//return *returnP;
		Future<coeffT> temp;
		return temp;
	}

	// 
	//	This Function is executed after running childrend
	//				  is for CompressOp from virtual function of PrimitiveOp
	//
	//	It should have a common parameter for data from children 
	//	such as "std::vector<Future<coeffT>>& v"
	//	This will be replaced with "CompressParameters" class with Future.
	//	CompressOp<T,NDIM>::compress_op(const keyT& key, const std::vector< Future<coeffT> >& v, bool nonstandard, bool redundant)
	//
	template <typename T, std::size_t NDIM>
	Future<GenTensor<T>>
	CompressOp<T,NDIM>::afterComputeB(const keyT& key, const std::vector<Future<coeffT>> &v)
	{
		// Copy child scaling coeffs into contiguous block
		tensorT d(this->_cdata.v2k);
		
		int i = 0;
		for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i)
		{
			d(child_patch(kit.key())) += v[i].get().full_tensor_copy();	
		}

		d = filter(d);

		typename dcT::accessor acc;
		_coeffs.find(acc, key);

		if (acc->second.has_coeff())
		{
			//printf (" stuff in compress_op");
		
			const tensorT c = acc->second.coeff().full_tensor_copy();
			if (c.dim(0) == _k)
				d(_cdata.s0) += c;
			else
				d += c;
		}

		// tighter thresh for intenal nodes
		TensorArgs targs2 = this->_targs;	// targs < Type of tensor to be used in the FunctionNode
		targs2.thresh*=0.1;
	
		// need the deep copy for contiguity
		coeffT ss = coeffT(copy(d(_cdata.s0)), targs2);

		if (key.level() > 0 && !_nonstandard)
			d(_cdata.s0) = 0.0;

		// insert either sum or difference coefficients
		if (this->_redundant)
		{
			coeffT dd = coeffT(d, targs2);
			acc->second.set_coeff(ss);
		}
		else
		{
			coeffT dd = coeffT(d, targs2);
			acc->second.set_coeff(dd);
		}

		//	Making Future form
		Future<GenTensor<T>> temp(ss);
		return temp;
	}

	template <typename T, std::size_t NDIM>
	std::vector<Slice>
	CompressOp<T,NDIM>::child_patch(const keyT& child) const 
	{
		std::vector<Slice> s(NDIM);
		const Vector<Translation,NDIM>& l = child.translation();
		for (std::size_t i = 0; i<NDIM; ++i)
		{
			s[i] = _cdata.s[l[i]&1];	// Lowest bit of translation
		}
		return s;
	}
	
	template <typename T, std::size_t NDIM>	
	GenTensor<T>
	CompressOp<T,NDIM>::filter(const coeffT& s) const
	{
		coeffT result = transform(s, _cdata.hgT);
		return result;
	}

	// calculate the wavelet coefficients using the sum coefficients of all child nodes
	// 
	// @param[in] key			this's key
	// @param[in] v				su coefficients of the child nodes
	// @param[in] nonstandard	keep the sum coefficients with the wavelet coeffients
	// @param[in] redundant		keep only the sum coefficients, discard the wavelet coefficients 
	// @return					the sum coefficients
	

	// isDone
    template<typename T, std::size_t NDIM>
	bool CompressOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
		if(!isE1) return isE1;
		bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf;
    }

}; /*fuset*/

#endif /* __fuset_CompressOp_h__ */
