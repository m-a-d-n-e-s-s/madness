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
#include "FuseTContainer.h"
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
		typedef Function<T,NDIM> KTREE;
		typedef FunctionNode<T,NDIM> KNODE;
		typedef Key<NDIM> keyT;
		typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>> dcT;
		typedef WorldObject<FunctionImpl<T,NDIM>> woT;
		typedef GenTensor<T> coeffT;
		typedef Tensor<T> tensorT;
	    
    public:
									CompressOp		(string opName, KTREE* output, const KTREE* i1);
		FuseTContainer<T>			compute			(const keyT& key, const FuseTContainer<T> &s);
		//Future<FuseTContainer<T>>	postCompute		(const keyT& key, const std::vector<Future<FuseTContainer<T>>> &s);

		bool notEmpty(map<int,bool>& notEmptyMap) const
		{
			unsigned long treeID = _i1->get_impl()->id().get_obj_id();
		    //cout<<"Checking for treeID : "<<treeID<<" and result is "<<notEmptyMap[treeID]<<endl;
		    return  notEmptyMap[treeID];
		}
		bool						isDone			(const keyT& key) const;
		bool						isPre			() const { return false; } // false does not work. but It should be false.
		bool						needsParameter	() const { return true; }
        void						reduce			(World& world){}
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
		dcT&								_coeffs_target;
		const FunctionCommonData<T,NDIM>&	_cdata; 
		TensorArgs							_targs;	

		int									_k;		// Wavelet order
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
		, _coeffs_target(output->get_impl()->get_coeffs())
		, _k(i1->get_impl()->get_k())
	{
	    // dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,true,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,true,false));

		woT(i1->world());
	}
	
	//
	//	it should hangle both a parent and a leaf node.
	//
	template <typename T, std::size_t NDIM>
	FuseTContainer<T>
	CompressOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s)
	{	
		// For Root
		if ( key == keyT(0) ) // need to be changed.
		{
			this->_result->get_impl()->compressed	= true;
			this->_result->get_impl()->nonstandard	= this->_nonstandard;
			this->_result->get_impl()->redundant	= this->_redundant;
		}
	
		// get fetches remote data (here actually local)
		FuseT_CoeffT<T>*		s_coeff;
		FuseT_VParameter<T>*	v_parameter;
		KNODE&					node = _coeffs.find(key).get()->second;
		if (node.has_children())
		{
			// 
			//	Intermediate Nodes
			//
			KNODE	temp;	// new parent which should be connected to its children.
			s_coeff				= new FuseT_CoeffT<T>();
			v_parameter			= new FuseT_VParameter<T>();
			v_parameter->value	= ((FuseT_VParameter<T>*)s.get())->value;
			tensorT d(this->_cdata.v2k);
			int i = 0;	

			// This could be a potential issue related to the children
			// What if a parent is executed before its all children are not done?
			// Besides, if I think the below way creates all possible children 
			// even if the original function has some of them.
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i)
			{
				d(child_patch(kit.key())) += ((FuseT_CoeffT<T>*)((v_parameter->value[i].get())))->value.full_tensor_copy();
				delete v_parameter->value[i].data;
				temp.set_has_children_recursive(_coeffs_target, kit.key());	 // children
			}


			d = filter(d);

			//this->_result->get_impl()->get_coeffs().replace(key,temp);
			_coeffs_target.replace(key, temp);

			// can we use MADNESS_ASSERT??
			typename dcT::accessor acc;
			typename dcT::accessor acc_t;
			_coeffs.find(acc, key);
			_coeffs_target.find(acc_t, key);

			// acc = parent	
			if (acc->second.has_coeff())
			{
				const tensorT c = acc->second.coeff().full_tensor_copy();
				if (c.dim(0) == _k)
					d(_cdata.s0) += c;
				else
					d += c;
			}

			// tighter thresh for intenal nodes
			TensorArgs targs2 = this->_targs;	
			targs2.thresh*=0.1;
		
			// need the deep copy for contiguity
			coeffT ss = coeffT(copy(d(_cdata.s0)), targs2);

			if (key.level() > 0 && !_nonstandard)
				d(_cdata.s0) = 0.0;

			// insert either sum or difference coefficients
			if (this->_redundant)
			{
				acc_t->second.set_coeff(ss);
			}
			else
			{
				coeffT dd = coeffT(d, targs2);
				acc_t->second.set_coeff(dd);
			}

			//	Making Future form
			s_coeff->value = ss.full_tensor_copy();
			FuseTContainer<T> result2(static_cast<Base<T>*>(s_coeff));
			return result2;
		}
		else
		{
			//
			//	Leaf Nodes
			//
			KNODE temp;
			s_coeff			= new FuseT_CoeffT<T>();
		
			if (!_i1->get_impl()->get_coeffs().probe(key))
			{
				print (node);
			}

			s_coeff->value	= node.coeff().full_tensor_copy();
			
			

	
			if (!_keepleaves)
			{	
				temp.clear_coeff();
				this->_result->get_impl()->get_coeffs().replace(key,temp);
			}
			FuseTContainer<T> result(static_cast<Base<T>*>(s_coeff));
			return result;
		}

		// Create a Return
		cout<<"This should not have happenned"<<endl; 
		FuseTContainer<T> resultV;
		return resultV;
	}
	
	//
	//	helper functions for compressOp
	//
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

	// isDone
    template<typename T, std::size_t NDIM>
	bool CompressOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
		if(!isE1) {
			printf ("possible?: %d\n", isE1);
			return isE1;
		}

		bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf;
    }


}; /*fuset*/

#endif /* __fuset_CompressOp_h__ */
