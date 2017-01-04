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
#ifndef __MADNESS_MRA_FUSET_INNER_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_INNER_OP__INCLUDED__

#include "PrimitiveOp.h"
#include "../mra.h"
#include "../function_common_data.h"
#include "../../world/MADworld.h"
#include "../../tensor/tensor.h"
#include "../../tensor/gentensor.h"

namespace madness 
{
	template<typename T, std::size_t NDIM>
	class InnerOp : public PrimitiveOp<T,NDIM> 
	{
		typedef Function<T,NDIM> KTREE;
		typedef FunctionNode<T,NDIM> KNODE;	// identical to nodeT
		typedef Key<NDIM> keyT;
		typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>> dcT;
		typedef WorldObject<FunctionImpl<T,NDIM>> woT;	///< Base class world object type
		typedef GenTensor<T> coeffT;	//Type of tensor used to hold coeffs
		typedef Tensor<T> tensorT;
	    
    public:
		InnerOp(string opName, KTREE* output, const KTREE* i1, const KTREE*i2);
		FuseTContainer<T> compute(const keyT& key, const FuseTContainer<T> &s);

		bool notEmpty(map<int,bool>& notEmptyMap) const{
		    unsigned long treeID = _i1->get_impl()->id().get_obj_id();
		    unsigned long treeID2 = _i2->get_impl()->id().get_obj_id();
		    //cout<<"Checking for treeID : "<<treeID<<" and result is "<<notEmptyMap[treeID]<<endl;
		    return  notEmptyMap[treeID2] && notEmptyMap[treeID];
		}

		//
		bool isDone(const keyT& key) const;
		bool isPre() const { return true; }
		bool needsParameter() const { return false; }
                void reduce(World& world);
	public:	// for Inner Product (specific)
		T _sum;

	private:
		//!Points to operand trees
		const KTREE* _i1, *_i2;
    
		//!Points to operand nodes of the tree
		KNODE *_t1, *_t2;

		//!Variables for CompressOp
		dcT& _coeffs;
		dcT& _coeffs2;
		const FunctionCommonData<T,NDIM>& _cdata; 
		TensorArgs _targs;	

		// Wavelet order	
		int _k;
    };

		
	// Constructor
	// World is needed for communication in the "compute" function
    template<typename T, std::size_t NDIM>
	InnerOp<T,NDIM>::InnerOp(string opName, KTREE* output, const KTREE* i1, const KTREE* i2)
	: PrimitiveOp<T,NDIM>(opName, output, false,true)
		, _i1(i1)
	    , _i2(i2)
		, _cdata(FunctionCommonData<T,NDIM>::get(i1->k()))
		, _targs(i1->get_impl()->get_tensor_args())
		, _coeffs(i1->get_impl()->get_coeffs())
		, _coeffs2(i2->get_impl()->get_coeffs())
		, _k(i1->get_impl()->get_k())
		, _sum(0.0)
	{

	    // output is itself.
	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i2,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));
	    
		woT(i1->world());
	}

	//
	template<typename T, std::size_t NDIM>
	FuseTContainer<T>
	InnerOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s)
	{
		KNODE& fnode = _coeffs.find(key).get()->second;

		if (fnode.has_coeff())
		{
			//if (_i2->get_impl()->get_coeffs().probe(key))
			if (_coeffs.probe(key))
			{
				//const KNODE& gnode = _i2->get_impl()->get_coeffs().find(key).get()->second;
				const KNODE& gnode = _coeffs2.find(key).get()->second;
				if (gnode.has_coeff())
				{
					if (gnode.coeff().dim(0) != fnode.coeff().dim(0))
					{
						cerr<<"functions have different k or compress/reconstruct error"<<endl;
					}
					this->_sum += fnode.coeff().trace_conj(gnode.coeff());	
				}
			}
		}

		//FuseTContainer<T> result(static_cast<Base<T>*>(new FuseT_Type<T>(this->_sum)));
		FuseTContainer<T> result;
		return result;
	}

    // isDone
    template<typename T, std::size_t NDIM>
	bool 
	InnerOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
		if(!isE1){ std::cout<<"ouch!"<<std::endl; return isE1;}

		bool isE2 = _i2->get_impl()->get_coeffs().probe(key);
		if(!isE2) { std::cout<<"oops!" <<std::endl; return isE2; }

		bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
		if (isLeaf)
		    return isLeaf;

		bool isLeaf2 = !_i2->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf2;
    }
    
    // Reduction
    template<typename T, std::size_t NDIM>
	void  
	InnerOp<T,NDIM>::reduce(World& world){
	    world.gop.sum(_sum);	
    }
}; /*fuset*/

#endif /* __fuset_InnerOp_h__ */
