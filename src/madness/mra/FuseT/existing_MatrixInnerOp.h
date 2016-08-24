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
#ifndef __MADNESS_MRA_FUSET_MATRIXINNER_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_MATRIXINNER_OP__INCLUDED__

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
	class MatrixInnerOp : public PrimitiveOp<T,NDIM> 
	{
		typedef Function<T,NDIM> KTREE;
		typedef FunctionNode<T,NDIM> KNODE;
		typedef Key<NDIM> keyT;
		typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>> dcT;
		typedef WorldObject<FunctionImpl<T,NDIM>> woT;
		typedef GenTensor<T> coeffT;
		typedef Tensor<T> tensorT;

    public:
									MatrixInnerOp	(string opName, KTREE* output, const std::vector<KTREE>& f, const std::vector<KTREE>& g, bool sym);
		FuseTContainer<T>			compute			(const keyT& key, const FuseTContainer<T> &s);

		bool notEmpty(map<int,bool>& notEmptyMap) const
		{
			unsigned long treeID = _i1->get_impl()->id().get_obj_id();
		    return  notEmptyMap[treeID];
		}

		bool						isDone			(const keyT& key) const;
		bool						isPre			() const { return true; } // false does not work. but It should be false.
		bool						needsParameter	() const { return true; }
        void						reduce			(World& world);

	public:	
		// MatrixInnerOpp
		Tensor<TENSOR_RESULT_TYPE(T, T)>* _r;

    private:
		//!Points to operand trees
		const KTREE*				_i1;
    
		//!Points to operand nodes of the tree
		KNODE						*_t1, *_t2;

		//!Variables for MatrixInnerOp
		std::vector<dcT>			_left_v_coeffs;
		std::vector<dcT>			_right_v_coeffs;
		//dcT&										_coeffs;
		//dcT&										_coeffs_target;

		bool										_sym;
        std::vector<const FunctionImpl<T,NDIM>* >	_left;
        std::vector<const FunctionImpl<T,NDIM>* >	_right;
	
		std::map<keyT, bool>	checkKeyDoneLeft;
		std::map<keyT, bool>	checkKeyDoneRight;

		bool										overallL;
		bool										overallR;

		int											_k;		// Wavelet order
    };
		
	// Constructor
	// World is needed for communication in the "compute" function
    template<typename T, std::size_t NDIM>
	MatrixInnerOp<T,NDIM>::MatrixInnerOp(string opName, KTREE* output, const std::vector<KTREE>& f, const std::vector<KTREE>& g, bool sym)
	: PrimitiveOp<T,NDIM>(opName, output, false, true)
	, _sym(sym)
	{
		this->_r = new Tensor<TENSOR_RESULT_TYPE(T,T)>(f.size(), g.size());

		for (unsigned int i=0; i<f.size(); i++)
			for (unsigned int j=0; j<g.size(); j++)
				(*this->_r)(i,j) = 0.0;

		for (unsigned int i=0; i<f.size(); i++) _left.push_back( f[i].get_impl().get() );
		for (unsigned int i=0; i<g.size(); i++) _right.push_back( g[i].get_impl().get() );

		for (unsigned int i=0; i<f.size(); i++) _left_v_coeffs.push_back( f[i].get_impl()->get_coeffs() );
		for (unsigned int j=0; j<g.size(); j++) _right_v_coeffs.push_back( g[j].get_impl()->get_coeffs() );

	    // dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();

		for (unsigned int i=0; i<f.size(); i++)
			this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(&f[i], true, true, false, false));

		for (unsigned int i=0; i<g.size(); i++)
			this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(&g[i], true, true, false, false));

	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,true,false,false));

		woT(f[0].world());
	}
	
	//
	//	it should hangle both a parent and a leaf node.
	//
	template <typename T, std::size_t NDIM>
	FuseTContainer<T>
	MatrixInnerOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s)
	{
		FuseT_VType<T>* possibleLists;
		
		// Processing for Paramter
		if (s.get() == 0)
		{
			possibleLists = new FuseT_VType<T>;
			for (int i=0; i<_left.size() + _right.size(); i++)
				possibleLists->value.push_back(1);
		}
		else
		{
			possibleLists  = new FuseT_VType<T>( ((FuseT_VType<T>*)s.get())->value );
		}		

		// Main --- O(M*N)
		bool tempLeft;
		bool tempRight;
		bool overallLeft = true;
		bool overallRight = true;
		FuseT_VType<T> whichNodesLeft;
		FuseT_VType<T> whichNodesRight;
		for (int i=0; i<_left.size(); i++)
		{
			if (possibleLists->value[i] != 0) 
			{
				const KNODE& fnode = _left_v_coeffs[i].find(key).get()->second;
				tempLeft	= _left_v_coeffs[i].find(key).get()->second.has_children();
				overallLeft = overallLeft && !tempLeft;

				if (fnode.has_coeff())
				{
					for (int j=0; j<_right.size(); j++)
					{
						if (possibleLists->value[_left.size() + j] != 0)
						{
							if (i == 0)
							{
								tempRight = _right_v_coeffs[j].find(key).get()->second.has_children();
								overallRight = overallRight && !tempRight;
							}

							const KNODE& gnode = _right_v_coeffs[j].find(key).get()->second;
							if (gnode.has_coeff()) //{} 
								(*this->_r)(i, j) += fnode.coeff().trace_conj(gnode.coeff());
						}
						else
						{
							if(i == 0)
							{
								tempRight = false;
								overallRight = overallRight && true;
							}
						}
						if (i == 0)
							whichNodesRight.value.push_back(tempRight);	
					}
				}
				else
				{
					if (i == 0) {
						for (int j=0; j<_right.size(); j++)
						{
							if (possibleLists->value[_left.size() + j] != 0)
							{
								tempRight = _right_v_coeffs[j].find(key).get()->second.has_children();
								overallRight = overallRight && !tempRight;	
							}
							else
							{	
								tempRight = false;
								overallRight = overallRight && true;
							}
							whichNodesRight.value.push_back(tempRight);	
						}
					}
				}
			}
			else
			{
				tempLeft = false;
				overallLeft = overallLeft && true;

				if (i == 0)
				{
					for (int j=0; j<_right.size(); j++)
					{
						if (possibleLists->value[_left.size() + j] != 0)
						{
							tempRight = _right_v_coeffs[j].find(key).get()->second.has_children();
							overallRight = overallRight && !tempRight;	
						}
						else
						{	
							tempRight = false;
							overallRight = overallRight && true;
						}
						whichNodesRight.value.push_back(tempRight);	
					}
				}
			}
			whichNodesLeft.value.push_back(tempLeft);
		}
		checkKeyDoneLeft.insert( std::pair<keyT,bool>(key,overallLeft) );
		checkKeyDoneRight.insert( std::pair<keyT,bool>(key,overallRight) );

		whichNodesLeft.value.insert(whichNodesLeft.value.end(), whichNodesRight.value.begin(), whichNodesRight.value.end());

		// 
		FuseT_VParameter<T> v_parameter;
		for (KeyChildIterator<NDIM> kit(key); kit; ++kit)
		{
			FuseTContainer<T> wrapper(static_cast<Base<T>*>(new FuseT_VType<T>(whichNodesLeft.value)));
			v_parameter.value.push_back(wrapper);
		}

		FuseTContainer<T> targets(static_cast<Base<T>*>(new FuseT_VParameter<T>(v_parameter.value)));
		return targets;
	}

	// isDone
    template<typename T, std::size_t NDIM>
	bool 
	MatrixInnerOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1;
		bool isE2;

		// O(M + N)
 		for (int i=0; i<_left.size(); i++)	
		{
			isE1 = _left[i]->get_coeffs().probe(key) || isE1;
		}
		if (!isE1) { std::cout<<key<<"!!!"<<std::endl; return isE1;}

		for (int i=0; i<_right.size(); i++)
		{	
			isE2 = _right[i]->get_coeffs().probe(key) || isE2;
		}
		if (!isE2) { std::cout<<key<<"???"<<std::endl;  return isE2;}

		if (checkKeyDoneLeft.find(key)->second)		return true;
		if (checkKeyDoneRight.find(key)->second)	return true;

		return false;
    }

    template<typename T, std::size_t NDIM>
	void  
	MatrixInnerOp<T,NDIM>::reduce(World& world){
        world.gop.sum(_r->ptr(),_left.size()*_right.size());
    }

}; /*fuset*/

#endif /* __fuset_MatrixInnerOp_h__ */
