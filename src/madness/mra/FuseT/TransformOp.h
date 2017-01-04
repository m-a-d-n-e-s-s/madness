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
#ifndef __MADNESS_MRA_FUSET_TRANSFORM_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_TRANSFORM_OP__INCLUDED__

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
	class TransformOp : public PrimitiveOp<T,NDIM> 
	{
		typedef Function<T,NDIM> KTREE;
		typedef FunctionNode<T,NDIM> KNODE;
		typedef Key<NDIM> keyT;
		typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>> dcT;
		typedef WorldObject<FunctionImpl<T,NDIM>> woT;
		typedef GenTensor<T> coeffT;
		typedef Tensor<T> tensorT;

    public:
									TransformOp	(string opName, std::vector<KTREE>& output, const std::vector<KTREE>& v, const DistributedMatrix<T>& g, bool sym);
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
		
		std::map<keyT, int>		candidatesLeft;
		std::map<keyT, int>		candidatesRight;	

		int						_k;		// Wavelet order
    };
		
	// Constructor
	// World is needed for communication in the "compute" function
    template<typename T, std::size_t NDIM>
	TransformOp<T,NDIM>::TransformOp(string opName, std::vector<KTREE>& output, const std::vector<KTREE>& v, const DistributedMatrix<T>& c, bool sym)
	: PrimitiveOp<T,NDIM>(opName, output, false, true)
	, _sym(sym)
	{
		long n = v.size();
		long m = c.rowdim();

		MADNESS_ASSERT(n == c.columndim());


		output = zero_functions_compressed<T,NDIM>(v[0].world(), m);

//		this->_r = new Tensor<TENSOR_RESULT_TYPE(T,T)>(f.size(), g.size());

//		for (unsigned int i=0; i<f.size(); i++)
//			for (unsigned int j=0; j<g.size(); j++)
//				(*this->_r)(i,j) = 0.0;

//		for (unsigned int i=0; i<f.size(); i++) _left.push_back( f[i].get_impl().get() );
	//	for (unsigned int i=0; i<g.size(); i++) _right.push_back( g[i].get_impl().get() );

//		for (unsigned int i=0; i<f.size(); i++) _left_v_coeffs.push_back( f[i].get_impl()->get_coeffs() );
	//	for (unsigned int j=0; j<g.size(); j++) _right_v_coeffs.push_back( g[j].get_impl()->get_coeffs() );

	    // dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	  //  this->_OpID = output->get_impl()->id().get_obj_id();

		for (unsigned int i=0; i<v.size(); i++)
			this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(&v[i], true, true, false, false));

	//	for (unsigned int i=0; i<g.size(); i++)
	//		this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(&g[i], true, true, false, false));

	    //this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,true,false,false));

		woT(v[0].world());
	}
	
	//
	//	it should hangle both a parent and a leaf node.
	//
	template <typename T, std::size_t NDIM>
	FuseTContainer<T>
	TransformOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s)
	{
/*
		FuseT_VParameter<T>*	inheritedWhole;
		FuseT_VType<T>*			inheritedLeft;
		FuseT_VType<T>*			inheritedRight;
		
		// Processing for Paramter
		if (s.get() == 0)
		{
			inheritedLeft	= new FuseT_VType<T>;
			inheritedRight	= new FuseT_VType<T>;

			for (unsigned int i=0; i<_left.size(); i++)
				inheritedLeft->value.push_back(i);
			for (unsigned int i=0; i<_right.size(); i++)
				inheritedRight->value.push_back(i);
		}
		else
		{
			inheritedWhole	= new FuseT_VParameter<T>( ((FuseT_VParameter<T>*)s.get())->value );

			inheritedLeft	= new FuseT_VType<T>(((FuseT_VType<T>*)(((inheritedWhole->value[0]).get())))->value);
			inheritedRight	= new FuseT_VType<T>(((FuseT_VType<T>*)(((inheritedWhole->value[1]).get())))->value);
		}		

		FuseT_VType<T> whichNodesLeft;		// value = std::vector<int>
		FuseT_VType<T> whichNodesRight;
	
		unsigned int indexLeft;
		unsigned int indexRight;
		unsigned int leftSize	= inheritedLeft->value.size();
		unsigned int rightSize	= inheritedRight->value.size();

		// The Pre-Computatio
		// Assumption: the size of coefficient --> 16*16*16 = 4096
		double* A = (double*)malloc(sizeof(double)*16*16*16*leftSize);
		double* B = (double*)malloc(sizeof(double)*16*16*16*rightSize);
		double* C = (double*)malloc(sizeof(double)*leftSize*rightSize);
		unsigned int k,l,m;

		//
		for (unsigned int i=0; i<leftSize; i++)
		{
			indexLeft = inheritedLeft->value[i];
			const KNODE& fnode = _left_v_coeffs[indexLeft].find(key).get()->second;
				
			if (_left_v_coeffs[indexLeft].find(key).get()->second.has_children())
				whichNodesLeft.value.push_back(indexLeft);

			// 3D array to 1D array with i for fnode and j for gnode
			if (fnode.has_coeff())
			{
				for (k=0; k<16; k++) 
					for (l=0; l<16; l++) 
						for (m=0; m<16; m++) 
							A[i*16*16*16 + k*16*16 + l*16 + m] = (fnode.coeff())(k,l,m);	
			}
			else
			{
				for (k=0; k<16; k++) 
					for (l=0; l<16; l++) 
						for (m=0; m<16; m++) 
							A[i*16*16*16 + k*16*16 + l*16 + m] = 0.0;
			}
		}

	
		//
		for (unsigned int i=0; i<rightSize; i++)
		{
			indexRight = inheritedRight->value[i];
			const KNODE& gnode = _right_v_coeffs[indexRight].find(key).get()->second;

			if (_right_v_coeffs[indexRight].find(key).get()->second.has_children())
				whichNodesRight.value.push_back(indexRight);

			// 3D array to 1D array with i for fnode and j for gnode
			if (gnode.has_coeff())
			{
				for (k=0; k<16; k++) 
					for (l=0; l<16; l++) 
						for (m=0; m<16; m++) 
							B[i*16*16*16 + k*16*16 + l*16 + m] = (gnode.coeff())(k,l,m);	
			}
			else
			{
				for (k=0; k<16; k++) 
					for (l=0; l<16; l++) 
						for (m=0; m<16; m++) 
							B[i*16*16*16 + k*16*16 + l*16 + m] = 0.0;
			}
		}

		// 
		for (k=0; k<leftSize; k++)
			for (l=0; l<rightSize; l++)
				C[k*rightSize + l] = 0.0;

		//	The Actual-Computation
		//	Return: left.size() * right.size();
		//
		//	A [leftSize * 4096] 
		//	B [rightSize * 4096] 
		//	C [leftSize*rightSize]
		//
		cblas::gemm(cblas::CBLAS_TRANSPOSE::Trans, cblas::CBLAS_TRANSPOSE::NoTrans, leftSize, rightSize, 16*16*16, 1, A, 16*16*16, B, 16*16*16, 1, C, leftSize);

		// The Post-Computation
		for (k=0; k<leftSize; k++)
		{	
			indexLeft = inheritedLeft->value[k];	
			for (l=0; l<rightSize; l++)
			{
				indexRight = inheritedRight->value[l];
				(*this->_r)(indexLeft, indexRight) += C[k + l*leftSize];	// k*rightSize + l --> row-major
			}
		}

		delete A;
		delete B;
		delete C;



		if (whichNodesLeft.value.size() == 0)
			checkKeyDoneLeft.insert(std::pair<keyT,bool>(key,true));
		else
			checkKeyDoneLeft.insert(std::pair<keyT,bool>(key,false));

		if (whichNodesRight.value.size() == 0)
			checkKeyDoneRight.insert(std::pair<keyT,bool>(key,true));
		else
			checkKeyDoneRight.insert(std::pair<keyT,bool>(key,false));


		// 
		FuseT_VParameter<T> v_parameter;
		FuseT_VParameter<T> inner_parameter;
	
		FuseTContainer<T>	candiParameter_L(static_cast<Base<T>*> (new FuseT_VType<T>(whichNodesLeft.value)));	
		FuseTContainer<T>	candiParameter_R(static_cast<Base<T>*> (new FuseT_VType<T>(whichNodesRight.value)));	
		inner_parameter.value.push_back(candiParameter_L);
		inner_parameter.value.push_back(candiParameter_R);

		for (KeyChildIterator<NDIM> kit(key); kit; ++kit)
		{
			FuseTContainer<T> wrapper(static_cast<Base<T>*>(new FuseT_VParameter<T>(inner_parameter.value)));
			v_parameter.value.push_back(wrapper);
		}

		// Return Parameters
		FuseTContainer<T> targets(static_cast<Base<T>*>(new FuseT_VParameter<T>(v_parameter.value)));
		return targets;
*/
	}

	// isDone
    template<typename T, std::size_t NDIM>
	bool 
	TransformOp<T,NDIM>::isDone(const keyT& key) const 
	{
/*
		bool isE1;
		bool isE2;

		// O(M + N)
 		for (unsigned int i=0; i<_left.size(); i++)	
		{
			isE1 = _left[i]->get_coeffs().probe(key) || isE1;
		}
		if (!isE1) { std::cout<<key<<"!!!"<<std::endl; return isE1;}

		for (unsigned int i=0; i<_right.size(); i++)
		{	
			isE2 = _right[i]->get_coeffs().probe(key) || isE2;
		}
		if (!isE2) { std::cout<<key<<"???"<<std::endl;  return isE2;}

		if (checkKeyDoneLeft.find(key)->second)		return true;
		if (checkKeyDoneRight.find(key)->second)	return true;
*/
		return false;
    }

    template<typename T, std::size_t NDIM>
	void  
	TransformOp<T,NDIM>::reduce(World& world){
  //      world.gop.sum(_r->ptr(),_left.size()*_right.size());
    }

}; /*fuset*/

#endif /* __fuset_TransformOp_h__ */
