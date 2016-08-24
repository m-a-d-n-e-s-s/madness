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
	MatrixInnerOp	(string opName, KTREE* output, const std::vector<KTREE>& f, const std::vector<KTREE>& g, bool sym, bool dgemm);
	FuseTContainer<T>			compute			(const keyT& key, const FuseTContainer<T> &s);

	bool notEmpty(map<int,bool>& notEmptyMap) const
	{
	    for(auto a:_left)
		{
			if(notEmptyMap[a->id().get_obj_id()]){
				for(auto b:_right){
					if(notEmptyMap[b->id().get_obj_id()]){
						return true;
					}
				}
			}	
	    }
	    return false;
	    //unsigned long treeID = _i1->get_impl()->id().get_obj_id();
	    //return  notEmptyMap[treeID];
	}

	bool						isDone			(const keyT& key) const;
	bool						isPre			() const { return true; } // false does not work. but It should be false.
	bool						needsParameter	() const { return false; }
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

	bool						_sym;
	bool						_dgemm;
	std::vector<const FunctionImpl<T,NDIM>* >	_left;
	std::vector<const FunctionImpl<T,NDIM>* >	_right;
	
		
	int						_k;		// Wavelet order
    };
		
    // Constructor
    // World is needed for communication in the "compute" function
    template<typename T, std::size_t NDIM>
	MatrixInnerOp<T,NDIM>::MatrixInnerOp(string opName, KTREE* output, const std::vector<KTREE>& f, const std::vector<KTREE>& g, bool sym, bool dgemm)
	: PrimitiveOp<T,NDIM>(opName, output, false, true)
	, _sym(sym)
	, _dgemm(dgemm)
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
			this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(&f[i], true, false, false, false));

		for (unsigned int i=0; i<g.size(); i++)
			this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(&g[i], true, false, false, false));

		this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));

		woT(f[0].world());
    }
	
    //
    //	it should hangle both a parent and a leaf node.
    //
    template <typename T, std::size_t NDIM>
	FuseTContainer<T>
	MatrixInnerOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s)
    {
	//	std::cout<<__func__<<", key: "<<key<<std::endl;

		FuseT_VType<T>* inheritedLeft;
		FuseT_VType<T>*	inheritedRight;

		inheritedLeft	= new FuseT_VType<T>;
		inheritedRight	= new FuseT_VType<T>;

		for (unsigned int i=0; i<_left.size(); i++) {
			if(_left[i]->get_coeffs().probe(key)) {
				inheritedLeft->value.push_back(i);
			}
		}

		for (unsigned int i=0; i<_right.size(); i++) {
			if(_right[i]->get_coeffs().probe(key))  {
				inheritedRight->value.push_back(i);
			}
		}

		if(inheritedLeft->value.empty())
		{
			return FuseTContainer<T>();
		}

		if(inheritedRight->value.empty()) 
		{
			return FuseTContainer<T>();
		}


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
	//
		for (unsigned int i=0; i<leftSize; i++)
		{
			indexLeft = inheritedLeft->value[i];
			const KNODE& fnode = _left_v_coeffs[indexLeft].find(key).get()->second;


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

		//
		//
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


		return FuseTContainer<T>();
	}

	// isDone
	template<typename T, std::size_t NDIM>
	bool 
	MatrixInnerOp<T,NDIM>::isDone(const keyT& key) const 
	{
	    bool isE1= false;
	    bool isE2 =false;


	    for (unsigned int i=0; i<_left.size(); i++) {
		if(_left[i]->get_coeffs().probe(key))
		    isE1 = isE1 || _left_v_coeffs[i].find(key).get()->second.has_children();
	    }
	    for (unsigned int i=0; i<_right.size(); i++) {
		if(_right[i]->get_coeffs().probe(key))
		    isE2 = isE2 || _right_v_coeffs[i].find(key).get()->second.has_children();
	    }
	    
	    return !(isE1 && isE2);

	}

	template<typename T, std::size_t NDIM>
	void  
	MatrixInnerOp<T,NDIM>::reduce(World& world){
		world.gop.sum(_r->ptr(),_left.size()*_right.size());
	}

}; /*fuset*/

#endif /* __fuset_MatrixInnerOp_h__ */
