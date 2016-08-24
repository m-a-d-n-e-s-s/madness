#ifndef __MADNESS_MRA_FUSET_NOTHING_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_NOTHING_OP__INCLUDED__

#include "PrimitiveOp.h"

namespace madness 
{
    template<typename T, std::size_t NDIM>
	class NothingOp : public PrimitiveOp<T,NDIM> {

	typedef Function<T,NDIM> KTREE;
	typedef FunctionNode<T,NDIM> KNODE;
	typedef Key<NDIM> keyT;
	typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;

///< Type of container holding the nodes
	    
    public:
	NothingOp(string opName, KTREE* output, const KTREE* i1);
    
		void						compute			(const keyT& key) { }
	FuseTContainer<T> compute(const keyT& key, const FuseTContainer<T> &s);

	//Future<GenTensor<T>> afterComputeB(const keyT& key, const std::vector<Future<GenTensor<T>>> &v) { }
	//Future<FuseTContainer<T>> postCompute(const keyT& key, const std::vector<Future<FuseTContainer<T>>> &s) { }

	bool isDone(const keyT& key) const;
	bool isPre() const { return true; }
	bool needsParameter() const { return false; }

    private:
	//!Points to operand trees
	const KTREE* _i1;
	dcT&								_coeffs;
    
	//!Points to operand nodes of the tree
	KNODE *_t1, *_t2;
    };

    template<typename T, std::size_t NDIM>
	NothingOp<T,NDIM>::NothingOp(string opName, KTREE* output, const KTREE* i1)
	: PrimitiveOp<T,NDIM>(opName, output, false),
	_i1(i1)
	, _coeffs(i1->get_impl()->get_coeffs())
	{
	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));
      
	}
   
	template<typename T, std::size_t NDIM> 
	FuseTContainer<T>
	NothingOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s)
	{
		std::cout<<"["<<__func__<<"] "<<key<<std::endl;
		typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();

		KNODE& node = _coeffs.find(key).get()->second;
		FuseTContainer<T> temp;
		return temp;
	}


    template<typename T, std::size_t NDIM>
	bool NothingOp<T,NDIM>::isDone(const keyT& key) const {
	
	bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
	if(!isE1) return isE1;
	bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
	return isLeaf;
	
    }

}; /*fuset*/

#endif /* __fuset_CopyOp_h__ */
