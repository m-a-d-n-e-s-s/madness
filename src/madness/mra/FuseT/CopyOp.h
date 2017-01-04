
#ifndef __MADNESS_MRA_FUSET_COPY_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_COPY_OP__INCLUDED__

#include "PrimitiveOp.h"

namespace madness 
{
    template<typename T, std::size_t NDIM>
	class CopyOp : public PrimitiveOp<T,NDIM> 
	{
		typedef Function<T,NDIM> KTREE;
		typedef FunctionNode<T,NDIM> KNODE;
		typedef Key<NDIM> keyT;
		typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;

		///< Type of container holding the nodes
	    
    public:
	CopyOp(string opName, KTREE* output, const KTREE* i1);

	//
	void compute(const keyT& key);

	//	Non-blocking
	FuseTContainer<T> compute(const keyT& key, const FuseTContainer<T> &s);
	
	// Blocking
	//Future<FuseTContainer<T>> postCompute(const keyT& key, const std::vector<Future<FuseTContainer<T>>> &s) { }

	bool isDone(const keyT& key) const;
	bool isPre() const { return true; }
	bool needsParameter() const { return false; }
	void reduce(World& world){}
    private:
	//!Points to operand trees
	const KTREE* _i1;
    
	//!Points to operand nodes of the tree
	KNODE *_t1, *_t2;
    };

    template<typename T, std::size_t NDIM>
	CopyOp<T,NDIM>::CopyOp(string opName, KTREE* output, const KTREE* i1)
	: PrimitiveOp<T,NDIM>(opName, output, false),
	_i1(i1)
	{
	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));
	}
    
	template<typename T, std::size_t NDIM>
	void		
	CopyOp<T,NDIM>::compute(const keyT& key)
	{
		typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();
        if (it == _i1->get_impl()->get_coeffs().end())	
		{
			cerr<<"This should not have happenned"<<endl;
			exit(0);
        }
        KNODE& node = it->second;
		this->_result->get_impl()->get_coeffs().replace(key,node);
	}

	template<typename T, std::size_t NDIM>
	FuseTContainer<T>
	CopyOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s)
	{
		typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();
        if (it == _i1->get_impl()->get_coeffs().end())	
		{
			cerr<<"This should not have happenned"<<endl;
			exit(0);
        }
	
        KNODE& node = it->second;
		this->_result->get_impl()->get_coeffs().replace(key,node);
	    
		FuseTContainer<T> temp;
		return temp;       
	}

    template<typename T, std::size_t NDIM>
	bool CopyOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
		if(!isE1) return isE1;
		bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf;
    }

}; /*fuset*/

#endif /* __fuset_CopyOp_h__ */
