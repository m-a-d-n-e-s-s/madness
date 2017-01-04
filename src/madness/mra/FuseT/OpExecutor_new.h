#ifndef __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__
#define __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__

#define DEBUG_1
#define DEBUG 0

#include <madness/world/MADworld.h>
#include <madness/world/world_object.h>
#include "PrimitiveOp.h"

namespace madness 
{
    /*!This class takes a primitive Op and computes the result tree*/
    template<typename T, std::size_t NDIM>
	class OpExecutor: public WorldObject<OpExecutor<T,NDIM>> 
	{
		typedef OpExecutor<T,NDIM> oeT;
		typedef Key<NDIM> keyT;
		typedef WorldObject<oeT> woT;
		typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>> dcT;
		typedef FunctionNode<T,NDIM> KNODE;
		typedef GenTensor<T> coeffT;

	public:	
		OpExecutor(World& world)
		: woT(world)
		,_world(world)
		{
			woT::process_pending();
		}

		Future<FuseTContainer<T>>	traverseTree	(const keyT key, const FuseTContainer<T> &s);
		FuseTContainer<T>			PostCompute		(const keyT key, const std::vector< Future < FuseTContainer<T> > > &v);

		//Starts the computation of pOp
		void				execute			(PrimitiveOp<T,NDIM>* pOp, bool isBlocking);

	private:
		World& _world;
		PrimitiveOp<T,NDIM>* _pOp;
		dcT* coeffs;
	
		//pointer to the resulting kary tree
		Function<T,NDIM>* _result;
    };
	
	//
	//
	template<typename T, std::size_t NDIM>
	Future< FuseTContainer<T> > 
	OpExecutor<T,NDIM>::traverseTree(const keyT key, const FuseTContainer<T> &s)
	{
/*
	    if(DEBUG){
		cout<<" Key : "<<key.level()<<" Translation : " ;
		for(auto a : key.translation())
		    cout<<a<<", ";
		cout<<endl;
	    }
*/	    
		FuseTContainer<T> temp;

		// Pre-Computation	
		if (_pOp->isPre())
		{
			temp = _pOp->compute(key, s);
		}

		//std::vector< Future<FuseTContainer<T> > > v;// = future_vector_factory<FuseTContainer<T> >(1<<NDIM);		
		std::vector< Future<FuseTContainer<T> > > v;// = future_vector_factory<FuseTContainer<T> >(1<<NDIM);		
		if (!_pOp->isDone(key))
		{
			if (!_pOp->isPre())
				v = future_vector_factory<FuseTContainer<T> >(1<<NDIM);		
			int i = 0;
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i)
			{
				const keyT& child = kit.key();
				
				if (_pOp->needsParameter())
				{
					if (_pOp->isPre())
					{
						woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child, ((FuseT_VParameter<T>*)(temp.get()))->value[i]);	
					} 
					else
					{
						v[i] = woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child, temp);
					}
				}
				else
				{
					woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child, temp);
				}
			}
		}

		// Post-Computation
		if (!_pOp->isPre())
		{
		    //std::cout<<"Running again"<<std::endl;
		    //Future<FuseTContainer<T>> sum(static_cast<Base<T>*>(new FuseT_VParameter<T>(k.value)));
			return woT::task(_world.rank(), &OpExecutor<T,NDIM>::PostCompute, key, v);	
			//	_world.gop.fence();
		}

		// for Pre-Computation (temporal)
		Future<FuseTContainer<T>> retValue(temp);
		return retValue;
	}

	//
	template<typename T, std::size_t NDIM>
	FuseTContainer<T>
    OpExecutor<T,NDIM>::PostCompute(const keyT key, const std::vector<Future<FuseTContainer<T>> > &v)
	{
	    vector<FuseTContainer <T> > temp;
	    for(auto a: v)
		temp.push_back(a.get());

	    return _pOp->compute(key, FuseTContainer<T>(static_cast<Base<T>*>(new FuseT_VParameter<T>(temp))));
	}

    //	execute with operators
    template<typename T, std::size_t NDIM>
	void 
	OpExecutor<T,NDIM>::execute(PrimitiveOp<T,NDIM>* pOp, bool hasParameters) 
    {
		_pOp	= pOp;
		_result = _pOp->_result;
		coeffs	= &_result->get_impl()->get_coeffs();

		if (_world.rank() == coeffs->owner(keyT(0)) )
		{
			FuseTContainer<T> initParameter;
			FuseTContainer<T> root = traverseTree(keyT(0), initParameter);
		}
		
		pOp->setComplete(true);
		_world.gop.fence();
		
		//for operators that require a reduction at the end
		if(pOp->_isReduce)
		    pOp->reduce(_world);
    }

}; /*fuset*/

#endif /*__fuset_OpExecutor_h__*/

