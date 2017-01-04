#ifndef __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__
#define __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__

#define DEBUG_1

#include <madness/world/MADworld.h>
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
		{ woT::process_pending();}

		FuseTContainer<T>	traverseTree	(const keyT key, FuseTContainer<T> &s);
		FuseTContainer<T>	PostCompute		(const keyT key, FuseTContainer<T> &v);

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
	FuseTContainer<T>
	OpExecutor<T,NDIM>::traverseTree(const keyT key, FuseTContainer<T> &s)
	{
		FuseTContainer<T> temp;

		// Pre-Computation	
		if (_pOp->isPre())
		{
			temp = _pOp->compute(key, s);
		}
		
		FuseT_VParameter<T>	k;
		if (!_pOp->isDone(key))
		{
			int i = 0;
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i)
			{
				const keyT& child = kit.key();
				if (temp.what() == WHAT_AM_I::FuseT_VCoeffT)
				{
					FuseTContainer<T> single(static_cast<Base<T>*>(new FuseT_CoeffT<T>( ((FuseT_VCoeffT<T>*)temp.get())->value[i]  )));
					woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child, single);	
				}
				else
				{	
					if (!_pOp->isPre())
						k.value.push_back(woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child, temp));
					else
						woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child, temp);
				}				
			}
		}

		// Post-Computation
		if (!_pOp->isPre())
		{
			Future<FuseTContainer<T>> sum(static_cast<Base<T>*>(new FuseT_VParameter<T>(k.value)));
			return woT::task(_world.rank(), &OpExecutor<T,NDIM>::PostCompute, key, sum);	
			_world.gop.fence();
		}

		// for Pre-Computation (temporal)
		return temp;
	}

	//
	template<typename T, std::size_t NDIM>
	FuseTContainer<T>
	OpExecutor<T,NDIM>::PostCompute(const keyT key, FuseTContainer<T> &v)
	{
		return _pOp->compute(key, v);
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
    }

}; /*fuset*/

#endif /*__fuset_OpExecutor_h__*/

