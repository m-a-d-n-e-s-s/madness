#ifndef __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__
#define __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__

#define DEBUG_1

#include "PrimitiveOp.h"
#include "BaseParameters.h"
#include <madness/world/MADworld.h>

namespace madness 
{
	//class BaseParameters;

    /*!This class takes a primitive Op and computes the result tree*/
    template<typename T, std::size_t NDIM>
	class OpExecutor: public WorldObject< OpExecutor<T,NDIM> > 
	{
		typedef OpExecutor<T,NDIM>									oeT;
		typedef Key<NDIM>											keyT;
		typedef WorldObject<oeT>									woT; ///< Base class world object type
		typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >	dcT; ///< Type of container holding the nodes
		typedef FunctionNode<T,NDIM>								KNODE;
		typedef GenTensor<T>										coeffT;	//Type of tensor used to hold coeffs
		dcT*														coeffs;

	public:	
		OpExecutor(World& world)
		: woT(world)
		,_world(world)
		{}

		//tree traversal for the primitive operator (non-blocking)
		void traverseTree(keyT key, const BaseParameters<T,NDIM> &s);

		// tree traversal for the primitive operator (blocking)
		Future<GenTensor<T>> traverseTreeB		(keyT key);
		Future<BaseParameters<T,NDIM>> 	traverseTreeC(keyT key);
		Future<GenTensor<T>> afterTraverseTreeB	(keyT key, const std::vector<Future<GenTensor<T>>>& v);
		Future<BaseParameters<T,NDIM>> afterTraverseTreeC(keyT key, const std::vector<Future<BaseParameters<T,NDIM>>>& v);

	
		//Starts the computation of pOp
		void execute(PrimitiveOp<T,NDIM>* pOp, bool isBlocking);
	
    private:	
		World&					_world;
		PrimitiveOp<T, NDIM>*	_pOp;
	
		//pointer to the resulting kary tree
		Function<T,NDIM>*		_result;
	};

	//
	//	based on Non-blocking
	//
    template<typename T, std::size_t NDIM>
	void 
	OpExecutor<T,NDIM>::traverseTree(keyT key, const BaseParameters<T,NDIM> &s)
    {
		BaseParameters<T,NDIM> temp;

		if(_pOp->isPre())
		{
			//std::cout << "pre: " << __func__ << ": key(" << key << std::endl;
			temp = _pOp->compute(key, s);	// it should return!!!
		}

		if(!_pOp->isDone(key)) 
		{
			//std::cout << __func__ << ": key(" << key << ")" << std::endl;
			int i = 0;
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i) 
			{	
				const keyT& child=kit.key();
				//std::cout << "\t child" << i << ": " << child << std::endl;

				BaseParameters<T,NDIM> spawn;
				if (temp._isVector == true)
				{
					spawn._coeff = temp._vCoeff[i];
				}
			
				woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child, spawn);
			}
		}

		if(!_pOp->isPre())
		{
			temp = _pOp->compute(key, s);
		}
	}

	//
	//	based on Blocking
	//
	template<typename T, std::size_t NDIM>
	Future<GenTensor<T>> 
	OpExecutor<T,NDIM>::traverseTreeB(keyT key)
	{
		Future<GenTensor<T>> fromComputeB;
		BaseParameters<T,NDIM> tempB;

		if (_pOp->isPre())
		{
			fromComputeB = _pOp->computeB(key, tempB);
		}

		if (!_pOp->isDone(key))
		{
			std::vector<Future<GenTensor<T>>> v = future_vector_factory<GenTensor<T>>(1<<NDIM);

			int i = 0;
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i)	
			{
				const keyT& child = kit.key();
				v[i] = woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTreeB, child);
			}
			// After receiving data from children
			fromComputeB = woT::task(this->_world.rank(), &OpExecutor<T,NDIM>::afterTraverseTreeB, key, v);
		}

		if (!_pOp->isPre())
		{
			fromComputeB = _pOp->computeB(key, tempB);
		}

		return fromComputeB;
	}

	// will replace B by C
	template<typename T, std::size_t NDIM>
	Future<BaseParameters<T,NDIM>> 
	OpExecutor<T,NDIM>::traverseTreeC(keyT key)
	{
		printf ("[%s] ", __func__);
		std::cout << key << std::endl;
		Future<BaseParameters<T,NDIM>>	fromTraverseTree;	// for a spatial tree (traverseTree level)
		BaseParameters<T,NDIM>			tempB;				// for a node in traverseTree.
	
		if (_pOp->isPre())
			fromTraverseTree = _pOp->computeC(key, tempB);

		if (!_pOp->isDone(key))
		{
			std::vector<Future<BaseParameters<T,NDIM>>> v = future_vector_factory<BaseParameters<T,NDIM>>(1<<NDIM);
		
			int i = 0;
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i)
			{
				const keyT& child = kit.key();
				v[i] = woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTreeC, child);
			}
			// Process After receiving data from children
			fromTraverseTree = woT::task(this->_world.rank(), &OpExecutor<T,NDIM>::afterTraverseTreeC, key, v);
		}

		if (!_pOp->isPre())
			fromTraverseTree = _pOp->computeC(key, tempB);

		return fromTraverseTree;
	}

	//
	//
	//
	template<typename T, std::size_t NDIM>
	Future<GenTensor<T>> 
	OpExecutor<T,NDIM>::afterTraverseTreeB(keyT key, const std::vector<Future<GenTensor<T>>>& v)
	{
		Future<GenTensor<T>> temp;
		temp = _pOp->afterComputeB(key, v);

		return temp;
	}

	template<typename T, std::size_t NDIM>
	Future<BaseParameters<T,NDIM>> 
	OpExecutor<T,NDIM>::afterTraverseTreeC(keyT key, const std::vector<Future<BaseParameters<T,NDIM>>>& v)
	{
		Future<BaseParameters<T,NDIM>> temp;
		temp = _pOp->afterComputeC(key, v);
		return temp;
	}
	
	//
	//	execute with operators
	//
	template<typename T, std::size_t NDIM>
    void 
    OpExecutor<T,NDIM>::execute(PrimitiveOp<T,NDIM>* pOp, bool isBlocking) 
	{
		//assign operator
		_pOp	= pOp;
		_result = _pOp->_result;
		
		//dcT* coeffs = &_result->get_impl()->get_coeffs();
		coeffs	= &_result->get_impl()->get_coeffs();

		if (_world.rank() == coeffs->owner(keyT(0)) )
		{
			if (isBlocking == false)
			{
				BaseParameters<T,NDIM> temp;
				temp._coeff = coeffT();
				traverseTree(keyT(0), temp);
			}
			else
			{
				//Future<GenTensor<T>> temp = traverseTreeC(keyT(0));
				Future<BaseParameters<T,NDIM>> temp = traverseTreeC(keyT(0));
			}
		}
 
		_world.gop.fence();
		pOp->setComplete(true);
	}

}; /*fuset*/

#endif /*__fuset_OpExecutor_h__*/

