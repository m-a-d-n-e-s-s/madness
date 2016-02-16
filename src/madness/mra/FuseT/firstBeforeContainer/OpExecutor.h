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
		void							traverseTree		(keyT key, const BaseParameters<T,NDIM> &s);

		// tree traversal for the primitive operator (blocking)
		Future<BaseParameters<T,NDIM>>	afterTraverseTreeB	(keyT key, const std::vector<Future<BaseParameters<T,NDIM>>>& v);
		Future<BaseParameters<T,NDIM>>	traverseTreeB		(keyT key, const BaseParameters<T,NDIM> &s);
	
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
	//	Through this function, we can pass the data from parent to childrun,
	//	then, the parent can accumulate the data from childrun.
	//
	template<typename T, std::size_t NDIM>
	Future<BaseParameters<T,NDIM>>
	OpExecutor<T,NDIM>::traverseTreeB(keyT key, const BaseParameters<T,NDIM> &s)
	{
		Future<BaseParameters<T,NDIM>>	fromTraverseTree;
		
		if (_pOp->isPre())
			fromTraverseTree = _pOp->computeC(key, s);
		
		if (!_pOp->isDone(key))
		{
			std::vector<Future<BaseParameters<T,NDIM>>> v = future_vector_factory<BaseParameters<T,NDIM>>(1<<NDIM);
		
			int i = 0;
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i)
			{
				const keyT& child = kit.key();
				v[i] = woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTreeB, child, s);
			}

			fromTraverseTree = woT::task(this->_world.rank(), &OpExecutor<T,NDIM>::afterTraverseTreeB, key, v);
		}
	
		if (!_pOp->isPre())
			fromTraverseTree = _pOp->computeC(key, s);

		return fromTraverseTree;	
	}

	// 
	//	
	//
	template<typename T, std::size_t NDIM>
	Future<BaseParameters<T,NDIM>> 
	OpExecutor<T,NDIM>::afterTraverseTreeB(keyT key, const std::vector<Future<BaseParameters<T,NDIM>>>& v)
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
		
		coeffs	= &_result->get_impl()->get_coeffs();

		if (_world.rank() == coeffs->owner(keyT(0)) )
		{
			BaseParameters<T,NDIM> temp;
			if (isBlocking == false)
			{
				std::cout << "Non-Blocking TraverseTree" << std::endl;
				temp._coeff = coeffT();
				traverseTree(keyT(0), temp);
			}
			else
			{
				std::cout << "Blocking TraverseTree" << std::endl;
				temp._coeff = coeffT();
				Future<BaseParameters<T,NDIM>> ghaly = traverseTreeB(keyT(0), temp);

				//
				//	IT IS POSSIBLE TO USE THE RETURN VALUE AFTER THE ROOT IS FINISHED.
				//
			}
		}
 
		_world.gop.fence();
		pOp->setComplete(true);
	}

}; /*fuset*/

#endif /*__fuset_OpExecutor_h__*/

