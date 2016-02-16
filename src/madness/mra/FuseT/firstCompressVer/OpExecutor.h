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
		typedef OpExecutor<T,NDIM>		oeT;
		typedef Key<NDIM>				keyT;
		typedef WorldObject<oeT>		woT; ///< Base class world object type
		typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT; ///< Type of container holding the nodes
		typedef FunctionNode<T,NDIM> KNODE;
		dcT* coeffs;
	public:	
		OpExecutor(World& world)
			:woT(world)
			,_world(world)
		{}

		//tree traversal for the primitive operator (non-blocking)
		void traverseTree(keyT key);

		// tree traversal for the primitive operator (blocking)
		Future<GenTensor<T>> traverseTreeB		(keyT key);
		Future<GenTensor<T>> afterTraverseTreeB	(keyT key, const std::vector<Future<GenTensor<T>>>& v);

	
		//Starts the computation of pOp
		void execute(PrimitiveOp<T,NDIM>* pOp, bool isBlocking);

	
    private:	
		World&					_world;
		PrimitiveOp<T, NDIM>*	_pOp;
	
		//pointer to the resulting kary tree
		Function<T,NDIM>* _result;
	};


    template<typename T, std::size_t NDIM>
	void OpExecutor<T,NDIM>::traverseTree(keyT key)
    {
		if(_pOp->isPre())
		{
			_pOp->compute(key);
		}

		if(!_pOp->isDone(key)) 
		{
			std::cout << __func__ << ": key(" << key << ")" << std::endl;
			for (KeyChildIterator<NDIM> kit(key); kit; ++kit) 
			{	
				const keyT& child=kit.key();
				std::cout << "\t child: " << child << std::endl;

				woT::task(coeffs->owner(child), &OpExecutor<T,NDIM>::traverseTree, child);
			}
		}

		if(!_pOp->isPre())
		{
			_pOp->compute(key);	
		}
	}

	//Future<BaseParameters<T,NDIM>>
	template<typename T, std::size_t NDIM>
	Future<GenTensor<T>>
	OpExecutor<T,NDIM>::traverseTreeB(keyT key)
	{
		Future<GenTensor<T>> fromComputeB;

		if (_pOp->isPre())
			fromComputeB = _pOp->computeB(key);

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
			fromComputeB = _pOp->computeB(key);

		return fromComputeB;
	}

	//
	//
	//
	template<typename T, std::size_t NDIM>
	Future<GenTensor<T>> 
	OpExecutor<T,NDIM>::afterTraverseTreeB(keyT key, const std::vector<Future<GenTensor<T>>>& v)
	{
		//std::cout << "\t[" << __func__ << "]: key (" << key << ")" << std::endl;

		Future<GenTensor<T>> temp;
		temp = _pOp->afterComputeB(key, v);

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
				traverseTree(keyT(0));
			}
			else
			{
				Future<GenTensor<T>> temp = traverseTreeB(keyT(0));
			}
		}
 
		_world.gop.fence();
		pOp->setComplete(true);
	}


//namespace archive {
//    /// Serialize an AST
//    template <class Archive, typename T, std::size_t NDIM>
//	struct ArchiveStoreImpl< Archive, shared_ptr<PrimitiveOp<T,NDIM> > > {
//	
//	typedef Function<T,NDIM> KTREE;
//	typedef FunctionImpl<T,NDIM> implT;
//	
//	static void store(const Archive& s, const shared_ptr<PrimitiveOp<T, NDIM> >& t) {
//	    s & t.get();
//	}
//    };
//    
//    /// Deserialize an AST
//    template <class Archive, typename T, std::size_t NDIM>
//	struct ArchiveLoadImpl< Archive, shared_ptr<PrimitiveOp<T,NDIM> > > {
//	
//	typedef Function<T,NDIM> KTREE;
//	typedef FunctionImpl<T,NDIM> implT;
//	
//	static void load(const Archive& s, shared_ptr<PrimitiveOp<T,NDIM> >& t) {		
//	    PrimitiveOp<T,NDIM>* temp = NULL;	    
//	    s & temp;
//	    t = temp;
//	}
//    };
//    
//}


}; /*fuset*/

#endif /*__fuset_OpExecutor_h__*/

