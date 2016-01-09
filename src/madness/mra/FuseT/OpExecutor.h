#ifndef __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__
#define __MADNESS_MRA_FUSET_OPEXECUTOR_H__INCLUDED__

#include "PrimitiveOp.h"
#include <madness/world/MADworld.h>
namespace madness {

    /*!This class takes a primitive Op and computes the result tree*/
    template<typename T, std::size_t NDIM>
	class OpExecutor: public WorldObject<OpExecutor<T,NDIM> > {

	typedef OpExecutor<T,NDIM> oeT;
	typedef Key<NDIM> keyT;
	typedef WorldObject<oeT> woT; ///< Base class world object type
        typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT; ///< Type of container holding the nodes
    public:	
	OpExecutor(World& world)
	    :woT(world)
	    ,_world(world)
	    
	{}

	//tree traversal for the primitive operator
	void traverseTree(keyT key);
	
	//Starts the computation of pOp
	void execute(PrimitiveOp<T,NDIM>* pOp);

	
    private:	
	World& _world;
	PrimitiveOp<T, NDIM>* _pOp;
	
	//pointer to the resulting kary tree
	Function<T,NDIM>* _result;
    };


    template<typename T, std::size_t NDIM>
	void 
	OpExecutor<T,NDIM>::traverseTree(keyT key)
    {
	if(_pOp->isPre())
	    _pOp->compute(key);
	if(!_pOp->isDone(key)) {
	    for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
		const keyT& child=kit.key();
		//Want to make a copy to be general
		//shared_ptr<PrimitiveOp<T,NDIM> > new_pOp = pOp;
		woT::task(_pOp->_result->get_impl()->get_coeffs().owner(child), &OpExecutor<T,NDIM>::traverseTree, child);
	    }
	}		
    if(!_pOp->isPre())
	_pOp->compute(key);	
}

template<typename T, std::size_t NDIM>
    void 
    OpExecutor<T,NDIM>::execute(PrimitiveOp<T,NDIM>* pOp) {
    //assign operator
    _pOp = pOp;
    _result = _pOp->_result;
    
    dcT* coeffs = &_result->get_impl()->get_coeffs();
    if (_world.rank() == coeffs->owner(keyT(0)) )
	traverseTree(keyT(0));
    
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

