
// Ghaly
//
// Reconstructs the function, transforming into scaling function basis. Possible non-blocking comm.
//
// By default fence=true meaning that this operation completes before returning,
// otherwise if fence=false it returns without fencing and the user must invoke
// world.gop.fence() to assure global eompletion before using the function
// for other purposes.
//
// Noop if already reconstructed or if no initialized.
// 
// Since reconstuction/compression do not discard information we define them
// as const ... "logical constness" not "bitwise constness"
//

/* <mra.h>
const Function<T,NDIM>& reconstruct(bool fence = true) const {
	PROFILE_MEMBER_FUNC(Function);
	if (!impl || !is_compressed()) return *this;
	const_cast<Function<T,NDIM>*>(this)->impl->reconstruct(fence);
	if (fence && VERIFY_TREE) verify_tree();	// Must be after in case nonstandard
	return *this;
} 
*/

// Verifies the tree data structure ... global sync implied
/*
void verify_tree() const {
	PROFILE_MEMBER_FUNC(Function);
	if (impl) impl->verify_tree();
}
*/

// Returns true if compressed, false otherwise. No communication.
//
// If the function is not initialized, returns false.
/*
bool is_compressed() const {
	PROFILE_MEMBER_FUNC(Function);
	if (imple)
		return impl->is_compressed();
	else
		return false;
}

*/

/* <funcimpl.h>
void reconstruct(bool fence);
*/

// Returns true if the function is compressed.
/*
bool is_compressed() const;
*/

#ifndef __MADNESS_MRA_FUSET_RECONSTRUCT_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_RECONSTRUCT_OP__INCLUDED__

#include "PrimitiveOp.h"

namespace madness {
    template<typename T, std::size_t NDIM>
	class ReconstructOp : public PrimitiveOp<T,NDIM> {

	typedef Function<T,NDIM> KTREE;
	typedef FunctionNode<T,NDIM> KNODE;
	typedef Key<NDIM> keyT;
	typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;

///< Type of container holding the nodes
	    
    public:
	ReconstructOp(string opName, KTREE* output, const KTREE* i1);
    
	void compute(const keyT& key);

	bool isDone(const keyT& key) const;

	bool isPre() const { return true; }

    private:
	//!Points to operand trees
	const KTREE* _i1;
    
	//!Points to operand nodes of the tree
	KNODE *_t1, *_t2;
    };

    template<typename T, std::size_t NDIM>
	ReconstructOp<T,NDIM>::ReconstructOp(string opName, KTREE* output, const KTREE* i1)
	: PrimitiveOp<T,NDIM>(opName, output, false),
	_i1(i1)
	{
	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));
      
	}
    
    template<typename T, std::size_t NDIM>
	void ReconstructOp<T,NDIM>::compute(const keyT& key) {

	typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();

        if (it == _i1->get_impl()->get_coeffs().end()) {
	    cerr<<"This should not have happenned"<<endl;
	    exit(0);
        }
	
        KNODE& node = it->second;
	this->_result->implP()->get_coeffs().replace(key,node);
	           
    }

    template<typename T, std::size_t NDIM>
	bool ReconstructOp<T,NDIM>::isDone(const keyT& key) const {
	
	bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
	if(!isE1) return isE1;
	bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
	return isLeaf;
	
    }

}; /*fuset*/

#endif /* __fuset_ReconstructOp_h__ */
