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
#ifndef __MADNESS_MRA_FUSET_INNER_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_INNER_OP__INCLUDED__

#include "PrimitiveOp.h"
#include "BaseParameters.h"
#include "../mra.h"
#include "../function_common_data.h"
#include "../../world/MADworld.h"
#include "../../tensor/tensor.h"
#include "../../tensor/gentensor.h"

namespace madness 
{
	template<typename T, std::size_t NDIM>
	class InnerOp : public PrimitiveOp<T,NDIM> 
	{
	typedef Function<T,NDIM>									KTREE;
	typedef FunctionNode<T,NDIM>								KNODE;	// identical to nodeT
	typedef Key<NDIM>											keyT;
	typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM>>	dcT;
	typedef WorldObject< FunctionImpl<T,NDIM> >					woT;	///< Base class world object type

	typedef GenTensor<T>										coeffT;	//Type of tensor used to hold coeffs
	typedef Tensor<T>											tensorT;

	///< Type of container holding the nodes
	    
    public:
										InnerOp			(string opName, KTREE* output, const KTREE* i1);
		BaseParameters<T,NDIM>			compute			(const keyT& key, const BaseParameters<T,NDIM> &s) { }

		//
		Future<BaseParameters<T,NDIM>>	computeC		(const keyT& key, const BaseParameters<T, NDIM> &s);
		Future<BaseParameters<T,NDIM>>	afterComputeC	(const keyT& key, const std::vector<Future<BaseParameters<T,NDIM>>> &v);
		
		bool							isDone			(const keyT& key) const;
		bool							isPre			() const { return true; }

	public:	// for CompressOp

    private:
		//!Points to operand trees
		const KTREE*						_i1;
    
		//!Points to operand nodes of the tree
		KNODE								*_t1, *_t2;

		//!Variables for CompressOp
		dcT&								_coeffs;
		const FunctionCommonData<T,NDIM>&	_cdata; 
		TensorArgs							_targs;	

		int									_k;				// Wavelet order	
    };

		
	// Constructor
	// World is needed for communication in the "compute" function
    template<typename T, std::size_t NDIM>
	InnerOp<T,NDIM>::InnerOp(string opName, KTREE* output, const KTREE* i1)
		: PrimitiveOp<T,NDIM>(opName, output, false)
		, _i1(i1)
		, _cdata(FunctionCommonData<T,NDIM>::get(i1->k()))
		, _targs(i1->get_impl()->get_tensor_args())
		, _coeffs(i1->get_impl()->get_coeffs())
		, _k(i1->get_impl()->get_k())
	{
		// output is itself.

	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));

		woT(i1->world());
	}


    // Called by result function to differentiate f
/*
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::diff(const DerivativeBase<T,NDIM>* D, const implT* f, bool fence) {
        typedef std::pair<keyT,coeffT> argT;
        typename dcT::const_iterator end = f->coeffs.end();
        for (typename dcT::const_iterator it=f->coeffs.begin(); it!=end; ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.has_coeff()) {
                Future<argT> left  = D->find_neighbor(f, key,-1);
                argT center(key,node.coeff());
                Future<argT> right = D->find_neighbor(f, key, 1);
                world.taskq.add(*this, &implT::do_diff1, D, f, key, left, center, right, TaskAttributes::hipri());
            }
            else {
                coeffs.replace(key,nodeT(coeffT(),true)); // Empty internal node
            }
        }
        if (fence) world.gop.fence();
    }
*/
	//
	//	
	//
	template <typename T, std::size_t NDIM>
	Future<BaseParameters<T,NDIM>> 
	InnerOp<T,NDIM>::computeC(const keyT& key, const BaseParameters<T, NDIM> &s) 
	{

 
	}

	//
	//
	//
	template <typename T, std::size_t NDIM>
	Future<BaseParameters<T,NDIM>>
	InnerOp<T,NDIM>::afterComputeC(const keyT& key, const std::vector<Future<BaseParameters<T,NDIM>>> &v)
	{


	}

	// isDone
    template<typename T, std::size_t NDIM>
	bool 
	InnerOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
		if(!isE1) return isE1;
		bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf;
    }

}; /*fuset*/


/*
        Key<NDIM> neighbor(const keyT& key, int step) const {
            Vector<Translation,NDIM> l = key.translation();
            l[axis] += step;
            if (!enforce_bc(bc(axis,0), bc(axis,1), key.level(), l[axis])) {
                return keyT::invalid();
            }
            else {
                return keyT(key.level(),l);
            }
        }

        Future<argT>
        find_neighbor(const implT* f, const Key<NDIM>& key, int step) const {
            keyT neigh = neighbor(key, step);
            if (neigh.is_invalid()) {
                return Future<argT>(argT(neigh,coeffT(vk,f->get_tensor_args()))); // Zero bc
            }
            else {
                Future<argT> result;
                f->task(f->get_coeffs().owner(neigh), &implT::sock_it_to_me, neigh, result.remote_ref(world), TaskAttributes::hipri());
                return result;
            }
        }
*/

#endif /* __fuset_InnerOp_h__ */
