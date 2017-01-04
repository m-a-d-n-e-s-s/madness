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
#ifndef __MADNESS_MRA_FUSET_DIFF_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_DIFF_OP__INCLUDED__

#include "PrimitiveOp.h"
#include "../mra.h"
#include "../function_common_data.h"
#include "../../world/MADworld.h"
#include "../../tensor/tensor.h"
#include "../../tensor/gentensor.h"

namespace madness 
{
	template<typename T, std::size_t NDIM>
	class DiffOp : public PrimitiveOp<T,NDIM> 
	{
		typedef Function<T,NDIM>									KTREE;
		typedef FunctionNode<T,NDIM>								KNODE;	// identical to nodeT
		typedef Key<NDIM>											keyT;
		typedef WorldContainer<Key<NDIM>, FunctionNode<T, NDIM>>	dcT;
		typedef WorldObject<FunctionImpl<T,NDIM>>					woT;	///< Base class world object type

		typedef GenTensor<T>										coeffT;	//Type of tensor used to hold coeffs
		typedef Tensor<T>											tensorT;
		typedef std::pair<keyT,coeffT>								argT;


	///< Type of container holding the nodes
	    
    public:
									DiffOp		(string opName, KTREE* output, const KTREE* i1);
		void						compute		(const keyT& key) { }
		FuseTContainer<T>			compute		(const keyT& key, const FuseTContainer<T> &s) { }

		//
		//Future<FuseTContainer<T>> postCompute(const keyT& key, const std::vector<Future<FuseTContainer<T>>> &s) { }
		bool						isDone			(const keyT& key) const;
		bool						isPre			() const { return true; }
		bool						needsParameter	() const { return false; }

	public:	// for DiffOp
	
		Key<NDIM>						neighbor		(const keyT& key, int step) const;
		argT							find_neighbor	(const Key<NDIM>& key, int step) const;

    private:
		//!Points to operand trees
		const KTREE*						_i1;
    
		//!Points to operand nodes of the tree
		KNODE								*_t1, *_t2;

		//!Variables for CompressOp
		dcT&								_coeffs;
		const FunctionCommonData<T,NDIM>&	_cdata; 
		TensorArgs							_targs;	
		bool								_temp;
		int									_num;

		int									_k;				// Wavelet order	
    };

		
	// Constructor
	// World is needed for communication in the "compute" function
    template<typename T, std::size_t NDIM>
	DiffOp<T,NDIM>::DiffOp(string opName, KTREE* output, const KTREE* i1)
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
	
		this->_num = 0;
		this->_temp = false;
		woT(i1->world());
	}

/*
    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::do_diff1(const DerivativeBase<T,NDIM>* D,
                                        const implT* f,
                                        const keyT& key,
                                        const std::pair<keyT,coeffT>& left,
                                        const std::pair<keyT,coeffT>& center,
                                        const std::pair<keyT,coeffT>& right) {
        D->do_diff1(f,this,key,left,center,right);
    }
*/
/*
        void do_diff1(const implT* f, implT* df, const keyT& key,
                      const argT& left,
                      const argT& center,
                      const argT& right) const {
            MADNESS_ASSERT(axis<NDIM);

            if ((!left.second.has_data()) || (!right.second.has_data())) {
                // One of the neighbors is below us in the tree ... recur down
                df->get_coeffs().replace(key,nodeT(coeffT(),true));
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    if ((child.translation()[axis]&1) == 0) {
                        // leftmost child automatically has right sibling
                        forward_do_diff1(f, df, child, left, center, center);
                    }
                    else {
                        // rightmost child automatically has left sibling
                        forward_do_diff1(f, df, child, center, center, right);
                    }
                }
            }
            else {
                forward_do_diff1(f, df, key, left, center, right);
            }
        }
*/
/*
        void forward_do_diff1(const implT* f, implT* df, const keyT& key,
                              const argT& left,
                              const argT& center,
                              const argT& right)  const {

            const dcT& coeffs = f->get_coeffs();
            ProcessID owner = coeffs.owner(key);

            if (owner == world.rank()) {
                if (!left.second.has_data()) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff1,
                            f, df, key, find_neighbor(f, key,-1), center, right,
                            TaskAttributes::hipri());
                }
                else if (!right.second.has_data()) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff1,
                            f, df, key, left, center, find_neighbor(f, key,1),
                            TaskAttributes::hipri());
                }
                // Boundary node
                else if (left.first.is_invalid() || right.first.is_invalid()) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff2b,
                            f, df, key, left, center, right);
                }
                // Interior node
                else {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff2i,
                            f, df, key, left, center, right);
                }
            }
            else {
                df->task(owner, &madness::FunctionImpl<T,NDIM>::forward_do_diff1,
                        this, f, key, left, center, right, TaskAttributes::hipri());
            }
        }
*/
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

	// neighbor from Derivative
	template<typename T, std::size_t NDIM>
	Key<NDIM>
	DiffOp<T,NDIM>::neighbor(const keyT& key, int step) const 
	{
		Vector<Translation,NDIM> l = key.translation();
		// l[axis] += step;
		// enforce_bc(??)


		return keyT::invalid();
	}

	// find_neighbor from Derivative
	template<typename T, std::size_t NDIM>
	std::pair<Key<NDIM>,GenTensor<T>>	
	DiffOp<T,NDIM>::find_neighbor(const Key<NDIM>& key, int step) const
	{
		keyT neigh = neighbor(key, step);
		
		if (neigh.is_invalid())
		{
			argT temp;
			return temp;
			//
		}
		else
		{
			argT result;
			// f->
			// sock_it_to_me!!!	
			return result;
		}
	}


	// isDone
    template<typename T, std::size_t NDIM>
	bool 
	DiffOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
		if(!isE1) return isE1;
		bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf;
    }

}; /*fuset*/

#endif /* __fuset_DiffOp_h__ */
