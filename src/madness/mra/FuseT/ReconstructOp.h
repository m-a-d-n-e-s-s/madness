
// Ghaly
//
// Reconstructs the function, transforming into scaling function basis. Possible non-blocking comm.
//
#ifndef __MADNESS_MRA_FUSET_RECONSTRUCT_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_RECONSTRUCT_OP__INCLUDED__

#include "PrimitiveOp.h"
#include "../function_common_data.h"

namespace madness 
{
    template<typename T, std::size_t NDIM>
	class ReconstructOp : public PrimitiveOp<T,NDIM>
	{
	public:
		typedef Function<T,NDIM>								KTREE;
		typedef FunctionNode<T,NDIM>							KNODE;
		typedef Key<NDIM>										keyT;
		typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>>	dcT;

		typedef GenTensor<T>									coeffT;
		typedef Tensor<T>										tensorT;
	    
    public:
							ReconstructOp	(string opName, KTREE* output, const KTREE* i1);
		FuseTContainer<T>	compute			(const keyT& key, const FuseTContainer<T> &s);

		bool				isDone			(const keyT& key) const;
		bool				isPre			() const { return true; }
		bool				needsParameter	() const { return true; }
        void reduce(World& world){}
	public:
		std::vector<Slice>			child_patch		(const keyT& child) const;

    private:
		//!Points to operand trees
		const KTREE* _i1;
		
		//!Points to operand nodes of the tree
		KNODE *_t1, *_t2;

		//!Variables for ReconstructOp
		const FunctionCommonData<T,NDIM>&	_cdata;
		TensorArgs							_targs;
    };

	// Constructor
    template<typename T, std::size_t NDIM>
	ReconstructOp<T,NDIM>::ReconstructOp(string opName, KTREE* output, const KTREE* i1)
		: PrimitiveOp<T,NDIM>(opName, output, false)
		, _i1(i1)
		, _cdata(FunctionCommonData<T,NDIM>::get(i1->k()))
		, _targs(i1->get_impl()->get_tensor_args())
	{
	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,true,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,true,false,false));
	}

	template<typename T, std::size_t NDIM>    
	FuseTContainer<T> 
	ReconstructOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s) 
	{
		FuseT_CoeffT<T>* s_coeff;
		if (s.get() == 0) 
		{
			// it should be called initially
            s_coeff = new FuseT_CoeffT<T>();
			s_coeff->value = coeffT();
		} 
		else 
		{
			s_coeff = new FuseT_CoeffT<T>();
			if (s.what() == WHAT_AM_I::FuseT_CoeffT) 
			{
				s_coeff->value = (((FuseT_CoeffT<T>*)s.get())->value).full_tensor_copy();
				delete s.data;
			}
			else
			{
				cerr<<"This should not have happenned"<<endl;
			}
		}

		// for Root
		if (key == keyT(0))
		{
			this->_result->get_impl()->compressed	= false;
			this->_result->get_impl()->nonstandard	= false;
			this->_result->get_impl()->redundant	= false;
		}	

		typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();
        if (it == _i1->get_impl()->get_coeffs().end()) 
		{
			_i1->get_impl()->get_coeffs().replace(key,KNODE(coeffT(), false));
			it = _i1->get_impl()->get_coeffs().find(key).get();
        }

		// The integral operator will correctly connect interior nodes
		// to children but may leave interior nodes without coefficients
		// ... but they still need to sum down so just give them zeros
		KNODE&		node = it->second;	// source 
		KNODE		targetNode;				// target (if target is initially copied from source, it should be perfect.)
		FuseTContainer<T>	returnP;
		if (node.has_children() && !node.has_coeff())
		{
			targetNode.set_coeff(coeffT(_cdata.v2k, _targs));
			this->_result->get_impl()->get_coeffs().replace(key,targetNode);
		}

		if (node.has_children() || node.has_coeff())
		{
			coeffT d = copy(node.coeff());
		
			if (!d.has_data()) {
				d = coeffT(_cdata.v2k, _targs);
			}
			if (key.level() > 0) 
			{
				d(_cdata.s0) += s_coeff->value; // this is the problem!!!
			}
			delete s_coeff;		

			if (d.dim(0) == 2*(_i1->get_impl()->get_k()))
			{
				FuseT_VParameter<T>		v_parameter;

				d = _i1->get_impl()->unfilter(d);
				targetNode.clear_coeff();
				targetNode.set_has_children(true);
				this->_result->get_impl()->get_coeffs().replace(key, targetNode);
	
				for (KeyChildIterator<NDIM> kit(key); kit; ++kit) 
				{
					FuseT_CoeffT<T>		s_coeff_s;

					const keyT& child = kit.key();
					coeffT ss = copy(d(child_patch(child)));
					ss.reduce_rank(_i1->thresh());
					s_coeff_s.value = ss;

					FuseTContainer<T> wrapper(static_cast<Base<T>*>(new FuseT_CoeffT<T>(s_coeff_s.value)));
					v_parameter.value.push_back(wrapper);
				}

				//	Wrapping ReconstructV<T> to FuseTContainer<T>
				FuseTContainer<T> temp(static_cast<Base<T>*>(new FuseT_VParameter<T>(v_parameter.value)));
				return temp;
			}
			else
			{
				MADNESS_ASSERT(node.is_leaf());	//???
				targetNode.coeff().reduce_rank(_targs.thresh);
				this->_result->get_impl()->get_coeffs().replace(key, targetNode);
			}
		}
		else
		{
			coeffT ss = s_coeff->value;
			if (s_coeff->value.has_no_data()) 
			{	
				ss = coeffT(_cdata.vk, _targs);
			}
			
			if (key.level())
			{
				targetNode.set_coeff(copy(ss));
				this->_result->get_impl()->get_coeffs().replace(key, targetNode);
			}
			else
			{
				// my example--ReconstructEx does visit this else statement.
				targetNode.set_coeff(ss);
				this->_result->get_impl()->get_coeffs().replace(key, targetNode);
			}
		}
		//
		//	 nothing... 
		delete s_coeff;
		return returnP;
	}

	template <typename T, std::size_t NDIM>
	std::vector<Slice>
	ReconstructOp<T,NDIM>::child_patch(const keyT& child) const 
	{
		std::vector<Slice> s(NDIM);
		const Vector<Translation,NDIM>& l = child.translation();
		for (std::size_t i = 0; i<NDIM; ++i)
		{
			s[i] = _cdata.s[l[i]&1];	// Lowest bit of translation
		}
		return s;
	}

    template<typename T, std::size_t NDIM>
	bool ReconstructOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = _i1->get_impl()->get_coeffs().probe(key);
		if(!isE1) return isE1;
		bool isLeaf = !_i1->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf;
    }

}; /*fuset*/

#endif /* __fuset_ReconstructOp_h__ */
