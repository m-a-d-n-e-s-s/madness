
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
		typedef Function<T,NDIM> KTREE;
		typedef FunctionNode<T,NDIM> KNODE;
		typedef Key<NDIM> keyT;
		typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;


	///< Type of container holding the nodes

		typedef GenTensor<T>										coeffT;	//Type of tensor used to hold coeffs
		typedef Tensor<T>											tensorT;
	    
    public:
		ReconstructOp(string opName, KTREE* output, const KTREE* i1);
    
		BaseParameters<T> compute(const keyT& key, const BaseParameters<T> &s);

		bool isDone(const keyT& key) const;
		bool isPre() const { return true; }

		//
		Future <BaseParameters<T>> computeC(const keyT& key, const BaseParameters<T> &s) { }
		Future <BaseParameters<T>> computeC	(const keyT& key, const std::vector<Future<BaseParameters<T>>> &v) { }

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



    template<typename T, std::size_t NDIM>
	ReconstructOp<T,NDIM>::ReconstructOp(string opName, KTREE* output, const KTREE* i1)
		: PrimitiveOp<T,NDIM>(opName, output, false)
		, _i1(i1)
		, _cdata(FunctionCommonData<T,NDIM>::get(i1->k()))
		, _targs(i1->get_impl()->get_tensor_args())
	{
	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,false,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,false,false,false));
      
	}
    
    template<typename T, std::size_t NDIM>
	BaseParameters<T>
	ReconstructOp<T,NDIM>::compute(const keyT& key, const BaseParameters<T> &s)
	{
		//printf ("ReconstructOp: %s\n", __func__);
		BaseParameters<T> returnP;

		// for Root
		if (key == _i1->world().rank())
		{
			if (_i1->get_impl()->compressed == false)
				exit(0);

			_i1->get_impl()->compressed		= false;
			_i1->get_impl()->nonstandard	= false;
			_i1->get_impl()->redundant		= false;
		}	

		typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();

        if (it == _i1->get_impl()->get_coeffs().end()) 
		{
			cout << __func__ << " key: " << key << "1st if" << endl;
			cerr<<"This should not have happenned"<<endl;
			exit(0);
        }

		// The integral operator will correctly connect interior nodes
		// to children but may leave interior nodes without coefficients
		// ... but they still need to sum down so just give them zeros
		KNODE& node = it->second;
		if (node.has_children() && !node.has_coeff())
		{
			//cout << __func__ << " key: " << key << "2nd if" << endl;
			node.set_coeff(coeffT(_cdata.v2k, _targs));
		}

		if (node.has_children() || node.has_coeff())
		{
			//cout << __func__ << " key: " << key << " [3rd if]" << endl;
			coeffT d = node.coeff();
		
			if (!d.has_data())
				d = coeffT(_cdata.v2k, _targs);

			if (key.level() > 0) 
				d(_cdata.s0) += s._coeff; // this is the problem!!!

			if (d.dim(0) == 2*(_i1->get_impl()->get_k()))
			{
				d = _i1->get_impl()->unfilter(d);
				node.clear_coeff();
				node.set_has_children(true);

				for (KeyChildIterator<NDIM> kit(key); kit; ++kit)
				{
					const keyT& child = kit.key();
					coeffT ss = copy(d(child_patch(child)));
					ss.reduce_rank(_i1->thresh());
			
					//cout << __func__ << " key: " << key << "ss: " << ss << endl;
					returnP._isVector = true;	
					returnP._vCoeff.push_back(ss);	
				}
			}
			else
			{
				node.coeff().reduce_rank(_targs.thresh);
			}
		
		}
		else
		{
			//cout << __func__ << " key: " << key << " [3rd else]" << endl;
			coeffT ss = s._coeff;
			if (s._coeff.has_no_data())
				ss = coeffT(_cdata.vk, _targs);

			if (key.level()) 
				node.set_coeff(copy(ss));
			else
				node.set_coeff(ss);

			//cout << "ss: " << ss << endl;
		}
	
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
