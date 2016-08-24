
// Ghaly
//
// Reconstructs the function, transforming into scaling function basis. Possible non-blocking comm.
//
#ifndef __MADNESS_MRA_FUSET_MULTIPLY_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_MULTIPLY_OP__INCLUDED__

#include "PrimitiveOp.h"
#include "../function_common_data.h"
#define DEBUG 0

namespace madness 
{
    
    template<typename T, std::size_t NDIM>
	class MultiplyOp : public PrimitiveOp<T,NDIM>
	{
	public:
		typedef Function<T,NDIM> KTREE;
		typedef FunctionNode<T,NDIM> KNODE;
		typedef Key<NDIM> keyT;
		typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>>	dcT;

		typedef GenTensor<T> coeffT;
		typedef Tensor<T> tensorT;
	    
    public:
		MultiplyOp(string opName, KTREE* output, const KTREE* i1, const KTREE* i2, double tol);
		FuseTContainer<T> compute(const keyT& key, const FuseTContainer<T> &s);

		bool isDone(const keyT& key) const;
		bool isPre() const { return true; }
		bool needsParameter() const { return true; }
        void reduce(World& world){}
	public:
		std::vector<Slice> child_patch(const keyT& child) const;
		void do_mul(const keyT& key, const Tensor<T>& left, const std::pair< keyT, Tensor<T> >& arg);
    
		double truncate_tol(double tol, const keyT& key) const;
		
    private:
		//!Points to operand trees
		const KTREE* _i1;
		const KTREE* _i2;
		double _tol;
		//!Variables for MultiplyOp
		const FunctionCommonData<T,NDIM>&  _cdata;
		TensorArgs _targs;
    };

	// Constructor
    template<typename T, std::size_t NDIM>
	MultiplyOp<T,NDIM>::MultiplyOp(string opName, KTREE* output, const KTREE* i1, const KTREE* i2, double tol)
		: PrimitiveOp<T,NDIM>(opName, output, false)
	        , _i1(i1)
	        , _i2(i2)
		, _cdata(FunctionCommonData<T,NDIM>::get(i1->k()))
        	, _targs(i1->get_impl()->get_tensor_args())
	        , _tol(tol)
        {
	    //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	    this->_OpID = output->get_impl()->id().get_obj_id();
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,true,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i2,true,true,false,false));
	    this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,true,false,false));
	}

	template<typename T, std::size_t NDIM>    
	FuseTContainer<T> 
	MultiplyOp<T,NDIM>::compute(const keyT& key, const FuseTContainer<T> &s) 
	{
	    typedef typename FunctionImpl<T,NDIM>::dcT::const_iterator literT;
	    typedef typename FunctionImpl<T,NDIM>::dcT::const_iterator riterT;
	    double lnorm=1e99, rnorm=1e99;
	    std::vector<coeffT> parameter;

	    if(s.get() == 0){
		parameter = std::vector<coeffT>(2);
		
	    }else if (s.what() == WHAT_AM_I::FuseT_VCoeffT) {
		parameter = ((FuseT_VCoeffT<T>*)s.get())->value;
		
		delete (s.data);
		
		//printf ("target: %p vs. parameter: %p\n", &parameter, &(((FuseT_VCoeffT<T>*)s.get())->value));

	    }else{
		cerr<<"This should not have happenned"<<endl;
	    }

	    //parmeter 0 is coeff from i1 and parameter 1 is coeff from i2
	    Tensor<T> lc = parameter[0];
	    if (lc.size() == 0) {
		literT it = _i1->get_impl()->get_coeffs().find(key).get();
		MADNESS_ASSERT(it != _i1->get_impl()->get_coeffs().end());
		lnorm = it->second.get_norm_tree();
		if (it->second.has_coeff())
		    lc = it->second.coeff().full_tensor_copy();
	    }

	    Tensor<T> rc = parameter[1];
	    if (rc.size() == 0) {
		riterT it = _i2->get_impl()->get_coeffs().find(key).get();
		MADNESS_ASSERT(it != _i2->get_impl()->get_coeffs().end());
		rnorm = it->second.get_norm_tree();
		if (it->second.has_coeff())
		    rc = it->second.coeff().full_tensor_copy();
	    }

	    // both nodes are leaf nodes: multiply and return
	    if (rc.size() && lc.size()) { // Yipee!
		do_mul(key, lc, std::make_pair(key,rc));
		return FuseTContainer<T>();
	    }

	    if (_tol) {
		if (lc.size())
		    lnorm = lc.normf(); // Otherwise got from norm tree above
		if (rc.size())
		    rnorm = rc.normf();
		if (lnorm*rnorm < truncate_tol(_tol, key)) {
		    this->_result->get_impl()->get_coeffs().replace(key, KNODE(coeffT(_cdata.vk,_targs),false)); // Zero leaf node
		    return FuseTContainer<T>();
		}
	    }

	    // Recur down
	    this->_result->get_impl()->get_coeffs().replace(key, KNODE(coeffT(),true)); // Interior node

	    Tensor<T> lss;
	    if (lc.size()) {
		Tensor<T> ld(_cdata.v2k);
		ld(_cdata.s0) = lc(___);
		lss = _i1->get_impl()->unfilter(ld);
	    }

	    Tensor<T> rss;
	    if (rc.size()) {
		Tensor<T> rd(_cdata.v2k);
		rd(_cdata.s0) = rc(___);
		rss = _i2->get_impl()->unfilter(rd);
	    }

			

	    FuseT_VParameter<T>*  v_parameter = new FuseT_VParameter<T>();

	    for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
		const keyT& child = kit.key();
		FuseT_VCoeffT<T>* tvec = new FuseT_VCoeffT<T>(2);
		if (lc.size())
		    tvec->value[0] = copy(lss(child_patch(child)));
		if (rc.size())
		    tvec->value[1] = copy(rss(child_patch(child)));
		
		FuseTContainer<T> wrapper(static_cast<Base<T>*>(tvec));
		v_parameter->value.push_back(wrapper);

	    }
	    FuseTContainer<T> temp(static_cast<Base<T>*>(v_parameter));
	    return temp;


        }

	
	template <typename T, std::size_t NDIM>
	std::vector<Slice>
	MultiplyOp<T,NDIM>::child_patch(const keyT& child) const 
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
	bool MultiplyOp<T,NDIM>::isDone(const keyT& key) const 
	{
		bool isE1 = this->_result->get_impl()->get_coeffs().probe(key);
		if(!isE1) return isE1;
		bool isLeaf = !this->_result->get_impl()->get_coeffs().find(key).get()->second.has_children();
		return isLeaf;
    }


    template<typename T, std::size_t NDIM>
	void MultiplyOp<T,NDIM>::do_mul(const keyT& key, const Tensor<T>& left, const std::pair< keyT, Tensor<T> >& arg) {
	PROFILE_MEMBER_FUNC(FunctionImpl);
	const keyT& rkey = arg.first;
	const Tensor<T>& rcoeff = arg.second;

	Tensor<T> rcube = this->_result->get_impl()->fcube_for_mul(key, rkey, rcoeff);
	//madness::print("do_mul: l", key, left.size());
	Tensor<T> lcube = this->_result->get_impl()->fcube_for_mul(key, key, left);


	Tensor<T> tcube(_cdata.vk,false);
	TERNARY_OPTIMIZED_ITERATOR(T, tcube, T, lcube, T, rcube, *_p0 = *_p1 * *_p2;);
	double scale = pow(0.5,0.5*NDIM*key.level())*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
	tcube = transform(tcube,_cdata.quad_phiw).scale(scale);
	if(DEBUG && key == keyT(0)){
	    std::cout<<"Computed Coeff at "<<key<<" Norm is :"<<tcube.normf()<<std::endl;
	}

	this->_result->get_impl()->get_coeffs().replace(key, KNODE(coeffT(tcube),false));
	return;
    }
    
    template<typename T, std::size_t NDIM>
	double MultiplyOp<T,NDIM>::truncate_tol(double tol, const keyT& key) const {
	const static double fac=1.0/std::pow(2,NDIM*0.5);
	tol*=fac;

	// RJH ... introduced max level here to avoid runaway
	// refinement due to truncation threshold going down to
	// intrinsic numerical error
	const int MAXLEVEL1 = 20; // 0.5**20 ~= 1e-6
	const int MAXLEVEL2 = 10; // 0.25**10 ~= 1e-6
	int truncate_mode = 1;
	if (truncate_mode == 0) {
	    return tol;
	}
	else if (truncate_mode == 1) {
	    double L = FunctionDefaults<NDIM>::get_cell_min_width();
	    return tol*std::min(1.0,pow(0.5,double(std::min(key.level(),MAXLEVEL1)))*L);
	}
	else if (truncate_mode == 2) {
	    double L = FunctionDefaults<NDIM>::get_cell_min_width();
	    return tol*std::min(1.0,pow(0.25,double(std::min(key.level(),MAXLEVEL2)))*L*L);
	}
	else {
	    MADNESS_EXCEPTION("truncate_mode invalid",truncate_mode);
	}
    }


    
}; /*fuset*/

#endif /* __fuset_MultiplyOp_h__ */
