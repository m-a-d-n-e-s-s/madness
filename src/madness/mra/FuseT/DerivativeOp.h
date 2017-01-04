
// Ghaly
//
// Derivatives the function, transforming into scaling function basis. Possible non-blocking comm.
//
#ifndef __MADNESS_MRA_FUSET_DERIVATIVE_OP__INCLUDED__
#define __MADNESS_MRA_FUSET_DERIVATIVE_OP__INCLUDED___OP__INCLUDED__

#include "PrimitiveOp.h"
#include "../function_common_data.h"
#define LEFT 0
#define CENTER 1
#define RIGHT 2
namespace madness 
{
    template<typename T, std::size_t NDIM>
	class DerivativeOp : public PrimitiveOp<T,NDIM>, public WorldObject<DerivativeOp<T,NDIM>> 
    {
    public:
	typedef Function<T,NDIM> KTREE;
	typedef FunctionNode<T,NDIM> KNODE;
	typedef GenTensor<T> coeffT;
	typedef Key<NDIM> keyT;
	typedef std::pair<keyT,coeffT>  argT;
	typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>>	dcT;
	typedef WorldContainer<keyT, Future<argT> > daT;
	typedef Tensor<T> tensorT;
	typedef DerivativeOp<T,NDIM> dOp;
	typedef WorldObject<dOp> woT;
	int pending_task, pending_from_left, pending_from_right;
    public:
	DerivativeOp(string opName, KTREE* output, const KTREE* i1, World& world, Derivative<T,NDIM>* dOp);
	Future<FuseTContainer<T> >  computeFuture(const keyT& key, const FuseTContainer<T> &s);
	bool isDone(const keyT& key) const;
	bool isPre() const { return true; }
	bool needsParameter() const { return true; }
	void reduce(World& world){}
	bool returnsFuture(){return true;}

    private:

	void createFuture(const keyT& key);
	void pushFromLeft(const keyT& key, const coeffT& coeff);
	void pushFromRight(const keyT& key, const coeffT& coeff);
	//FuseTContainer<T>  computeDerivative(const keyT& key, Future<argT>& left, argT& center, Future<argT>& right); 
	FuseTContainer<T>  computeDerivative(const keyT& key, argT& left, argT& center, argT& right); 

	//!Points to operand trees
	const KTREE* _i1;
	const Derivative<T,NDIM>* _D;
 
	//!Points to operand nodes of the tree
	KNODE *_t1, *_t2;

	daT _fromLeft, _fromRight;
		
	World& _world;
	//!Variables for DerivativeOp
	const FunctionCommonData<T,NDIM>& _cdata;

    };

    // Constructor
    template<typename T, std::size_t NDIM>
	DerivativeOp<T,NDIM>::DerivativeOp(string opName, KTREE* output, const KTREE* i1, World& world, Derivative<T,NDIM>* dOp)
	: PrimitiveOp<T,NDIM>(opName, output, false)
	, _i1(i1)
	, _D(dOp)
	, _cdata(FunctionCommonData<T,NDIM>::get(i1->k()))
	, woT(world)
        , _world(world)
	
    {
	pending_from_left = 0;
	pending_from_right = 0;
	pending_task = 0;
	woT::process_pending();
	const std::shared_ptr<WorldDCPmapInterface<Key<NDIM> > > _pmap(FunctionDefaults<NDIM>::get_pmap());		    
	_fromLeft = daT(world,_pmap,false);
	_fromRight = daT(world,_pmap,false);
	_fromLeft.process_pending();
	_fromRight.process_pending();

	//dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
	this->_OpID = output->get_impl()->id().get_obj_id();
	this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(i1,true,true,false,false));
	this->_dInfoVec.push_back(DependencyInfo<T,NDIM>(output,true,true,false,false));
    }


    
    template<typename T, std::size_t NDIM>    
	Future<FuseTContainer<T> >
	DerivativeOp<T,NDIM>::computeFuture(const keyT& key, const FuseTContainer<T> &s) 
    {
	//cout<<key<<" Pending task "<<pending_task<<" Pending from Left "<<pending_from_left<<" Pending from Right "<<pending_from_right<<endl;
	argT left, right, center;

	/********Unpack the paramters left, center, right************/

	//create left, center and right
	if(s.get() == 0){
	    left  = std::make_pair(key,coeffT());
	    right = std::make_pair(key,coeffT());
	    center = std::make_pair(key,coeffT());
	}else{
	    if(s.what() == WHAT_AM_I::FuseT_VArgT){
		left = (((FuseT_VArgT<T>*)s.get())->value)[LEFT];
		center = (((FuseT_VArgT<T>*)s.get())->value)[CENTER];
		right = (((FuseT_VArgT<T>*)s.get())->value)[RIGHT];
	    }else{
		cerr<<"This should not have happenned"<<endl;
	    }

	}

	/***********************************************************
	    if center does not exist, fetch center. Then push it to left
            and right keys. 
	************************************************************/
	keyT leftKey, rightKey;

	if(!center.second.has_data()){
	    typename dcT::iterator it= _i1->get_impl()->get_coeffs().find(key).get();
	    KNODE&	node = it->second;	// source 
	    if(!node.has_children()){
		//cout<<key<<"has Children: False"<<endl;
		center = std::make_pair(key,node.coeff());
	    }else{
		//cout<<key<<"has Children: True"<<endl;
	    }
	    leftKey = _D->neighbor(key,-1);
	    rightKey = _D->neighbor(key,+1);
		
	    if(!leftKey.is_invalid()){
		pending_from_right++;
		woT::task(_fromRight.owner(leftKey),&DerivativeOp<T,NDIM>::pushFromRight, leftKey, node.coeff());

	    }	    else
		left.first = leftKey;

	    if(!rightKey.is_invalid()){
		pending_from_left++;
		woT::task(_fromLeft.owner(rightKey),&DerivativeOp<T,NDIM>::pushFromLeft, rightKey, node.coeff());

	    }	    else
		right.first = rightKey;

	}


	//creates futures for this key from left and right
	//incase it is not created by fromLeft or fromRight yet
	createFuture(key);

	/***************************************************************
             Set up the futures that the actual computation depends upon
	***************************************************************/

	Future<argT> leftPara;
	Future<argT> rightPara;
	    
	if(left.second.has_data() || left.first.is_invalid())
	    leftPara = Future<argT>(left);
	else{
	    typename daT::accessor acc;
	    if(_fromLeft.find(acc,key)) 
		leftPara = acc->second;  
	    else   
		cerr<<"This should not have happenned"<<endl;
	}


	if(right.second.has_data() || right.first.is_invalid())
	    rightPara = Future<argT>(right);
	else{
	    typename daT::accessor acc;
	    if(_fromRight.find(acc,key)) 
		rightPara = acc->second;  
	    else   
		cerr<<"This should not have happenned"<<endl;
	}

	/********************************************************************
	    Finally schedule task for actual computation
	********************************************************************/
	Future<FuseTContainer<T> > returnParameter;
	pending_task++;
	returnParameter = woT::task(_world.rank(),&DerivativeOp<T,NDIM>::computeDerivative, key, leftPara, center, rightPara);
	return returnParameter;
    }

    template<typename T, std::size_t NDIM>    
	FuseTContainer<T> 
	//DerivativeOp<T,NDIM>::computeDerivative(const keyT& key, Future<argT>& left, argT& center, Future<argT>& right) 
	DerivativeOp<T,NDIM>::computeDerivative(const keyT& key, argT& left, argT& center, argT& right) 
    {
	
	//std:://cout<<key<<endl;
	bool c = center.second.has_data();
	bool l = left.second.has_data();
	bool r = right.second.has_data();
	bool lInv = left.first.is_invalid();
	bool rInv = right.first.is_invalid();
/*	bool l = left.get().second.has_data();
	bool r = right.get().second.has_data();
	bool lInv = left.get().first.is_invalid();
	bool rInv = right.get().first.is_invalid();*/
	bool isInterior = false;

	//all Coefficients are available, we can compute and interior leaf node
	if(l && c && r ){
	    //cout<<key<<" Leaf : True Computing Leaf "<<endl;
	    //_D->do_diff2i(_i1->get_impl().get(), this->_result->get_impl().get(), key, left.get(), center, right.get());
	    _D->do_diff2i(_i1->get_impl().get(), this->_result->get_impl().get(), key, left, center, right);

	    //All Coefficients are available, but the key is a boundary, so use boundary calculations
	}else if (c && ( (rInv && lInv) || (l&&rInv) || (r&&lInv) ) ){
	    //_D->do_diff2b(_i1->get_impl().get(), this->_result->get_impl().get(), key, left.get(), center, right.get());
	    //cout<<key<<" Border : True.  Center Has Data :"<<c<<endl;
	    _D->do_diff2b(_i1->get_impl().get(), this->_result->get_impl().get(), key, left, center, right);
		
	    //Some coefficients are not available, so this is an empy interior node
	}else{
	    //cout<<key<<" Interior : True "<<endl;
	    isInterior = true;
	    this->_result->get_impl()->get_coeffs().replace(key,KNODE(coeffT(),true)); // Empty internal node
	}
	pending_task--;
	if(!isInterior){
	    return FuseTContainer<T>();
	}

	//Construct Parameters if recursive call is required
	FuseT_VParameter<T>*  v_parameter = new FuseT_VParameter<T>();
	for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
	    const keyT& child = kit.key();
	    FuseT_VArgT<T>* tvec = new FuseT_VArgT<T>();
	    if ((child.translation()[_D->get_axis()]&1) == 0) {
		//tvec->value.push_back(left.get());
		tvec->value.push_back(left);
		tvec->value.push_back(center);
		tvec->value.push_back(center);
	    }else{
		tvec->value.push_back(center);
		tvec->value.push_back(center);
		//tvec->value.push_back(right.get());		
		tvec->value.push_back(right);		
	    }		
	    FuseTContainer<T> wrapper(static_cast<Base<T>*>(tvec));
	    v_parameter->value.push_back(wrapper);
	}

	FuseTContainer<T> temp(static_cast<Base<T>*>(v_parameter));
	return temp;
    }
	
    template<typename T, std::size_t NDIM>    
	void
	DerivativeOp<T,NDIM>::createFuture(const keyT& key){
	typename daT::accessor acc;
	if(!_fromLeft.find(acc,key)) 
	    _fromLeft.insert(acc,key);
	
	if(!_fromRight.find(acc,key)) 
	    _fromRight.insert(acc,key);
    }
	    

    template<typename T, std::size_t NDIM>    
	void
	DerivativeOp<T,NDIM>::pushFromLeft(const keyT& key, const coeffT& coeff){
	if(_fromLeft.owner(key) == _world.rank()){
	    typename daT::accessor acc;
	    if(!_fromLeft.find(acc,key)) 
		_fromLeft.insert(acc,key);
	    acc->second = Future<argT>(std::make_pair(_D->neighbor(key,-1),coeff));
	}else{
	    cerr<<"Should not have happened. "<<endl;
	}
	pending_from_left--;
    } 

    template<typename T, std::size_t NDIM>    
	void
	DerivativeOp<T,NDIM>::pushFromRight(const keyT& key, const coeffT& coeff){
	if(_fromLeft.owner(key) == _world.rank()){
	    typename daT::accessor acc;
	    if(!_fromRight.find(acc,key)) 
		_fromRight.insert(acc,key);
	    acc->second = Future<argT>(std::make_pair(_D->neighbor(key,+1),coeff));
	}else{
	    cerr<<"Should not have happened. "<<endl;
	}
	pending_from_right--;
    } 



    template<typename T, std::size_t NDIM>
	bool DerivativeOp<T,NDIM>::isDone(const keyT& key) const 
    {

	bool isLeaf = !this->_result->get_impl()->get_coeffs().find(key).get()->second.has_children();
	//cout<<key<<" Is Done : "<<isLeaf<<endl;
	return isLeaf;
    }

}; /*fuset*/

#endif /* __fuset_DerivativeOp_h__ */
