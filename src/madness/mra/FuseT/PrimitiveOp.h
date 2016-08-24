#ifndef __MADNESS_MRA_FUSET_PRIMITIVE_OP_H__INCLUDED__
#define __MADNESS_MRA_FUSET_PRIMITIVE_OP_H__INCLUDED__

#include "madness/mra/mra.h"
#include "madness/world/archive.h"
#include <string>
#include "FuseTContainer.h"

namespace madness 
{
    using namespace std;
/*	!spatial relation, Same Node Relation (PSI), Some Sibling Relation (SIGMA), Ancestor Relation (ALPHA), 
 *	Descendent Relation (DELTA). Dependency Info contains information about type of dependency to parameters*/
    template<typename T, std::size_t NDIM>
    struct DependencyInfo
	{
		const Function<T,NDIM>* _producerTree;
		string _treeName;
		unsigned long _treeID;

		bool _psi;
		bool _alpha;
		bool _delta;
		bool _sigma;

/*		!this is used by FuseT. _producerIndex stores the index of operator that produces this tree/ 
 *		_consumerIndex stores the index of the operator that consumes this tree*/
		int _producerIndex;
		int _consumerIndex;

		DependencyInfo(const Function<T,NDIM>* producerTree, bool psi, bool alpha, bool delta, bool sigma)
		{
			_psi = psi;
			_alpha = alpha;
			_delta = delta;
			_sigma = sigma;
			
			_producerTree = producerTree;
			_treeID	= _producerTree->get_impl()->id().get_obj_id();
			_treeName = _producerTree->_treeName;
		}
    };

/*The _OpID is used by the fusion compiler to analyze dependencies relating to this operation. 
  It "MUST" be set equal to the _treeID of the result tree. 
  In MADNESS, the tree ID of a Function is given by func->get_impl()->id().get_obj_id()*/
    template<typename T, std::size_t NDIM> class PrimitiveOp { 
	typedef Function<T,NDIM> KTREE; 
	typedef FunctionNode<T,NDIM> KNODE; 
	typedef Key<NDIM> keyT; 
	typedef WorldContainer<Key<NDIM>, FunctionNode<T,NDIM> > dcT; 
///< Type of container holding the nodes

    public:
    PrimitiveOp(string opName="Unknown",KTREE* result=NULL, bool isComplete=false, bool isReduce=false)  
	:_opName(opName),
	 _result(result),
	 _isComplete(isComplete),
     _isReduce(isReduce)
    {}

	virtual ~PrimitiveOp() {}
   
	virtual FuseTContainer<T> compute (const keyT& key, const FuseTContainer<T> &s)	{ return FuseTContainer<T>();}
	virtual Future< FuseTContainer<T> > computeFuture (const keyT& key, const FuseTContainer<T> &s) { return Future<FuseTContainer<T> >();} 

	virtual bool isDone(const keyT& key) const = 0;
	virtual bool isPre() const = 0;
	virtual bool needsParameter() const = 0;
	virtual void reduce(World& world) = 0;
	virtual bool returnsFuture(){ return false;}

	//!used for postCompute ops to see if it needs to be pushed to the compute stack or not
	virtual bool notEmpty(map<int,bool>& notEmptyMap) const{return true;}
	void setComplete(bool isComplete) { _isComplete=isComplete; }

	void PrintDependencyInfo(){

	    std::cout<<std::endl<<"Dependency Info for OpName : "<<_opName<<"  OpID : "<<_OpID<<std::endl;
	    for(unsigned int i = 0; i < _dInfoVec.size(); i++){
		std::cout<<"Type : ";
		if(_dInfoVec[i]._psi)
		    std::cout<<"PSI  ";
		if(_dInfoVec[i]._alpha)
		    std::cout<<"ALPHA";
		if(_dInfoVec[i]._delta)
		    std::cout<<"DELTA";
		if(_dInfoVec[i]._sigma)
		    std::cout<<"SIGMA";		
		std::cout<<"   On : "<<_dInfoVec[i]._treeName<<"     Tree ID : "<<_dInfoVec[i]._treeID<<"."<<std::endl;
	    }
	}

    public:
	T _sum;
	string _opName;
	unsigned long _OpID;
	KTREE* 	_result;
	bool _isComplete;
	std::vector< DependencyInfo <T, NDIM> > _dInfoVec;
	bool _isReduce;

    };

} /*fuset*/

#endif /* __MADNESS_MRA_PRIMITIVE_OP_H__INCLUDED__ */
