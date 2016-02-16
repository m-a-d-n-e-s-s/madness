#ifndef __MADNESS_MRA_FUSET_PRIMITIVE_OP_H__INCLUDED__
#define __MADNESS_MRA_FUSET_PRIMITIVE_OP_H__INCLUDED__

#include "madness/mra/mra.h"
#include "madness/world/archive.h"
#include <string>
#include "BaseParameters.h"

namespace madness {

    using namespace std;
    /*!spatial relation, Same Node Relation (PSI), Some Sibling Relation (SIGMA), Ancestor Relation (ALPHA), Descendent Relation (DELTA). Dependency Info contains information about type of dependency to parameters*/
    template<typename T, std::size_t NDIM>
    struct DependencyInfo{
	const Function<T,NDIM>* _producerTree;
	string _treeName;
        unsigned long _treeID;

	bool _psi;
	bool _alpha;
	bool _delta;
	bool _sigma;


	/*!this is used by FuseT. _producerIndex stores the index of operator that produces this tree/ _consumerIndex stores the index of the operator that consumes this tree*/
	int _producerIndex;
	int _consumerIndex;
	DependencyInfo(const Function<T,NDIM>* producerTree, bool psi, bool alpha, bool delta, bool sigma){
	    _psi =psi;
	    _alpha = alpha;
	    _delta=delta;
	    _sigma=sigma;
	    
	    _producerTree = producerTree;
	    _treeID=_producerTree->get_impl()->id().get_obj_id();
	    _treeName = _producerTree->_treeName;
	   
	}
    };

//    struct DependencyInfo<double>;
    
    template<typename T, std::size_t NDIM>
	class PrimitiveOp{
	typedef Function<T,NDIM> KTREE;
	typedef FunctionNode<T,NDIM> KNODE;
	typedef Key<NDIM> keyT;
	typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT; ///< Type of container holding the nodes


    public:
    PrimitiveOp(string opName="Unknown",KTREE* result=NULL, bool isComplete=false)  
	: _opName(opName),
	_result(result),
	    _isComplete(isComplete)
	    {}

	virtual ~PrimitiveOp() {}
    
	virtual void compute(const keyT& key)		= 0;
	//virtual BaseParameters<T,NDIM>& computeB(const keyT& key) = 0;
	virtual Future<GenTensor<T>> computeB(const keyT& key)		= 0;
	virtual Future<GenTensor<T>> afterComputeB(const keyT& key, const std::vector<Future<GenTensor<T>>> &v) = 0;

	virtual bool isDone(const keyT& key) const	= 0;
	virtual bool isPre() const					= 0;
 
	//!used for postCompute ops to see if it needs to be pushed to the compute stack or not
	virtual bool notEmpty(map<int,bool>& emptyMap) const{return true;}
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
	string _opName;
	unsigned long _OpID;
	KTREE* _result;
	bool _isComplete;
	std::vector< DependencyInfo <T, NDIM> > _dInfoVec;
    };


//    namespace archive {
//	/// Serialize an AST
//	template <class Archive, typename T, std::size_t NDIM>
//	    struct ArchiveStoreImpl< Archive, PrimitiveOp<T,NDIM> > {
//
//	    typedef Function<T,NDIM> KTREE;
//	    typedef FunctionImpl<T,NDIM> implT;
//	    
//	    static void store(const Archive& s, const PrimitiveOp<T, NDIM>& t) {
//		s & t._opName;
//		s & t._result.get();
//	    }
//	};
//
//	/// Deserialize an AST
//	template <class Archive, typename T, std::size_t NDIM>
//	    struct ArchiveLoadImpl< Archive, PrimitiveOp<T,NDIM> > {
//
//	    typedef Function<T,NDIM> KTREE;
//	    typedef FunctionImpl<T,NDIM> implT;
//	    
//	    static void load(const Archive& s, PrimitiveOp<T,NDIM>& t) {
//
//		string opName = "";
//		implT* impl = NULL;
//		shared_ptr<KTREE> result;
//		bool isComplete=false;		
//
//		s & opName;
//		s & impl & isComplete;
//		result = impl;
//		PrimitiveOp<T,NDIM> temp(opName,result,isComplete);
//		t = temp;
//	    }
//	};
//	
//    }
//    
} /*fuset*/

#endif /* __MADNESS_MRA_PRIMITIVE_OP_H__INCLUDED__ */
