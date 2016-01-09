#ifndef __fuset_FusedExecutor_h__
#define __fuset_FusedExecutor_h__

//#include "PrimitiveOp.h"
#include "FuseT.h"
#include <queue>
//#define DEBUG 0
namespace fuset {

    //!this is local information passed on to each recursive call.
    //!It can be thought of as the data structure holding information
    //!about what needs to computed and what does not
    template<typename Type>
	struct LocalFuseInfo{
	    
	    map<int,bool> _notEmpty; 
	    //input operands for postCompute Ops
	    vector<int> _postOperands;
	
	    vector<int> _preCompute;
	    vector<int> _postCompute;
	};

    template<typename Type>
	class FusedExecutor{
    private:
	void fusedTraversal(const LocalFuseInfo<Type>& lInfo, int n, TRANS l);
	FusedOpSequence<Type>* _fOps;

    
    public:
	FusedExecutor(FusedOpSequence<Type>* fOp);
	void execute();    
    };

    template<typename Type>
	FusedExecutor<Type>::FusedExecutor(FusedOpSequence<Type>* fOp)
    {
	_fOps = fOp;
    }
    
    

    template<typename Type>
	void FusedExecutor<Type>::execute(){

	LocalFuseInfo<Type> lInfo;
	for(int i = 0; i< _fOps->_validOpDags.size(); i++){
	  
	    for(int j : _fOps->_validOpDags[i]._preOps)
		lInfo._preCompute.push_back(j);	  
	  
	    for(int j : _fOps->_validOpDags[i]._postOps)
		lInfo._postCompute.push_back(j);
	  
	    for(int j : _fOps->_validOpDags[i]._postOperands)
		lInfo._postOperands.push_back(j);
	  
	    fusedTraversal(lInfo,0,0);	  
	}
    }
  
    template<typename Type>
	void FusedExecutor<Type>::fusedTraversal(const LocalFuseInfo<Type>& lInfo, int n, TRANS l){
	if(DEBUG) cout<<"n : "<<n<<" l : "<<l<<endl;
	LocalFuseInfo<Type> newlInfo;

	if(DEBUG) cout<<"Before PreCompute"<<endl;
	for(int i: lInfo._preCompute){
	    PrimitiveOp<Type> *  preOp = _fOps->_sequence[i];

	    //Real Work
	    preOp->compute(n,l);

	    //Control Code for identifying if this operator needs
	    //to be executed in the next executive call	  
	    if(!preOp->isDone(n,l))
		newlInfo._preCompute.push_back(i);

	}

	if(DEBUG) cout<<"After PreCompute"<<endl;
	
	//Control Code check if any of the producers
	//for the postOps are already empty. Some false positive might
	//be registered in this loop as some of the producers that are
	//postOps might not be completed. These will be rectified in the
	//next loop      
	for(int i : lInfo._postOperands){
	    
	    const KaryTree<Type>* source = _fOps->_trees[i];       
	    if(DEBUG) cout<<source->_treeName<<endl;
	    bool isE;	    
	    Node<Type>* tnode = source->getNode(n,l,isE);
	    if(isE || tnode->getType() == LEAF)
		newlInfo._notEmpty[source->_treeID] = false;
	    else
		newlInfo._notEmpty[source->_treeID] = true;
	}

	if(DEBUG) cout<<"After notEmpty"<<endl;

	//identifying if a postCompute op needs to be computed at child locations
	for(int i : lInfo._postCompute){
	    PrimitiveOp<Type>* postOp = _fOps->_sequence[i];
	    if(postOp->notEmpty(newlInfo._notEmpty)){
		if(DEBUG) cout<<"Not empty postOp : "<<postOp->_result->_treeName<<endl;
		newlInfo._notEmpty[postOp->_OpID]  = true;
		newlInfo._postCompute.push_back(i);
	    }
	}

	if(DEBUG) cout<<"After Finding Post Copmpute Operators"<<endl;

	//identifying postOperands for postCompute
	for(int i : newlInfo._postCompute)
	    for(int j : _fOps->_operandTrees[i])
		newlInfo._postOperands.push_back(j);

	if(DEBUG) cout<<"After Finding Post Operands"<<endl;
	
	if(!newlInfo._preCompute.empty() || !newlInfo._postCompute.empty())
	    for(int childId = 0; childId < _fOps->_trees[0]->getNumChild(); childId++)
		fusedTraversal(newlInfo,n+1,_fOps->_trees[0]->childTranslation(l,childId));


	for(int i : lInfo._postCompute){
	    PrimitiveOp<Type>* postOp = _fOps->_sequence[i];
	    postOp->compute(n,l);

	}

    }

    
}; /*fuset*/

#endif /*__fuset_FusedExecutor_h__*/

