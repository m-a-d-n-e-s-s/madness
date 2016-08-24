#ifndef __fuset_FuseT_h__
#define __fuset_FuseT_h__
#define DEBUG 0

#include "PrimitiveOp.h"
#include <algorithm>
namespace madness {
    using namespace std;

        struct ValidOpDag{
	vector<int> _preOps;
	vector<int> _postOps;
	vector<int> _postOperands;
    };

    /*! A Node in the OpDAG for _sequence*/
	template<typename T, std::size_t NDIM>
    struct FusedOpSequence{
	vector<PrimitiveOp<T,NDIM>*> _sequence;
	vector<const Function<T,NDIM>*> _trees;
	vector<ValidOpDag> _validOpDags;
	vector<vector<int> > _operandTrees;	

	FusedOpSequence(vector<PrimitiveOp<T,NDIM>*> sequence, vector<const Function<T,NDIM>*> trees){
	    _sequence = sequence;
	    _trees = trees;
	}
    };


	//vector<DependencyInfo> _consumerList;
	//
	////helper variables
	//bool _isCompleted, _isVisitd, _isPre, _isPost;       
	//};

	template<typename T, std::size_t NDIM>
	class FuseT{
	typedef DependencyInfo<T,NDIM> DInfo;
    private:
	//!creates an opDAG that represents the sequence
	//!of operators in _sequence
	void generateOpDAG();

	//!generates a sequence of valid OpDAGs from the opDAG
	//!for sequence of operators in _sequence
	void getValidOpDAGSequence();
	
	/*! This method bins operators in each validOpDAG 
	  into pre or post compute ops*/
	void getPrePostOps();

	//! Once a preOp is identified, this methods adds
	//!all the producers ops in the chain in preOps
	void addPreOps(int opIndex, int validOpDagID);

	//! Once a postOp is identified, this method adds
	//!all the consumers ops in this chain in postOps
	void addPostOps(int opIndex, int validOpDagID);

	//!Helper Method for getPrePostOps.
	//!Identifies if an oprator should be in preOps 
	inline bool isPostOp(int opIndex, int validOpDagID)
	{
	    if(DEBUG) cout<<endl<<"Is Post "<<opIndex<<endl;
	    for (DInfo pn : _producerList[opIndex]){
		//if cn has ancestor dependency on op
		int index = pn._producerIndex;
		if(DEBUG) cout<<"Producer "<<index<<endl;
		//(index < sequence.size()) for all the operators
		//it is only false for source trees
		if((index < _sequence.size()) && pn._delta &&
		   _idValidOpDAG[index]==validOpDagID){
		    if(DEBUG) cout<<"Returning true"<<endl<<endl;
		    return true;
		}
	    }
	    if(DEBUG) cout<<"Returning false"<<endl<<endl;
	    return false;

	}
	
	//!Helper Method for getPrePostOps.
	//!Identifies if an oprator should be in postOps 
	bool isPreOp(int opIndex, int validOpDagID){
	    if(DEBUG) cout<<endl<<"Is Pre "<<opIndex<<endl;
	    for (DInfo cn : _consumerList[opIndex]){
		if(DEBUG) cout<<"Consumer "<<cn._consumerIndex<<endl;
		//if cn has ancestor dependency on op
		if(cn._alpha &&
		   _idValidOpDAG[cn._consumerIndex]==validOpDagID){
		    if(DEBUG) cout<<"Returning true"<<endl<<endl;
		    return true;
		}
	    }
	    if(DEBUG) cout<<"Returning false"<<endl<<endl;
	    return false;
	}

	/*! IsSource is used to identify if an operator in a valid opDAG is a source operator or not. It returns true if all the input operands have already been computed by the time the excution gets to the current valid opDAG*/
	inline bool isSource(int opIndex) const{
	    bool source = true;
	    for(DInfo di : _producerList[opIndex])
		if(di.producerIndex < _numOps)
		    source = source && _isCompleted[di._producerIndex];
	    return source;
	}
	
	
    private: 
	//number of operators
	int _numOps;
	
	/*! Pointer to the sequence of operators in the program. 
	  Each element of this vector corresponds to a tree traversal operator, in the order that it is meant to be executed */
	vector<PrimitiveOp<T,NDIM>*> _sequence;

	/*! Pointer to all the trees in the program. _sequence[i]->result 
 	  is the same pointer as _trees[i]. If _sequence[i] is does not exist, then          the tree is not created by any operators in the opDAG */
	
	vector<const Function<T,NDIM>*> _trees;	
	
	/*!maps the operator/tree id to the index value in _sequence. 
	  For example, if (*_sequence)[10]->_OpID = i, then _opIdToIndex[i] = 10*/
	std::map<int,int> _opIdToIndex;
	
	/*!This will hold the output of the the algorithm for generating seqeuence of valid OpDAGs from _sequence. 
	  This provides a sequence of OpDAGs that can be fused together*/
	vector<vector<int>> _validOpDAGs;
	
	/*!stores pre Operators and post Operators for each of the valid opDAGs*/       
	vector<vector<int> > _preOps;
	vector<vector<int> > _postOps;
	
	/* The following set of variables describe attributes 
	  for each operator or opNode in the OpDag.
	  /////////////////////////////////////////////////////////////////////
	  ///////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////
	
	! for each operator this stores that validOpDagID that it belongs to
	  this is redundant information, that can be obtained from _validOpDags, 
	  stored just for convienence in isPreOp and isPostOp methods*/
	vector<int> _idValidOpDAG;
	
	/*! This is a directed acyclic graph where each node contains a vector containing all its producers, and the type of dependency it has with that producer. _producerList[3]={0,1} implies that Operator (*_sequence)[3] depends on operators (*_sequence)[0] (*_sequence)[1] */ 
	vector<vector<DInfo > > _producerList;       

	/*! Similar to _producerList but contains all the consumers*/
	vector<vector<DInfo> > _consumerList;

	//!during generation of preOps and postOps, isCompleted is used to indentify if a operator has already been completed
	vector<bool> _isCompleted;

	//!for computing preOPs and postOps
	vector<bool> _isVisited;

	//!stores if an operator is in pre or post compute
	vector<bool> _isPre;
	vector<bool> _isPost;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//true after processing is complete. Used by print to
	//figure out if the analysis is complete
	bool _isSequenceProcessed;

    public:

	FuseT(vector<PrimitiveOp<T,NDIM>*> sequence){
	    _sequence = sequence;
	    _numOps = _sequence.size();
	    _isSequenceProcessed = false;
	}
    
	void processSequence(){

	    generateOpDAG();
	    getValidOpDAGSequence();
	    getPrePostOps();
	    _isSequenceProcessed = true;
	    
	}
	
	void printValidSequences();

	void printOpsAndTrees();
	
	FusedOpSequence<T,NDIM> getFusedOpSequence();


    };


    template<typename T, std::size_t NDIM>
	FusedOpSequence<T,NDIM> FuseT<T,NDIM>::getFusedOpSequence(){

	//set _sequence and _trees
	FusedOpSequence<T,NDIM> temp(_sequence,_trees);

	//set _validOpDags
	int numValidOpDags = _validOpDAGs.size();	
	for(int i = 0; i< numValidOpDags; i++){
	    ValidOpDag vdag;
	    vdag._preOps = _preOps[i];
	    vdag._postOps = _postOps[i];

	    for(int op : vdag._postOps)
		for(auto di : _producerList[op])
		    vdag._postOperands.push_back(di._producerIndex);

	    //sort and eleminate redundancy
	    auto *pOps = &vdag._postOperands;
	    sort(pOps->begin(), pOps->end());
	    auto it = unique(pOps->begin(),pOps->end());
	    pOps->resize(distance(pOps->begin(),it));

	    //adding valid opDAGs
	    temp._validOpDags.push_back(vdag);
	    	    
	}

	//set OperandTrees
	for(int i =0; i<_producerList.size();i++){
	    temp._operandTrees.push_back(vector<int> ());
	    for(auto di: _producerList[i])
		temp._operandTrees[i].push_back(di._producerIndex);
	}
	
	return temp;
	
    }


    template<typename T, std::size_t NDIM>
	void FuseT<T,NDIM>::addPreOps(int opIndex, int validOpDagID){
	if (DEBUG) cout<<"Checking for adding to Pre ID : "<<opIndex<<endl;
	if(!_isVisited[opIndex]) {
	    _isVisited[opIndex] = true;
	    _isPre[opIndex] = true;
	    if(DEBUG) cout<<"Adding to Pre ID : "<<opIndex<<endl;
	    _preOps[validOpDagID].push_back(opIndex);
	    for (DInfo di: _producerList[opIndex]){ 
		int index = di._producerIndex;
		if((index<_numOps) &&_idValidOpDAG[index] == validOpDagID){
		    if(DEBUG) cout<<"Producer ID of Pre : "<<index<<endl;
		    addPreOps(index,validOpDagID);
		}
	    }
	}
    }
    
    template<typename T, std::size_t NDIM>
	void FuseT<T,NDIM>::addPostOps(int opIndex, int validOpDagID){
	if (DEBUG) cout<<"Checking for adding to Post ID : "<<opIndex<<endl;
	if(!_isVisited[opIndex]) {
	    _isVisited[opIndex] = true;
	    _isPost[opIndex] = true;
	    if(DEBUG) cout<<"Adding to Post ID : "<<opIndex<<endl;
	    _postOps[validOpDagID].push_back(opIndex);
	    for (DInfo di: _consumerList[opIndex]) 
		if(_idValidOpDAG[di._consumerIndex] == validOpDagID)
		{
		    if(DEBUG) cout<<"Consumer ID of Post : "<<di._consumerIndex<<endl;
		    addPostOps(di._consumerIndex,validOpDagID);

		}
	}
    }
    
    /*! fills up _preOps and _postOps for each valid OpDag*/
    template<typename T, std::size_t NDIM>
	void FuseT<T,NDIM>::getPrePostOps(){

	_isCompleted.assign(_numOps,false);
	_isVisited.assign(_numOps,false);
	_isPre.assign(_numOps,false);
	_isPost.assign(_numOps,false);

	for(int i =0; i< _validOpDAGs.size(); i++)
	{
	    vector<int> validOpDAG = _validOpDAGs[i];
//	    vector<int> source;	    
//
//	    for(int op : validOpDAG)
//		if (isSource(op))
//		    source.push_back(op);

	    for (int op : validOpDAG)
		if (!_isVisited[op] && isPreOp(op,i))
		    addPreOps(op,i);

	    for (int op : validOpDAG)
		_isVisited[op] = false;

	    for (int op : validOpDAG)
		if (!_isVisited[op] && isPostOp(op,i))
		    addPostOps(op,i);

	    for (int op : validOpDAG) 
		if (_isPost[op] && _isPre[op])
		    std::cerr<<"Should not be happening"<<std::endl;

	    for(int op : validOpDAG)
		if(!_isPost[op] && !_isPre[op])
		    _preOps[i].push_back(op);

	    std::sort(_preOps[i].begin(),_preOps[i].end());
	    std::sort(_postOps[i].begin(),_postOps[i].end());
	    
	    
	    for(int op : validOpDAG)
		_isCompleted[op] = true;
	}
    }
	
    /*! fills up _opIdToIndex, _producerList, and _consumerList*/
    template<typename T, std::size_t NDIM>
	void FuseT<T,NDIM>::generateOpDAG(){

	//set up _opIdToIndex
	for(int i =0;i<_sequence.size();i++){
	    _opIdToIndex[_sequence[i]->_OpID] = i;
	    _producerList.push_back(vector<DInfo >());
	    _consumerList.push_back(vector<DInfo >());
	    _trees.push_back(_sequence[i]->_result);
	}

	for(int i =0; i<_sequence.size(); i++){

	    int numProducers = _sequence[i]->_dInfoVec.size();

	    for(int j = 0; j<numProducers;j++){

		DInfo treeInfo = _sequence[i]->_dInfoVec[j];
		
		/*!treeID and opID are the same for the op that produces the tree
		  if the treeID does not exist in the sequence of operators 
		  then it must mean that the tree was not produced by any 
		  operator in the program. It might be a true input.*/
		if(_opIdToIndex.find(treeInfo._treeID) == _opIdToIndex.end()){
		    _trees.push_back(treeInfo._producerTree);
		    _opIdToIndex[treeInfo._treeID]  = _trees.size() - 1 ;
		}

		treeInfo._producerIndex = _opIdToIndex[treeInfo._treeID];
		treeInfo._consumerIndex = i;


		//pure inputs are also inserted as part of the producerList.
		//This is used by postOps for identifying which ops are done
		//and which are not
		_producerList[i].push_back(treeInfo);
		//cout<<treeInfo._producerIndex<<endl;

		//these are not operators, but trees source
		//input that have already been computed. No need to
		//find the consumers for these trees 
		if(treeInfo._producerIndex < _sequence.size())
		   _consumerList[treeInfo._producerIndex].push_back(treeInfo);
		
	    }
	    
	}
    }
    
    //this is incoplete right now. It just returns the entire sequence for now
    template<typename T, std::size_t NDIM>
	void FuseT<T,NDIM>::getValidOpDAGSequence(){

	_validOpDAGs.push_back(vector<int>());
	_idValidOpDAG.assign(_sequence.size(), -1);
	for(int i = 0; i< _sequence.size(); i++){
	    _validOpDAGs[0].push_back(i);
	    _idValidOpDAG[i] = 0;	   
	}

	//initialize _preOps and _postOps
	for(int i =0; i< _validOpDAGs.size(); i++)
	{
	    _preOps.push_back(vector<int>());
	    _postOps.push_back(vector<int>());
	    
	}

    }

    template<typename T, std::size_t NDIM>
	void FuseT<T,NDIM>::printValidSequences(){
	if(!_isSequenceProcessed){
	    cout<<"Sequence not processed!"<<endl;
	    return;
	}

	cout<<endl<<"Original Sequence "<<endl<<endl;
	for(int i =0; i<_numOps; i++)
	    cout<<i<<" "<<_sequence[i]->_opName<<" OpID : "<<_sequence[i]->_OpID<< endl;
	cout<<"________________________"<<endl;
	cout<<"________________________"<<endl;

	
	cout<<endl<<endl<<"Fused Sequences"<<endl<<endl;
	for(int i =0; i< _validOpDAGs.size(); i++)
	    {
		cout<<"Valid OpDAG Sequence : "<<i<<endl;
		cout<<"______Pre Ops_______ : "<<endl;
		for(int id : _preOps[i])
		    cout<<_sequence[id]->_opName<<endl;
		cout<<"_____Post Ops_______ : "<<endl;
		for(int id : _postOps[i])
		    cout<<_sequence[id]->_opName<<endl;
		cout<<"____________________"<<endl;

	    }	
    }
    

        template<typename T, std::size_t NDIM>
	void FuseT<T,NDIM>::printOpsAndTrees(){
	    if(!_isSequenceProcessed){
		cout<<"Sequence not processed!"<<endl;
		return;
	    }

	    cout<<endl<<"Ops and Trees "<<endl<<endl;
	    for(int i =0; i<_numOps; i++)
		cout<<i<<" "<<_sequence[i]->_opName<<" OpID : "<<_sequence[i]->_OpID<<" TreeID :"<<_trees[i]->get_impl()->id().get_obj_id()<<" TreeName : "<<_trees[i]->_treeName<< endl;
	    cout<<"________________________"<<endl;
	    cout<<"________________________"<<endl;
	    cout<<endl<<"Source Trees "<<endl<<endl;
	    for(int i =_numOps; i<_trees.size(); i++)
		cout<<i<<" "<<" TreeID :"<<_trees[i]->get_impl()->id().get_obj_id()<<" TreeName : "<<_trees[i]->_treeName<< endl;
	    cout<<"________________________"<<endl;
	    cout<<"________________________"<<endl;
	    
	}

    
};
#endif
