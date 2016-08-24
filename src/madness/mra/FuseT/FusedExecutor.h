
#ifndef __fuset_FusedExecutor_h__
#define __fuset_FusedExecutor_h__

//#include "PrimitiveOp.h"
#include <madness/world/MADworld.h>
#include "FuseT.h"
#include <queue>
#define DEBUG 0
#define DEBUG1 0
namespace madness {

    //!this is local information passed on to each recursive call.
    //!It can be thought of as the data structure holding information
    //!about what needs to computed and what does not
    template<typename T, std::size_t NDIM>
	struct LocalFuseInfo{
	    
	    map<int,bool> _notEmpty; 
	    //input operands for postCompute Ops
	    vector<int> _postOperands;
	
	    vector<int> _preCompute;
	    vector<int> _postCompute;
	    
	    template<typename Archive>
	    void serialize(Archive& ar) { ar & _notEmpty & _postOperands & _preCompute & _postCompute; }
	};

    /*The FusedExecutor takes a fusedOpSequence and computes it by
     doing the appropriate traversal. The traversal is done using
     fusedTraversal method. The fusedTraversal method is a wrapper
     for continueTraversal that actually does the
     work. ContinueTraversal consists of two parts :1) Self
     Recusion : If an operator returns a future, then continue
     traversal calls itself to wait for the future to be
     available. It then continues from where it left off. 2) Calls
     to FusedTraversal: After it computes all the preOps, it calls
     the fusedTraversal at the children nodes.
    */
    template<typename T, std::size_t NDIM>
	class FusedExecutor:public WorldObject<FusedExecutor<T,NDIM> >{
    public:
	typedef FusedExecutor<T,NDIM> feT;
	typedef Key<NDIM>keyT;
	typedef WorldObject<feT> woT;
	typedef WorldContainer<Key<NDIM>,FunctionNode<T,NDIM>>dcT;
	typedef FunctionNode<T,NDIM>KNODE;
	typedef GenTensor<T>coeffT;
	typedef map<int, FuseTContainer<T> > paraMap;
	typedef vector<Future<paraMap > > postParameter;
	
    private:
	World& _world;
	Future <paraMap> fusedTraversal(keyT key, const LocalFuseInfo<T,NDIM>& lInfo, paraMap& pmap); 
	Future <paraMap> continueTraversal(keyT key, const LocalFuseInfo<T,NDIM>& lInfo, paraMap& pMap, LocalFuseInfo<T,NDIM>& newlInfo, vector<paraMap> &pMapVec, int lastIdx, FuseTContainer<T>& lastReturn);

	FusedOpSequence<T,NDIM>* _fOps;
	map<int,FuseTContainer<T> >  fusedPostCompute(keyT Key, const vector<int> postComputeOps, vector<Future<paraMap> >& v);
	dcT _coeffs; 

    public: 
	/*passes world, fOp and coeffs. The last parameters is used to figure out 
	  the distribution of keys amoung multiple MPI ranks*/
    FusedExecutor(World& world, FusedOpSequence<T,NDIM>* fOp): woT(world), _world(world){

	    woT::process_pending();
	const std::shared_ptr<WorldDCPmapInterface<Key<NDIM> > > _pmap(FunctionDefaults<NDIM>::get_pmap());		    
	_coeffs = dcT(world,_pmap,false);
	_coeffs.process_pending();
	_fOps = fOp;

    }
	void execute(); 
    };


    template<typename T, std::size_t NDIM>
	void FusedExecutor<T,NDIM>::execute()
    {
	for(int i = 0; i< _fOps->_validOpDags.size(); i++)
	{
	    //cout<<"Beginning"<<endl<<fflush;
		LocalFuseInfo<T,NDIM> lInfo;
		for(int j : _fOps->_validOpDags[i]._preOps)
		    lInfo._preCompute.push_back(j);	  
		//cout<<"Middle"<<endl<<fflush;
		for(int j : _fOps->_validOpDags[i]._postOps)
		    lInfo._postCompute.push_back(j);
		//cout<<"Second Middle"<<endl<<fflush;
		for(int j : _fOps->_validOpDags[i]._postOperands)
		    lInfo._postOperands.push_back(j);
		//contains a mapping of parameters for each operation
		//the key is the order in which the operations appear in _fOps->sequence
		//cout<<"Second Execution Time"<<endl<<fflush;
		if(_world.rank() == _coeffs.owner(keyT(0))){  
		    paraMap pMap;
		    fusedTraversal(keyT(0),lInfo,pMap);
		}	    
		if(DEBUG) cout<<"Exited from traversal before fence"<<endl;		
		_world.gop.fence();
		if(DEBUG) cout<<"Exited from traversal after fence"<<endl;		
		//check if any reduction is required
		for(int j : _fOps->_validOpDags[i]._preOps)
		    if(_fOps->_sequence[j]->_isReduce)
			_fOps->_sequence[j]->reduce(_world);
		    
		for(int j : _fOps->_validOpDags[i]._postOps)
		    if(_fOps->_sequence[j]->_isReduce)
			_fOps->_sequence[j]->reduce(_world);		
		
	}
    }

    //typedef map<int, FuseTContainer<T> > paraMap;
    template<typename T, std::size_t NDIM>
	Future <map<int,FuseTContainer<T> > > 
	FusedExecutor<T,NDIM>::continueTraversal(keyT key, const LocalFuseInfo<T,NDIM>& lInfo, paraMap& pMap, LocalFuseInfo<T,NDIM> &newlInfo, vector<paraMap> &pMapVec, int lastIdx, FuseTContainer<T>& lastReturn){

	//process the preCompute Op that was computed before this call
	if(lastIdx != -1){
	    int i = lInfo._preCompute[lastIdx];	   

	    PrimitiveOp<T,NDIM> *  preOp = _fOps->_sequence[i];
	    //Control Code for identifying if this operator needs
	    //to be executed in the next executive call	  
	    if(!preOp->isDone(key)){
		newlInfo._preCompute.push_back(i);

		//creating parameters for the recursive call
		//creates pMapVec[j][i] represents the parameters for operator i
		//anf the jth child node
		int j =0;
		for (KeyChildIterator<NDIM> kit(key); kit; ++kit,j++)				
		{
		    //parameter must be passed so the return type of preOp must ne FuseT_VParameter
		    if (preOp->needsParameter())
			pMapVec[j][i] = ((FuseT_VParameter<T>*)lastReturn.get())->value[j];
		    //parameter doesnt have to be passed so create an empty Container		
		    else
			pMapVec[j][i] = FuseTContainer<T>();
		}
	    }
	}

	    //process the remainder of the preCompute Operators
	for(int idx = lastIdx+1; idx<lInfo._preCompute.size(); idx++)
	{
	    int i = lInfo._preCompute[idx];
	    
	    PrimitiveOp<T,NDIM> *  preOp = _fOps->_sequence[i];

	    if(preOp->returnsFuture()){

		Future<FuseTContainer <T> > temp = preOp->computeFuture(key,pMap[i]);
		return woT::task(_world.rank(), &feT::continueTraversal, key, lInfo, pMap, newlInfo, pMapVec, idx, temp);
	    }

	    //Real Work
	    FuseTContainer<T> temp;
	    temp = preOp->compute(key,pMap[i]);

	    //Control Code for identifying if this operator needs
	    //to be executed in the next executive call	  
	    if(!preOp->isDone(key)){
		newlInfo._preCompute.push_back(i);

		//creating parameters for the recursive call
		//creates pMapVec[j][i] represents the parameters for operator i
		//anf the jth child node
		int j =0;
		for (KeyChildIterator<NDIM> kit(key); kit; ++kit,j++)				
		{
		    //parameter must be passed so the return type of preOp must ne FuseT_VParameter
		    if (preOp->needsParameter())
			pMapVec[j][i] = ((FuseT_VParameter<T>*)temp.get())->value[j];
		    //parameter doesnt have to be passed so create an empty Container		
		    else
			pMapVec[j][i] = FuseTContainer<T>();
		}
	    }
	}

	if(DEBUG) cout<<"After PreCompute"<<endl;
		
	//Control Code check if any of the producers
	//for the postOps are already empty. Some false positive might
	//be registered in this loop as some of the producers that are
	//postOps might not be completed. These will be rectified in the
	//next loop      
	for(int i : lInfo._postOperands){
			
	    const Function<T,NDIM>* source = _fOps->_trees[i];       
	    unsigned long treeID = source->get_impl()->id().get_obj_id();



	    bool exists, notEmpty;
	    exists = source->get_impl()->get_coeffs().probe(key);
            if(exists){
		notEmpty = source->get_impl()->get_coeffs().find(key).get()->second.has_children();
	    }else
		notEmpty = false;



	    newlInfo._notEmpty[treeID] = notEmpty;
	    if(DEBUG) cout<<"Source Tree : "<<source->_treeName<<" Tree ID "<<treeID<<" Not Empty : "<<newlInfo._notEmpty[treeID]<<endl;

	}

	if(DEBUG) cout<<endl<<"After notEmpty"<<endl;

	//identifying if a postCompute op needs to be computed at child locations
	for(int i : lInfo._postCompute){
	    PrimitiveOp<T,NDIM>* postOp = _fOps->_sequence[i];
	    if(postOp->notEmpty(newlInfo._notEmpty)){
		if(DEBUG1) cout<<endl<<"Not empty postOp : "<<postOp->_OpID<<endl;
		newlInfo._notEmpty[postOp->_OpID]  = true;
		newlInfo._postCompute.push_back(i);
	    }
	}

	if(DEBUG) cout<<endl<<"After Finding Post Compute Operators"<<endl;

	//identifying postOperands for postCompute
	for(int i : newlInfo._postCompute)
	    for(int j : _fOps->_operandTrees[i])
		newlInfo._postOperands.push_back(j);

	if(DEBUG) cout<<"After Finding Post Operands"<<endl;
	
	postParameter v;	// FuseTConatiner has coeffT (not vector)
	if(!newlInfo._preCompute.empty() || !newlInfo._postCompute.empty()){
	    if(DEBUG) cout<<"Recursive Call"<<endl;
	    // Main-Computation
	    v = future_vector_factory<paraMap>(1<<NDIM);
	    int i =0;
	    for (KeyChildIterator<NDIM> kit(key); kit; ++kit, i++) {
		const keyT& child = kit.key();
		v[i] = woT::task(_coeffs.owner(child), &feT::fusedTraversal, child, newlInfo, pMapVec[i]);
	    }
	}

	if(DEBUG) cout<<"After Recursive Call"<<endl;
	if(lInfo._postCompute.size() > 0){
	    if(DEBUG) cout<<"Calling Post Compute"<<endl;
	    return woT::task(_world.rank(), &feT::fusedPostCompute, key, lInfo._postCompute, v);
	}

	if(DEBUG) cout<<"No Post Compute . Exit"<<endl;
	paraMap temp;
	return Future<paraMap>(temp);	

    }


    //typedef map<int, FuseTContainer<T> > paraMap;
    template<typename T, std::size_t NDIM>
	Future <map<int,FuseTContainer<T> > > 
	FusedExecutor<T,NDIM>::fusedTraversal(keyT key, const LocalFuseInfo<T,NDIM>& lInfo, paraMap& pMap)
    {
	if(DEBUG1){ cout<<" Key : "<<key.level()<<" Translation : " ;
	for(auto a : key.translation())
	    cout<<a<<", ";
	cout<<endl;
	cout<<"PreCompute Size : "<<lInfo._preCompute.size()<<endl;
	cout<<"postCompute Size : "<<lInfo._postCompute.size()<<endl;	
	}
	LocalFuseInfo<T,NDIM> newlInfo;
	//will hold new parameters
	vector<paraMap> pMapVec = vector<paraMap>(1<<NDIM);
	FuseTContainer<T> temp = FuseTContainer<T>();
	//cout<<"About to enter continue traversal"<<endl;
	return continueTraversal(key,lInfo,pMap,newlInfo,pMapVec, -1, temp);
    }

        template<typename T, std::size_t NDIM>
	    map<int,FuseTContainer<T> >  FusedExecutor<T,NDIM>::fusedPostCompute(keyT key, const vector<int> postComputeOps, vector<Future<paraMap> >& v){

	if(DEBUG) cout<<"Entering Post Compute"<<endl;
	paraMap returnPara;
	for(int i : postComputeOps){
	    
	    //vector of return parameters form each child for this operator
	    vector<FuseTContainer<T> > temp;
	    //the first condition makes sure that v is not empty
	    //the second condition makes sure that operator i produces coeffs at a child node
	    if(v.size() > 0 && v[0].get().find(i) != v[0].get().end()){
		int j = 0;
		for (KeyChildIterator<NDIM> kit(key); kit; ++kit, j++) {
		    FuseTContainer<T> t = v[j].get()[i];
		    temp.push_back(t);
		}
	    }
	    FuseTContainer<T> para(static_cast<Base<T>*>(new FuseT_VParameter<T>(temp)));
	    PrimitiveOp<T,NDIM>* postOp = _fOps->_sequence[i];
	    if(DEBUG) cout<<"Computing in  Post Compute"<<endl;
	    returnPara[i] = postOp->compute(key,para);
	}
	if(DEBUG) cout<<"Exiting Post Compute"<<endl;
	return returnPara;

       

    }

    

    
}; /*fuset*/

#endif /*__fuset_FusedExecutor_h__*/

