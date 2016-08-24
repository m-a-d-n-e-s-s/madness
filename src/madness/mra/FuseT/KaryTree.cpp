#include "KaryTree.h"
#include <cmath>
using namespace std;
template<typename DTYPE> 
KaryTree<DTYPE>::~KaryTree(){
    if(isInitialized){
	for(int i =0; i< _maxDepth; i++){
	    for(typename TreeDS::iterator it = _mapKT[i].begin(); it!=_mapKT[i].end(); it++)
	    {
		it->second->DeleteData();
	    }
	    _mapKT[i].clear();
	}
	delete[] _mapKT;
    }
}


template<typename DTYPE> 
void KaryTree<DTYPE>::DeleteTreeOnly(){
    if(isInitialized){
	for(int i =0; i< _maxDepth; i++)
	    _mapKT[i].clear();
	delete[] _mapKT;
    }
}

template<typename DTYPE> 
KaryTree<DTYPE>::KaryTree(int maxDepth,int numChild, int dataLength){
    _maxDepth = maxDepth;
    _mapKT = new TreeDS[_maxDepth];
    _dataLength = dataLength;
    _numChild = numChild;
    isInitialized = true;
    _rseed = 0;
}


//KaryTree::KaryTree(KaryTree<DTYPE> kt);
template<typename DTYPE> 
void KaryTree<DTYPE>::GenerateRandom(GeneratedData gd, DTYPE d){
    //delete the tree if it exists
    //this->~KaryTree();
    //recreate a clean tree
    //KaryTree(_maxDepth,_numChild,_dataLength);
    int seed= std::time(NULL)%100 + KaryTree<DTYPE>::_rseed;
    srand(seed);
    std::cout<<"Seed : "<<seed<<" _rseed"<<_rseed<<std::endl;
    KaryTree<DTYPE>::_rseed+=8;
    RecurseRandom(0,0,gd,d);
}


template<typename DTYPE> 
void KaryTree<DTYPE>::RecurseRandom(int n, int l, GeneratedData gd, DTYPE d){
    Node<DTYPE>* newNode = CreateNode(n,l);
    if(gd == ZERO)
	newNode->SetZero();
    if(gd == SOME_NUM)
	newNode->SetNum(d);
    else
	newNode->SetRandom();
    
    double r = (double)(rand() % 1000);

    //probabilty is high at low depth and low at high depth
    double prob = 1.0-pow((double)n/double(_maxDepth),0.5);
    //std::cout<<"("<<r<<"<"<<prob*1000<<")"<<std::endl;
    //creates children nodes with decreasing probability
    if(r<prob*1000.0 && n+1<_maxDepth)
	for(int childId = 0; childId<_numChild; childId++)
	    RecurseRandom(n+1,this->ChildTranslation(l,childId), gd,d);
    else
	newNode->SetType(LEAF);
}



template<typename DTYPE> 
void KaryTree<DTYPE>::PrintTree(int n, int l){
    Node<DTYPE>* node = _mapKT[n][l];
    
    for(int i =0; i<n; i++)
	std::cout<<"|       ";
    
    std::cout<<"|---->"<<"Node( "<<n<<", "<<l<<" ) : "
	     <<node->GetType()<<" Data : ("
	     <<node->GetData()[0]<<", "
	     <<node->GetData()[1]<<",... )"<<endl;

    if(node->GetType()==INTERIOR)
	for(int childId = 0; childId < _numChild; childId++)
	    PrintTree(n+1,this->ChildTranslation(l,childId));
    

}


template class KaryTree<double>;
//template int KaryTree<double>::ChildTranslation(int a, int b);


