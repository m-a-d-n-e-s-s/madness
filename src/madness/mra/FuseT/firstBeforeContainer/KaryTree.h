#ifndef __fuset_KaryTree_h__
#define __fuset_KaryTree_h__
#include <math.h>
#include "Node.h"
namespace fuset {

  //#define DEBUG 1
    using namespace std;
    enum GeneratedData {ZERO, SOME_NUM, RANDOM};
    typedef unsigned long long TRANS;
  template<typename Type>
  class KaryTree {
    typedef std::map<TRANS, Node<Type>*> TreeDS;
  private:    

    TreeDS* _mapKT;
    int _dataLength;
    int _maxDepth;
    int _numChild;
    bool isInitialized;
    void DeleteTreeOnly();
    //delete data at a tree level from all nodes
    void DeleteData(int level);
    void recurseRandom(int n, TRANS l, GeneratedData gd, Type d);

  public:
    string _treeName;
    static int _rseed;
    static int _numTrees;

    int _treeID;

    inline TRANS childTranslation(TRANS l, int childId) const{
      return l*_numChild + childId;
    }

    inline TRANS parentTranslation(TRANS l) const {
      return l/_numChild;
    }
    
    //!Constructor
    KaryTree(string name, int maxDepth,int numChild, int dataLength);
    
    //!Delete Tree including the data
    ~KaryTree();

    //!generates a random tree with the parameters defined by the constructor
    void generateRandom(GeneratedData gd, Type d);

    //!print out the child tree starting at depth n and translation l
    void printTree(int n, TRANS l);

    Type getNorm();
    
    inline bool hasNode(int depth, TRANS translation){
      return (_mapKT[depth].find(translation) != _mapKT[depth].end()) ;
    }
    inline Node<Type>* getNode(int depth, TRANS translation, bool& isEmpty) const{
      isEmpty = (_mapKT[depth].find(translation) == _mapKT[depth].end());
      if(!isEmpty)
        return (*this)(depth, translation);

      Node<Type>* temp = new Node<Type>();
      return temp;
    
    }

  
    inline Node<Type>* setNode(int depth, TRANS translation, Node<Type>* n){
      _mapKT[depth][translation] = n;
      return _mapKT[depth][translation];
    }
    
    //creates a new node and allocates memory for data at that node
    inline Node<Type>* createNode(int depth, TRANS translation){
      Node<Type>* temp = new Node<Type>(INTERIOR,_dataLength);	
      _mapKT[depth][translation] = temp;
      return _mapKT[depth][translation];
    }

    //!returns the node at given depth and translation
    inline Node<Type>* operator()(int depth, TRANS translation) const{
      bool isEmpty = (_mapKT[depth].find(translation) == _mapKT[depth].end());
      if(!isEmpty)
        return _mapKT[depth][translation];
      std::cerr<<"The node is empty"<<std::endl;
      return new Node<Type>();
    }

    inline int getNumChild() const {return _numChild;}
  };


  template<typename Type> 
  KaryTree<Type>::~KaryTree(){
    if(isInitialized){
      for(int i =0; i< _maxDepth; i++){
        for(typename TreeDS::iterator it = _mapKT[i].begin(); it!=_mapKT[i].end(); it++)
          {
            it->second->deleteData();
          }
        _mapKT[i].clear();
      }
      delete[] _mapKT;
    }
  }


  template<typename Type> 
  void KaryTree<Type>::DeleteTreeOnly(){
    if(isInitialized){
      for(int i =0; i< _maxDepth; i++)
        _mapKT[i].clear();
      delete[] _mapKT;
    }
  }

  template<typename Type> 
      KaryTree<Type>::KaryTree(string name,int maxDepth,int numChild, int dataLength){
      _treeName = name;
    _maxDepth = maxDepth;
    _mapKT = new TreeDS[_maxDepth];
    _dataLength = dataLength;
    _numChild = numChild;
    isInitialized = true;
    _rseed = 0;
    _treeID = KaryTree<Type>::_numTrees++;
  }


  //KaryTree::KaryTree(KaryTree<Type> kt);
  template<typename Type> 
  void KaryTree<Type>::generateRandom(GeneratedData gd, Type d){
    //delete the tree if it exists
    //this->~KaryTree();
    //recreate a clean tree
    //KaryTree(_maxDepth,_numChild,_dataLength);
    //int seed= std::time(NULL)%100 + KaryTree<Type>::_rseed;
    int seed= KaryTree<Type>::_rseed;
    srand(seed);
    std::cout<<"Seed : "<<seed<<std::endl;
    KaryTree<Type>::_rseed+=8;
    recurseRandom(0,0,gd,d);
  }


  template<typename Type> 
  void KaryTree<Type>::recurseRandom(int n, TRANS l, GeneratedData gd, Type d){
    Node<Type>* newNode = createNode(n,l);
    if(gd == ZERO)
      newNode->setZero();
    if(gd == SOME_NUM)
      newNode->setVal(d);
    else
      newNode->setRandom();
    
    double r = (double)(rand() % 1000);

    //probabilty is high at low depth and low at high depth
    double prob = 1.0-pow((double)n/double(_maxDepth),0.5);
    //std::cout<<"("<<r<<"<"<<prob*1000<<")"<<std::endl;
    //creates children nodes with decreasing probability
    if(r<prob*1000.0 && n+1<_maxDepth)
      for(int childId = 0; childId<_numChild; childId++)
        recurseRandom(n+1,this->childTranslation(l,childId), gd,d);
    else
      newNode->setType(LEAF);
  }

  template<typename Type> 
  Type KaryTree<Type>::getNorm(){
      Type retValue = 0;
      for(int d = 0; d<_maxDepth; d++)
	  for(auto p : _mapKT[d])
	      retValue+=p.second->innerLocal();	      
      return pow(retValue,0.5);
  }

  template<typename Type> 
  void KaryTree<Type>::printTree(int n, TRANS l){
    Node<Type>* node = _mapKT[n][l];
    
    for(int i =0; i<n; i++)
      std::cout<<"|       ";
    
    std::cout<<"|---->"<<"Node( "<<n<<", "<<l<<" ) : "
	     <<node->getType()<<" Data : ("
	     <<node->getData()[0]<<", "
	     <<node->getData()[1]<<",... )"<<std::endl;

    if(node->getType()==INTERIOR)
      for(int childId = 0; childId < _numChild; childId++)
        printTree(n+1,this->childTranslation(l,childId));
  }


  template <> int KaryTree<double>::_rseed = 0;
  template <> int KaryTree<double>::_numTrees = 0;

}; /*fuset*/

#endif /*__fuset_KarrayTree_h__*/
