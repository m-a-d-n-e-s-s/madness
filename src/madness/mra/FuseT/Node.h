#include <iostream>
#include <cstdlib>
#include <map>
#include <ctime>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>


namespace fuset {

    
  //#define DEBUG 1
  using namespace std;
  enum NodeType {EMPTY, INTERIOR, LEAF};

  
  template<typename Type>
  class Node {
  public:
    Node();
    
    Node(const Node<Type>& node);
    
    Node(NodeType nt, const std::vector<Type> &data);

    Node(NodeType nt, int length);

    //FIXME: not the right way to delete. object maybe stack allocated
    void deleteEmptyNode();

    //!sets the data at this node to zero
    void setZero();

    void setVal(Type val);

    //!set random data
    void setRandom();


    inline NodeType getType() { return _nt; }
    // const std::vector<Type> &data();
    // //Type* operator()(){return _data;}
    // //int getLength(){return _dataLength;}
    // size_t size();
    
    // void deleteData();   
    // void setType(NodeType nt);
    // void setData(const std::vector<Type>& data);
    // void setNode(NodeType nt, const std::vector<Type>& data);
    // Node<Type>& operator=(Node<Type>& node);
    
    inline const std::vector<Type>& getData() const { return _data;}
    inline std::vector<Type>& getData() { return _data;}
    //Type* operator()(){return _data;}
    inline int getLength(){return _data.size();}
    
    inline void deleteData() { _data.clear(); }
    
    inline void setType(NodeType nt) {_nt = nt;}
    inline void setData(const std::vector<Type>& data) {_data = data;}
    inline void setNode(NodeType nt, const std::vector<Type>& data){_nt = nt; _data = data;}
    inline Node<Type>& operator=(Node<Type>& node) { swap(node); return *this; }
    
    inline Type innerLocal(){
	Type retValue = (Type)0;
	if(_nt != EMPTY)
	    for(auto d :  _data)
		retValue+=d*d;
	else
	    cerr<<"Error at inner Local"<<endl;
	return retValue;
    }
    

    
  private:
    NodeType _nt;
    std::vector<Type> _data;
  };


  template<typename Type>
  Node<Type>::Node() 
    : _nt(EMPTY) {}
    
  template<typename Type>
  Node<Type>::Node(const Node<Type>& node) 
    : _nt(node.nt),
      _data(node.data) {}
    
  template<typename Type>
  Node<Type>::Node(NodeType nt, const std::vector<Type> &data)
    : _nt(nt),
      _data(data) {}
  
  /*!creates a node and allocates memory for data
    also initialized the data to 0.0*/
  template<typename Type>
  Node<Type>::Node(NodeType nt, int length)
    : _nt(nt),
      _data(length) {}

  //not the right way to delete. object maybe stack allocated
  template<typename Type>
  void
  Node<Type>::deleteEmptyNode() {
    if(_nt == EMPTY)
      delete this;
  }

  //!sets the data at this node to zero
  template<typename Type>
  void
  Node<Type>::setZero() {
    setVal(Type(0));
  }

  template<typename Type>
  void
  Node<Type>::setVal(Type val){
    if(_data.empty())
	std::cerr<<"Writing to NULL data in setZero"<<std::endl;    
    std::fill(_data.begin(), _data.end(), val);
  }

  int RandomNumber () { return (std::rand()%100); }

  //!set random data
  template<typename Type>
  void
  Node<Type>::setRandom() {
    std::generate(_data.begin(), _data.end(), RandomNumber);
    if(_data.empty())
      std::cerr<<"Writing to NULL data in SetRandom"<<std::endl;
  }
    
};
