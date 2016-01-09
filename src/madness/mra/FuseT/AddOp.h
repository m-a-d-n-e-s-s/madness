#ifndef __fuset_AddOp_h__
#define __fuset_AddOp_h__

#include "PrimitiveOp.h"

namespace fuset {

  template<typename Type>
  class AddOp : public PrimitiveOp<Type>{
    typedef KaryTree<Type> KTREE;
    typedef Node<Type> KNODE;

  public:
    /*!Adds trees i1 and i2 and creates result tree*/
    AddOp(string name, KTREE *result, const KTREE *i1, const KTREE *i2);    
    
    /*!Adds two nodes from i1 and i2. 
      If one exists and another is empty, the it just replaces the values at the result node with data from the node that exists*/
    void compute(int n, TRANS l);

    /*!Checks to see if the add operation needs to recurse further down*/
    bool isDone(int n, TRANS l) const;

    bool isPre() const { return true; }

    bool notEmpty(map<int,bool> &emptyMap) const{
	bool retValue = (emptyMap[_i1->_treeID] || emptyMap[_i2->_treeID]);
	//cout<<"Add Not Empty : "<<retValue<<endl;
	return retValue;
    }
  private:
    //!Points to operand trees
    const KTREE *_i1, *_i2;
    
    //!Points to operand nodes of the tree
    KNODE *_t1, *_t2;
  };

  template<typename Type>
      AddOp<Type>::AddOp(string name, KTREE *output, const KTREE *i1, const KTREE *i2)
      : PrimitiveOp<Type>(name, output,false),
      _i1(i1),
      _i2(i2)
  {   
      //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
      this->_OpID = output->_treeID;
      this->_dInfoVec.push_back(DependencyInfo<Type>(i1,true,false,false,false));
      this->_dInfoVec.push_back(DependencyInfo<Type>(i2,true,false,false,false));
  }
    
  template<typename Type>
  void 
  AddOp<Type>::compute(int n, TRANS l) {
    //    std::cout<<"Entering n : "<<n<<" l : "<<l<<std::endl;
    bool isE1, isE2;
    _t1 = _i1->getNode(n,l,isE1);
    _t2 = _i2->getNode(n,l,isE2);
    //    std::cout<<"Before if n : "<<n<<" l : "<<l<<std::endl;
    if(!isE1 || !isE2){
      KNODE* t3 = this->_result->createNode(n,l);
      t3->setZero();
      if(!isE1)
        for(int i =0; i< _t1->getLength(); i++)
          t3->getData()[i] += _t1->getData()[i];
      //	std::cout<<"Before isE2 n : "<<n<<" l : "<<l<<std::endl;
      if(!isE2)
        for(int i =0; i< _t2->getLength(); i++)
          t3->getData()[i] += _t2->getData()[i];
      //	std::cout<<"Before isDone n : "<<n<<" l : "<<l<<std::endl;
      if(isDone(n,l))
        t3->setType(LEAF);
      //	std::cout<<"After isDone n : "<<n<<" l : "<<l<<std::endl;
    }
    //    std::cout<<"Exiting n : "<<n<<" l : "<<l<<std::endl;
    _t1->deleteEmptyNode();
    _t2->deleteEmptyNode();
  }

  template<typename Type>
  bool
  AddOp<Type>::isDone(int n, TRANS l) const {
    if(_t1->getType() !=INTERIOR
       && _t2->getType() != INTERIOR)
      return true;	
    return false;
  }

}; /*fuset*/

#endif /* __fuset_AddOp_h__ */
