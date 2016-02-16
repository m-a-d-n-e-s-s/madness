#ifndef __fuset_CcOp_h__
#define __fuset_CcOp_h__

#include "PrimitiveOp.h"

namespace fuset {

  template<typename Type>
  class CcOp : public PrimitiveOp<Type> {
    typedef KaryTree<Type> KTREE;
    typedef Node<Type> KNODE;
  private:
    //!Points to operand trees
    const KTREE *_i1;

    int _numChild;
    //!Points to operand nodes of the tree
    KNODE *_t1, *_t2;

  public:
    CcOp(string name, KTREE *output, const KTREE *i1);
    
    void compute(int n, TRANS l);

    bool isDone(int n, TRANS l) const;

    bool isPre() const { return false; }

    bool notEmpty(map<int,bool>& emptyMap) const{
	//cout<<"Not Empty i1 : "<<emptyMap[_i1->_treeID]<<endl;
	return  emptyMap[_i1->_treeID];
    }

	
  };

  template<typename Type>
      CcOp<Type>::CcOp(string name, KTREE *output, const KTREE *i1)
      : PrimitiveOp<Type>(name, output,false),
      _i1(i1),
      _numChild(output->getNumChild())
  {
      //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
      this->_OpID = output->_treeID;
      this->_dInfoVec.push_back(DependencyInfo<Type>(i1,true,false,false,false));
      this->_dInfoVec.push_back(DependencyInfo<Type>(output,false,false,true,false));


  }

  template<typename Type>
  void CcOp<Type>::compute(int n, TRANS l) {
    bool isE1;
    _t1 = _i1->getNode(n,l,isE1);
    KNODE* rNode = NULL;
    
    if(!isE1){
      rNode = this->_result->createNode(n,l);
      rNode->setZero();

      for(int i =0; i< _t1->getLength(); i++)
        rNode->getData()[i] += _t1->getData()[i];
    } else {
      std::cerr<<"Not supposed to be empty "<<n<<", "<<l<<std::endl;
      exit(0);
    }

    if(!isDone(n,l)){

      int n_child = n+1;

      for(int childID = 0; childID<_numChild; childID++){

        TRANS l_child = this->_result->childTranslation(l,childID);

        KNODE* temp = this->_result->getNode(n_child,l_child,isE1);

        if(isE1){
          std::cerr<<"Not supposed to be empty"<<n<<", "<<l<<std::endl;
          exit(0);
        }
	    
        for(int i =0; i< _t1->getLength(); i++)
          rNode->getData()[i] += temp->getData()[i];
      }
	
    }else{
      rNode->setType(LEAF);
    }
  }

  template<typename Type>
  bool CcOp<Type>::isDone(int n, TRANS l) const {
      //cout<<"Is done "<<n<<", "<<l<<endl;
      bool isE1;
    Node<Type>* temp = _i1->getNode(n,l,isE1);
    if(temp->getType() == LEAF)
      return true;	
    return false;
  }

}; /*fuset*/

#endif /* __fuset_CcOp_h__ */

