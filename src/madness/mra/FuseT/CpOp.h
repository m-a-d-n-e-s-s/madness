#ifndef __fuset_CpOp_h__
#define __fuset_CpOp_h__

#include "PrimitiveOp.h"

namespace fuset {
  template<typename Type>
  class CpOp : public PrimitiveOp<Type> {
    typedef KaryTree<Type> KTREE;
    typedef Node<Type> KNODE;
  public:
    CpOp(string opName, KTREE *output, const KTREE *i1);
    
    void compute(int n, TRANS l);

    bool isDone(int n, TRANS l) const;

    bool isPre() const { return true; }

  private:
    //!Points to operand trees
    const KTREE *_i1;
    
    //!Points to operand nodes of the tree
    KNODE *_t1, *_t2;
  };

  template<typename Type>
      CpOp<Type>::CpOp(string opName, KTREE *output, const KTREE *i1)
      : PrimitiveOp<Type>(opName, output, false),
      _i1(i1)
  {
      //dependnecy Info PSI, ALPHA, DELTA,SIGMA, ID
      this->_OpID = output->_treeID;
      
      this->_dInfoVec.push_back(DependencyInfo<Type>(i1,true,false,false,false));
      this->_dInfoVec.push_back(DependencyInfo<Type>(output,false,true,false,false));
      
  }
    
  template<typename Type>
  void CpOp<Type>::compute(int n, TRANS l) {
    bool isE1;
    _t1 = _i1->getNode(n,l,isE1);

    if(!isE1){

      KNODE* rNode = this->_result->createNode(n,l);
      rNode->setZero();

      if(n !=0){
        int n_parent = n-1;
        TRANS l_parent = this->_result->parentTranslation(l);
        bool isEmpty;
        _t2 = this->_result->getNode(n_parent, l_parent, isEmpty);
        if(!isEmpty) {
          rNode->setData(_t2->getData());
          /* for(int i =0; i< _t2->getLength(); i++) */
          /*   rNode->getData()[i] = _t2->getData()[i]; */
        }
        else{
	    
	    std::cerr<<"Not supposed to be empty (n,l) :"<<n<<", "<<l<<std::endl;
	    
	  exit(0);
	}
      }
      
      //rNode->setData(_t1->getData());
      for(int i =0; i< _t1->getLength(); i++) 
	  rNode->getData()[i] += _t1->getData()[i]; 
      
      if(isDone(n,l))
	  rNode->setType(LEAF);
    }	           
  }

  template<typename Type>
  bool CpOp<Type>::isDone(int n, TRANS l) const {
    bool isE1;
    Node<Type>* temp = _i1->getNode(n,l,isE1);
    if(temp->getType() != INTERIOR)
      return true;	
    return false;
  }

}; /*fuset*/

#endif /* __fuset_CpOp_h__ */
