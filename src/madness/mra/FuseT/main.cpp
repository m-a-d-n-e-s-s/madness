#include "KaryTree.h"
#include "CpOp.h"
#include "AddOp.h"
#include "CcOp.h"
#include "OpExecutor.h"
#include "FuseT.h"
#include "FusedExecutor.h"
#include "helper.h"
#define RESULTTREE 0
#define MAXDEPTH 12
#define NUMCHILD 8
#define DATALENGTH 1000
//template class KaryTree<double>;

using namespace fuset;
using namespace std;

typedef double dtype;

int main(int argc, char* argv[]) {
    int isFused = 0;
    if ( argc != 2 ){
	cout<<"Incorrect command line parameters"<<endl;
    	cout<<"Enter as : executable  isFused"<<endl;
	cout<<"Enter 1 for fused and 0 for unfused"<<endl;
	exit(0);
    }else {
	isFused=atoi(argv[1]);
    }
    bool printInfo = true;
    
  //The operands are maxDepth, numChild, dataLength
  //inputs
    KaryTree<dtype> t1("t1",MAXDEPTH,NUMCHILD,DATALENGTH), t2("t2",MAXDEPTH,NUMCHILD,DATALENGTH), t3("t3",MAXDEPTH,NUMCHILD,DATALENGTH), t4("t4",MAXDEPTH,NUMCHILD,DATALENGTH);

  /*!the seed for srand is the time(NULL), this adding the sleep
   * to generate different trees*/
  t1.generateRandom(SOME_NUM,1.0);
  t2.generateRandom(SOME_NUM,2.0);
  t3.generateRandom(SOME_NUM,1.5);
  t4.generateRandom(SOME_NUM,2.0);

  //results
  KaryTree<dtype> rCp0("rCp0",MAXDEPTH,NUMCHILD,DATALENGTH);
  KaryTree<dtype> rAdd0("rAdd0",MAXDEPTH,NUMCHILD,DATALENGTH);
  KaryTree<dtype> rCc0("rCc0",MAXDEPTH,NUMCHILD,DATALENGTH);
  KaryTree<dtype> rCp1("rCp1",MAXDEPTH,NUMCHILD,DATALENGTH);
  KaryTree<dtype> rAdd1("rAdd1",MAXDEPTH,NUMCHILD,DATALENGTH);
  KaryTree<dtype> rCc1("rCc1",MAXDEPTH,NUMCHILD,DATALENGTH);

  //initializing operators

  AddOp<dtype> opAdd0("opAdd0", &rAdd0, &t2, &t3);
  CpOp<dtype> opCp0("opCp0", &rCp0, &rAdd0);
  CcOp<dtype> opCc0("opCc0", &rCc0,&t4);
  CpOp<dtype> opCp1("opCp1", &rCp1, &t1);
  AddOp<dtype> opAdd1("opAdd1", &rAdd1, &rCc0, &t3);
  CcOp<dtype> opCc1("opCc1", &rCc1,&rAdd1);

  opAdd0.PrintDependencyInfo();  
  opCp0.PrintDependencyInfo();
  opCc0.PrintDependencyInfo();
  opCp1.PrintDependencyInfo();
  opAdd1.PrintDependencyInfo();  
  opCc1.PrintDependencyInfo();

  std::vector<PrimitiveOp<dtype>*> sequence;
  sequence.push_back(&opAdd0);
  sequence.push_back(&opCp0);
  sequence.push_back(&opCc0);
//  sequence.push_back(&opCp1);
//  sequence.push_back(&opAdd1);
//  sequence.push_back(&opCc1);


  double start,end;

  if(isFused){
      cout<<"***********************"<<endl;  
      cout<<"***********************"<<endl;  
      cout<<"Fused Tree Traversal "<<endl;
      cout<<"***********************"<<endl;
      cout<<"***********************"<<endl;
      
      FuseT<dtype> odag(sequence);
      odag.processSequence();

      if(printInfo){
	  odag.printOpsAndTrees();
	  odag.printValidSequences();
      }
      FusedOpSequence<dtype> fsequence = odag.getFusedOpSequence();
      FusedExecutor<dtype> fexecuter(&fsequence);
      start = rtclock();
      fexecuter.execute();
      end = rtclock();
  }else{
      cout<<"***********************"<<endl;  
      cout<<"***********************"<<endl;     
      cout<<"Unfused Tree Traversal "<<endl;
      cout<<"***********************"<<endl;
      cout<<"***********************"<<endl;  
    
      //initializing executor
      OpExecutor<dtype> opx;
      start = rtclock();
      for(auto op : sequence)
	  opx.execute(op, false);
      end =rtclock();
  }

  cout<<"Total Run time is "<<end-start<<endl;
  if(RESULTTREE){  
      std::cout<<std::endl<<"Result Tree AddOp0"<<std::endl;
      rAdd0.printTree(0,0);
      std::cout<<std::endl<<"Result Tree CpOp0"<<std::endl;
      rCp0.printTree(0,0);
      std::cout<<std::endl<<"Result Tree CcOp0"<<std::endl;
      rCc0.printTree(0,0);
      std::cout<<std::endl<<"Result Tree CpOp1"<<std::endl;
      rCp1.printTree(0,0);
      std::cout<<std::endl<<"Result Tree AddOp1"<<std::endl;
      rAdd1.printTree(0,0);
      std::cout<<std::endl<<"Result Tree CcOp1"<<std::endl;
      rCc1.printTree(0,0);
  }
  double sum = 0.0;
  for(auto a : sequence){
      double norm = a->_result->getNorm();
      cout<<"Norm is "<<norm<<endl;
      sum+=norm;
  }

  cout<<"Sum of Norm is "<<sum<<endl;




  return 0;
}
