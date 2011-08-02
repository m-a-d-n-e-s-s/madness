#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <world/world.h>
#include <world/worlddc.h>

using namespace std;
using namespace madness;

/*
template <typename keyI, typename keyO, typename valueI, typename valueO>
class SqF : public MapFunctor<keyI, keyO, valueI, valueO>{
  public:

    SqF() : MapFunctor<keyI, keyO, valueI, valueO>() {}

    virtual std::pair<keyO,valueO> call(const keyI &k, const valueI &vIn){
      print("<", k," ",vIn, ">");
      return std::pair<keyO,valueO>(k, vIn*vIn);
    }

    virtual ~SqF(){}
};
*/

class SquareMapFunctor : public MapFunctor<int, int, int, int>{
  public:

    //SqF() : MapFunctor<keyI, keyO, valueI, valueO>() {}

    virtual std::pair<int,int> call(const int &k, const int &vIn){
      //print("<", k," ",vIn, ">");
      return std::pair<int,int>(4*k, vIn*vIn);
    }

    //virtual ~SqF(){}
};

class SqThF : public MapToManyFunctor<int, int, int, int>{
  public:

    //SqF() : MapFunctor<keyI, keyO, valueI, valueO>() {}

    virtual std::vector<std::pair<int,int>> call(const int &k, const int &vIn){
      //print("<", k," ",vIn, ">");
      std::vector<std::pair<int, int>> out;
      out.push_back(std::pair<int,int>(4*k, vIn*vIn));
      out.push_back(std::pair<int,int>(6*k, vIn*vIn*vIn));
      return out;
    }

    //virtual ~SqF(){}
};

class IdJoin : public JoinOp<int,int,int>{
  public:

    virtual int join (const int& k, const int& vIn){
      return k;
    }

};

class SqJoin : public JoinOp<int,int,int>{
  public:

    virtual int join (const int& k, const int& vIn){
      return k*k;
    }

};

class NilJoin : public JoinOp<int,int,int>{
  public:

    virtual int join (const int& k, const int& vIn){
      return 0;
    }

};

class IntUpdateOp : public UpdateOp<int, int, int, int>{
  public:

    IntUpdateOp() : UpdateOp<int, int, int, int>() {}

    virtual int call(const int& kI, const int& vI, const int& kO, 
                      const int* vO, const bool& present){
      if (present){
        return *vO;
      }
      else{
        return vI;
      }
    }

    virtual ~IntUpdateOp(){}
};

class SumUpdateOp : public UpdateOp<int, int, int, int>{
  public:

    SumUpdateOp() : UpdateOp<int, int, int, int>() {}

    virtual int call(const int& kI, const int& vI, const int& kO, 
                      const int* vO, const bool& present){
      if (present){
        return *vO + vI;
      }
      else{
        return vI;
      }
    }

    virtual ~SumUpdateOp(){}
};

class SumLocalReduce : public LocalReduce<int, int, int>{
  public:

    void merge(const int& kI, const int& vI, int * redInOut){
      *redInOut = *redInOut + kI + vI;
    }

};

class SumParallelReduce : public ParallelReduce<int>{
  public:

    void merge(const int * redIn, int * redInOut){
      *redInOut = *redInOut + *redIn;
    }

};

int main(int argc, char** argv){
  initialize(argc, argv);
  madness::World world(MPI::COMM_WORLD);

  WorldContainer<int, int> dfc(world);

  dfc.replace(0,0);
  dfc.replace(1,1);
  dfc.replace(2,2);
  dfc.replace(3,3);
  dfc.replace(4,4);
  dfc.replace(5,5);
  dfc.replace(6,6);
  dfc.replace(7,7);
  dfc.fence();

  cout << "Original container:" << endl;
  cout << dfc;

  cout << "--" ;
  //SqF<int,int,int,int> * sqf = new SqF<int,int,int,int>;

  //sqf->call(1,1);

  SquareMapFunctor * sqf = new SquareMapFunctor;

  WorldContainer<int,int> dfc1 = dfc.map(sqf, dfc.get_pmap());
  dfc1.fence();

  cout << "container 1 (map)" << endl;
  cout << dfc1;

  WorldContainer<int,int> dfc2 = dfc1.combine(dfc);
  dfc2.fence();

  cout << "container 2 (combine)" << endl;
  cout << dfc2;

  WorldContainer<int,int> dfc3 = dfc.combine(dfc1);
  dfc3.fence();

  cout << "container 3 (combine)" << endl;
  cout << dfc3;

  IdJoin * idj = new IdJoin;  

  WorldContainer<int,std::pair<std::pair<int,int>,std::pair<int,bool>>> dfc4 = dfc3.join2(dfc1, idj);
  dfc4.fence();

  cout << "container 4" << endl;
  cout << dfc4;

  WorldContainer<int,std::pair<std::pair<int,int>,std::pair<int,bool>>> dfc5 = dfc1.join2<int,int>(dfc3, idj);
  dfc5.fence();

  cout << "container 5" << endl;
  cout << dfc5;

  SqThF * thF = new SqThF;

  WorldContainer<int,int> dfc6 = dfc.mapToMany<int,int>(thF, dfc.get_pmap());
  dfc6.fence();

  cout << "container 6" << endl;
  cout << dfc6;

  SqJoin * sj = new SqJoin;
  NilJoin * nlj = new NilJoin;

  IntUpdateOp * iuo = new IntUpdateOp;
  SumUpdateOp * suo = new SumUpdateOp;
  WorldContainer<int,int> dfc7 = dfc.mapUpdate(dfc, idj, iuo);
  dfc7.fence();

  cout << "container 7" << endl;
  cout << dfc7;

  WorldContainer<int,int> dfc8 = dfc.mapUpdate(dfc1, sj, iuo);
  dfc8.fence();

  cout << "container 8" << endl;
  cout << dfc8;

  WorldContainer<int,int> dfc9 = dfc.mapUpdate(dfc, nlj, suo);
  dfc9.fence();

  cout << "container 9" << endl;
  cout << dfc9;

  SumLocalReduce * sl = new SumLocalReduce;
  SumParallelReduce * pl = new SumParallelReduce;
  int * sum = new int;
  *sum = 0;
  dfc.reduce(sl, pl, sum);

  cout << "Key-value-sum: " << *sum << endl;

  finalize();
};

