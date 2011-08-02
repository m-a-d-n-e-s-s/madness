//// \file mapreduce.h
//// \brief Provides MapReduce-like functor base abstract classes for MADNESS

//// \defgroup MapReduce
//// @{
namespace madness{

//// A predicate functor
template <typename keyT, typename valueT>
class Pred{
  public:
    Pred(){}

    //// pure virtual predicate function
    //// @param[in] k the key of an entry in a WorldContainer
    //// @param[in] v the value of an entry in a WorldContainer
    //// @return the value of the predicate applied to the entry 
    virtual bool call(const keyT& k, const valueT& v) = 0;

    virtual ~Pred(){}
};

template <typename keyI, typename keyO, typename valueI, typename valueO>
class MapFunctor{
  public:
    MapFunctor(){}

    virtual std::pair<keyO,valueO> call(const keyI &k, const valueI &vIn) = 0;
//      print("BASE");
//    }

    virtual ~MapFunctor(){}
};

template <typename keyI, typename keyO, typename valueI, typename valueO>
class MapToManyFunctor{
  public:
    MapToManyFunctor(){}

    virtual std::vector<std::pair<keyO,valueO>> call(const keyI &k, const valueI &vIn) = 0;
//      print("BASE");
//    }

    virtual ~MapToManyFunctor(){}
};

template<typename keyO, typename keyI, typename valueI>
class JoinOp{
  public:
    JoinOp(){}

    virtual keyO join(const keyI& k, const valueI& val) = 0;

    virtual ~JoinOp(){}
};

template<typename keyI, typename valueI, typename keyO, typename valueO>
class UpdateOp{
  public:
    UpdateOp(){}

    virtual valueO call(const keyI& kI, const valueI& vI, const keyO& kO, const valueO* vO, const bool& present) = 0;
   
    virtual ~UpdateOp(){}
};

template<typename keyI, typename valueI, typename reduceO>
class LocalReduce{
  public:
    LocalReduce(){}

    virtual void merge(const keyI& kI, const valueI& vI, reduceO * redInOut) = 0;

    virtual ~LocalReduce(){}
};

template<typename reduceO>
class ParallelReduce{
  public:
    ParallelReduce(){}

    virtual void merge(const reduceO * redIn, reduceO * redInOut) = 0;

    virtual ~ParallelReduce(){}
}; 

}

//// @}
