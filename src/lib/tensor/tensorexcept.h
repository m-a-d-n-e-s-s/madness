#ifndef TENSOREXCEPT_H
#define TENSOREXCEPT_H

/// \file tensorexcept.h
/// \brief Declares and implements TensorException

namespace madness {
/// Tensor is intended to throws only TensorExceptions
class TensorException {
public:
  const char* msg;
  const char* assertion;
  const int value;
  const BaseTensor* t;
  const int line;
  const char *function;
  const char *filename;

  // Capturing the line/function/filename info is best done with the macros below
  TensorException(const char* s, const char *a, int err, const BaseTensor* tp,
		  int lin, const char *func, const char *file) 
    : msg(s)
    , assertion(a)
    , value(err)
    , t(tp)
    , line(lin)
    , function(func)
    , filename(file) 
    {};
};

// implemented in tensor.cc
std::ostream& operator <<(std::ostream& out, const TensorException& e);

#define TENSOR_EXCEPTION(msg,value,t) \
throw TensorException(msg,0,value,t,__LINE__,__FUNCTION__,__FILE__)

#define TENSOR_ASSERT(condition,msg,value,t) \
do {if (!(condition)) \
    throw TensorException(msg,#condition,value,t,__LINE__,__FUNCTION__,__FILE__); \
   } while (0)

}
  
#endif
