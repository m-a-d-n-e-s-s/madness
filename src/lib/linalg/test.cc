
/// \file ../linalg/test.cc
/// \brief Test code for LAPACK, Tensor+LAPACK, etc.

#include <iostream>
#include <linalg/tensor_lapack.h>

using namespace madness;

int
main(int argc, char* argv[])
{
  bool testok = test_tensor_lapack();
  std::cout << "Test " << (testok ? "passed" : "did not pass") << std::endl;
  return 0;
}

