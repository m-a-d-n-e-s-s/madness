#include <iostream>
#include <fstream>
#include <string>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/print.h>



using namespace madness;

int main()
{

  std::ifstream file;
  file.open("integrals.dat");

  // create needed variables  
  int nbf;
  int nocc;
  double enrep;
  double ** S;//Overlap Matrix
  double ** KE;//Kinetic Energy 
  double ** mux;//
  double ** muy;// 
  double ** muz;//
  double ** Electron;
  

  return 0; 
}

