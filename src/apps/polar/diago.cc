#include <iostream>


#include <tensor/tensor.h>

#include <linalg/tensor_lapack.h>
#include <linalg/clapack.h>

using namespace madness;
using namespace std;


#ifndef SCALAR
// #define SCALAR std::complex<float>
#define SCALAR double
#endif
#ifndef REALX
// #define SCALAR std::complex<float>
#define REALX real8
#endif

#ifndef STATIC
#define STATIC static
#endif

typedef double RealX;

namespace madness {
    template <typename T>
    void test_syev1(int n) {
        Tensor<T> a(n,n), V;
        Tensor< typename Tensor<T>::scalar_type > e;
        a.fillrandom();
//        a += madness::my_conj_transpose(a);
        a += conj_transpose(a);
//        a += transpose(a); //fool
        syev(a,V,e);
    cout << "a syev \n" << a; 
    cout << "v syev \n" << V; 
    cout << "e syev \n" << e; 
        double err = 0.0;
        for (int i=0; i<n; ++i) {
            cout << "hola \n";
          err = max(err,(double) (inner(a,V(_,i)) - V(_,i)*e(i)).normf());
        }
        cout <<  err << endl;
    }
}

int main(int argc, char* argv[]) {

   long n = 3;   


   test_syev1<std::complex<double>>(n);
   
   return 0;
}
