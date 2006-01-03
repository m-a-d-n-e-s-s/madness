#include "legendre.h"

using namespace madness;
using namespace std;

double f(double x) {
  return 1.0/(1.0+x*x);
}

int main() {

    double x[100], w[100], err[100];

    for(int n=1; n<=32; n++) {
       if (!gauss_legendre(n, -1.0, 1.0, x, w)) throw "bad stuff";
       double sum = 0.0;
       for (int i=0; i<n; i++) sum += f(x[i])*w[i];
       err[n] = sum-1.5707963267948966192;
    }

    for (int n=2; n<=32; n++) 
       cout << n << " " << err[n] << "   " << err[n-1]/err[n] << endl;

    return 0;
}



