#include <iostream>
#include <cstdio>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "tensor.h"
using madness::Tensor;


int main() {

    int k = 10;
    int twok = 2*k;
    double ops = 2*3*twok*twok*twok*twok;
    long times = 10000;
    double million = 1e6;
    double used;
    double mops;
    double start;

    Tensor<double> x = Tensor<double>(2*k,2*k,2*k);
    Tensor<double> r = Tensor<double>(2*k,2*k,2*k);
    Tensor<double> w = Tensor<double>(2*k,2*k,2*k);
    Tensor<double> c = Tensor<double>(2*k,2*k);

    start = std::clock();
    for (long i=0; i<times; i++) {
        r = transform(x,c);
    }
    used = (double)(std::clock()-start)/(double)CLOCKS_PER_SEC;
    mops = ((double)times*ops)/(used*million);
    std::cout << "TRANSFORM MOPS=" << mops << "   "
    << used << "   "<< times << "   "<< ops << "   "<< million << std::endl;


    start = std::clock();
    for (long i=0; i<times; i++) {
        fast_transform(x,c,r,w);
    }
    used = (double)(std::clock()-start)/(double)CLOCKS_PER_SEC;
    mops = ((double)times*ops)/(used*million);
    std::cout << "TRANSFORM MOPS=" << mops << "   "
    << used << "   "<< times << "   "<< ops << "   "<< million << std::endl;

    start = std::clock();
    for (long i=0; i<times; i++) {
        r = transform3d(x,c);
    }
    used = (double)(std::clock()-start)/(double)CLOCKS_PER_SEC;
    mops = ((double)times*ops)/(used*million);
    std::cout << "TRANSFORM MOPS=" << mops << "   "
    << used << "   "<< times << "   "<< ops << "   "<< million << std::endl;

}
