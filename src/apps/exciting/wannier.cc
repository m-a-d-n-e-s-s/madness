#include <mra/mra.h>
#include <iostream>
#include <cmath>

extern "C" void readinput_();
extern "C" void wann_init1_();
extern "C" void wann_unk_(int* n,double* vpl,double* vrc,double* val);

//template<typename T, int NDIM>
//class Wannier: public FunctionFunctorInterface<T, NDIM>
//{
//public:
//  typedef Vector<double, NDIM> coordT;
//  const int _n;
//  const coordT _center;
//  const Vector<int,3> _nk;
//
//  Wannier(const int& n, const coordT& center, const Vector<int,3>& nk) :
//    _n(n), _center(center), _nk(nk) {}
//
//  T operator()(const coordT& x) const
//  {
//    // initialize value to zero
//    double[2] val;
//    val[0] = 0.0; val[1] = 0.0;
//    // number of k-points
//    int nkpt = _nk[0] * _nk[1] * _nkp[2];
//    //compute \sum_{k}u_{nk}(r)
//    for (int j1 = 0; j1 < _nk[0]; j1++)
//    {
//      for (int j2 = 0; j2 < _nk[1]; j2++)
//      {
//        for (int j3 = 0; j3 < _nk[2]; j3++)
//        {
//          //k-point in lattice coordinates
//          vkl[0] = 1.0 * j1 / _nk[0];
//          vkl[1] = 1.0 * j2 / _nk[1];
//          vkl[2] = 1.0 * j3 / _nk[2];
//          //get value of u_{nk}(r)
//          wann_unk_(&n, vkl, vr, val0);
//          val[0] += val0[0];
//          val[1] += val0[1];
//        }
//      }
//    }
//    val[0]=val[0]/nkpt;
//    val[1]=val[1]/nkpt;
//
//    return sqrt(val[0]*val[0] + val[1]*val[1]);
//  }
//};



int main(int argn, char** argv) {
    //k-mesh division
//    Vector<int,3> nk(4);
//    //cener point of the box
//    Vector<double,3> center(0.0);
//    //index of Wannier function
//    int n=3;

    readinput_();
    wann_init1_();

    return 0;
}
