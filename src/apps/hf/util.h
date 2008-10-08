#ifndef UTIL_H_
#define UTIL_H_

#include <mra/mra.h>
#include <world/world.h>

namespace madness {
//  void printfunc(const World& world, Function<double,3> f, int npts)
//  {
//    Tensor<double> LL = FunctionDefaults<3>::get_cell_width();
//    double L = LL[0];
//    double bstep = L / npts;
//    f.reconstruct();
//    for (int i = 0; i <= npts; i++)
//    {
//      Vector<double,3> p(-L/2 + i * bstep);
//      if (world.rank() == 0) printf("%.2f\t\t%.8f\n", p[0], f(p));
//    }
//    if (world.rank() == 0) printf("\n");
//  }

//  void printfunc(const World& world, Function<double,3> f1, Function<double,3> f2, int npts)
//  {
//    Tensor<double> LL = FunctionDefaults<3>::get_cell_width();
//    double L = LL[0];
//    double bstep = L / npts;
//    f1.reconstruct();
//    f2.reconstruct();
//    for (int i = 0; i <= npts; i++)
//    {
//      Vector<double,3> p(-L/2 + i * bstep);
//      if (world.rank() == 0) printf("%.2f\t\t%.8f\t%.8f\n", p[0], f1(p), f2(p));
//    }
//    if (world.rank() == 0) printf("\n");
//  }
}
//
//#include <mra/mra.h>
//#include <world/world.h>
//#include <vector>
//
//namespace madness
//{
//  class OnesFunctor :
//  public FunctionFunctorInterface<double,3>
//  {
//  private:
//
//  public:
//    //*************************************************************************
//    OnesFunctor()
//    {
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    virtual ~OnesFunctor() {}
//    //*************************************************************************
//
//    //*************************************************************************
//    double operator()(const coordT& x) const
//    {
//      return 1.0;
//    }
//    //*************************************************************************
//  };
//
//  //***************************************************************************
//  class ZerosFunctor :
//  public FunctionFunctorInterface<double,3>
//  {
//  private:
//
//  public:
//    //*************************************************************************
//    ZerosFunctor()
//    {
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    virtual ~ZerosFunctor() {}
//    //*************************************************************************
//
//    //*************************************************************************
//    double operator()(const coordT& x) const
//    {
//      return 0.0;
//    }
//    //*************************************************************************
//  };
//  //***************************************************************************
//}

#endif
