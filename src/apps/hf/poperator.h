#include "mra/operator.h"

#define WST_PI 3.1415926535897932384

namespace madness
{
  const double acut1e_6 = 0.24; //0.6450626287524907;
  /*    template <typename Q, int NDIM>
   Q laplacian(const Q& f) {
   Q lapf = diff(diff(f,0),0);
   for (int i=1; i<NDIM; ++i) lapf += diff(diff(f,i),i);
   return lapf;
   }
   */
  template <typename Q> class PBSHFunctor
  {
  public:
    int kmax;
    double coeff, expnt;
    PBSHFunctor(int kmax, double coeff, double expnt) :
      kmax(kmax), coeff(coeff), expnt(expnt)
    {
    }
    Q operator()(double x) const
    { //x \in [-1,1] as checked
      Q sx0 = 0;
      for (int k=-kmax; k<=kmax; ++k)
        sx0 += exp(-expnt*(x-k)*(x-k));
      return sx0*coeff;
    }
  };

  double foo(const Tensor<double>& coeff, const Tensor<double>& expnt, double r)
  {
    double rval = 0.0;
    for (int i = 0; i < coeff.dim[0]; i++)
    {
      rval += coeff[i] * exp(-expnt[i]*r*r);
    }
    return rval;
  }

  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM> PBSHOperator(World& world, double mu, long k,
      double lo, double eps)
  {
    const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
    //    Tensor<double> cell_width(FunctionDefaults<NDIM>::cell(_, 1)-FunctionDefaults<
    //        NDIM>::cell(_, 0));
    //double hi = sqrt(double(NDIM))*cell_width.normf(); // Diagonal width of cell
    //const double pi = 3.14159265358979323846264338328;
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(mu, lo, 10.0, eps, &coeff, &expnt, true); //eps /(4.0*pi)
    //
    if (mu == 0.0)
    {
      // WSTHORNTON
      if (world.rank() == 0)
        printf("MU IS ZERO!!\n\n");
      coeff.scale(4.0 * WST_PI);
    }

    // WSTHORNTON
    double pt1 = 1e-1;
    double pt2 = 1e-2;
    double pt3 = 1e-3;
    double pt4 = 2;
    double f1 = foo(coeff, expnt, pt1);
    double f2 = foo(coeff, expnt, pt2);
    double f3 = foo(coeff, expnt, pt3);
    double f4 = foo(coeff, expnt, pt4);
    if (world.rank() == 0)
      printf("mu = %.3e\t\tlo = %.3e\n\n", mu, lo);
    if (world.rank() == 0)
      printf("pt1 = %.5e\t\tf1 = %.5e\n", pt1, f1);
    if (world.rank() == 0)
      printf("pt2 = %.5e\t\tf2 = %.5e\n", pt2, f2);
    if (world.rank() == 0)
      printf("pt3 = %.5e\t\tf3 = %.5e\n", pt3, f3);
    if (world.rank() == 0)
      printf("pt4 = %.5e\t\tf4 = %.5e\n", pt4, f4);

    int i;
    for (i = 0; i < coeff.dim[0]; ++i)
    {
      if (expnt[i] <= acut1e_6)
        break;
      coeff[i] = pow(coeff[i], 1.0 / double(NDIM));
      // WSTHORNTON
      if (world.rank() == 0)
        printf("expt[%d] = %.8f\t\tcoeff[%d] = %.8f\n", i, expnt[i], i, coeff[i]);
    }

    std::vector< SharedPtr< Convolution1D<Q> > > ops(i);
    double dum = log(2.0 / eps);
    for (int j = 0; j < i; ++j)
    {
      int kmax = ceil(sqrt(dum / expnt[j]));
      if (world.rank() == 0)
        printf("j = %d\t\tkmax = %d\n\n", j, kmax);
      ops[j] = SharedPtr<Convolution1D<Q> > (new GenericConvolution1D<Q,
          PBSHFunctor<Q> > (k, PBSHFunctor<Q> (kmax, coeff[j], expnt[j])));
    }
    printf("op effective rank: %d, max shift %d\n", i, (int) ceil(sqrt(dum
        / expnt[i - 1])));
    return SeparatedConvolution<Q, NDIM> (world, k, ops);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename Q, int NDIM>
  class wstPeriodicCoulombOp
  {
  public:
    wstPeriodicCoulombOp(World& world, double mu, long k, double lo, double eps)
    {
      // bsh_fit generates representation for 1/4Pir but we want 1/r
      // so have to scale eps by 1/4Pi
      Tensor<double> coeff, expnt;
      //if (mu==0) eps /= 4.0*pi;
      bsh_fit(mu, lo, 10.0, eps, &coeff, &expnt, true); //eps /(4.0*pi)
      //
      if (mu==0.0) coeff.scale(4.0*WST_PI);

      // WSTHORNTON
      // test that it's really 1/r
      double pt1 = 1e-1;
      double pt2 = 1e-2;
      double pt3 = 1e-3;
      double pt4 = 2;
      double f1 = foo(coeff, expnt, pt1);
      double f2 = foo(coeff, expnt, pt2);
      double f3 = foo(coeff, expnt, pt3);
      double f4 = foo(coeff, expnt, pt4);
      if (world.rank() == 0) printf("mu = %.3e\t\tlo = %.3e\n\n", mu, lo);
      if (world.rank() == 0) printf("pt1 = %.5e\t\tf1 = %.5e\n", pt1, f1);
      if (world.rank() == 0) printf("pt2 = %.5e\t\tf2 = %.5e\n", pt2, f2);
      if (world.rank() == 0) printf("pt3 = %.5e\t\tf3 = %.5e\n", pt3, f3);
      if (world.rank() == 0) printf("pt4 = %.5e\t\tf4 = %.5e\n", pt4, f4);

      // Scale coefficients according to the dimensionality
      int i, min_i;
      for (i=0; i < coeff.dim[0] && expnt[i] <= acut1e_6; ++i)
      {
        if (expnt[i] >= 60000) min_i = i;
        coeff[i]=pow(coeff[i], 1.0/double(NDIM));
//        if (world.rank() == 0) printf("expt[%d] = %.8f\t\tcoeff[%d] = %.8f\n", i, expnt[i], i, coeff[i]);
      }


      std::vector< SharedPtr< Convolution1D<Q> > > ops(i);
      double dum = log(2.0/eps);
      for (int j=0; j<i; ++j)
      {
        // kmax is the maximum number of lattice translations in a given dimension
        int kmax = ceil(sqrt(dum/expnt[j]));
        if (world.rank() == 0) printf("j = %d\t\tkmax = %d\n\n", j, kmax);
        ops[j]
            = SharedPtr< Convolution1D<Q> >(new GenericConvolution1D< Q,PBSHFunctor<Q> >(k,PBSHFunctor<Q>(kmax, coeff[j], expnt[j])));
      }
      printf("op effective rank: %d, max shift %d\n", i, (int)ceil(sqrt(dum
          /expnt[i-1])));
      _sepconv = new SeparatedConvolution<Q, NDIM>(world, k, ops);
    }
    //***************************************************************************

    //***************************************************************************
    ~wstPeriodicCoulombOp()
    {
      delete _sepconv;
    }
    //***************************************************************************

    //***************************************************************************
    void apply() {}
    //***************************************************************************

  private:
    //***************************************************************************
    SeparatedConvolution<Q, NDIM>* _sepconv;
    //***************************************************************************
  };
};
