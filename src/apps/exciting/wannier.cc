#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
#include <iostream>
#include <cmath>
#include <complex>

extern "C" void exciting_init_();
extern "C" void f_veff_(double* vrc,double* val);
extern "C" void f_rho_(double* vrc,double* val);
extern "C" void f_wann_(int* n,int* ispn,double* d0,int* itr,double* vrc,double* val);
extern "C" void f_lvec_(double(*)[3], double(*)[3], double(*)[3], double(*)[3]);
extern "C" void f_nkpts_(int*);
extern "C" void f_kpts_(int*, double*, double*, double*, double*);
extern "C" void f_ngridk_(int*, int*, int*);
extern "C" void f_wann_c_(int*, int*, int*, double*);
extern "C" void f_nwann_(int*);
extern "C" void f_wann_center_(int*,double*);

using namespace madness;

typedef SharedPtr< FunctionFunctorInterface< std::complex<double> ,3> > cfunctorT;
typedef FunctionFactory<std::complex<double>,3> cfactoryT;
typedef Function<std::complex<double>,3> cfunctionT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef FunctionFactory<double,3> factoryT;
typedef Function<double,3> functionT;
typedef Vector<double,3> coordT;
typedef Vector<double,3> vec3dT;

////***************************************************************************
//double abs(double x) {return x;}
////***************************************************************************

//***************************************************************************
double real(double x) {return x;}
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct abs_square_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim, t.dim);
    BINARY_OPTIMIZED_ITERATOR(Q, t, resultT, result, resultT d = abs(*_p0); *_p1 = d*d);
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> abs_square(const Function<Q,NDIM>& func)
{
  return unary_op(func, abs_square_op<Q,NDIM>());
}
//***************************************************************************

std::complex<double> func(int n, int t1, int t2, int t3, double xx, double yy, double zz)
{
  // initialize value to zero
  double val[2];
  val[0] = 0.0; val[1] = 0.0;
  double vr[3];
  vr[0] = xx; vr[1] = yy; vr[2] = zz;
  double d1=10.0;
  int ispn=1;
  int itr[3];
  itr[0]=t1; itr[1]=t2; itr[2]=t3;

  f_wann_(&n, &ispn, &d1, &itr[0], &vr[0], &val[0]);

  return std::complex<double>(val[0], val[1]);
}

template<typename T, int NDIM>
class Wannier: public FunctionFunctorInterface<T, NDIM>
{
public:
  typedef Vector<double, NDIM> coordT;
  int _n;
  coordT _center;
  Vector<int,3> _itr;

  Wannier(int n, Vector<int,3>& itr, const coordT& center) :
    _n(n), _center(center), _itr(itr) {}

  T operator()(const coordT& x) const
  {
    return func(_n, _itr[0], _itr[1], _itr[2], x[0], x[1], x[2]);
  }
};

class ComplexExpFunction: public FunctionFunctorInterface<std::complex<double>, 3>
{
public:
  double _kx;
  double _ky;
  double _kz;

  ComplexExpFunction(double kx, double ky, double kz)
  {
    _kx = kx;
    _ky = ky;
    _kz = kz;
  }

  std::complex<double> operator()(const coordT& r) const
  {
    std::complex<double> arg = std::complex<double>(0.0, _kx*r[0] + _ky*r[2] + _kz*r[2]);
    return arg;
  }
};

class ExcitingApp
{

public:

  struct wannierf
  {
    int n;
    int rindx;
    cfunctionT f;

    wannierf(int n, int rindx, const cfunctionT& g) :
      n(n), rindx(rindx), f(g) {}
  };

  // The world
  World& _world;

  // reciprocal lattice vectors in cartesian coords
  std::vector<vec3dT> _gvecs;
  std::vector<int> _gshells;
  // lattice vectors, lengths and shells
  std::vector<vec3dT> _lvecs;
  std::vector<double> _rlens;
  std::vector<int> _rshells;
  // matrices to convert from cartesian to
  // lattice coordinates
  Tensor<double> _avec;
  Tensor<double> _bvec;
  Tensor<double> _ainv;
  Tensor<double> _binv;
  // k-mesh grid size
  int _ngridk0;
  int _ngridk1;
  int _ngridk2;
  // k-points in lattice coordinates
  std::vector<vec3dT> _vkl;
  // k-points in cartesian coordinates
  std::vector<vec3dT> _vkc;
  // k-point weights
  std::vector<vec3dT> _wkpt;
  // wannier functions
  std::vector<wannierf> _wfs;

  ExcitingApp(World& world) : _world(world)
  {
    // initialize exciting
    exciting_init_();
    // read lattice from exciting
    read_lattice();
    // generate g-vectors
    generate_gvectors();
    // generate l-vectors
    generate_lvectors();
    // get k-points
    read_kpoints();
    // read wannier functions
//    read_wannier_functions();

//    // compute response
//    vec3dT qv;
//    qv[0] = 1.0; qv[1] = 0.0; qv[2] = 0.0;
//    response(qv);
  }

  vec3dT r3mv(const Tensor<double>& m, const vec3dT& v)
  {
    vec3dT r;
    r[0] = m(0,0)*v[0] + m(0,1)*v[1] + m(0,2)*v[2];
    r[1] = m(1,0)*v[0] + m(1,1)*v[1] + m(1,2)*v[2];
    r[2] = m(2,0)*v[0] + m(2,1)*v[1] + m(2,2)*v[2];
    return r;
  }

  int findkpt(vec3dT kvecl)
  {
    // is k-point in first BZ
    double eps = 1e-8;
    if ((kvecl[0] < (-eps)) || (kvecl[0] > (1.0 - eps))) return -2;
    if ((kvecl[1] < (-eps)) || (kvecl[1] > (1.0 - eps))) return -2;
    if ((kvecl[2] < (-eps)) || (kvecl[2] > (1.0 - eps))) return -2;

    for (unsigned int idx = 0; idx < _vkl.size(); idx++)
    {
      vec3dT tkvl = _vkl[idx];
      if ((abs(kvecl[0] - tkvl[0]) < eps) &&
          (abs(kvecl[1] - tkvl[1]) < eps) &&
          (abs(kvecl[2] - tkvl[2]) < eps))
      {
        return idx;
      }
    }
    return -1;
  }

//  Tensor<std::complex<double> > compute_A_coeffs(vec3dT qvec)
//  {
//    // Maximum number of r-shells, g-shells, number of wannier functions
//    int nrshells = 2;
//    int ngshells = 3;
//    int nwann = 0;
//    f_nwann_(&nwann);
//
//    // create return value
//    Tensor<std::complex<double> > acoeffs()
//
//    // create exponential functions
//    std::vector<cfunctionT> gfuns;
//    int igsh = 0;
//    for (int ig = 0; igsh < ngshells; ig++)
//    {
//      igsh = _gshells[ig];
//      if (igsh < ngshells)
//      {
//        vec3dT g = _gvecs[ig];
//        double qpg0 = g[0] + qvec[0];
//        double qpg1 = g[1] + qvec[1];
//        double qpg2 = g[2] + qvec[2];
//        cfunctionT gf = cfactoryT(_world).functor(cfunctorT(new ComplexExpFunction(qpg0, qpg1, qpg2)));
//        gfuns.push_back(gf);
//      }
//    }
//
//    // Loop through i wannier functions
//    for (unsigned int iwf = 0; iwf < _wfs.size(); iwf++)
//    {
//      // search for wavefunction jwf and with translation ir == 0
//      cfunctionT wannf_i;
//      bool ifound = false;
//      for (unsigned ifn = 0; ifn < _wfs.size() && !ifound; ifn++)
//      {
//        int wfn = _wfs[ifn].n;
//        unsigned int wfir = _wfs[ifn].rindx;
//        if ((iwf == wfn) && (wfir == 0))
//        {
//          wannf_i = _wfs[ifn].f;
//          ifound = true;
//        }
//      }
//      if (!ifound) MADNESS_EXCEPTION("Did not find iwf wavefunction!", 0);
//
//      // Loop through j wannier functions
//      for (unsigned int jwf = 0; jwf < iwf; jwf++)
//      {
//        // Loop over lattice translations less than nrshells
//        int irsh = 0;
//        for (int ir = 0; irsh < nrshells; ir++)
//        {
//          irsh = _rshells[ir];
//          if (irsh < nrshells)
//          {
//            // search for wavefunction jwf and with translation ir
//            cfunctionT wannf_j;
//            bool jfound = false;
//            for (unsigned ifn = 0; ifn < _wfs.size() && !jfound; ifn++)
//            {
//              int wfn = _wfs[ifn].n;
//              unsigned int wfir = _wfs[ifn].rindx;
//              if ((jwf == wfn) && (wfir == ir))
//              {
//                wannf_j = _wfs[ifn].f;
//                jfound = true;
//              }
//            }
//            if (!jfound) MADNESS_EXCEPTION("Did not find jwf wavefunction!", 0);
//          }
//        }
//
//        // Loop through exponential functions of q+G
//        for (unsigned int iqpg = 0; iqpg < gfuns.size(); iqpg++)
//        {
//          cfunctionT expqpg = gfuns[iqpg];
//          std::complex<double> t1 = inner(wannf_i,expqpg*wannf_j);
//
//        }
//      }
//    }
//
//  }

  void read_wannier_functions()
  {
    // get number of wannier functions
    int nwann = 0;
    f_nwann_(&nwann);

    for (int iwf = 1; iwf <= nwann; iwf++)
    {
      // get center of wave function
      double tcenter[3];
      f_wann_center_(&iwf, &tcenter[0]);
      Vector<double,3> center;
      center[0] = tcenter[0]; center[1] = tcenter[1]; center[2] = tcenter[2];
      //printf("%1d%25.5f%15.5f%15.5f\n", iwf, tcenter[0], tcenter[1], tcenter[2]);
      int ishell = 0;
      for (int ir = 0; ishell < 2; ir++)
      {
        // set shell
        ishell = _rshells[ir];
        if (ishell < 2)
        {
          // get lattice vector and convert to lattice coordinates
          vec3dT rvecc = _lvecs[ir];
          vec3dT rvecl = r3mv(_ainv, rvecc);
          Vector<int,3> itr;
          itr[0] = (int) floor(rvecl[0]); itr[1] = (int) floor(rvecl[1]); itr[2] = (int) floor(rvecl[2]);
          printf("%2d    %1d%20.5f%12.5f%12.5f%20.3f%6.3f%6.3f\n", ir, ishell, rvecc[0], rvecc[0], rvecc[0],
              rvecl[0], rvecl[1], rvecl[2]);
          // convert
          cfunctionT w = cfactoryT(_world).functor(cfunctorT(new Wannier<std::complex<double>,3>(iwf, itr, center)));
          w.reconstruct();
          _wfs.push_back(wannierf(iwf-1, ir, w));
        }
      }
    }
  }

  void read_lattice()
  {
    // initialize structures
    _avec = Tensor<double>(3,3);
    _bvec = Tensor<double>(3,3);
    _ainv = Tensor<double>(3,3);
    _binv = Tensor<double>(3,3);
    // create temporary arrays for exciting
    double avec[3][3];
    double ainv[3][3];
    double bvec[3][3];
    double binv[3][3];
    // read lattice info into temp arrays
    f_lvec_(&avec[0], &ainv[0], &bvec[0], &binv[0]);
    // convert (switch indicies) and store in tensors
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        _avec(i,j) = avec[j][i];
        _bvec(i,j) = bvec[j][i];
        _ainv(i,j) = ainv[j][i];
        _binv(i,j) = binv[j][i];
      }
    }
  }

  struct vectorLengthFunctor : public binary_function<vec3dT, vec3dT, bool>
  {
      bool operator()( vec3dT lhs, vec3dT rhs)
      {
        double llen = sqrt(lhs[0]*lhs[0] + lhs[1]*lhs[1] + lhs[2]*lhs[2]);
        double rlen = sqrt(rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);
        return (llen < rlen);
      }
  };

  void generate_gvectors()
  {
    double gmaxlen = 20.0;
    int ghi = 31;
    int glo = -30;

    for (int i0 = glo; i0 < ghi; i0++)
    {
      for (int i1 = glo; i1 < ghi; i1++)
      {
        for (int i2 = glo; i2 < ghi; i2++)
        {
          vec3dT t1;
          t1[0] = i0; t1[1] = i1; t1[2] = i2;
          vec3dT gvec = r3mv(_bvec, t1);
          _gvecs.push_back(gvec);
        }
      }
    }
    // sort g-vectors according to length
    std::sort(_gvecs.begin(), _gvecs.end(), vectorLengthFunctor());
    // lenths of g-vectors in cartesian
    std::vector<double> glens(_gvecs.size());
    for (unsigned int i = 0; i < _gvecs.size(); i++)
    {
      vec3dT gvec = _gvecs[i];
      double len = sqrt(gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2]);
      glens[i] = len;
    }
//    for (unsigned int i = 0; i < gvecs.size(); i++)
//    {
//      vec3dT gvec = gvecs[i];
//      double glen = glens[i];
//      printf("%12.7f%12.7f%12.7f%25.7f\n", gvec[0], gvec[1], gvec[2], glen);
//    }
  }

  void generate_lvectors()
  {
    double rmaxlen = 20.0;
    int hi = 21;
    int lo = -20;

    bool bstop = false;
    for (int i0 = lo; i0 < hi && !bstop; i0++)
    {
      for (int i1 = lo; i1 < hi && !bstop; i1++)
      {
        for (int i2 = lo; i2 < hi && !bstop; i2++)
        {
          vec3dT t1;
          t1[0] = i0; t1[1] = i1; t1[2] = i2;
          vec3dT lvec = r3mv(_avec, t1);
          double lsize = sqrt(lvec[0]*lvec[0] + lvec[1]*lvec[1] + lvec[2]*lvec[2]);
          if (lsize < rmaxlen)
            _lvecs.push_back(lvec);
        }
      }
    }
    // sort g-vectors according to length
    std::sort(_lvecs.begin(), _lvecs.end(), vectorLengthFunctor());
    // lenths of g-vectors in cartesian
    _rlens = std::vector<double>(_lvecs.size());
    for (unsigned int i = 0; i < _rlens.size(); i++)
    {
      vec3dT lvec = _lvecs[i];
      double len = sqrt(lvec[0]*lvec[0] + lvec[1]*lvec[1] + lvec[2]*lvec[2]);
      _rlens[i] = len;
    }

    // make vector of shell indicies
    double eps = 1e-8;
    _rshells = std::vector<int>(_lvecs.size());
    int ishell = 0;
    for (unsigned int il = 0; il < _rlens.size(); il++)
    {
      if (il < (_rlens.size()-1))
      {
        _rshells[il] = ishell;
        if (abs(_rlens[il+1] - _rlens[il]) > eps)
        {
          ishell++;
        }
      }
      else
      {
        if (abs(_rlens[il] - _rlens[il-1]) > eps)
        {
          ishell++;
        }
        _rshells[il] = ishell;
      }
    }

    // debug print stuff
    for (unsigned int i = 0; i < _lvecs.size(); i++)
    {
      vec3dT rvec = _lvecs[i];
      double rlen = _rlens[i];
      int ishell = _rshells[i];
      printf("%12.7f%12.7f%12.7f%25.7f%25d\n", rvec[0], rvec[1], rvec[2], rlen, ishell);
    }
  }

  void read_kpoints()
  {
    // get size of k-mesh
    f_ngridk_(&_ngridk0, &_ngridk1, &_ngridk2);
    // number of k-points
    int nkpts = 0;
    f_nkpts_(&nkpts);
    // create vectors associated with k-points
    _vkl = std::vector<vec3dT>(nkpts);
    _vkc = std::vector<vec3dT>(nkpts);
    _wkpt = std::vector<vec3dT>(nkpts);
    for (int ik = 1; ik <= nkpts; ik++)
    {
      // get k-point in lattice coordinates
      double t0, t1, t2, t3;
      f_kpts_(&ik, &t0, &t1, &t2, &t3);
      vec3dT kvecl;
      kvecl[0] = t0; kvecl[1] = t1; kvecl[2] = t2;
      _vkl[ik-1] = kvecl;
      // convert to cartesian
      vec3dT kvecc = r3mv(_bvec, kvecl);
      _vkc[ik-1] = kvecc;
      // weight
      _wkpt[ik-1] = t3;
      //printf("%12.7f%12.7f%12.7f%25.7f\n", kvecl[0], kvecl[1], kvecl[2], t3);
    }
  }

  void response(vec3dT qvec)
  {
    // convert q to lattice coordinates
    vec3dT qvl = qvec;
    qvl[0] = qvl[0] / _ngridk0 + 1e-10;
    qvl[1] = qvl[1] / _ngridk1 + 1e-10;
    qvl[2] = qvl[2] / _ngridk2 + 1e-10;
    // find the G vector that brings the q vector back to the first BZ
    vec3dT gvec;
    for (int i = 0; i < 3; i++)
    {
      gvec[i] = floor(qvl[i]);
      qvl[i] -= gvec[i];
    }
    //printf("\n\n%12.7f%12.7f%12.7f\n\n", qvl[0], qvl[1], qvl[2]);

    // create vector of k+q vector indicies for response
    double eps = 1e-8;
    std::vector<unsigned int> indxkpq(_vkl.size());
    for (unsigned int ik = 0; ik < _vkl.size(); ik++)
    {
      vec3dT tv1;
      vec3dT vkl = _vkl[ik];
      for (int i = 0; i < 3; i++)
      {
        tv1[i] = qvl[i] + vkl[i];
        // make exception for any of the components that equal 1.0
        if (abs(tv1[i] - 1.0) < eps)
          tv1[i] -= 1.0;
      }
      int indx = findkpt(tv1);
      // abort if error
      if (indx < 0) abort();
      indxkpq[ik] = indx;

//      // debug output
//      vec3dT kv2 = _vkl[indx];
//      printf("%12.7f%12.7f%12.7f%30.7f%12.7f%12.7f\n", vkl[0], vkl[1], vkl[2], kv2[0], kv2[1], kv2[2]);
    }

    // get number of wannier functions
    int nwann = 0;
    f_nwann_(&nwann);
    // read in wannier coeffs
    // this looks like the following
    // wann_c('wf index', 'band index', 'k point')
    Tensor< std::complex<double> > wann_c(nwann, nwann, _vkl.size());
    for (int iwf = 1; iwf <= nwann; iwf++)
    {
      for (int iband = 1; iband <= nwann; iband++)
      {
        for (int ik = 1; ik <= _vkl.size(); ik++)
        {
          double val[2];
          f_wann_c_(&iwf, &iband, &ik, &val[0]);
          // need the hermitian conjugate
          wann_c(iband-1, iwf-1, ik-1) = std::complex<double>(val[0],-val[1]);
        }
      }
    }

    // Compute the product of coeff's (boy, this is nasty)
    Tensor< std::complex<double> > c1kc2kq(nwann, nwann, nwann, nwann, _vkl.size());
    for (int iwf = 0; iwf < nwann; iwf++)
    {
      for (int jwf = 0; jwf < nwann; jwf++)
      {
        for (int iband = 0; iband < nwann; iband++)
        {
          for (int jband = 0; jband < nwann; jband++)
          {
            for (int ik = 0; ik < _vkl.size(); ik++)
            {
              std::complex<double> first = wann_c(iwf, iband, ik);
              std::complex<double> second = wann_c(jwf, jband, indxkpq[ik]);
              c1kc2kq(iwf, iband, jwf, jband, ik) = conj(first) * second;
            }
          }
        }
      }
    }
  }
};

void test_wannier(World& world)
{
    //cener point of the box
    Vector<double,3> center(0.0);
    //index of Wannier function
    int n=2;

    exciting_init_();

    // Function defaults
    int funck = 5;
    double thresh = 1e-3;
    double bsize = 20.0;
    FunctionDefaults<3>::set_k(funck);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);
    Vector<int,3> itr1;
    Vector<int,3> itr2;

    itr1[0]=0;itr1[1]=0;itr1[2]=0;
    itr2[0]=1;itr2[1]=0;itr2[2]=0;


    cfunctionT w = cfactoryT(world).functor(cfunctorT(new Wannier<std::complex<double>,3>(n, itr1, center)));
    w.reconstruct();
    cfunctionT w2 = cfactoryT(world).functor(cfunctorT(new Wannier<std::complex<double>,3>(n, itr2, center)));
    w2.reconstruct();
    double wnorm = w.norm2();
    double wnorm2 = w2.norm2();
    double tmp = real(inner(w,w2));
    if (world.rank() == 0) {
      printf("Normalization of wannier function is: %14.8f\n\n", wnorm);
      printf("Normalization of translated wannier function is: %14.8f\n\n", wnorm2);
      printf("Inner product is: %14.8f\n\n", tmp);
    }
//    // Plot to OpenDX
//    vector<long> npt(3,101);
//    plotdx(w, "wannier2.dx", FunctionDefaults<3>::get_cell(), npt);

}

void test_wannier2(World& world)
{
    //cener point of the box
    Vector<double,3> center(0.0);
    //index of Wannier function

    ExcitingApp app(world);

}

//void test_wannier3(World& world)
//{
//    //k-mesh division
//    Vector<int,3> nk(4);
//    //cener point of the box
//    Vector<double,3> center(0.0);
//
//    readinput_();
//    wann_init1_();
//
//    // Function defaults
//    int funck = 6;
//    double thresh = 1e-4;
//    double bsize = 6.0;
//    FunctionDefaults<3>::set_k(funck);
//    FunctionDefaults<3>::set_thresh(thresh);
//    FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);
//
//    std::vector<cfunctionT> w = zero_functions<std::complex<double>,3>(world,5);
//
//    for (int n = 1; n <= 5; n++)
//    {
//      w[n-1] = cfactoryT(world).functor(cfunctorT(new Wannier<std::complex<double>,3>(n, center, nk)));
//      double wnorm = w[n-1].norm2();
//      if (world.rank() == 0)
//        printf("Normalization of wannier function #%d is: %14.8f\n", n, wnorm);
//    }
//    if (world.rank() == 0) printf("\n\n");
//
//    for (int i = 1; i <= 5; i++)
//    {
//      for (int j = 1; j < i; j++)
//      {
//        double tmp = real(inner(w[i],w[j]));
//        if (world.rank() == 0)
//          printf("Inner product (%d,%d) = %15.8f\n", i, j, tmp);
//      }
//    }
//}
//
//void compute_U(World& world)
//{
//    //k-mesh division
//    Vector<int,3> nk(4);
//    //cener point of the box
//    Vector<double,3> center(0.0);
//
//    readinput_();
//    wann_init1_();
//
//    // Function defaults
//    int funck = 7;
//    double thresh = 1e-5;
//    double bsize = 6.0;
//    FunctionDefaults<3>::set_k(funck);
//    FunctionDefaults<3>::set_thresh(thresh);
//    FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);
//
//    SeparatedConvolution<double,3> cop = CoulombOperator<double>(world,
//        FunctionDefaults<3>::get_k(), 1e-8, thresh * 0.1);
//
//    // Vector of cfunctionT's
//    std::vector<cfunctionT> w;
//    // Compress wannier functions and add to vector
//    for (int n = 1; n <= 5; n++)
//    {
//      cfunctionT w_ = cfactoryT(world).functor(cfunctorT(new Wannier<std::complex<double>,3>(n, center, nk)));
//      w_.truncate();
//      w.push_back(w_);
//    }
//
//    for (int n = 1; n <= 5; n++)
//    {
//      functionT tmp_n = abs_square(w[n-1]);
//      tmp_n.truncate();
//      functionT tmp_n_apply = apply(cop, tmp_n);
//      tmp_n_apply.truncate();
//      for (int p = 1; p <= n; p++)
//      {
//        functionT tmp_p = abs_square(w[p-1]);
//        tmp_p.truncate();
//        // Compute inner product
//        double ip = abs(inner(w[n-1],w[p-1]));
//        // Compute U
//        double U = inner(tmp_p, tmp_n_apply);
//        if (world.rank() == 0) printf("%4.1d%4.1d%15.8f%15.8f\n", n, p, ip, U);
//      }
//    }
//    if (world.rank() == 0) printf("\n\n");
//
//}
//
//void test_wannier2(World& world)
//{
//    //k-mesh division
//    Vector<int,3> nk(4);
//    //cener point of the box
//    Vector<double,3> center(0.0);
//    //index of Wannier function
//    int n=3;
//
//    exciting_init_();
//
//
//    if (world.rank() == 0)  printf("\n");
//    double bsize = 6.0;
//    int npts = 30000;
//
//    for (int k = 0; k < npts; k++)
//    {
//      double z = (k+0.5) * (2.0*bsize/npts) - bsize;
//      double fval = abs(func(n, nk, 0.0, 0.0, z));
//      if (world.rank() == 0)
//        printf("%20.12f%20.12f\n", z, fval);
//    }
//}

#define TO_STRING(s) TO_STRING2(s)
#define TO_STRING2(s) #s

//*****************************************************************************
int main(int argc, char** argv)
{
  initialize(argc, argv);
  World world(MPI::COMM_WORLD);
  if (world.rank() == 0)
  {
    print("");
    print("--------------------------------------------");
    print("   MADNESS", " multiresolution testsuite");
    print("--------------------------------------------");
    print("");
    print("   number of processors ...", world.size());
    print("    processor frequency ...", cpu_frequency());
    print("            host system ...", TO_STRING(HOST_SYSTEM));
    print("             byte order ...", TO_STRING(MADNESS_BYTE_ORDER));
    print("          configured by ...", MADNESS_CONFIGURATION_USER);
    print("          configured on ...", MADNESS_CONFIGURATION_HOST);
    print("          configured at ...", MADNESS_CONFIGURATION_DATE);
    print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
    print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
#ifdef WORLD_WATCHDOG
    print("               watchdog ...", WATCHDOG_BARK_INTERVAL,
        WATCHDOG_TIMEOUT);
#endif
#ifdef OPTERON_TUNE
    print("             tuning for ...", "opteron");
#elif defined(CORE_DUO_TUNE)
    print("             tuning for ...", "core duo");
#else
    print("             tuning for ...", "core2");
#endif
#ifdef BOUNDS_CHECKING
    print(" tensor bounds checking ...", "enabled");
#endif
#ifdef TENSOR_INSTANCE_COUNT
    print("  tensor instance count ...", "enabled");
#endif
    print(" ");
  }

  try
  {
    printf("WSTHORNTON: Starting up the world ... \n");

    startup(world,argc,argv);
    if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
    test_wannier2(world);
  }
  catch (const MPI::Exception& e)
  {
    //        print(e);
    error("caught an MPI exception");
  }
  catch (const madness::MadnessException& e)
  {
    print(e);
    error("caught a MADNESS exception");
  }
  catch (const madness::TensorException& e)
  {
    print(e);
    error("caught a Tensor exception");
  }
  catch (const char* s)
  {
    print(s);
    error("caught a string exception");
  }
  catch (const std::string& s)
  {
    print(s);
    error("caught a string (class) exception");
  }
  catch (const std::exception& e)
  {
    print(e.what());
    error("caught an STL exception");
  }
  catch (...)
  {
    error("caught unhandled exception");
  }

  if (world.rank() == 0)
    print("entering final fence");
  world.gop.fence();
  if (world.rank() == 0)
    print("done with final fence");
  if (world.rank() == 0)
    print("Final tensor instance count", BaseTensor::get_instance_count());

  finalize();
  return 0;
}

