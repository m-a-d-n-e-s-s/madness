/*
 * electronicstructureapp.h
 *
 *  Created on: Nov 5, 2008
 *      Author: eh7
 */

#ifndef ELECTRONICSTRUCTUREAPP_H_
#define ELECTRONICSTRUCTUREAPP_H_

#include <mra/mra.h>
#include <misc/ran.h>
#include "electronicstructureparams.h"
#include "poperator.h"
#include "libxc.h"
#include "complexfun.h"

struct KOrbital;

typedef SharedPtr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<std::complex<double>,3> > functorT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > rfunctorT;
typedef Function<std::complex<double>,3> functionT;
typedef Function<double,3> rfunctionT;
typedef vector<functionT> vecfuncT;
typedef vector<rfunctionT> rvecfuncT;
typedef vector<KOrbital> kvecfuncT;
typedef Tensor< std::complex<double> > tensorT;
typedef Tensor<double> rtensorT;
typedef FunctionFactory<std::complex<double>,3> factoryT;
typedef FunctionFactory<double,3> rfactoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef SharedPtr<operatorT> poperatorT;

void print_cube(World& world, const Function<double,3>& f, int npts)
{
  f.reconstruct();
  if (world.rank() == 0)  printf("\n");
  Tensor<double> csize = FunctionDefaults<3>::get_cell_width();

  for (int i = 0; i < npts; i++)
  {
    for (int j = 0; j < npts; j++)
    {
      for (int k = 0; k < npts; k++)
      {
        double x = (i+0.5) * (csize[0]/npts) - csize[0]/2;
        double y = (j+0.5) * (csize[1]/npts) - csize[1]/2;
        double z = (k+0.5) * (csize[2]/npts) - csize[2]/2;
        coordT p(0.0);
        p[0] = x; p[1] = y; p[2] = z;
        if (world.rank() == 0)
          printf("%10.2f%10.2f%10.2f%15.8f\n", x, y, z, f(p));
      }
    }
  }
}

void print_cube(World& world, const Function<double,3>& f1, const Function<double,3>& f2, int npts)
{
  f1.reconstruct();
  f2.reconstruct();
  if (world.rank() == 0)  printf("\n");
  Tensor<double> csize = FunctionDefaults<3>::get_cell_width();

  for (int i = 0; i < npts; i++)
  {
    for (int j = 0; j < npts; j++)
    {
      for (int k = 0; k < npts; k++)
      {
        double x = (i+0.5) * (csize[0]/npts) - csize[0]/2;
        double y = (j+0.5) * (csize[1]/npts) - csize[1]/2;
        double z = (k+0.5) * (csize[2]/npts) - csize[2]/2;
        coordT p(0.0);
        p[0] = x; p[1] = y; p[2] = z;
        if (world.rank() == 0)
          printf("%10.2f%10.2f%10.2f%15.8f%15.8f\n", x, y, z, f1(p), f2(p));
      }
    }
  }
}

void print_cube(World& world, const Function<double,3>& f1, const Function<double,3>& f2,
    const Function<double,3>& f3, int npts)
{
  f1.reconstruct();
  f2.reconstruct();
  f3.reconstruct();
  if (world.rank() == 0)  printf("\n");
  Tensor<double> csize = FunctionDefaults<3>::get_cell_width();

  for (int i = 0; i < npts; i++)
  {
    for (int j = 0; j < npts; j++)
    {
      for (int k = 0; k < npts; k++)
      {
        double x = (i+0.5) * (csize[0]/npts) - csize[0]/2;
        double y = (j+0.5) * (csize[1]/npts) - csize[1]/2;
        double z = (k+0.5) * (csize[2]/npts) - csize[2]/2;
        coordT p(0.0);
        p[0] = x; p[1] = y; p[2] = z;
        if (world.rank() == 0)
          printf("%10.2f%10.2f%10.2f%15.8f%15.8f%15.8f\n", x, y, z, f1(p), f2(p), f3(p));
      }
    }
  }
}

class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};

    LevelPmap(World& world) : nproc(world.nproc()) {}

    /// Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;
        if (n <= 3 || (n&0x1)) hash = key.hash();
        else hash = key.parent().hash();
        //hashT hash = key.hash();
        return hash%nproc;
    }
};

class MolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
public:
    MolecularPotentialFunctor(const MolecularEntity& mentity)
        : _mentity(mentity)
    {}

    double operator()(const coordT& x) const
    {
      return _mentity.nuclear_attraction_potential(x[0], x[1], x[2]);
    }
};

class MolecularNuclearChargeDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
    const double R;
    const bool periodic;
public:
    MolecularNuclearChargeDensityFunctor(const MolecularEntity& mentity, const double& R,
        const bool& periodic) : _mentity(mentity), R(R), periodic(periodic) {
    }

    double operator()(const coordT& x) const
    {
        double big = 0.5*R + 6.0*_mentity.smallest_length_scale();
        // Only one contribution at any point due to the short
        // range of the nuclear charge density
        if (periodic)
        {
            for (int xr = -1; xr <= 1; xr += 1)
            {
                double xx = x[0] + xr*R;
                if (xx < big && xx > -big)
                {
                    for (int yr = -1; yr <= 1; yr += 1)
                    {
                        double yy = x[1] + yr*R;
                        if (yy < big && yy > -big)
                        {
                            for (int zr = -1; zr <= 1; zr += 1)
                            {
                                double zz = x[2] + zr*R;
                                if (zz < big && zz > -big)
                                {
                                    return _mentity.nuclear_charge_density(xx, yy, zz);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            return _mentity.nuclear_charge_density(x[0], x[1], x[2]);
        }
        return 0.0;
    }
};

class AtomicBasisFunctor : public FunctionFunctorInterface<double,3> {
private:
  const AtomicBasisFunction aofunc;
  const double R;
  const bool periodic;
public:
  AtomicBasisFunctor(const AtomicBasisFunction& aofunc, double R, bool periodic) : aofunc(aofunc),
    R(R), periodic(periodic) {}

  double operator()(const coordT& x) const
  {
    double value = 0.0;
    if (periodic)
    {
      for (int xr = -1; xr <= 1; xr += 1)
      {
        for (int yr = -1; yr <= 1; yr += 1)
        {
          for (int zr = -1; zr <= 1; zr += 1)
          {
            value += aofunc(x[0]+xr*R, x[1]+yr*R, x[2]+zr*R);
          }
        }
      }
    }
    else
    {
      value = aofunc(x[0], x[1], x[2]);
    }
    return value;
  }
};

struct KPoint
{
  coordT k;
  double weight;

  KPoint()
  {
    k[0] = 0.0; k[1] = 0.0; k[2] = 0.0;
    weight = 0.0;
  }

  KPoint(const coordT& k, const double& weight)
   : k(k), weight(weight) {}

  template <typename Archive>
  void serialize(Archive& ar) {
      ar & k & weight;
  }

};

std::istream& operator >> (std::istream& is, KPoint& kpt)
{
  for (int i = 0; i < kpt.k.size(); i++)
    is >> kpt.k[i];
  is >> kpt.weight;
  return is;
}

struct KOrbital
{
  coordT k;
  double weight;
  functionT orbital;

  KOrbital(const coordT& k, const double& weight, const functionT& orbital)
   : k(k), weight(weight), orbital(orbital) {}
};

double rsquared(const coordT& r) {
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

template <typename Q, int NDIM>
Function<Q,NDIM> pdiff(const Function<Q,NDIM>& f, int axis, bool fence = true)
{
  Function<Q,NDIM>& g = const_cast< Function<Q,NDIM>& >(f);
  // Check for periodic boundary conditions
  Tensor<int> oldbc = g.get_bc();
  Tensor<int> bc(NDIM,2);
  bc(___) = 1;
  g.set_bc(bc);
  // Do calculation
  Function<Q,NDIM> rf = diff(g,axis,fence);
  // Restore previous boundary conditions
  g.set_bc(oldbc);
  return rf;
}

class ElectronicStructureApp
{
public:
  ElectronicStructureApp(World& world, const std::string& filename)
   : _world(world)
  {
    init(filename);
  }

  std::vector<KPoint> read_kpoints(const std::string& filename)
  {
    ifstream indata;
    indata.open(filename.c_str());
    if(!indata)
    {
      // file couldn't be opened
      MADNESS_EXCEPTION("Error: file could not be opened", 0);
    }

    std::vector<KPoint> kpoints;
    char tmpstr[120];
    indata.getline(tmpstr, 120, '\n');
    while(!indata.eof())
    {
      int num;
      KPoint kpt;
      int junk;
      indata >> num;
      indata >> kpt;
      indata >> junk;
      kpoints.push_back(kpt);
    }
    // Delete last entry (duplicate)
    kpoints.pop_back();

    return kpoints;
  }

  void init(const std::string& filename)
  {
    // params
    if (_world.rank() == 0)
    {
      _params.read_file(filename);
      //_params.fractional = false;
    }
    _world.gop.broadcast_serializable(_params, 0);
    if (_params.fractional)
      FunctionDefaults<3>::set_cubic_cell(0,_params.L);
    else
      FunctionDefaults<3>::set_cubic_cell(-_params.L/2,_params.L/2);
    FunctionDefaults<3>::set_thresh(_params.thresh);
    FunctionDefaults<3>::set_k(_params.waveorder);

    // mentity and aobasis
    if (_world.rank() == 0)
    {
      _aobasis.read_file("sto-3g");
      _mentity.read_file(filename, _params.fractional);
      _mentity.center();
      if (_params.periodic && _params.kpoints)
        _kpoints = read_kpoints("KPOINTS.OUT");
      else
        _kpoints.push_back(KPoint(coordT(0.0), 1.0));
    }
    _world.gop.broadcast_serializable(_mentity, 0);
    _world.gop.broadcast_serializable(_aobasis, 0);
    _world.gop.broadcast_serializable(_kpoints, 0);

  }

  void make_nuclear_potential()
  {
    if (_world.rank() == 0) print("Making nuclear potential ..\n\n");
    Tensor<double> csize = FunctionDefaults<3>::get_cell();
    if (_world.rank() == 0)
    {
      print("cell(x) is ",csize(0,0), csize(0,1));
      print("cell(y) is ",csize(1,0), csize(1,1));
      print("cell(z) is ",csize(2,0), csize(2,1));
    }
    if (_params.ispotential) // potential
    {
      _vnucrhon = rfactoryT(_world).functor(rfunctorT(new MolecularPotentialFunctor(_mentity))).thresh(_params.thresh * 0.1).truncate_on_project();
      _vnuc = copy(_vnucrhon);
      _vnuc.reconstruct();
    }
    else // charge density
    {
      std::vector<coordT> specialpts;
      for (int i = 0; i < _mentity.natom(); i++)
      {
        coordT pt(0.0);
        Atom atom = _mentity.get_atom(i);
        pt[0] = atom.x; pt[1] = atom.y; pt[2] = atom.z;
        specialpts.push_back(pt);
        print("Special point: ", pt);
      }
      double now = wall_time();
      _vnucrhon = rfactoryT(_world).functor(
          rfunctorT(new MolecularNuclearChargeDensityFunctor(_mentity, _params.L, _params.periodic))).
          thresh(_params.thresh).initial_level(6).truncate_on_project().specialpts(specialpts);
      
      if (_world.rank() == 0) printf("%f\n", wall_time() - now);
      if (_world.rank() == 0) print("calculating trace of rhon ..\n\n");
      double rtrace = _vnucrhon.trace();
      if (_world.rank() == 0) print("rhon trace = ", rtrace);
      SeparatedConvolution<double,3>* op = 0;
      if (_params.periodic) // periodic
      {
        Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
        op = PeriodicCoulombOpPtr<double,3>(_world, _params.waveorder,_params.lo, _params.thresh, cellsize);
      }
      else // not periodic
      {
        op = CoulombOperatorPtr<double>(_world, _params.waveorder,_params.lo, _params.thresh);
      }
      now = wall_time();
      _vnucrhon.truncate();
      _vnuc = apply(*op, _vnucrhon);
      if (_world.rank() == 0) printf("%f\n", wall_time() - now);
      if (_world.rank() == 0) print("Done creating nuclear potential ..\n");
      delete op;
    }
    
    vector<long> npt(3,101);
    plotdx(_vnuc, "vnuc.dx", FunctionDefaults<3>::get_cell(), npt);
  }

  struct GuessDensity : public FunctionFunctorInterface<double,3> {
      const MolecularEntity& mentity;
      const AtomicBasisSet& aobasis;
      double R;
      const bool periodic;

      double operator()(const coordT& x) const
      {
        double value = 0.0;
        if (periodic)
        {
          for (int xr = -1; xr <= 1; xr += 1)
          {
            for (int yr = -1; yr <= 1; yr += 1)
            {
              for (int zr = -1; zr <= 1; zr += 1)
              {
                value += aobasis.eval_guess_density(mentity,
                    x[0]+xr*R, x[1]+yr*R, x[2]+zr*R);
              }
            }
          }
        }
        else
        {
          value = aobasis.eval_guess_density(mentity, x[0], x[1], x[2]);
        }
        return value;
      }

      GuessDensity(const MolecularEntity& mentity, const AtomicBasisSet& aobasis,
          const double& R, const bool& periodic)
      : mentity(mentity), aobasis(aobasis), R(R), periodic(periodic) {}
  };

  static double munge(double r) {
      if (r < 1e-12) r = 1e-12;
      return r;
  }

  rfunctionT
  make_lda_potential(World& world,
                     const rfunctionT& arho,
                     const rfunctionT& brho,
                     const rfunctionT& adelrhosq,
                     const rfunctionT& bdelrhosq)
  {
//      MADNESS_ASSERT(!_params.spinpol);
      rfunctionT vlda = copy(arho);
      vlda.unaryop(&::libxc_ldaop);
      return vlda;
  }

  rvecfuncT project_ao_basis(World& world) {
      rvecfuncT ao(_aobasis.nbf(_mentity));

      Level initial_level = 3;
      for (int i=0; i < _aobasis.nbf(_mentity); i++) {
          rfunctorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
              _params.L, _params.periodic));
          ao[i] = rfactoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
      }
      world.gop.fence();

      vector<double> norms;

      while (1) {
          norms = norm2(world, ao);
          initial_level += 2;
          if (initial_level >= 11) throw "project_ao_basis: projection failed?";
          int nredone = 0;
          for (int i=0; i<_aobasis.nbf(_mentity); i++) {
              if (norms[i] < 0.5) {
                  nredone++;
                  if (world.rank() == 0) print("re-projecting ao basis function", i,"at level",initial_level);
                  rfunctorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
                      _params.L, _params.periodic));
                  ao[i] = rfactoryT(world).functor(aofunc).initial_level(6).truncate_on_project().nofence();
              }
          }
          world.gop.fence();
          if (nredone == 0) break;
      }

      norms = norm2(world, ao);

      for (int i=0; i<_aobasis.nbf(_mentity); i++) {
          if (world.rank() == 0 && fabs(norms[i]-1.0)>1e-3) print(i," bad ao norm?", norms[i]);
          norms[i] = 1.0/norms[i];
      }

      scale(world, ao, norms);

      norms = norm2(world, ao);

      for (int i=0; i<_aobasis.nbf(_mentity); i++) {
          if (world.rank() == 0 && fabs(norms[i]-1.0)>1e-3) print(i," bad ao norm?", norms[i]);
          norms[i] = 1.0/norms[i];
      }

      scale(world, ao, norms);

      return ao;
  }

  tensorT kinetic_energy_matrix(World& world, const rvecfuncT& v, const KPoint k = KPoint(coordT(0.0),0.0))
  {
      reconstruct(world, v);
      int n = v.size();
      tensorT c(n,n);
      const std::complex<double> I = std::complex<double>(0.0,1.0);
      double k0 = k.k[0]; double k1 = k.k[1]; double k2 = k.k[2];
      if (_params.periodic)
      {
        for (int i = 0; i < n; i++)
        {
          functionT dv_i_0 = function_real2complex(pdiff(v[i],0)) - I*k0*v[i];
          functionT dv_i_1 = function_real2complex(pdiff(v[i],1)) - I*k1*v[i];
          functionT dv_i_2 = function_real2complex(pdiff(v[i],2)) - I*k2*v[i];
          for (int j = 0; j < n; j++)
          {
            functionT dv_j_0 = function_real2complex(pdiff(v[j],0)) + I*k0*v[j];
            functionT dv_j_1 = function_real2complex(pdiff(v[j],1)) + I*k1*v[j];
            functionT dv_j_2 = function_real2complex(pdiff(v[j],2)) + I*k2*v[j];
            c(i,j) = inner(dv_i_0,dv_j_0) + inner(dv_i_1,dv_j_1) + inner(dv_i_2,dv_j_2);
          }
        }
      }
      else
      {
        rtensorT r(n,n);
        for (int axis=0; axis<3; axis++)
        {
            rvecfuncT dv = diff(world,v,axis);
            r += matrix_inner(world, dv, dv, true);
            dv.clear(); // Allow function memory to be freed
        }
        c = tensor_real2complex(r);
      }

//      // DEBUG
//      rtensorT r(n,n);
//      for (int axis=0; axis<3; axis++) {
//          rvecfuncT dv = diff(world,v,axis);
//          r += matrix_inner(world, dv, dv, true);
//          dv.clear(); // Allow function memory to be freed
//      }
//
//      for (int i = 0; i < c.dim[0]; i++)
//      {
//        for (int j = 0; j < c.dim[1]; j++)
//        {
//          printf("%10.5f", imag(c(i,j)));
//        }
//        printf("\n");
//      }
//      printf("\n");
//      printf("\n");
//      for (int i = 0; i < r.dim[0]; i++)
//      {
//        for (int j = 0; j < r.dim[1]; j++)
//        {
//          printf("%10.5f", r(i,j));
//        }
//        printf("\n");
//      }
      return c.scale(0.5);
  }


  /// Initializes alpha and beta mos, occupation numbers, eigenvalues
  void initial_guess()
  {
    // Get initial guess for the electronic density
    if (_world.rank() == 0) print("Guessing rho ...\n\n");
    rfunctionT rho = rfactoryT(_world).functor(rfunctorT(
        new GuessDensity(_mentity, _aobasis, _params.L, _params.periodic)));
    rho.scale(_params.nelec/rho.trace());

    print(rho.trace());

    vector<long> npt(3,101);
    plotdx(rho, "rho_initial.dx", FunctionDefaults<3>::get_cell(), npt);

    rfunctionT vlocal;
    // Is this a many-body system?
    if (_params.nelec > 1)
    {
      if (_world.rank() == 0) print("Creating Coulomb op ...\n\n");
      SeparatedConvolution<double, 3>* op = 0;
      // Is this system periodic?
      if (_params.periodic)
      {
        Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
        op = PeriodicCoulombOpPtr<double, 3> (_world, _params.waveorder,
            _params.lo, _params.thresh * 0.1, cellsize);
      }
      else
      {
        op = CoulombOperatorPtr<double> (_world, _params.waveorder,
            _params.lo, _params.thresh * 0.1);
      }
      if (_world.rank() == 0) print("Building effective potential ...\n\n");
      rfunctionT vc = apply(*op, rho);
      //      vlocal = _vnuc + apply(*op, rho); //.scale(1.0-1.0/nel); // Reduce coulomb to increase binding
      vlocal = _vnuc + vc; //.scale(1.0-1.0/nel); // Reduce coulomb to increase binding
      rho.scale(0.5);
      // Do the LDA
      rfunctionT vlda = make_lda_potential(_world, rho, rho, rfunctionT(), rfunctionT());
      //vlocal = vlocal + make_lda_potential(_world, rho, rho, rfunctionT(), rfunctionT());
      vlocal = vlocal + vlda;
      delete op;
      vector<long> npt(3,101);
      plotdx(vc, "vc.dx", FunctionDefaults<3>::get_cell(), npt);
//      print_cube(_world, rho, vc, vlda, 5);
    }
    else
    {
      vlocal = _vnuc;
    }

    // Clear these functions
    rho.clear();
    vlocal.reconstruct();

    // These are our initial basis functions
    if (_world.rank() == 0) print("Projecting atomic orbitals ...\n\n");
    rvecfuncT ao = project_ao_basis(_world);

    for (unsigned int ai = 0; ai < ao.size(); ai++)
    {
      std::ostringstream strm;
      strm << "aod" << ai << ".dx";
      std::string fname = strm.str();
      vector<long> npt(3,101);
      plotdx(ao[ai], fname.c_str(), FunctionDefaults<3>::get_cell(), npt);
    }

//    {
//      rfunctionT dao0 = pdiff(ao[0], 0);
//      rfunctionT dao1 = pdiff(ao[0], 1);
//      rfunctionT dao2 = pdiff(ao[0], 2);
//      if (_world.rank() == 0)  printf("\n");
//      double L = _params.L;
//      double bstep = L / 100.0;
//      dao0.reconstruct();
//      dao1.reconstruct();
//      dao2.reconstruct();
//      ao[0].reconstruct();
//      for (int i = 0; i < 101; i++)
//      {
//        coordT p(-L / 2 + i * bstep);
//        if (_world.rank() == 0)
//          printf("%.2f\t\t%.8f\t%.8f\t%.8f\n", p[0], dao0(p), dao1(p), dao2(p));
////        if (_world.rank() == 0)
////          printf("%.2f\t\t%.8f\n", p[0], ao[0](p));
//      }
//      if (_world.rank() == 0) printf("\n");
//    }

    // Get size information from k-points and ao_basis so that we can correctly size
    // the _orbitals data structure and the eigs tensor
    int nao = ao.size();
    int nkpts = _kpoints.size();
    int norbs = nao * nkpts;
    //_orbitals = std::vector<functionT>(norbs, factoryT(_world));
    _eigs = rtensorT(norbs);
    _occs = rtensorT(norbs);
    for (int i = 0; i < norbs; i++) _occs[i] = _params.maxocc;
    print("norbs = ", norbs);

    // Build the overlap matrix
    if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
    rtensorT roverlap = matrix_inner(_world, ao, ao, true);
    // Convert to a complex tensor
    tensorT overlap = tensor_real2complex<double>(roverlap);

    // Build the potential matrix
    reconstruct(_world, ao);
    if (_world.rank() == 0) print("Building potential energy matrix ...\n\n");
    rvecfuncT vpsi = mul_sparse(_world, vlocal, ao, _params.thresh);
    _world.gop.fence();
    _world.gop.fence();
    compress(_world, vpsi);
    truncate(_world, vpsi);
    compress(_world, ao);

    rtensorT rpotential = matrix_inner(_world, vpsi, ao, true);
    // Convert to a complex tensor
    tensorT potential = tensor_real2complex<double>(rpotential);
    _world.gop.fence();
    vpsi.clear();
    _world.gop.fence();

    int kp = 0;
    if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
    // Need to do kinetic piece for every k-point
    for (int ki = 0; ki < nkpts; ki++)
    {
      // Get k-point from list
      KPoint kpt = _kpoints[ki];

        if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
      tensorT kinetic = kinetic_energy_matrix(_world, ao, kpt);

        if (_world.rank() == 0) print("Constructing Fock matrix ...\n\n");
      tensorT fock = kinetic + potential;
      fock = 0.5 * (fock + transpose(fock));

      tensorT c; rtensorT e;
        if (_world.rank() == 0) print("Diagonlizing Fock matrix ...\n\n");
      sygv(fock, overlap, 1, &c, &e);

      compress(_world, ao);
      _world.gop.fence();
      vecfuncT tmp_orbitals = transform(_world, ao, c(_, Slice(0, nao - 1)));
      _world.gop.fence();
      truncate(_world, tmp_orbitals);
      normalize(_world, tmp_orbitals);
      rtensorT tmp_eigs = e(Slice(0, nao - 1));

      if (_world.rank() == 0) printf("(%8.4f,%8.4f,%8.4f)\n",kpt.k[0], kpt.k[1], kpt.k[2]);
      if (_world.rank() == 0) print(tmp_eigs);
      if (_world.rank() == 0) print("\n");

      // DEBUG
      for (int i = 0; i < kinetic.dim[0]; i++)
      {
        for (int j = 0; j < kinetic.dim[1]; j++)
        {
          if (_world.rank() == 0) printf("%10.5f", real(kinetic(i,j)));
        }
        if (_world.rank() == 0) printf("\n");
      }
      if (_world.rank() == 0) printf("\n");
      if (_world.rank() == 0) printf("\n");

      for (int i = 0; i < potential.dim[0]; i++)
      {
        for (int j = 0; j < potential.dim[1]; j++)
        {
          if (_world.rank() == 0) printf("%10.5f", real(potential(i,j)));
        }
        if (_world.rank() == 0) printf("\n");
      }
      if (_world.rank() == 0) printf("\n");
      if (_world.rank() == 0) printf("\n");

      for (int i = 0; i < fock.dim[0]; i++)
      {
        for (int j = 0; j < fock.dim[1]; j++)
        {
          if (_world.rank() == 0) printf("%10.5f", real(fock(i,j)));
        }
        if (_world.rank() == 0) printf("\n");
      }
      if (_world.rank() == 0) printf("\n");
      if (_world.rank() == 0) printf("\n");

      // Fill in orbitals and eigenvalues
      int kend = kp + nao;
      for (int oi = kp, ti = 0; oi < kend; oi++, ti++)
      {

//        {
//          if (_world.rank() == 0)  printf("\n");
//          double L = _params.L;
//          double bstep = L / 100.0;
//          print("ti = ", ti);
//          tmp_orbitals[ti].reconstruct();
//          for (int i = 0; i < 101; i++)
//          {
//            coordT p(-L / 2 + i * bstep);
//            if (_world.rank() == 0)
//              printf("%5.2f%15.8f\n", p[0], (tmp_orbitals[ti])(p));
//          }
//          if (_world.rank() == 0) printf("\n");
//        }

//        _orbitals[oi] = tmp_orbitals[ti];
//        _eigs[oi] = tmp_eigs[ti];
        if (_world.rank() == 0) print(oi, ti, kp, kend);
        _orbitals.push_back(tmp_orbitals[ti]);
        _eigs[oi] = tmp_eigs[ti];
      }

      kp += nao;
    }
  }

  vecfuncT orbitals()
  {
    return _orbitals;
  }

  rtensorT eigs()
  {
    return _eigs;
  }

  ElectronicStructureParams params()
  {
    return _params;
  }

  MolecularEntity entity()
  {
    return _mentity;
  }

  rfunctionT vnucrhon()
  {
    return _vnucrhon;
  }

private:
  World& _world;
  MolecularEntity _mentity;
  AtomicBasisSet _aobasis;
  ElectronicStructureParams _params;
  rfunctionT _vnuc;
  rfunctionT _vnucrhon;
  vecfuncT _orbitals;
  rtensorT _eigs;
  rtensorT _occs;
  std::vector<KPoint> _kpoints;
};

#endif /* ELECTRONICSTRUCTUREAPP_H_ */
