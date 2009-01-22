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

typedef double valueT;
typedef SharedPtr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef vector<functionT> vecfuncT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef SharedPtr<operatorT> poperatorT;

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

    double operator()(const coordT& x) const {
        return _mentity.nuclear_attraction_potential(x[0], x[1], x[2]);
    }
};

class MolecularNuclearChargeDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
public:
    MolecularNuclearChargeDensityFunctor(const MolecularEntity& mentity)
        : _mentity(mentity)
    {}

    double operator()(const coordT& x) const {
        return _mentity.nuclear_charge_density(x[0], x[1], x[2]);
    }
};

class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
    const AtomicBasisSet& _aobasis;
public:
    MolecularGuessDensityFunctor(const MolecularEntity& mentity, const AtomicBasisSet& aobasis)
        : _mentity(mentity), _aobasis(aobasis)
    {}

    double operator()(const coordT& x) const {
        return _aobasis.eval_guess_density(_mentity, x[0], x[1], x[2]);
    }
};


class AtomicBasisFunctor : public FunctionFunctorInterface<double,3> {
private:
  const AtomicBasisFunction aofunc;
  const double R;
  const bool periodic;
public:
  AtomicBasisFunctor(const AtomicBasisFunction& aofunc, double R, bool periodic) : aofunc(aofunc),
    R(R), periodic(periodic)
  {}

  double operator()(const coordT& x) const
  {
    double value = 0.0;
    if (periodic)
    {
      for (int xr = -R; xr <= R; xr += R)
      {
        for (int yr = -R; yr <= R; yr += R)
        {
          for (int zr = -R; zr <= R; zr += R)
          {
            value += aofunc(x[0]+xr, x[1]+yr, x[2]+zr);
//            double diff = fabs(aofunc(x[0], x[1], x[2])-aofunc(x[0]+xr, x[1]+yr, x[2]+zr));
//            if (diff > 1e-3) printf("%10.8e%16.8e%16.8e\n", diff, aofunc(x[0], x[1], x[2]), aofunc(x[0]+xr, x[1]+yr, x[2]+zr));
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

/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const coordT& x) const {
        return x[axis];
    }
};

double rsquared(const coordT& r) {
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

/// A MADNESS functor to compute the cartesian moment x^i * y^j * z^k (i, j, k integer and >= 0)
class MomentFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int i, j, k;
public:
    MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
    double operator()(const coordT& r) const {
        double xi=1.0, yj=1.0, zk=1.0;
        for (int p=0; p<i; p++) xi *= r[0];
        for (int p=0; p<j; p++) yj *= r[1];
        for (int p=0; p<k; p++) zk *= r[2];
        return xi*yj*zk;
    }
};

tensorT sqrt(const tensorT& s, double tol=1e-8) {
    int n=s.dim[0], m=s.dim[1];
    MADNESS_ASSERT(n==m);
    tensorT c, e;
    //s.gaxpy(0.5,transpose(s),0.5); // Ensure exact symmetry
    syev(s, &c, &e);
    for (int i=0; i<n; i++) {
        if (e(i) < -tol) {
            MADNESS_EXCEPTION("Matrix square root: negative eigenvalue",i);
        }
        else if (e(i) < tol) { // Ugh ..
            print("Matrix square root: Warning: small eigenvalue ", i, e(i));
            e(i) = tol;
        }
        e(i) = 1.0/sqrt(e(i));
    }
    for (int j=0; j<n; j++) {
        for (int i=0; i<n; i++) {
            c(j,i) *= e(i);
        }
    }
    return c;
}

tensorT energy_weighted_orthog(const tensorT& s, const tensorT eps) {
    int n=s.dim[0], m=s.dim[1];
    MADNESS_ASSERT(n==m);
    tensorT d(n,n);
    for (int i=0; i<n; i++) d(i,i) = eps(i);
    tensorT c, e;
    sygv(d, s, 1, &c, &e);
    return c;
}


template <typename T, int NDIM>
Cost lbcost(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) {
  return 1;
}

class ElectronicStructureApp
{
public:
  ElectronicStructureApp(World& world, const std::string& filename)
   : _world(world)
  {
    init(filename);
  }

  void init(const std::string& filename)
  {
    if (_world.rank() == 0)
    {
      _params.read_file(filename);
      _aobasis.read_file("sto-3g");
      _mentity.read_file(filename);
      _mentity.center();
    }
    _world.gop.broadcast_serializable(_mentity, 0);
    _world.gop.broadcast_serializable(_params, 0);
    _world.gop.broadcast_serializable(_aobasis, 0);

    FunctionDefaults<3>::set_cubic_cell(-_params.L,_params.L);
    FunctionDefaults<3>::set_thresh(_params.thresh);
    FunctionDefaults<3>::set_k(_params.waveorder);
  }

  void make_nuclear_potential()
  {
    if (_world.rank() == 0) print("Making nuclear potential ..\n\n");
    if (_params.ispotential)
    {
      _vnucrhon = factoryT(_world).functor(functorT(new MolecularPotentialFunctor(_mentity))).thresh(_params.thresh * 0.1).truncate_on_project();
      _vnuc = copy(_vnucrhon);
      _vnuc.reconstruct();
    }
    else
    {
      _vnucrhon = factoryT(_world).functor(
          functorT(new MolecularNuclearChargeDensityFunctor(_mentity))).
          thresh(_params.thresh * 0.1).initial_level(6).truncate_on_project();
      if (_world.rank() == 0) print("calculating trace of rhon ..\n\n");
      double rtrace = _vnucrhon.trace();
      if (_world.rank() == 0) print("rhon trace = ", rtrace);
      SeparatedConvolution<double,3>* op = 0;
      if (_params.periodic)
      {
        Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
        op = PeriodicCoulombOpPtr<double,3>(_world, _params.waveorder,_params.lo, _params.thresh, cellsize);
      }
      else
      {
        op = CoulombOperatorPtr<double,3>(_world, _params.waveorder,_params.lo, _params.thresh);
      }
      _vnuc = apply(*op, _vnucrhon);
      delete op;
    }
  }

  struct GuessDensity : public FunctionFunctorInterface<double,3> {
      const MolecularEntity& mentity;
      const AtomicBasisSet& aobasis;
      double operator()(const coordT& x) const {
          return aobasis.eval_guess_density(mentity, x[0], x[1], x[2]);
      }
      GuessDensity(const MolecularEntity& mentity, const AtomicBasisSet& aobasis)
          : mentity(mentity), aobasis(aobasis) {}
  };

  static double munge(double r) {
      if (r < 1e-12) r = 1e-12;
      return r;
  }

//  static void ldaop(const Key<3>& key, tensorT& t) {
//      UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;x_rks_s__(&r, &q, &dq1);c_rks_vwn5__(&r, &q, &dq2); *_p0 = dq1+dq2);
//  }
//
//  static void ldaeop(const Key<3>& key, tensorT& t) {
//      UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q1; double q2; double dq;x_rks_s__(&r, &q1, &dq);c_rks_vwn5__(&r, &q2, &dq); *_p0 = q1+q2);
//  }


  functionT
  make_lda_potential(World& world,
                     const functionT& arho,
                     const functionT& brho,
                     const functionT& adelrhosq,
                     const functionT& bdelrhosq)
  {
//      MADNESS_ASSERT(!_params.spinpol);
      functionT vlda = copy(arho);
      vlda.unaryop(&::libxc_ldaop);
      return vlda;
  }


//  functionT make_lda_potential(const functionT& rho)
//  {
//    functionT V_rho = 0.5 * copy(rho);
//    V_rho.reconstruct();
//    V_rho.unaryop(&dft_xc_lda_V<3>);
//    return V_rho;
//  }
//
  vecfuncT project_ao_basis(World& world) {
      vecfuncT ao(_aobasis.nbf(_mentity));

      Level initial_level = 3;
      for (int i=0; i < _aobasis.nbf(_mentity); i++) {
          functorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
              _params.L, _params.periodic));
//             if (world.rank() == 0) {
//                 aobasis.get_atomic_basis_function(molecule,i).print_me(cout);
//             }
          ao[i] = factoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
      }
      world.gop.fence();
      //truncate(world, ao);

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
                  functorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
                      _params.L, _params.periodic));
                  ao[i] = factoryT(world).functor(aofunc).initial_level(6).truncate_on_project().nofence();
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

      return ao;
  }

  tensorT kinetic_energy_matrix(World& world, const vecfuncT& v) {
      reconstruct(world, v);
      int n = v.size();
      tensorT r(n,n);
      for (int axis=0; axis<3; axis++)
      {
//        vecfuncT dv = wst_diff(world,v,axis,_params.periodic);
        vecfuncT dv = diff(world,v,axis);
        r += matrix_inner(world, dv, dv, true);
        dv.clear(); world.gop.fence(); // Allow function memory to be freed
      }

      return r.scale(0.5);
  }


  /// Initializes alpha and beta mos, occupation numbers, eigenvalues
  void initial_guess()
  {
    if (_world.rank() == 0) print("Guessing rho ...\n\n");
    functionT rho = factoryT(_world).functor(functorT(
        new GuessDensity(_mentity, _aobasis)));

//    {
//      rho.reconstruct();
//      if (_world.rank() == 0)  printf("\n");
//      double L = _params.L;
//      double bstep = L / 100.0;
//      for (int i = 0; i < 101; i++)
//      {
//        coordT p(-L / 2 + i * bstep);
//        if (_world.rank() == 0)
//          printf("%.2f\t\t%.8f\n", p[0], rho(p));
//      }
//    }

    functionT vlocal;
    if (_params.nelec > 1)
    {
      if (_world.rank() == 0) print("Creating Coulomb op ...\n\n");
      SeparatedConvolution<double, 3>* op = 0;
      if (_params.periodic)
      {
        Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
        op = PeriodicCoulombOpPtr<double, 3> (_world, _params.waveorder,
            _params.lo, _params.thresh, cellsize);
      }
      else
      {
        op = CoulombOperatorPtr<double, 3> (_world, _params.waveorder,
            _params.lo, _params.thresh);
      }
      if (_world.rank() == 0) print("Building effective potential ...\n\n");
      vlocal = _vnuc + apply(*op, rho); //.scale(1.0-1.0/nel); // Reduce coulomb to increase binding
      rho.scale(0.5);
      vlocal = vlocal + make_lda_potential(_world, rho, rho, functionT(), functionT());
      delete op;
    }
    else
    {
      vlocal = _vnuc;
    }

//    {
//      //functionT Vxc = make_lda_potential(_world, rho, rho, functionT(), functionT());
//      if (_world.rank() == 0)  printf("\n");
//      double L = _params.L;
//      double bstep = L / 100.0;
//      rho.reconstruct();
//      _vnuc.reconstruct();
//      for (int i = 0; i < 101; i++)
//      {
//        coordT p(-L / 2 + i * bstep);
//        if (_world.rank() == 0)
//          printf("%.2f\t\t%.8f\t%.8f\n", p[0], rho(p), _vnuc(p));
//      }
//      if (_world.rank() == 0) printf("\n");
//    }

    rho.clear();
    vlocal.reconstruct();

    if (_world.rank() == 0) print("Projecting atomic orbitals ...\n\n");
    vecfuncT ao = project_ao_basis(_world);
//    if (_params.periodic)
//    {
//
//    }
//    else
//    {
      if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
      tensorT overlap = matrix_inner(_world, ao, ao, true);
      if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
      tensorT kinetic = kinetic_energy_matrix(_world, ao);

      reconstruct(_world, ao);
      //      vecfuncT vpsi = mul(world, vlocal, ao);
      if (_world.rank() == 0) print("Building potential energy matrix ...\n\n");
      vecfuncT vpsi = mul_sparse(_world, vlocal, ao, _params.thresh);
      _world.gop.fence();
      _world.gop.fence();
      compress(_world, vpsi);
      truncate(_world, vpsi);
      compress(_world, ao);

      tensorT potential = matrix_inner(_world, vpsi, ao, true);
      _world.gop.fence();
      vpsi.clear();
      _world.gop.fence();

      if (_world.rank() == 0) print("Constructing Fock matrix ...\n\n");
      tensorT fock = kinetic + potential;
      fock = 0.5 * (fock + transpose(fock));

      for (int fi = 0; fi < fock.dim[0]; fi++)
      {
        for (int fj = 0; fj < fock.dim[1]; fj++)
        {
          printf("%10.5f", fock(fi,fj));
        }
        printf("\n");
      }

      tensorT c, e;
      if (_world.rank() == 0) print("Diagonlizing Fock matrix ...\n\n");
      sygv(fock, overlap, 1, &c, &e);

      if (_world.rank() == 0)
      {
        print("initial eigenvalues");
        print(e);
      }

      compress(_world, ao);
      _world.gop.fence();
      _orbitals = transform(_world, ao, c(_, Slice(0, _params.nbands - 1)));
      _world.gop.fence();
      truncate(_world, _orbitals);
      normalize(_world, _orbitals);
      if (_world.rank() == 0)
        print("Analysis of initial alpha MO vectors");
      //      analyze_vectors(world, amo);

      _eigs = e(Slice(0, _params.nbands - 1));

      _occs = tensorT(_params.nbands);
      for (int i = 0; i < _params.nbands; i++)
        _occs[i] = _params.maxocc;
//    }
  }

  vecfuncT orbitals()
  {
    return _orbitals;
  }

  tensorT eigs()
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

  functionT vnucrhon()
  {
    return _vnucrhon;
  }

private:
  World& _world;
  MolecularEntity _mentity;
  AtomicBasisSet _aobasis;
  ElectronicStructureParams _params;
  functionT _vnuc;
  functionT _vnucrhon;
  vecfuncT _orbitals;
  tensorT _eigs;
  tensorT _occs;
};

#endif /* ELECTRONICSTRUCTUREAPP_H_ */
