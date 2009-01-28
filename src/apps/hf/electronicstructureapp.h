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
    if (_params.ispotential) // potential
    {
      _vnucrhon = rfactoryT(_world).functor(rfunctorT(new MolecularPotentialFunctor(_mentity))).thresh(_params.thresh * 0.1).truncate_on_project();
      _vnuc = copy(_vnucrhon);
      _vnuc.reconstruct();
    }
    else // charge density
    {
      _vnucrhon = rfactoryT(_world).functor(
          rfunctorT(new MolecularNuclearChargeDensityFunctor(_mentity))).
          thresh(_params.thresh * 0.1).initial_level(6).truncate_on_project();
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
              _params.L, false));
//             if (world.rank() == 0) {
//                 aobasis.get_atomic_basis_function(molecule,i).print_me(cout);
//             }
          ao[i] = rfactoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
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
                  rfunctorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
                      _params.L, false));
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

      return ao;
  }

  tensorT kinetic_energy_matrix(World& world, const rvecfuncT& v, const coordT k = coordT(0.0))
  {
      reconstruct(world, v);
      int n = v.size();
      rtensorT r(n,n);
      for (int axis=0; axis<3; axis++)
      {
//        rvecfuncT dv = wst_diff(world,v,axis,_params.periodic);
        rvecfuncT dv = diff(world,v,axis);
        r += matrix_inner(world, dv, dv, true);
        dv.clear(); world.gop.fence(); // Allow function memory to be freed
      }
      tensorT c(n,n);
      tensor_real2complex<double>(r.scale(0.5),c);
      return c;
  }


  /// Initializes alpha and beta mos, occupation numbers, eigenvalues
  void initial_guess()
  {
    if (_world.rank() == 0) print("Guessing rho ...\n\n");
    rfunctionT rho = rfactoryT(_world).functor(rfunctorT(
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

    rfunctionT vlocal;
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
      vlocal = vlocal + make_lda_potential(_world, rho, rho, rfunctionT(), rfunctionT());
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
    rvecfuncT ao = project_ao_basis(_world);
    if (_params.periodic)
    {

    }
    else
    {
      if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
      rtensorT roverlap = matrix_inner(_world, ao, ao, true);
      tensorT overlap(roverlap.dim[0], roverlap.dim[1]);
      tensor_real2complex<double>(roverlap,overlap);

      if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
      tensorT kinetic = kinetic_energy_matrix(_world, ao);

      reconstruct(_world, ao);
      if (_world.rank() == 0) print("Building potential energy matrix ...\n\n");
      rvecfuncT vpsi = mul_sparse(_world, vlocal, ao, _params.thresh);
      _world.gop.fence();
      _world.gop.fence();
      compress(_world, vpsi);
      truncate(_world, vpsi);
      compress(_world, ao);

      rtensorT rpotential = matrix_inner(_world, vpsi, ao, true);
      tensorT potential(rpotential.dim[0], rpotential.dim[1]);
      tensor_real2complex<double>(rpotential,potential);
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
          printf("%10.5f", abs(fock(fi,fj)));
        }
        printf("\n");
      }

      tensorT c, e;
      if (_world.rank() == 0) print("Diagonlizing Fock matrix ...\n\n");
//      sygv(fock, overlap, 1, &c, &e);

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

      _eigs = tensor_real(e(Slice(0, _params.nbands - 1)));

      _occs = rtensorT(_params.nbands);
      for (int i = 0; i < _params.nbands; i++)
        _occs[i] = _params.maxocc;
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
  std::vector<coordT> _kpoints;
};

#endif /* ELECTRONICSTRUCTUREAPP_H_ */
