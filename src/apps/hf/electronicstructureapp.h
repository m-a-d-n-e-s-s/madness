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
#include "esolver.h"

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
    const std::vector<coordT> _specialpts;
public:
    MolecularNuclearChargeDensityFunctor(const MolecularEntity& mentity, const double& R,
        const bool& periodic, const std::vector<coordT>& specialpts)
      : _mentity(mentity), R(R), periodic(periodic), _specialpts(specialpts) {
    }

    virtual std::vector<coordT> special_points() const
    {
      return _specialpts;
    }

    virtual Level special_level()
    {
      return 10;
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
  std::vector<coordT> _specialpts;
public:
  AtomicBasisFunctor(const AtomicBasisFunction& aofunc, double R, 
                     bool periodic)
  : aofunc(aofunc), R(R), periodic(periodic) 
  {
    double x, y, z;
    aofunc.get_coords(x,y,z);
    coordT r;
    r[0]=x; r[1]=y; r[2]=z;
    _specialpts=vector<coordT>(1,r);
  }

  virtual std::vector<coordT> special_points() const
  {
    return _specialpts;
  }

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
      //FunctionDefaults<3>::set_cubic_cell(0,_params.L);
      FunctionDefaults<3>::set_cubic_cell(-_params.L/2,_params.L/2);
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
      // set number of electrons to the total nuclear charge of the mentity
      _params.nelec = _mentity.total_nuclear_charge();
      // total number of bands include empty
      _params.nbands = (_params.nelec/2) + _params.nempty;
      if ((_params.nelec % 2) == 1) _params.nelec++;

//      if (_params.periodic && _params.kpoints)
//        _kpoints = read_kpoints("KPOINTS.OUT");
//      else
//        _kpoints.push_back(KPoint(coordT(0.0), 1.0));
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
          rfunctorT(new MolecularNuclearChargeDensityFunctor(_mentity, _params.L, _params.periodic, specialpts))).
          thresh(_params.thresh).initial_level(6).truncate_on_project();
      
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

      Level initial_level = 2;
      for (int i=0; i < _aobasis.nbf(_mentity); i++) {
          rfunctorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
              _params.L, _params.periodic));
          ao[i] = rfactoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
      }
      world.gop.fence();

      vector<double> norms;

//      while (1) {
//          norms = norm2(world, ao);
//          initial_level += 2;
//          if (initial_level >= 11) throw "project_ao_basis: projection failed?";
//          int nredone = 0;
//          for (int i=0; i<_aobasis.nbf(_mentity); i++) {
//              if (norms[i] < 0.5) {
//                  nredone++;
//                  if (world.rank() == 0) print("re-projecting ao basis function", i,"at level",initial_level);
//                  rfunctorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
//                      _params.L, _params.periodic));
//                  ao[i] = rfactoryT(world).functor(aofunc).initial_level(6).truncate_on_project().nofence();
//              }
//          }
//          world.gop.fence();
//          if (nredone == 0) break;
//      }

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

  /// Initializes alpha and beta mos, occupation numbers, eigenvalues
  void initial_guess()
  {
    // Get initial guess for the electronic density
    if (_world.rank() == 0) print("Guessing rho ...\n\n");
    rfunctionT rho = rfactoryT(_world).functor(rfunctorT(
        new GuessDensity(_mentity, _aobasis, _params.L, _params.periodic)));
    rho.scale(_params.nelec/rho.trace());

    print(rho.trace());

//    vector<long> npt(3,101);
//    plotdx(rho, "rho_initial.dx", FunctionDefaults<3>::get_cell(), npt);

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
      vlocal = _vnuc + vc; //.scale(1.0-1.0/nel); // Reduce coulomb to increase binding
      rho.scale(0.5);
      // Do the LDA
      rfunctionT vlda = make_lda_potential(_world, rho, rho, rfunctionT(), rfunctionT());
      vlocal = vlocal + vlda;
      delete op;
//      vector<long> npt(3,101);
//      plotdx(vc, "vc.dx", FunctionDefaults<3>::get_cell(), npt);
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

//    for (unsigned int ai = 0; ai < ao.size(); ai++)
//    {
//      std::ostringstream strm;
//      strm << "aod" << ai << ".dx";
//      std::string fname = strm.str();
//      vector<long> npt(3,101);
//      plotdx(ao[ai], fname.c_str(), FunctionDefaults<3>::get_cell(), npt);
//    }

    // Get size information from k-points and ao_basis so that we can correctly size
    // the _orbitals data structure and the eigs tensor
    // number of orbitals in the basis set
    int nao = ao.size();
    // number of kpoints
    int nkpts = _kpoints.size();
    // total number of orbitals to be processed (no symmetry)
    int norbs = _params.nbands * nkpts;
    // Check to see if the basis set can accomodate the number of bands
    if (_params.nbands > nao)
      MADNESS_EXCEPTION("Error: basis not large enough to accomodate number of bands", 0);
    // set the number of orbitals
    _eigs = std::vector<double>(norbs, 0.0);
    _occs = std::vector<double>(norbs, 0.0);
    // Build the overlap matrix
    if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
    rtensorT roverlap = matrix_inner(_world, ao, ao, true);
    // Convert to a complex tensor
    ctensorT overlap = tensor_real2complex<double>(roverlap);
    // Build the potential matrix
    reconstruct(_world, ao);
    if (_world.rank() == 0) print("Building potential energy matrix ...\n\n");
    rvecfuncT vpsi = mul_sparse(_world, vlocal, ao, _params.thresh);
    // I don't know why fence is called twice here
    _world.gop.fence();
    _world.gop.fence();
    compress(_world, vpsi);
    truncate(_world, vpsi);
    compress(_world, ao);
    // Build the potential matrix
    rtensorT rpotential = matrix_inner(_world, vpsi, ao, true);
    // Convert to a complex tensor
    ctensorT potential = tensor_real2complex<double>(rpotential);
    _world.gop.fence();
    // free memory
    vpsi.clear();
    _world.gop.fence();

    int kp = 0;
    if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
    // Need to do kinetic piece for every k-point
    for (int ki = 0; ki < nkpts; ki++)
    {
      // Set occupation numbers
      if (_params.spinpol)
      {
        MADNESS_EXCEPTION("spin polarized not implemented", 0);
      }
      else
      {
        int filledbands = _params.nelec / 2;
        int occstart = kp;
        int occend = kp + filledbands;
        for (int i = occstart; i < occend; i++) _occs[i] = _params.maxocc;
        if ((_params.nelec % 2) == 1)
          _occs[occend] = 1.0;
      }
      // Get k-point from list
      KPoint kpt = _kpoints[ki];
      // Build kinetic matrx
      ctensorT kinetic = ::kinetic_energy_matrix(_world, ao, _params.periodic, kpt);
      // Construct and diagonlize Fock matrix
      ctensorT fock = potential + kinetic;
      fock = 0.5 * (fock + transpose(fock));
      ctensorT c; rtensorT e;
      sygv(fock, overlap, 1, &c, &e);

      compress(_world, ao);
      _world.gop.fence();
      // Take linear combinations of the gaussian basis orbitals as the starting
      // orbitals for solver
      vecfuncT tmp_orbitals = transform(_world, ao, c(_, Slice(0, nao - 1)));
      _world.gop.fence();
      truncate(_world, tmp_orbitals);
      normalize(_world, tmp_orbitals);
      rtensorT tmp_eigs = e(Slice(0, nao - 1));

      if (_world.rank() == 0) printf("(%8.4f,%8.4f,%8.4f)\n",kpt.k[0], kpt.k[1], kpt.k[2]);
      if (_world.rank() == 0) print(tmp_eigs);
      if (_world.rank() == 0) print("\n");

      // DEBUG
//      for (int i = 0; i < kinetic.dim[0]; i++)
//      {
//        for (int j = 0; j < kinetic.dim[1]; j++)
//        {
//          if (_world.rank() == 0) printf("%10.5f", real(kinetic(i,j)));
//        }
//        if (_world.rank() == 0) printf("\n");
//      }
//      if (_world.rank() == 0) printf("\n");
//      if (_world.rank() == 0) printf("\n");
//
//      for (int i = 0; i < potential.dim[0]; i++)
//      {
//        for (int j = 0; j < potential.dim[1]; j++)
//        {
//          if (_world.rank() == 0) printf("%10.5f", real(potential(i,j)));
//        }
//        if (_world.rank() == 0) printf("\n");
//      }
//      if (_world.rank() == 0) printf("\n");
//      if (_world.rank() == 0) printf("\n");
//
//      for (int i = 0; i < fock.dim[0]; i++)
//      {
//        for (int j = 0; j < fock.dim[1]; j++)
//        {
//          if (_world.rank() == 0) printf("%10.5f", real(fock(i,j)));
//        }
//        if (_world.rank() == 0) printf("\n");
//      }
//      if (_world.rank() == 0) printf("\n");
//      if (_world.rank() == 0) printf("\n");

      // Fill in orbitals and eigenvalues
      int kend = kp + _params.nbands;
      _kpoints[kp].begin = kp;
      _kpoints[kp].end = kend;
      for (int oi = kp, ti = 0; oi < kend; oi++, ti++)
      {
        if (_world.rank() == 0) print(oi, ti, kp, kend);
        _orbitals.push_back(tmp_orbitals[ti]);
        _eigs[oi] = tmp_eigs[ti];
      }

      kp += _params.nbands;
    }
  }

  vecfuncT orbitals()
  {
    return _orbitals;
  }

  std::vector<double> eigs()
  {
    return _eigs;
  }

  std::vector<KPoint> kpoints()
  {
    return _kpoints;
  }

  std::vector<double> occs()
  {
    return _occs;
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
  std::vector<double> _eigs;
  std::vector<double> _occs;
  std::vector<KPoint> _kpoints;
};

#endif /* ELECTRONICSTRUCTUREAPP_H_ */
