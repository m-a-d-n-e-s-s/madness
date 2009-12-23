#include <mra/mra.h>
#include <world/world.h>
#include <linalg/solvers.h>
#include <vector>
#include <fortran_ctypes.h>
#include <cmath>

#include "poperator.h"
#include "libxc.h"
#include "electronicstructureparams.h"
#include "complexfun.h"
#include "esolver.h"

#ifndef SOLVER_H_

/*! 
  \ingroup applications
  @{
  \defgroup periodic_solver Periodic Solver
  \brief The Periodic Solver group is a group that contains the software 
  objects that are needed to solve a periodic Kohn-Sham hamiltonian.
  @{
*/

//*****************************************************************************
static double onesfunc(const coordT& x)
{
  return 1.0;
}
//*****************************************************************************

namespace madness
{
  /*!
  
   \brief The SubspaceK class is a container class holding previous orbitals 
   and residuals.
   \par
   The Solver class uses the Krylov Accelerated Inexact Newton Solver (KAIN) 
   accelerate the convergence a given calculation. The KAIN solver needs to 
   store a subspace of previous orbitals and residuals. In the case is this 
   implementation, the orbitals are store according to which k-point to which 
   they belong.
  */
  
  //***************************************************************************
  template <typename T, int NDIM>
  class SubspaceK
  {
    // Typedef's
    typedef std::complex<T> valueT;
    typedef Function<valueT,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef std::vector<pairvecfuncT> subspaceT;
    typedef Tensor<valueT> tensorT;
    typedef std::vector<tensorT> vectensorT;
    typedef std::vector<subspaceT> vecsubspaceT;

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    vectensorT _Q;
    //*************************************************************************

    //*************************************************************************
    vecsubspaceT _subspace;
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> _kpoints;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

    //*************************************************************************
    double _residual;
    //*************************************************************************

  public:

    //*************************************************************************
    SubspaceK(World& world, const ElectronicStructureParams& params,
            const std::vector<KPoint>& kpoints) : _world(world), _kpoints(kpoints),
            _params(params)
    {
      _residual = 1e6;
      for (unsigned int kp = 0; kp < _kpoints.size(); kp++)
      {
        _Q.push_back(tensorT(1,1));
        _subspace.push_back(subspaceT());
      }
    }
    //*************************************************************************

    //*************************************************************************
    void update_subspace(vecfuncT& awfs_new,
                         vecfuncT& bwfs_new,
                         const vecfuncT& awfs_old,
                         const vecfuncT& bwfs_old)
    {
      // compute residuals (total)
      vecfuncT t_rm = sub(_world, awfs_old, awfs_new);
      if (_params.spinpol)
      {
        vecfuncT br = sub(_world, bwfs_old, bwfs_new);
        t_rm.insert(t_rm.end(), br.begin(), br.end());
      }
      std::vector<double> rnvec = norm2<valueT,NDIM>(_world, t_rm);
      if (_world.rank() == 0)
      {
        double rnorm = 0.0;
        for (unsigned int i = 0; i < rnvec.size(); i++) rnorm += rnvec[i];
        if (_world.rank() == 0) print("residual = ", rnorm);
        _residual = rnorm;
      }
      _world.gop.broadcast(_residual, 0);

      for (unsigned int kp = 0; kp < _kpoints.size(); kp++)
      {
        KPoint kpoint = _kpoints[kp];

        vecfuncT k_phisa(awfs_old.begin() + kpoint.begin, awfs_old.begin() + kpoint.end);
        vecfuncT k_phisb(bwfs_old.begin() + kpoint.begin, bwfs_old.begin() + kpoint.end);
        vecfuncT k_awfs(awfs_new.begin() + kpoint.begin, awfs_new.begin() + kpoint.end);
        vecfuncT k_bwfs(bwfs_new.begin() + kpoint.begin, bwfs_new.begin() + kpoint.end);

        // compute residuals for k-point
        // concatentate up and down spins
        vecfuncT k_rm = sub(_world, k_phisa, k_awfs);
        vecfuncT k_vm(k_phisa);
        if (_params.spinpol)
        {
          vecfuncT k_br = sub(_world, k_phisb, k_bwfs);
          k_rm.insert(k_rm.end(), k_br.begin(), k_br.end());
          k_vm.insert(k_vm.end(), k_phisb.begin(), k_phisb.end());
        }

        // Update subspace and matrix Q
        compress(_world, k_vm, false);
        compress(_world, k_rm, false);
        _world.gop.fence();
        subspaceT k_subspace = _subspace[kp];
        k_subspace.push_back(pairvecfuncT(k_vm,k_rm));

        int m = k_subspace.size();
        tensorT ms(m);
        tensorT sm(m);
        for (int s = 0; s < m; s++)
        {
            const vecfuncT& k_vs = k_subspace[s].first;
            const vecfuncT& k_rs = k_subspace[s].second;
            for (unsigned int i = 0; i < k_vm.size(); i++)
            {
                ms[s] += k_vm[i].inner_local(k_rs[i]);
                sm[s] += k_vs[i].inner_local(k_rm[i]);
            }
        }
        _world.gop.sum(ms.ptr(),m);
        _world.gop.sum(sm.ptr(),m);

        tensorT newQ(m,m);
        if (m > 1) newQ(Slice(0,-2),Slice(0,-2)) = _Q[kp];
        newQ(m-1,_) = ms;
        newQ(_,m-1) = sm;

        _Q[kp] = newQ;
        if (_world.rank() == 0) print(_Q[kp]);

        // Solve the subspace equations
        tensorT c;
        if (_world.rank() == 0) {
            double rcond = 1e-12;
            while (1) {
                c = KAIN(_Q[kp],rcond);
                if (abs(c[m-1]) < 3.0) {
                    break;
                }
                else if (rcond < 0.01) {
                    if (_world.rank() == 0)
                      print("Increasing subspace singular value threshold ", c[m-1], rcond);
                    rcond *= 100;
                }
                else {
                    if (_world.rank() == 0)
                      print("Forcing full step due to subspace malfunction");
                    c = 0.0;
                    c[m-1] = 1.0;
                    break;
                }
            }
        }

        _world.gop.broadcast_serializable(c, 0);
        if (_world.rank() == 0) {
            //print("Subspace matrix");
            //print(Q);
            print("Subspace solution", c);
        }

        // Form linear combination for new solution
        vecfuncT k_phisa_new = zero_functions<valueT,NDIM>(_world, k_phisa.size());
        vecfuncT k_phisb_new = zero_functions<valueT,NDIM>(_world, k_phisb.size());
        compress(_world, k_phisa_new, false);
        compress(_world, k_phisb_new, false);
        _world.gop.fence();
        std::complex<double> one = std::complex<double>(1.0,0.0);
        unsigned int norbitals = awfs_old.size() / _kpoints.size();
        for (unsigned int m = 0; m < k_subspace.size(); m++)
        {
            const vecfuncT& k_vm = k_subspace[m].first;
            const vecfuncT& k_rm = k_subspace[m].second;
            // WSTHORNTON Stopped here!
            const vecfuncT  vma(k_vm.begin(),k_vm.begin() + norbitals);
            const vecfuncT  rma(k_rm.begin(),k_rm.begin() + norbitals);
            const vecfuncT  vmb(k_vm.end() - norbitals, k_vm.end());
            const vecfuncT  rmb(k_rm.end() - norbitals, k_rm.end());

            gaxpy(_world, one, k_phisa_new, c(m), vma, false);
            gaxpy(_world, one, k_phisa_new,-c(m), rma, false);
            gaxpy(_world, one, k_phisb_new, c(m), vmb, false);
            gaxpy(_world, one, k_phisb_new,-c(m), rmb, false);
        }
        _world.gop.fence();

        if (_params.maxsub <= 1) {
            // Clear subspace if it is not being used
            k_subspace.clear();
        }
        else if (k_subspace.size() == _params.maxsub) {
            // Truncate subspace in preparation for next iteration
            k_subspace.erase(k_subspace.begin());
            _Q[kp] = _Q[kp](Slice(1,-1),Slice(1,-1));
        }
        // Save subspace
        _subspace[kp] = k_subspace;

        for (unsigned int wi = kpoint.begin, fi = 0; wi < kpoint.end;
          wi++, fi++)
        {
          awfs_new[wi] = k_phisa_new[fi];
          bwfs_new[wi] = k_phisb_new[fi];
        }
      }
    }
    //*************************************************************************

  };

  /*!
   \ingroup periodic_solver
  
   \brief The SubspaceK class is a container class holding previous orbitals 
   and residuals.
   \par
   The Solver class uses the Krylov Accelerated Inexact Newton Solver (KAIN) 
   accelerate the convergence a given calculation. The KAIN solver needs to 
   store a subspace of previous orbitals and residuals. 
  */

  //***************************************************************************
  template <typename T, int NDIM>
  class Subspace
  {
    // Typedef's
    typedef std::complex<T> valueT;
    typedef Function<valueT,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef std::vector<pairvecfuncT> subspaceT;
    typedef Tensor<valueT> tensorT;

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    tensorT _Q;
    //*************************************************************************

    //*************************************************************************
    subspaceT _subspace;
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> _kpoints;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

    //*************************************************************************
    double _residual;
    //*************************************************************************

  public:

    //*************************************************************************
    Subspace(World& world, const ElectronicStructureParams& params)
      : _world(world), _params(params)
    {
      _residual = 1e6;
    }
    //*************************************************************************

    //*************************************************************************
    void update_subspace(vecfuncT& awfs_new,
                         vecfuncT& bwfs_new,
                         const vecfuncT& awfs_old,
                         const vecfuncT& bwfs_old,
                         const vecfuncT& rm)
    {
      // concatentate up and down spins
      vecfuncT vm = awfs_old;
      if (_params.spinpol)
      {
        vm.insert(vm.end(), bwfs_old.begin(), bwfs_old.end());
      }

      // Update subspace and matrix Q
      compress(_world, vm, false);
      compress(_world, rm, false);
      _world.gop.fence();
      _subspace.push_back(pairvecfuncT(vm,rm));

      int m = _subspace.size();
      tensorT ms(m);
      tensorT sm(m);
      for (int s=0; s<m; s++)
      {
          const vecfuncT& vs = _subspace[s].first;
          const vecfuncT& rs = _subspace[s].second;
          for (unsigned int i=0; i<vm.size(); i++)
          {
              ms[s] += vm[i].inner_local(rs[i]);
              sm[s] += vs[i].inner_local(rm[i]);
          }
      }
      _world.gop.sum(ms.ptr(),m);
      _world.gop.sum(sm.ptr(),m);

      tensorT newQ(m,m);
      if (m > 1) newQ(Slice(0,-2),Slice(0,-2)) = _Q;
      newQ(m-1,_) = ms;
      newQ(_,m-1) = sm;

      _Q = newQ;
      if (_world.rank() == 0) print(_Q);

      // Solve the subspace equations
      tensorT c;
      if (_world.rank() == 0) {
          double rcond = 1e-12;
          while (1) {
              c = KAIN(_Q,rcond);
              if (abs(c[m-1]) < 3.0) {
                  break;
              }
              else if (rcond < 0.01) {
                  if (_world.rank() == 0)
                    print("Increasing subspace singular value threshold ", c[m-1], rcond);
                  rcond *= 100;
              }
              else {
                  if (_world.rank() == 0)
                    print("Forcing full step due to subspace malfunction");
                  c = 0.0;
                  c[m-1] = 1.0;
                  break;
              }
          }
      }

      _world.gop.broadcast_serializable(c, 0);
      if (_world.rank() == 0) {
          //print("Subspace matrix");
          //print(Q);
          print("Subspace solution", c);
      }

      // Form linear combination for new solution
      vecfuncT phisa_new = zero_functions<valueT,NDIM>(_world, awfs_old.size());
      vecfuncT phisb_new = zero_functions<valueT,NDIM>(_world, bwfs_old.size());
      compress(_world, phisa_new, false);
      compress(_world, phisb_new, false);
      _world.gop.fence();
      std::complex<double> one = std::complex<double>(1.0,0.0);
      for (unsigned int m=0; m<_subspace.size(); m++) {
          const vecfuncT& vm = _subspace[m].first;
          const vecfuncT& rm = _subspace[m].second;
          const vecfuncT  vma(vm.begin(),vm.begin()+awfs_old.size());
          const vecfuncT  rma(rm.begin(),rm.begin()+awfs_old.size());
          const vecfuncT  vmb(vm.end()-bwfs_old.size(), vm.end());
          const vecfuncT  rmb(rm.end()-bwfs_old.size(), rm.end());

          gaxpy(_world, one, phisa_new, c(m), vma, false);
          gaxpy(_world, one, phisa_new,-c(m), rma, false);
          gaxpy(_world, one, phisb_new, c(m), vmb, false);
          gaxpy(_world, one, phisb_new,-c(m), rmb, false);
      }
      _world.gop.fence();

      if (_params.maxsub <= 1) {
          // Clear subspace if it is not being used
          _subspace.clear();
      }
      else if (_subspace.size() == _params.maxsub) {
          // Truncate subspace in preparation for next iteration
          _subspace.erase(_subspace.begin());
          _Q = _Q(Slice(1,-1),Slice(1,-1));
      }
      awfs_new = phisa_new;
      bwfs_new = phisb_new;
    }
    //*************************************************************************
  };
  //***************************************************************************

  /*! \ingroup periodic_solver
      \brief The main class of the periodic DFT solver
      \f[
      z = frac{x}{1 - y^2}
      \f]
  */

  //***************************************************************************
  template <typename T, int NDIM>
  class Solver
  {
    // Typedef's
    typedef std::complex<T> valueT;
    typedef Function<T,NDIM> rfuntionT;
    typedef FunctionFactory<T,NDIM> rfactoryT;
    typedef Function<valueT,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef FunctionFactory<valueT,NDIM> factoryT;
    typedef Vector<double,NDIM> kvecT;
    typedef SeparatedConvolution<T,3> operatorT;
    typedef SharedPtr<operatorT> poperatorT;
    typedef Tensor<double> rtensorT;
    typedef Tensor<std::complex<double> > ctensorT;
    typedef Tensor<valueT> tensorT;
    typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef std::vector<pairvecfuncT> subspaceT;
    typedef std::vector<tensorT> vectensorT;
    typedef std::vector<subspaceT> vecsubspaceT;

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    // This variable could either be a nuclear potiential or a nuclear charge
    // density depending on the "ispotential" boolean variable in the
    // ElectronicStructureParams class.
    rfuntionT _vnucrhon;
    //*************************************************************************

    //*************************************************************************
    vecfuncT _phisa;
    //*************************************************************************

    //*************************************************************************
    vecfuncT _phisb;
    //*************************************************************************

    //*************************************************************************
    std::vector<T> _eigsa;
    //*************************************************************************

    //*************************************************************************
    std::vector<T> _eigsb;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> _kpoints;
    //*************************************************************************

    //*************************************************************************
    std::vector<double> _occs;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _rhoa;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _rhob;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _rho;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _vnuc;
    //*************************************************************************

    //*************************************************************************
    SeparatedConvolution<T,NDIM>* _cop;
    //*************************************************************************

//    //*************************************************************************
//    vecsubspaceT _subspace;
//    //*************************************************************************
//
//    //*************************************************************************
//    vectensorT _Q;
//    //*************************************************************************

    //*************************************************************************
    Subspace<T,NDIM>* _subspace;
    //*************************************************************************

    //*************************************************************************
    MolecularEntity _mentity;
    //*************************************************************************

    //*************************************************************************
    double _residual;
    //*************************************************************************

    //*************************************************************************
    AtomicBasisSet _aobasis;
    //*************************************************************************

    //*************************************************************************
    double _maxthresh;
    //*************************************************************************
  public:

    //*************************************************************************
    double ttt, sss;
    void START_TIMER(World& world) {
        world.gop.fence(); ttt=wall_time(); sss=cpu_time();
    }

    void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
    }
    //*************************************************************************
   
    //*************************************************************************
    Solver(World& world, const std::string& filename) : _world(world)
    {
      init(filename);
      _residual = 1e5;
      make_nuclear_potential();
      initial_guess();
//      for (unsigned int kp = 0; kp < _kpoints.size(); kp++)
//      {
//        _Q.push_back(tensorT(1,1));
//        _subspace.push_back(subspaceT());
//      }
      _subspace = new Subspace<T,NDIM>(world, _params);
    }
    //*************************************************************************

    //*************************************************************************
    void init(const std::string& filename)
    {
      // params
      if (_world.rank() == 0)
      {
        _params.read_file(filename);
        //_params.fractional = false;
      }
      // Send params
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
        _aobasis.read_file(_params.basis);
        _mentity.read_file(filename, _params.fractional);
        _mentity.center();
      }
      // Send mentity and aobasis
      _world.gop.broadcast_serializable(_mentity, 0);
      _world.gop.broadcast_serializable(_aobasis, 0);
      // set number of electrons to the total nuclear charge of the mentity
      _params.nelec = _mentity.total_nuclear_charge();
      // total number of bands include empty
      _params.nbands = (_params.nelec/2) + _params.nempty;
      if ((_params.nelec % 2) == 1) _params.nelec++;

      // kmesh
      if (_params.restart == 0)
      {
        if (_params.periodic) // PERIODIC
        {
          // GAMMA POINT
          if ((_params.ngridk0 == 1) && (_params.ngridk1 == 1) && (_params.ngridk2 == 1))
          {
            _kpoints.push_back(KPoint(coordT(0.0), 1.0));
          }
          if ((_params.ngridk0 == 0) && (_params.ngridk1 == 0) && (_params.ngridk2 == 0))
          {
            _kpoints.push_back(KPoint(coordT(0.0), 0.5));
            _kpoints.push_back(KPoint(coordT(0.0), 0.5));
          }
          else // NORMAL BANDSTRUCTURE
          {
            _kpoints = genkmesh(_params.ngridk0, _params.ngridk1,
                                _params.ngridk2, _params.koffset0,
                                _params.koffset1, _params.koffset2,
                                _params.L);
          }
        }
        else // NOT-PERIODIC
        {
          _kpoints.push_back(KPoint(coordT(0.0), 1.0));
        }
        if (_world.rank() == 0) print("kmesh:");
        for (unsigned int i = 0; i < _kpoints.size() && _world.rank() == 0; i++)
        {
          KPoint kpoint = _kpoints[i];
          print(kpoint.k[0], kpoint.k[1], kpoint.k[2], kpoint.weight);
        }
      }
      else
      {
        load_orbitals();
      }
    }
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> genkmesh(unsigned int ngridk0, unsigned ngridk1, unsigned int ngridk2, 
                                 double koffset0, double koffset1, double koffset2, double R)
    {
      std::vector<KPoint> kmesh;
      double step0 = 1.0/ngridk0;
      double step1 = 1.0/ngridk1;
      double step2 = 1.0/ngridk2;
      double weight = 1.0/(ngridk0*ngridk1*ngridk2);
      double TWO_PI = 2.0 * madness::constants::pi;
      for (unsigned int i = 0; i < ngridk0; i++)
      {
        for (unsigned int j = 0; j < ngridk1; j++)
        {
          for (unsigned int k = 0; k < ngridk2; k++)
          {
            //double k0 = (i*step0 - step0/2) * TWO_PI/R;
            //double k1 = (j*step1 - step1/2) * TWO_PI/R;
            //double k2 = (k*step2 - step2/2) * TWO_PI/R;
            double k0 = ((i*step0)+koffset0) * TWO_PI/R;
            double k1 = ((j*step1)+koffset1) * TWO_PI/R;
            double k2 = ((k*step2)+koffset2) * TWO_PI/R;
            KPoint kpoint(k0, k1, k2, weight);
            kmesh.push_back(kpoint);
          }
        }
      }
      return kmesh;
    }
    //*************************************************************************

    //*************************************************************************
    void save_orbitals()
    {
      ParallelOutputArchive ar(_world, "orbitals", _params.nio);
      ar & _params.spinpol;
      ar & (unsigned int)(_kpoints.size());
      for (unsigned int i = 0; i < _kpoints.size(); i++) ar & _kpoints[i];
      ar & (unsigned int)(_occs.size());
      for (unsigned int i = 0; i < _occs.size(); i++) ar & _occs[i];
      ar & (unsigned int)(_phisa.size());
      for (unsigned int i = 0; i < _phisa.size(); i++) ar & _phisa[i];
      if (_params.spinpol)
      {
        MADNESS_EXCEPTION("spin polarized is not implemented", 0);
      }
    }
    //*************************************************************************

    //*************************************************************************
    void load_orbitals()
    {
      const double thresh = FunctionDefaults<3>::get_thresh();
      const int k = FunctionDefaults<3>::get_k();

      ParallelInputArchive ar(_world, "orbitals");

      // spin-polarized
      bool spinrest;
      ar & spinrest;
      // kpoints
      unsigned int nkpts;
      ar & nkpts;
      _kpoints.clear();
      for (unsigned int i = 0; i < nkpts; i++)
      {
        KPoint tkpt;
        ar & tkpt;
        _kpoints.push_back(tkpt);
      }
      // occs
      unsigned int noccs;
      ar & noccs;
      _occs.clear();
      for (unsigned int i = 0; i < noccs; i++)
      {
        double tocc;
        ar & tocc;
        _occs.push_back(tocc);
      }
      // orbitals
      unsigned int norbs;
      ar & norbs;
      _phisa.clear();
      _eigsa.clear();
      for (unsigned int i = 0; i < norbs; i++)
      {
        functionT tfunc;
        ar & tfunc;
        _phisa.push_back(tfunc);
        _eigsa.push_back(-0.1);
      }
      // check for k mismatch
      if (_phisa[0].k() != k)
      {
        reconstruct(_world,_phisa);
        for (unsigned int i = 0; i < _phisa.size(); i++)
          _phisa[i] = madness::project(_phisa[i], k, thresh, false);
        _world.gop.fence();
      }
      // orthonormalize
      for (unsigned int i = 0; i < _kpoints.size(); i++)
        gram_schmidt(_phisa, _kpoints[i]);
      if (_params.spinpol)
      {
        MADNESS_EXCEPTION("spin-polarized not implemented", 0);
      }
    }
    //*************************************************************************

    //*************************************************************************
    void make_nuclear_potential()
    {

      if (_params.periodic) // periodic
      {
        Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
        _cop = PeriodicCoulombOpPtr<double,3>(_world, _params.waveorder,_params.lo, _params.thresh * 0.1, cellsize);
      }
      else // not periodic
      {
        _cop = CoulombOperatorPtr<double>(_world, _params.waveorder,_params.lo, _params.thresh * 0.1);
      }

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
        double vnuctrace = _vnuc.trace();
        //if (_world.rank == 0) print("vnuc trace:     ", vnuctrace);
        _vnuc.reconstruct();
        vnuctrace = _vnuc.trace();
        //if (_world.rank == 0) print("vnuc trace:     ", vnuctrace);
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
          if (_world.rank() == 0) print("Special point: ", pt);
        }
        double now = wall_time();
        _vnucrhon = rfactoryT(_world).functor(
            rfunctorT(new MolecularNuclearChargeDensityFunctor(_mentity, _params.L, _params.periodic, specialpts))).
            thresh(_params.thresh).initial_level(6).truncate_on_project();

        if (_world.rank() == 0) printf("%f\n", wall_time() - now);
        if (_world.rank() == 0) print("calculating trace of rhon ..\n\n");
        double rtrace = _vnucrhon.trace();
        if (_world.rank() == 0) print("rhon trace = ", rtrace);
        now = wall_time();
        _vnucrhon.truncate();
        _vnuc = apply(*_cop, _vnucrhon);
        if (_world.rank() == 0) printf("%f\n", wall_time() - now);
        if (_world.rank() == 0) print("Done creating nuclear potential ..\n");
      }

      std::vector<long> npt(3,101);
      plotdx(_vnuc, "vnuc.dx", FunctionDefaults<3>::get_cell(), npt);
    }
    //*************************************************************************

    //*************************************************************************
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
    //*************************************************************************

    //*************************************************************************
    rfunctionT
    make_lda_potential(World& world,
                       const rfunctionT& arho,
                       const rfunctionT& brho,
                       const rfunctionT& adelrhosq,
                       const rfunctionT& bdelrhosq)
    {
  //      MADNESS_ASSERT(!_params.spinpol);
        rfunctionT vlda = copy(arho);
        vlda.unaryop(&::libxc_ldaop<double>);
        return vlda;
    }

    vecfuncT project_ao_basis(World& world, KPoint kpt) {
        vecfuncT ao(_aobasis.nbf(_mentity));

        Level initial_level = 3;
        for (int i=0; i < _aobasis.nbf(_mentity); i++) {
            functorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
                _params.L, _params.periodic, kpt));
            ao[i] = factoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
        }
        world.gop.fence();

        std::vector<double> norms;

        norms = norm2s(world, ao);

        for (int i=0; i<_aobasis.nbf(_mentity); i++) {
            if (world.rank() == 0 && fabs(norms[i]-1.0)>1e-3) print(i," bad ao norm?", norms[i]);
            norms[i] = 1.0/norms[i];
        }

        scale(world, ao, norms);

        norms = norm2s(world, ao);

        for (int i=0; i<_aobasis.nbf(_mentity); i++) {
            if (world.rank() == 0 && fabs(norms[i]-1.0)>1e-3) print(i," bad ao norm?", norms[i]);
            norms[i] = 1.0/norms[i];
        }

        scale(world, ao, norms);

        return ao;
    }
    //*************************************************************************

    //*************************************************************************
    // Constructor
    Solver(World& world, 
           rfuntionT vnucrhon,
           vecfuncT phisa,
           vecfuncT phisb,
           std::vector<T> eigsa,
           std::vector<T> eigsb,
           ElectronicStructureParams params,
           MolecularEntity mentity)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phisa), _phisb(phisb),
       _eigsa(eigsa), _eigsb(eigsb), _params(params), _mentity(mentity)
    {
      _residual = 1e5;
      if (params.periodic)
      {
        Tensor<double> box = FunctionDefaults<NDIM>::get_cell_width();
        _cop = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1, box);
      }
      else
      {
        _cop = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1);
      }

      if (params.ispotential)
      {
        _vnuc = copy(_vnucrhon);
      }
      else
      {
        _vnuc = apply(*_cop, _vnucrhon);
      }
    }
    //*************************************************************************

    /// Initializes alpha and beta mos, occupation numbers, eigenvalues
    //*************************************************************************
    void initial_guess()
    {
      // Get initial guess for the electronic density
      if (_world.rank() == 0) print("Guessing rho ...\n\n");
      rfunctionT rho = (_params.restart == 0) ? rfactoryT(_world).functor(rfunctorT(
          new GuessDensity(_mentity, _aobasis, _params.L, _params.periodic))).initial_level(3) :
            compute_rho(_phisa, _kpoints);
      if (_params.restart == 1)
      {
        if (_params.spinpol)
          rho += compute_rho(_phisb, _kpoints);
        else
          rho.scale(2.0);
      }

      // This is a cheat
      double rtrace = rho.trace();
      if (_world.rank() == 0) print("trace of rho = ", rtrace);
      rho.scale(_params.nelec/rho.trace());

      if (_params.restart != 1)
      {
        // build effective potential
        rfunctionT vlocal;
        // Is this a many-body system?
        int rank = _world.rank();
        if (_world.rank() == 0) print("rank ", rank, "nelec ", _params.nelec);
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
          double vctrace = vc.trace();
          double vldatrace = vlda.trace();
          double vlocaltrace = vlocal.trace();
          if (_world.rank() == 0) print("vctrace:     ", vctrace);
          if (_world.rank() == 0) print("vldatrace:     ", vldatrace);
          if (_world.rank() == 0) print("vlocaltrace:     ", vlocaltrace);
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

        // Get size information from k-points and ao_basis so that we can correctly size
        // the _orbitals data structure and the eigs tensor
        // number of orbitals in the basis set
        int nao = _aobasis.nbf(_mentity);
        // number of kpoints
        int nkpts = _kpoints.size();
        // total number of orbitals to be processed (no symmetry)
        int norbs = _params.nbands * nkpts;
        // Check to see if the basis set can accomodate the number of bands
        if (_params.nbands > nao)
          MADNESS_EXCEPTION("Error: basis not large enough to accomodate number of bands", 0);
        // set the number of orbitals
        _eigsa = std::vector<double>(norbs, 0.0);
        _eigsb = std::vector<double>(norbs, 0.0);
        _occs = std::vector<double>(norbs, 0.0);
        int kp = 0;
        if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
        // Need to do kinetic piece for every k-point
        for (int ki = 0; ki < nkpts; ki++)
        {
          // These are our initial basis functions
          if (_world.rank() == 0) print("Projecting atomic orbitals ...\n\n");
          cvecfuncT ao = project_ao_basis(_world, _kpoints[ki]);

      //    for (unsigned int ai = 0; ai < ao.size(); ai++)
      //    {
      //      std::ostringstream strm;
      //      strm << "aod" << ai << ".dx";
      //      std::string fname = strm.str();
      //      vector<long> npt(3,101);
      //      plotdx(ao[ai], fname.c_str(), FunctionDefaults<3>::get_cell(), npt);
      //    }

          // Get k-point from list
          KPoint& kpt = _kpoints[ki];
          // Build kinetic matrx
          ctensorT kinetic = ::kinetic_energy_matrix(_world, ao, _params.periodic, kpt);
          // Build the overlap matrix
          if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
          ctensorT overlap = matrix_inner(_world, ao, ao, true);
          // Build the potential matrix
          reconstruct(_world, ao);
          if (_world.rank() == 0) print("Building potential energy matrix ...\n\n");
          cvecfuncT vpsi = mul_sparse(_world, vlocal, ao, _params.thresh);
          // I don't know why fence is called twice here
          _world.gop.fence();
          _world.gop.fence();
          compress(_world, vpsi);
          truncate(_world, vpsi);
          compress(_world, ao);
          // Build the potential matrix
          ctensorT potential = matrix_inner(_world, vpsi, ao, true);
          _world.gop.fence();
          // free memory
          vpsi.clear();
          _world.gop.fence();

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
          // Construct and diagonlize Fock matrix
          ctensorT fock = potential + kinetic;
          fock = 0.5 * (fock + transpose(fock));
          ctensorT c; rtensorT e;
          sygv(fock, overlap, 1, &c, &e);

          ctensorT ck; rtensorT ek;
          sygv(kinetic, overlap, 1, &ck, &ek);

          ctensorT cp; rtensorT ep;
          sygv(potential, overlap, 1, &cp, &ep);

          ctensorT co; rtensorT eo;
          syev(overlap, &co, &eo);

          if (_world.rank() == 0)
          {
            print("kinetic eigenvalues");
            print(ek);
          }

          if (_world.rank() == 0)
          {
            print("potential eigenvalues");
            print(ep);
          }

          if (_world.rank() == 0)
          {
            print("overlap eigenvalues");
            print(eo);
          }

          compress(_world, ao);
          _world.gop.fence();
          // Take linear combinations of the gaussian basis orbitals as the starting
          // orbitals for solver
          vecfuncT tmp_orbitals = transform(_world, ao, c(_, Slice(0, nao - 1)));
          _world.gop.fence();
          truncate(_world, tmp_orbitals);
          normalize(_world, tmp_orbitals);
          // Build the overlap matrix
          if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
          ctensorT overlap2 = matrix_inner(_world, tmp_orbitals, tmp_orbitals, true);

          rtensorT tmp_eigs = e(Slice(0, nao - 1));

          if (_world.rank() == 0) printf("(%8.4f,%8.4f,%8.4f)\n",kpt.k[0], kpt.k[1], kpt.k[2]);
          if (_world.rank() == 0) print(tmp_eigs);
          if (_world.rank() == 0) print("\n");

          //if (_world.rank() == 0) print("kinetic energy for kp = ", kp);
          //if (_world.rank() == 0) print(kinetic);
          //if (_world.rank() == 0) print("\n");

          // DEBUG
          if (_world.rank() == 0) printf("Overlap: \n");
          for (int i = 0; i < kinetic.dim(0); i++)
          {
            for (int j = 0; j < kinetic.dim(1); j++)
            {
              if (_world.rank() == 0) printf("%10.5f", real(overlap(i,j)));
            }
            if (_world.rank() == 0) printf("\n");
          }
          if (_world.rank() == 0) printf("\n");
          if (_world.rank() == 0) printf("\n");

          if (_world.rank() == 0) printf("Kinetic: \n");
          for (int i = 0; i < kinetic.dim(0); i++)
          {
            for (int j = 0; j < kinetic.dim(1); j++)
            {
              if (_world.rank() == 0) printf("%10.5f", real(kinetic(i,j)));
            }
            if (_world.rank() == 0) printf("\n");
          }
          if (_world.rank() == 0) printf("\n");
          if (_world.rank() == 0) printf("\n");

          if (_world.rank() == 0) printf("V: \n");
          for (int i = 0; i < potential.dim(0); i++)
          {
            for (int j = 0; j < potential.dim(1); j++)
            {
              if (_world.rank() == 0) printf("%10.5f", real(potential(i,j)));
            }
            if (_world.rank() == 0) printf("\n");
          }
          if (_world.rank() == 0) printf("\n");
          if (_world.rank() == 0) printf("\n");

          if (_world.rank() == 0) printf("Fock: \n");
          for (int i = 0; i < fock.dim(0); i++)
          {
            for (int j = 0; j < fock.dim(1); j++)
            {
              if (_world.rank() == 0) printf("%10.5f", real(fock(i,j)));
            }
            if (_world.rank() == 0) printf("\n");
          }
          if (_world.rank() == 0) printf("\n");
          if (_world.rank() == 0) printf("\n");

          if (_world.rank() == 0) printf("New overlap: \n");
          for (int i = 0; i < overlap2.dim(0); i++)
          {
            for (int j = 0; j < overlap2.dim(1); j++)
            {
              if (_world.rank() == 0) printf("%10.5f", real(overlap2(i,j)));
            }
            if (_world.rank() == 0) printf("\n");
          }
          if (_world.rank() == 0) printf("\n");
          if (_world.rank() == 0) printf("\n");

          // Fill in orbitals and eigenvalues
          int kend = kp + _params.nbands;
          kpt.begin = kp;
          kpt.end = kend;
          for (int oi = kp, ti = 0; oi < kend; oi++, ti++)
          {
            //if (_world.rank() == 0) print(oi, ti, kpt.begin, kpt.end);
            // normalize the orbitals
            //tmp_orbitals[ti].scale(1.0/tmp_orbitals[ti].norm2());
            _phisa.push_back(tmp_orbitals[ti]);
            _phisb.push_back(tmp_orbitals[ti]);
            _eigsa[oi] = tmp_eigs[ti];
            _eigsb[oi] = tmp_eigs[ti];
          }

          kp += _params.nbands;
        }

      }
    }
    //*************************************************************************

    //*************************************************************************
    // Constructor
    Solver(World& world, 
           const rfuntionT& vnucrhon,
           const vecfuncT& phis,
           const std::vector<T>& eigs,
           const ElectronicStructureParams& params,
           MolecularEntity mentity)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phis), _phisb(phis),
       _eigsa(eigs), _eigsb(eigs), _params(params), _mentity(mentity)
    {
      _residual = 1e5;
      if (params.periodic)
      {
        Tensor<double> box = FunctionDefaults<NDIM>::get_cell_width();
        _cop = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1, box);
      }
      else
      {
        _cop = CoulombOperatorPtr<T>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1);
      }

      if (params.ispotential)
      {
        _vnuc = copy(_vnucrhon);
      }
      else
      {
        _vnuc = apply(*_cop, _vnucrhon);
      }
    }
    //*************************************************************************

    //*************************************************************************
    // Constructor
    Solver(World& world, 
           rfuntionT vnucrhon,
           vecfuncT phis,
           std::vector<T> eigs,
           std::vector<KPoint> kpoints,
           std::vector<double> occs,
           ElectronicStructureParams params,
           MolecularEntity mentity)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phis), _phisb(phis),
         _eigsa(eigs), _eigsb(eigs), _params(params),
         _kpoints(kpoints), _occs(occs), _mentity(mentity)
    {
      _residual = 1e5;
      if (params.periodic)
      {
        Tensor<double> box = FunctionDefaults<NDIM>::get_cell_width();
        _cop = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1, box);
      }
      else
      {
        _cop = CoulombOperatorPtr<T>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1);
      }

      if (params.ispotential)
      {
        _vnuc = copy(_vnucrhon);
      }
      else
      {
        _vnuc = apply(*_cop, _vnucrhon);
      }
    }
    //*************************************************************************

    //*************************************************************************
    virtual ~Solver()
    {
      delete _subspace;
    }
    //*************************************************************************

    //***************************************************************************
    rfuntionT compute_rho(const vecfuncT& phis, std::vector<KPoint> kpoints)
    {
      // Electron density
      rfuntionT rho = rfactoryT(_world);
      _world.gop.fence();
      if (_world.rank() == 0) print("computing rho ...");
      // Loop over k-points
      for (unsigned int kp = 0; kp < kpoints.size(); kp++)
      {
        // get k-point
        KPoint kpoint = kpoints[kp];
        // loop through bands
        for (unsigned int j = kpoint.begin; j < kpoint.end; j++)
        {
          // Get phi(j) from iterator
          const functionT& phij = phis[j];
          // Compute the j-th density
          rfuntionT prod = abs_square(phij);
          double rnrm = prod.trace();
          prod.scale(1/rnrm);
          rho += 0.5 * _occs[j] * kpoint.weight * prod;
        }
      }
      rho.truncate();
      return rho;
    }
    //***************************************************************************

    //***************************************************************************
    std::vector<poperatorT> make_bsh_operators(const std::vector<T>& eigs)
    {
      // Make BSH vector
      std::vector<poperatorT> bops;
      // Get defaults
      int k = FunctionDefaults<NDIM>::get_k();
      double tol = FunctionDefaults<NDIM>::get_thresh();
      // Loop through eigenvalues, adding a BSH operator to bops
      // for each eigenvalue
      int sz = eigs.size();
      for (int i = 0; i < sz; i++)
      {
          T eps = eigs[i];
          if (eps > 0)
          {
              if (_world.rank() == 0)
              {
                  std::cout << "bsh: warning: positive eigenvalue" << i << eps << std::endl;
              }
              eps = -0.1;
          }
          if (_params.periodic)
          {
            Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
            bops.push_back(poperatorT(PeriodicBSHOpPtr<T,NDIM>(_world, sqrt(-2.0*eps), k, _params.lo, tol * 0.1,
                cellsize)));
          }
          else
          {
            bops.push_back(poperatorT(BSHOperatorPtr3D<T>(_world, sqrt(-2.0*eps), k, _params.lo, tol * 0.1)));
          }
      }
      return bops;
    }
    //*************************************************************************

    //*************************************************************************
    double calculate_kinetic_energy()
    {
      double ke = 0.0;
      if (!_params.periodic)
      {
        for (unsigned int i = 0; i < _phisa.size(); i++)
        {
          for (int axis = 0; axis < 3; axis++)
          {
            functionT dpsi = diff(_phisa[i], axis);
            ke += 0.5 * real(inner(dpsi, dpsi));
          }
        }
        if (_params.spinpol)
        {
          for (unsigned int i = 0; i < _phisb.size(); i++)
          {
            for (int axis = 0; axis < 3; axis++)
            {
              functionT dpsi = diff(_phisb[i], axis);
              ke += 0.5 * real(inner(dpsi, dpsi));
            }
          }
        }
        else
        {
          ke *= 2.0;
        }
      }
      return ke;
    }
    //*************************************************************************

    //*************************************************************************
    void apply_potential(int iter, vecfuncT& pfuncsa,
        vecfuncT& pfuncsb, const vecfuncT& phisa,
        const vecfuncT& phisb, const rfuntionT& rhoa, const rfuntionT& rhob,
        const rfuntionT& rho)
    {
      // Nuclear and coulomb potentials
      rfuntionT vc = apply(*_cop, rho);
      rfuntionT vlocal = _vnuc + vc;
      // WSTHORNTON
      std::ostringstream strm;
      strm << "vlocal" << iter << ".out";
      std::string fname = strm.str();
      //plot_line(fname.c_str(), 100, coordT(-5.0), coordT(5.0), _vnuc, vc);
      // Calculate energies for Coulomb and nuclear
      double ce = 0.5*inner(vc,rho);
      double pe = inner(_vnuc,rho);
      double xc = 0.0;
      double ke = calculate_kinetic_energy();
      // Exchange
      if (_params.functional == 1)
      {
        // LDA, is calculation spin-polarized?
        if (_params.spinpol)
        {
          // potential
          rfuntionT vxca = binary_op(rhoa, rhob, &::libxc_ldaop_sp<double>);
          rfuntionT vxcb = binary_op(rhob, rhoa, &::libxc_ldaop_sp<double>);
          pfuncsa = mul_sparse(_world, vlocal + vxca, phisa, _params.thresh * 0.1);
          pfuncsb = mul_sparse(_world, vlocal + vxcb, phisb, _params.thresh * 0.1);
          // energy
          rfuntionT fca = binary_op(rhoa, rhob, &::libxc_ldaeop_sp<double>);
          rfuntionT fcb = binary_op(rhob, rhoa, &::libxc_ldaeop_sp<double>);
          xc = fca.trace() + fcb.trace();
        }
        else
        {
          // potential
          rfuntionT vxc = copy(rhoa);
          vxc.unaryop(&::libxc_ldaop<double>);
          rfuntionT vxc2 = binary_op(rhoa, rhoa, &::libxc_ldaop_sp<double>);
          pfuncsa = mul_sparse(_world, vlocal + vxc2, phisa, _params.thresh * 0.1);
          // energy
          rfuntionT fc = copy(rhoa);
          fc.unaryop(&::ldaeop<double>);
          xc = fc.trace();
        }
      }
      std::cout.precision(8);
      if (_world.rank() == 0)
      {
        print("Energies:");
        print("Kinetic energy:\t\t ", ke);
        print("Potential energy:\t ", pe);
        print("Coulomb energy:\t\t ", ce);
        print("Exchage energy:\t\t ", xc, "\n");
        print("Total energy:\t\t ", ke + pe + ce + xc, "\n\n");
      }
    }
    //*************************************************************************

    //*************************************************************************
    void reproject()
    {
      double tol = 1e-9;
      if (_params.thresh > (_maxthresh - tol))
      {
       if (_world.rank() == 0)
       {
         printf("Reprojecting:\n");
         printf("_params.thresh = %.5e\t_maxthresh = %.5e\n", _params.thresh, _maxthresh);
        _params.thresh *= 10;
       }
      }
    }
    //*************************************************************************

    //*************************************************************************
    void solve()
    {
      
      // WSTHORNTON
      // multiply gamma-point WF by random phase
      double_complex t1 = RandomValue<double_complex>();
      double_complex t2 = exp(t1);
      std::vector<double_complex> phase(_phisa.size(), t2);
      scale(_world, _phisa, phase);
       
       // set maximum thresh from _params
      //_maxthresh = _params.thresh;
      //_params.thresh *= 1e2;

      for (int it = 0; it < _params.maxits; it++)
      {
        //if ((it > 0) && ((it % 4) == 0))
        //{
        //  reproject();
        //}

        if (_world.rank() == 0) print("it = ", it);
       
        // Compute density
        _rhoa = compute_rho(_phisa, _kpoints);
        _rhob = (_params.spinpol) ? compute_rho(_phisb, _kpoints) : _rhoa;
        _rho = _rhoa + _rhob;
        double rtrace = _rho.trace();
        if (_world.rank() == 0) print("trace of rho", rtrace);

        std::vector<functionT> pfuncsa =
                zero_functions<valueT,NDIM>(_world, _phisa.size());
        std::vector<functionT> pfuncsb =
                zero_functions<valueT,NDIM>(_world, _phisb.size());

        // Apply the potentials to the orbitals
        if (_world.rank() == 0) print("applying potential ...\n");
        START_TIMER(_world);
        apply_potential(it, pfuncsa, pfuncsb, _phisa, _phisb, _rhoa, _rhob, _rho);
        END_TIMER(_world,"apply potential");

        // Do right hand side for all k-points
        do_rhs(_phisa, pfuncsa, _eigsa, _kpoints);
        
        // Make BSH Green's function
        std::vector<poperatorT> bopsa = make_bsh_operators(_eigsa);
        std::vector<T> sfactor(pfuncsa.size(), -2.0);
        scale(_world, pfuncsa, sfactor);

        // Apply Green's function to orbitals
        if (_world.rank() == 0) print("applying BSH operator ...\n");
        truncate<valueT,NDIM>(_world, pfuncsa);
        START_TIMER(_world);
        std::vector<functionT> tmpa = apply(_world, bopsa, pfuncsa);
        END_TIMER(_world,"apply BSH");
        bopsa.clear();

        // Do other spin
        vecfuncT tmpb = zero_functions<valueT,NDIM>(_world, _phisb.size());
        if (_params.spinpol)
        {
          do_rhs(_phisb, pfuncsb, _eigsb, _kpoints);
          std::vector<poperatorT> bopsb = make_bsh_operators(_eigsb);
          scale(_world, pfuncsb, sfactor);
          truncate<valueT,NDIM>(_world, pfuncsb);
          tmpb = apply(_world, bopsb, pfuncsb);
          bopsb.clear();
        }

        // Update orbitals
        update_orbitals(tmpa, tmpb, _kpoints);

//        std::cout.precision(8);
//        if (_world.rank() == 0)
//        {
//          print("Iteration: ", it, "\nEigenvalues for alpha spin: \n");
//          for (unsigned int i = 0; i < _eigsa.size(); i++)
//          {
//            print(_eigsa[i]);
//          }
//          print("\n\n");
//        }
//        if (_params.spinpol)
//        {
//          if (_world.rank() == 0)
//          {
//            print("Eigenvalues for beta spin: \n");
//            for (unsigned int i = 0; i < _eigsb.size(); i++)
//            {
//              print(_eigsb[i]);
//            }
//            print("\n\n");
//          }
//        }
      }

      save_orbitals();
    }
    //*************************************************************************
    
    //*************************************************************************
    ctensorT matrix_exponential(const ctensorT& A) {
        const double tol = 1e-13;
        MADNESS_ASSERT(A.dim(0) == A.dim(1));

        // Scale A by a power of 2 until it is "small"
        double anorm = A.normf();
        int n = 0;
        double scale = 1.0;
        while (anorm*scale > 0.1)
        {
            n++;
            scale *= 0.5;
        }
        tensorT B = scale*A;    // B = A*2^-n

        // Compute exp(B) using Taylor series
        ctensorT expB = ctensorT(2, B.dims());
        for (int i = 0; i < expB.dim(0); i++) expB(i,i) = std::complex<T>(1.0,0.0);

        int k = 1;
        ctensorT term = B;
        while (term.normf() > tol)
        {
            expB += term;
            term = inner(term,B);
            k++;
            term.scale(1.0/k);
        }

        // Repeatedly square to recover exp(A)
        while (n--)
        {
            expB = inner(expB,expB);
        }

        return expB;
    }
    //*************************************************************************

    //*************************************************************************
    void do_rhs(vecfuncT& wf,
                vecfuncT& vwf,
                std::vector<T>& alpha,
                std::vector<KPoint> kpoints)
    {
      // tolerance
      double trantol = 0.1*_params.thresh/std::min(30.0,double(wf.size()));
      double thresh = 1e-4;

      for (unsigned int kp = 0; kp < kpoints.size(); kp++)
      {
        // Get k-point and orbitals for this k-point
        KPoint kpoint = kpoints[kp];
//        double k0 = 2.0 * madness::constants::pi * kpoint.k[0];
//        double k1 = 2.0 * madness::constants::pi * kpoint.k[1];
//        double k2 = 2.0 * madness::constants::pi * kpoint.k[2];
        double k0 = kpoint.k[0];
        double k1 = kpoint.k[1];
        double k2 = kpoint.k[2];
        // WSTHORNTON
        // Extract the relevant portion of the list of orbitals and the list of the
        // V times the orbitals
        vecfuncT k_wf(wf.begin() + kpoint.begin, wf.begin() + kpoint.end);
        vecfuncT k_vwf(vwf.begin() + kpoint.begin, vwf.begin() + kpoint.end);

        // Build fock matrix
        tensorT fock = build_fock_matrix(k_wf, k_vwf, kpoint);
        tensorT overlap = matrix_inner(_world, k_wf, k_wf, true);

        // Do right hand side stuff for kpoint
        bool isgamma = (is_equal(k0,0.0,1e-5) && is_equal(k1,0.0,1e-5) && is_equal(k2,0.0,1e-5));
        if (_params.periodic && !isgamma) // Non-zero k-point
        {
          // Do the gradient term and k^2/2
          vecfuncT d_wf = zero_functions<valueT,NDIM>(_world, k_wf.size());
          for (unsigned int i = 0; i < k_wf.size(); i++)
          {
            // gradient
            functionT dx_wf = pdiff(k_wf[i], 0, true);
            functionT dy_wf = pdiff(k_wf[i], 1, true);
            functionT dz_wf = pdiff(k_wf[i], 2, true);
            d_wf[i] = std::complex<T>(0.0,k0)*dx_wf + 
                      std::complex<T>(0.0,k1)*dy_wf + 
                      std::complex<T>(0.0,k2)*dz_wf;
            // k^/2
            double ksq = k0*k0 + k1*k1 + k2*k2;
            k_vwf[i] += 0.5 * ksq * k_wf[i];
          }
          gaxpy(_world, 1.0, k_vwf, -1.0, d_wf);
        }

        if (_params.canon) // canonical orbitals
        {

          // Debug output
          if (_world.rank() == 0)
          { print("Overlap matrix:");
            print(overlap);
          }
          if (_world.rank() == 0)
          { print("Fock matrix:");
            print(fock);
          }

          ctensorT U; rtensorT e;
          sygv(fock, overlap, 1, &U, &e);

          unsigned int nmo = k_wf.size();
          // Fix phases.
          long imax;
          for (long j = 0; j < nmo; j++)
          {
              // Get index of largest value in column
              U(_,j).absmax(&imax);
              T ang = arg(U(imax,j));
              std::complex<T> phase = exp(std::complex<T>(0.0,-ang));
              // Loop through the rest of the column and divide by the phase
              for (long i = 0; i < nmo; i++)
              {
                U(i,j) *= phase;
              }
          }

          // Within blocks with the same occupation number attempt to
          // keep orbitals in the same order (to avoid confusing the
          // non-linear solver).  Have to run the reordering multiple
          // times to handle multiple degeneracies.
          for (int pass = 0; pass < 5; pass++)
          {
              long j;
              for (long i = 0; i < nmo; i++)
              {
                U(_, i).absmax(&j);
                if (i != j)
                {
                  tensorT tmp = copy(U(_, i));
                  U(_, i) = U(_, j);
                  U(_, j) = tmp;
                  //swap(e[i], e[j]);
                  T ti = e[i];
                  T tj = e[j];
                  e[i] = tj; e[j] = ti;
                }
              }
          }

          // Rotations between effectively degenerate states confound
          // the non-linear equation solver ... undo these rotations
          long ilo = 0; // first element of cluster
          while (ilo < nmo-1) {
              long ihi = ilo;
              while (fabs(real(e[ilo]-e[ihi+1])) < thresh*10.0*max(fabs(real(e[ilo])),1.0)) {
                  ihi++;
                  if (ihi == nmo-1) break;
              }
              long nclus = ihi - ilo + 1;
              if (nclus > 1) {
                  if (_world.rank() == 0) print("   found cluster", ilo, ihi);
                  tensorT q = copy(U(Slice(ilo,ihi),Slice(ilo,ihi)));
                  //print(q);
                  // Special code just for nclus=2
                  // double c = 0.5*(q(0,0) + q(1,1));
                  // double s = 0.5*(q(0,1) - q(1,0));
                  // double r = sqrt(c*c + s*s);
                  // c /= r;
                  // s /= r;
                  // q(0,0) = q(1,1) = c;
                  // q(0,1) = -s;
                  // q(1,0) = s;

                  // Iteratively construct unitary rotation by
                  // exponentiating the antisymmetric part of the matrix
                  // ... is quadratically convergent so just do 3
                  // iterations
                  ctensorT rot = matrix_exponential(-0.5*(q - conj_transpose(q)));
                  q = inner(q,rot);
                  ctensorT rot2 = matrix_exponential(-0.5*(q - conj_transpose(q)));
                  q = inner(q,rot2);
                  ctensorT rot3 = matrix_exponential(-0.5*(q - conj_transpose(q)));
                  q = inner(rot,inner(rot2,rot3));
                  U(_,Slice(ilo,ihi)) = inner(U(_,Slice(ilo,ihi)),q);
              }
              ilo = ihi+1;
          }

          // Debug output
          if (_world.rank() == 0)
          { print("Overlap matrix:");
            print(overlap);
          }
          if (_world.rank() == 0)
          { print("Fock matrix:");
            print(fock);
          }
          if (_world.rank() == 0)
          { print("U matrix: (eigenvectors)");
            print(U);
          }

          if (_params.solver == 0)
          {
            // transform orbitals and V * (orbitals)
            k_vwf = transform(_world, k_vwf, U, 1e-5 / std::min(30.0, double(k_wf.size())), false);
            k_wf = transform(_world, k_wf, U, FunctionDefaults<3>::get_thresh() / std::min(30.0, double(k_wf.size())), true);
          }

          for (unsigned int ei = kpoint.begin, fi = 0; ei < kpoint.end;
            ei++, fi++)
          {
            valueT t1 = (_params.solver == 0) ? e(fi,fi) : fock(fi,fi);
            //if (real(e(fi,fi)) > -0.1)
            if (real(t1) > -0.1)
            {
              alpha[ei] = -0.5;
              k_vwf[fi] += (alpha[ei]-real(t1))*k_wf[fi];
            }
            else
            {
              alpha[ei] = real(t1);
            }
          }
          for (unsigned int ei = 0; ei < e.dim(0); ei++)
          {
            if (_world.rank() == 0)
              print("kpoint ", kp, "ei ", ei, "eps ", real(e(ei,ei)));
          }
        }
        else // non-canonical orbitals
        {
          // diagonlize just to print eigenvalues
          tensorT overlap = matrix_inner(_world, k_wf, k_wf, true);
          ctensorT c; rtensorT e;
          sygv(fock, overlap, 1, &c, &e);
          for (unsigned int ei = 0; ei < e.dim(0); ei++)
          {
            double diffe = (ei == 0) ? 0.0 : real(e(ei,ei))-real(e(ei-1,ei-1));
            if (_world.rank() == 0)
              print("kpoint ", kp, "ei ", ei, "eps ", real(e(ei,ei)), "\tdiff\t", diffe);
          }

          for (unsigned int ei = kpoint.begin, fi = 0; 
            ei < kpoint.end; ei++, fi++)
          {
            alpha[ei] = std::min(-0.1, real(fock(fi,fi)));
            fock(fi,fi) -= std::complex<T>(alpha[ei], 0.0);
          }

          std::vector<functionT> fwf = transform(_world, k_wf, fock, trantol);
          gaxpy(_world, 1.0, k_vwf, -1.0, fwf);
          fwf.clear();
        }
        for (unsigned int wi = kpoint.begin, fi = 0; wi < kpoint.end;
          wi++, fi++)
        {
          wf[wi] = k_wf[fi];
          vwf[wi] = k_vwf[fi];
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    tensorT build_fock_matrix(vecfuncT& psi,
                              vecfuncT& vpsi,
                              KPoint kpoint)
    {
      // Build the potential matrix
      START_TIMER(_world);
      tensorT potential = matrix_inner(_world, psi, vpsi, true);
      _world.gop.fence();
      END_TIMER(_world,"kinetic energy matrix");

      START_TIMER(_world);
      if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
        tensorT kinetic = ::kinetic_energy_matrix(_world, psi, 
                                                  _params.periodic,
                                                  kpoint);
      END_TIMER(_world,"potential energy matrix");

      if (_world.rank() == 0) print("Constructing Fock matrix ...\n\n");
      tensorT fock = potential + kinetic;
      fock = 0.5 * (fock + transpose(fock));
      _world.gop.fence();

      return fock;
    }
    //*************************************************************************

    //*************************************************************************
    void gram_schmidt(vecfuncT& f, KPoint kpoint)
    {
      for (unsigned int fi = kpoint.begin; fi < kpoint.end; ++fi)
      {
        // Project out the lower states
        for (unsigned int fj = kpoint.begin; fj < fi; ++fj)
        {
          valueT overlap = inner(f[fj], f[fi]);
          f[fi] -= overlap*f[fj];
        }
        f[fi].scale(1.0/f[fi].norm2());
      }
    }
    //*************************************************************************

    vecfuncT compute_residual(const vecfuncT& awfs,
                              const vecfuncT& bwfs)
    {
      // vector of residual functions
      vecfuncT rm = sub(_world, _phisa, awfs);
      // if spin-polarized insert beta spin-orbital residual functions
      if (_params.spinpol)
      {
        vecfuncT br = sub(_world, _phisb, bwfs);
        rm.insert(rm.end(), br.begin(), br.end());
      }
      // scalar residual
      std::vector<double> rnvec = norm2s<valueT,NDIM>(_world, rm);
      double rnorm = 0.0;
      for (unsigned int i = 0; i < rnvec.size(); i++) rnorm += rnvec[i];
      // renormalize and print
      _residual = rnorm / rnvec.size();
      if (_world.rank() == 0) print("residual = ", _residual);
      if (_world.rank() == 0)
      {
        printf("\n");
        for (unsigned int i = 0; i < rnvec.size(); i++)
        {
          print("residual", i, "\t\t", rnvec[i]);
        }
        printf("\n");
      }
      
      return rm;
    }

    //*************************************************************************
    void update_orbitals(vecfuncT& awfs,
                         vecfuncT& bwfs,
                         std::vector<KPoint> kpoints)
    {
      // truncate before we do anyting
      truncate<valueT,NDIM> (_world, awfs);
      truncate<valueT,NDIM> (_world, _phisa);
      if (_params.spinpol)
      {
        truncate<valueT,NDIM> (_world, bwfs);
        truncate<valueT,NDIM> (_world, _phisb);
      }
      // compute residual
      vecfuncT rm = compute_residual(awfs, bwfs);
      if (_params.solver > 0 && _params.maxsub > 1)
      {
        // nonlinear solver
        _subspace->update_subspace(awfs, bwfs, _phisa, _phisb, rm);
      }

      // do step restriction
      step_restriction(_phisa, awfs, 0);
      if (_params.spinpol)
      {
        step_restriction(_phisb, bwfs, 1);
      }
      // do gram-schmidt
      for (unsigned int kp = 0; kp < kpoints.size(); kp++)
      {
        gram_schmidt(awfs, kpoints[kp]);
        if (_params.spinpol)
        {
          gram_schmidt(bwfs, kpoints[kp]);
        }
      }
      // update alpha and beta orbitals
      truncate<valueT,NDIM>(_world, awfs);
      for (unsigned int ai = 0; ai < awfs.size(); ai++) {
          _phisa[ai] = awfs[ai].scale(1.0 / awfs[ai].norm2());
      }
      if (_params.spinpol)
      {
        truncate<valueT,NDIM>(_world, bwfs);
        for (unsigned int bi = 0; bi < bwfs.size(); bi++) {
            _phisb[bi] = bwfs[bi].scale(1.0 / bwfs[bi].norm2());
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    void step_restriction(vecfuncT& owfs,
                          vecfuncT& nwfs,
                          int aorb)
    {
      std::vector<double> rnorm = norm2s(_world, sub(_world, owfs, nwfs));
      // Step restriction
      int nres = 0;
      for (unsigned int i = 0; i < owfs.size(); i++)
      {
        if (rnorm[i] > _params.maxrotn)
        {
          double s = _params.maxrotn / rnorm[i];
          nres++;
          if (_world.rank() == 0)
          {
            if (!aorb && nres == 1) printf("  restricting step for alpha orbitals:");
            if (aorb && nres == 1) printf("  restricting step for beta orbitals:");
            printf(" %d", i);
          }
          nwfs[i].gaxpy(s, owfs[i], 1.0 - s, false);
        }
      }
      if (nres > 0 && _world.rank() == 0) printf("\n");
      _world.gop.fence();
    }
    //*************************************************************************

    //*************************************************************************
    void update_eigenvalues(const vecfuncT& wavefs,
        const vecfuncT& pfuncs, const vecfuncT& phis,
        std::vector<T>& eigs)
    {
      // Update e
      if (_world.rank() == 0) printf("Updating e ...\n\n");
      for (unsigned int ei = 0; ei < eigs.size(); ei++)
      {
        functionT r = wavefs[ei] - phis[ei];
        double tnorm = wavefs[ei].norm2();
        // Compute correction to the eigenvalues
        T ecorrection = -0.5*real(inner(pfuncs[ei], r)) / (tnorm*tnorm);
        T eps_old = eigs[ei];
        T eps_new = eps_old + ecorrection;
//        if (_world.rank() == 0) printf("ecorrection = %.8f\n\n", ecorrection);
//        if (_world.rank() == 0) printf("eps_old = %.8f eps_new = %.8f\n\n", eps_old, eps_new);
        // Sometimes eps_new can go positive, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
        // I bounce the new eigenvalue back into the negative side of the real axis. I
        // keep doing this until it's good or I've already done it 10 times.
        int counter = 50;
        while (eps_new >= 0.0 && counter < 20)
        {
          // Split the difference between the new and old estimates of the
          // pi-th eigenvalue.
          eps_new = eps_old + 0.5 * (eps_new - eps_old);
          counter++;
        }
        // Still no go, forget about it. (1$ to Donnie Brasco)
        if (eps_new >= 0.0)
        {
          if (_world.rank() == 0) printf("FAILURE OF WST: exiting!!\n\n");
          _exit(0);
        }
        // Set new eigenvalue
        eigs[ei] = eps_new;
      }
    }
    //*************************************************************************

//    //*************************************************************************
//    double get_eig(int indx)
//    {
//      return _solver->get_eig(indx);
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    functionT get_phi(int indx)
//    {
//      return _solver->get_phi(indx);
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    const std::vector<double>& eigs()
//    {
//      return _solver->eigs();
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    const vecfuncT& phis()
//    {
//      return _solver->phis();
//    }
//    //*************************************************************************

  };
  //***************************************************************************

}
#define SOLVER_H_
/*
  @}
  @}
 */

#endif /* SOLVER_H_ */
