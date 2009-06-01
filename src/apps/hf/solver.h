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
//#include "eigsolver.h"

#ifndef SOLVER_H_

//*****************************************************************************
static double onesfunc(const coordT& x)
{
  return 1.0;
}
//*****************************************************************************

namespace madness
{
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
    typedef pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef vector<pairvecfuncT> subspaceT;

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    // This variable could either be a nuclear potiential or a nuclear charge
    // density depending on the "ispotential" variable in the
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

    //*************************************************************************
    bool newscheme() {return true;}
    //*************************************************************************

    //*************************************************************************
    subspaceT _subspace;
    //*************************************************************************

    //*************************************************************************
    tensorT _Q;
    //*************************************************************************

  public:
    //*************************************************************************
    // Constructor
    Solver(World& world, 
           rfuntionT vnucrhon,
           vecfuncT phisa,
           vecfuncT phisb,
           std::vector<T> eigsa,
           std::vector<T> eigsb,
           ElectronicStructureParams params)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phisa), _phisb(phisb),
       _eigsa(eigsa), _eigsb(eigsb), _params(params)
    {
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

    //*************************************************************************
    // Constructor
    Solver(World& world, 
           const rfuntionT& vnucrhon,
           const vecfuncT& phis,
           const std::vector<T>& eigs,
           const ElectronicStructureParams& params)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phis), _phisb(phis),
       _eigsa(eigs), _eigsb(eigs), _params(params)
    {
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
           ElectronicStructureParams params)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phis), _phisb(phis),
         _eigsa(eigs), _eigsb(eigs), _params(params),
         _kpoints(kpoints), _occs(occs)
    {
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

    //***************************************************************************
    rfuntionT compute_rho(const vecfuncT& phis)
    {
      // Electron density
      rfuntionT rho = rfactoryT(_world);
      _world.gop.fence();
      // Loop over all wavefunctions to compute density
      for (unsigned int j = 0; j < phis.size(); j++)
      {
        // Get phi(j) from iterator
        const functionT& phij = phis[j];
        // Compute the j-th density
        //functionT prod = square(phij);
        rfuntionT prod = abs_square(phij);
        rho += 0.5*_occs[j]*prod;
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
          rfuntionT vxca = binary_op(rhoa, rhob, &::libxc_ldaop_sp);
          rfuntionT vxcb = binary_op(rhob, rhoa, &::libxc_ldaop_sp);
          pfuncsa = mul_sparse(_world, vlocal + vxca, phisa, _params.thresh * 0.1);
          pfuncsb = mul_sparse(_world, vlocal + vxcb, phisb, _params.thresh * 0.1);
          // energy
          rfuntionT fca = binary_op(rhoa, rhob, &::libxc_ldaeop_sp);
          rfuntionT fcb = binary_op(rhob, rhoa, &::libxc_ldaeop_sp);
          xc = fca.trace() + fcb.trace();
        }
        else
        {
          // potential
          rfuntionT vxc = copy(rhoa);
          vxc.unaryop(&::libxc_ldaop);
          rfuntionT vxc2 = binary_op(rhoa, rhoa, &::libxc_ldaop_sp);
          pfuncsa = mul_sparse(_world, vlocal + vxc2, phisa, _params.thresh * 0.1);
          // energy
          rfuntionT fc = copy(rhoa);
          fc.unaryop(&::ldaeop);
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
    virtual ~Solver() {}
    //*************************************************************************

    //*************************************************************************
    void solve()
    {
      for (int it = 0; it < _params.maxits; it++)
      {
        if (_world.rank() == 0) print("it = ", it);
        // Compute density
        _rhoa = compute_rho(_phisa);
        _rhob = (_params.spinpol) ? compute_rho(_phisb) : _rhoa;
        _rho = _rhoa + _rhob;

        vector<functionT> pfuncsa =
                zero_functions<valueT,NDIM>(_world, _phisa.size());
        vector<functionT> pfuncsb =
                zero_functions<valueT,NDIM>(_world, _phisb.size());

        // Apply the potentials to the orbitals
        if (_world.rank() == 0) print("applying potential ...\n");
        apply_potential(it, pfuncsa, pfuncsb, _phisa, _phisb, _rhoa, _rhob, _rho);

        // Do right hand side for all k-points
        do_rhs(_phisa, pfuncsa, _eigsa, _kpoints);
        
        // Make BSH Green's function
        std::vector<poperatorT> bopsa = make_bsh_operators(_eigsa);
        vector<T> sfactor(pfuncsa.size(), -2.0);
        scale(_world, pfuncsa, sfactor);

        // Apply Green's function to orbitals
        if (_world.rank() == 0) print("applying BSH operator ...\n");
        truncate<valueT,NDIM>(_world, pfuncsa);
        vector<functionT> tmpa = apply(_world, bopsa, pfuncsa);
        bopsa.clear();

        plot_line("pfunc", 100, coordT(-_params.L/2), coordT(_params.L/2), pfuncsa[0], pfuncsa[1], pfuncsa[2]);
        plot_line("tmpa", 100, coordT(-_params.L/2), coordT(_params.L/2), tmpa[0], tmpa[1], tmpa[2]);
        exit(0);

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

        std::cout.precision(8);
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
    }
    //*************************************************************************
    
    //*************************************************************************
    void do_rhs(vecfuncT& wf,
                vecfuncT& vwf,
                std::vector<T>& eps,
                std::vector<KPoint> kpoints)
    {
      // tolerance
      double trantol = 0.1*_params.thresh/min(30.0,double(wf.size()));


      for (unsigned int kp = 0; kp < kpoints.size(); kp++)
      {
        // Get k-point and orbitals for this k-point
        KPoint kpoint = kpoints[kp];
        // WSTHORNTON
        print("kpoint info:");
        print(wf.size(), vwf.size(), kpoint.begin, kpoint.end);
        vecfuncT k_wf(wf.begin() + kpoint.begin, wf.begin() + kpoint.end);
        vecfuncT k_vwf(vwf.begin() + kpoint.begin, vwf.begin() + kpoint.end);

        // Build fock matrix
        tensorT fock = build_fock_matrix(k_wf, k_vwf, kpoint);

        if (_params.canon)
        {
          tensorT overlap = matrix_inner(_world, k_wf, k_wf, true);
          ctensorT c; rtensorT e;
          sygv(fock, overlap, 1, &c, &e);
          for (unsigned int ei = kpoint.begin, fi = 0; ei < kpoint.end;
            ei++, fi++)
          {
            // WSTHORNTON
            print("\n");
            print(ei, fi);
            if (real(e(fi,fi)) > -0.1)
            {
              eps[ei] = -0.1;
              vwf[ei] -= (real(e(fi,fi))-eps[ei])*wf[ei];
            }
            else
            {
              eps[ei] = e(fi,fi);
            }
            eps[ei] = std::min(-0.1, real(e(fi,fi)));
          }
          // WSTHORNTON
          // this will work if there is only 1 k-point
          if (_world.rank() == 0) print("eigenvalues:\n");
          for (unsigned int ei = 0; ei < eps.size(); ei++)
          {
            if (_world.rank() == 0) print(real(e(ei,ei)), eps[ei]);
          }
        }
        else
        {
          for (unsigned int ei = kpoint.begin, fi = 0; 
            ei < kpoint.end; ei++, fi++)
          {
            eps[ei] = std::min(-0.1, real(fock(fi,fi)));
            fock(fi,fi) -= std::complex<T>(eps[ei], 0.0);
          }

          vector<functionT> fwf = transform(_world, k_wf, fock, trantol);
          gaxpy(_world, 1.0, k_vwf, -1.0, fwf);
          fwf.clear();
          for (unsigned int wi = kpoint.begin, fi = 0; wi < kpoint.end;
            wi++, fi++)
          {
            vwf[wi] = k_vwf[fi];
          }
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
      tensorT potential = matrix_inner(_world, psi, vpsi, true);
      _world.gop.fence();

      if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
        tensorT kinetic = ::kinetic_energy_matrix(_world, psi, 
                                                  _params.periodic,
                                                  kpoint);

      if (_world.rank() == 0) print("Constructing Fock matrix ...\n\n");
      tensorT fock = potential + kinetic;
      fock = 0.5 * (fock + transpose(fock));
      _world.gop.fence();

      print(kinetic);
      print(potential);
      print(fock);
      
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
      if (_params.maxsub > 1)
      {
        // nonlinear solver
        update_subspace(awfs, bwfs);
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
    void update_subspace(vecfuncT& awfs,
                         vecfuncT& bwfs)
    {
      // compute residuals
      vecfuncT rm = sub(_world, _phisa, awfs);
      if (_params.spinpol)
      {
        vecfuncT br = sub(_world, _phisb, bwfs);
        rm.insert(rm.end(), br.begin(), br.end());
      }
      std::vector<double> rnvec = norm2<valueT,NDIM>(_world, rm);
      if (_world.rank() == 0)
      {
        double rnorm = 0.0;
        for (unsigned int i = 0; i < rnvec.size(); i++) rnorm += rnvec[i];
        print("residual = ", rnorm);
      }
      // concatentate up and down spins
      vecfuncT vm = _phisa;
      if (_params.spinpol)
      {
        vm.insert(vm.end(), _phisb.begin(), _phisb.end());
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
      print(_Q);

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
                  print("Increasing subspace singular value threshold ", c[m-1], rcond);
                  rcond *= 100;
              }
              else {
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
      vecfuncT phisa_new = zero_functions<valueT,NDIM>(_world, _phisa.size());
      vecfuncT phisb_new = zero_functions<valueT,NDIM>(_world, _phisb.size());
      compress(_world, phisa_new, false);
      compress(_world, phisb_new, false);
      _world.gop.fence();
      std::complex<double> one = std::complex<double>(1.0,0.0);
      for (unsigned int m=0; m<_subspace.size(); m++) {
          const vecfuncT& vm = _subspace[m].first;
          const vecfuncT& rm = _subspace[m].second;
          const vecfuncT  vma(vm.begin(),vm.begin()+_phisa.size());
          const vecfuncT  rma(rm.begin(),rm.begin()+_phisa.size());
          const vecfuncT  vmb(vm.end()-_phisb.size(), vm.end());
          const vecfuncT  rmb(rm.end()-_phisb.size(), rm.end());

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
      awfs = phisa_new;
      bwfs = phisb_new;
    }
    //*************************************************************************

    //*************************************************************************
    void step_restriction(vecfuncT& owfs,
                          vecfuncT& nwfs,
                          int aorb)
    {
      vector<double> rnorm = norm2(_world, sub(_world, owfs, nwfs));
      // Step restriction
      double maxrotn = 0.1;
      int nres = 0;
      for (unsigned int i = 0; i < owfs.size(); i++)
      {
        if (rnorm[i] > maxrotn)
        {
          double s = maxrotn / rnorm[i];
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


#endif /* SOLVER_H_ */
