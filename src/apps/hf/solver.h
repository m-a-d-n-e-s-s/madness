#include <mra/mra.h>
#include <world/world.h>
#include <vector>

#include "poperator.h"
#include "libxc.h"
#include "electronicstructureparams.h"

#ifndef SOLVER_H_

namespace madness
{
  //***************************************************************************
  template <typename T, int NDIM>
  class Solver
  {
    // Typedef's
    typedef Function<T,NDIM> funcT;
    typedef Vector<double,NDIM> kvecT;
    typedef SeparatedConvolution<double,3> operatorT;
    typedef SharedPtr<operatorT> poperatorT;

    //*************************************************************************
    World _world;
    //*************************************************************************

    //*************************************************************************
    // This variable could either be a nuclear potiential or a nuclear charge
    // density depending on the "ispotential" variable in the
    // ElectronicStructureParams class.
    funcT _vnucrhon;
    //*************************************************************************

    //*************************************************************************
    std::vector<funcT> _phisa;
    //*************************************************************************

    //*************************************************************************
    std::vector<funcT> _phisb;
    //*************************************************************************

    //*************************************************************************
    std::vector<double> _eigsa;
    //*************************************************************************

    //*************************************************************************
    std::vector<double> _eigsb;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

    //*************************************************************************
    funcT _rhoa;
    //*************************************************************************

    //*************************************************************************
    funcT _rhob;
    //*************************************************************************

    //*************************************************************************
    funcT _rho;
    //*************************************************************************

    //*************************************************************************
    funcT _vnuc;
    //*************************************************************************

    //*************************************************************************
    SeparatedConvolution<T,NDIM>* _cop;
    //*************************************************************************


  public:
    //*************************************************************************
    // Constructor
    Solver(World& world, funcT vnucrhon, std::vector<funcT> phisa,
      std::vector<funcT> phisb, std::vector<double> eigsa, std::vector<double> eigsb,
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
      Solver(World& world, funcT vnucrhon, std::vector<funcT> phis,
          std::vector<double> eigs, ElectronicStructureParams params)
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

    //***************************************************************************
    funcT compute_rho(std::vector<funcT> phis)
    {
      // Electron density
      funcT rho = FunctionFactory<double,NDIM>(_world);
      // Loop over all wavefunctions to compute density
      for (unsigned int j = 0; j < phis.size(); j++)
      {
        // Get phi(j) from iterator
        const funcT& phij = phis[j];
        // Compute the j-th density
        funcT prod = square(phij);
        //rho += occs[j]*prod;
        rho += prod;
      }
      rho.truncate();
      return rho;
    }
    //***************************************************************************

    //***************************************************************************
    std::vector<poperatorT> make_bsh_operators(const std::vector<double>& eigs)
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
          double eps = eigs[i];
          if (eps > 0)
          {
              if (_world.rank() == 0)
              {
                  std::cout << "bsh: warning: positive eigenvalue" << i << eps << std::endl;
              }
              eps = -0.1;
          }
          bops.push_back(poperatorT(BSHOperatorPtr<double,NDIM>(_world, sqrt(-2.0*eps), k, _params.lo, tol * 0.1)));
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
            funcT dpsi = diff(_phisa[i], axis);
            ke += 0.5 * inner(dpsi, dpsi);
          }
        }
        if (_params.spinpol)
        {
          for (unsigned int i = 0; i < _phisb.size(); i++)
          {
            for (int axis = 0; axis < 3; axis++)
            {
              funcT dpsi = diff(_phisb[i], axis);
              ke += 0.5 * inner(dpsi, dpsi);
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
    void apply_potential(std::vector<funcT>& pfuncsa,
        std::vector<funcT>& pfuncsb, const std::vector<funcT>& phisa,
        const std::vector<funcT>& phisb, const funcT& rhoa, const funcT& rhob,
        const funcT& rho)
    {
      // Nuclear and coulomb potentials
      funcT vc = apply(*_cop, rho);
      funcT vlocal = _vnuc + vc;
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
        }
        else
        {
          // potential
          funcT vxc = copy(rhoa);
          vxc.unaryop(&::libxc_ldaop);
          pfuncsa = mul_sparse(_world, vlocal + vxc, phisa, _params.thresh * 0.1);
          // energy
          funcT fc = copy(rhoa);
          fc.unaryop(&::ldaeop);
          xc = fc.trace();
        }
      }
      std::cout.precision(8);
      print("Energies:");
      print("Kinetic energy:\t\t ", ke);
      print("Potential energy:\t ", pe);
      print("Coulomb energy:\t\t ", ce);
      print("Exchage energy:\t\t ", xc, "\n");
      print("Total energy:\t\t ", ke + pe + ce + xc, "\n\n");
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
        // Compute density
        _rhoa = compute_rho(_phisa);
        _rhob = (_params.spinpol) ? compute_rho(_phisb) : _rhoa;
        _rho = _rhoa + _rhob;

        vector<funcT> pfuncsa(_phisa.size()), pfuncsb(_phisb.size());
        for (unsigned int pi = 0; pi < _phisa.size(); pi++)
          pfuncsa[pi] = FunctionFactory<T, NDIM>(_world);
        for (unsigned int pi = 0; pi < _phisb.size(); pi++)
          pfuncsb[pi] = FunctionFactory<T, NDIM>(_world);

        // Apply the potentials to the orbitals
        apply_potential(pfuncsa, pfuncsb, _phisa, _phisb, _rhoa, _rhob, _rho);

        // Make BSH Green's function
        std::vector<poperatorT> bopsa = make_bsh_operators(_eigsa);
        vector<double> sfactor(pfuncsa.size());
        for (unsigned int si = 0; si < sfactor.size(); si++) sfactor[si] = -2.0;
        scale(_world, pfuncsa, sfactor);

        // Apply Green's function to orbitals
        vector<funcT> tmpa = apply(_world, bopsa, pfuncsa);

//        {
//          if (_world.rank() == 0) printf("\n");
//          Tensor<double> boxsize = FunctionDefaults<NDIM>::get_cell_width();
//          double L = boxsize[0];
//          double bstep = L / 100.0;
//          _phisa[0].reconstruct();
//          pfuncsa[0].reconstruct();
//          tmpa[0].reconstruct();
//          for (int i = 0; i < 101; i++)
//          {
//           Vector<T,NDIM> p(-L / 2 + i * bstep);
//           if (_world.rank() == 0)
//             printf("%.2f\t\t%.8f\t%.8f\t%.8f\n", p[0], _phisa[0](p), pfuncsa[0](p), tmpa[0](p));
//          }
//          if (_world.rank() == 0) printf("\n");
//        }

        // Gram-Schmidt
        gram_schmidt(tmpa, _phisa);

        // Update eigenvalues
        update_eigenvalues(tmpa, pfuncsa, _phisa, _eigsa);

        // Update orbitals
        truncate(_world, tmpa);
        for (unsigned int ti = 0; ti < tmpa.size(); ti++)
        {
          _phisa[ti] = tmpa[ti].scale(1.0/tmpa[ti].norm2());
        }

        // Do other spin
        if (_params.spinpol)
        {
          std::vector<poperatorT> bopsb = make_bsh_operators(_eigsb);
          scale(_world, pfuncsb, sfactor);
          vector<funcT> tmpb = apply(_world, bopsb, pfuncsb);
          gram_schmidt(tmpb, _phisb);
          // Update orbitals
          truncate(_world, tmpb);
          update_eigenvalues(tmpb, pfuncsb, _phisb, _eigsb);
          for (unsigned int ti = 0; ti < tmpb.size(); ti++)
          {
            _phisb[ti] = tmpb[ti].scale(1.0/tmpb[ti].norm2());
          }
        }

        std::cout.precision(8);
        if (_world.rank() == 0)
        {
          print("Iteration: ", it, "\nEigenvalues for alpha spin: \n");
          for (unsigned int i = 0; i < _eigsa.size(); i++)
          {
            print(_eigsa[i]);
          }
          print("\n\n");
        }
        if (_params.spinpol)
        {
          if (_world.rank() == 0)
          {
            print("Eigenvalues for beta spin: \n");
            for (unsigned int i = 0; i < _eigsb.size(); i++)
            {
              print(_eigsb[i]);
            }
            print("\n\n");
          }
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    void gram_schmidt(std::vector<funcT>& a, const std::vector<funcT>& b)
    {
      // Do Gram-Schmidt
      if (_world.rank() == 0) printf("Gram-Schmidt ...\n\n");
      for (unsigned int ai = 0; ai < a.size(); ++ai)
      {
        // Project out the lower states
        for (unsigned int bj = 0; bj < ai; ++bj)
        {
          double overlap = inner(a[ai], b[bj]);
          a[ai] -= overlap*b[bj];
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    void update_eigenvalues(const std::vector<funcT>& tmp,
        const std::vector<funcT>& pfuncs, const std::vector<funcT>& phis,
        std::vector<double>& eigs)
    {
      // Update e
      if (_world.rank() == 0) printf("Updating e ...\n\n");
      for (unsigned int ei = 0; ei < eigs.size(); ei++)
      {
        funcT r = tmp[ei] - phis[ei];
        double tnorm = tmp[ei].norm2();
        // Compute correction to the eigenvalues
        double ecorrection = -0.5*inner(pfuncs[ei], r) / (tnorm*tnorm);
        double eps_old = eigs[ei];
        double eps_new = eps_old + ecorrection;
//        if (_world.rank() == 0) printf("ecorrection = %.8f\n\n", ecorrection);
//        if (_world.rank() == 0) printf("eps_old = %.8f eps_new = %.8f\n\n", eps_old, eps_new);
        // Sometimes eps_new can go positive, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
        // I bounce the new eigenvalue back into the negative side of the real axis. I
        // keep doing this until it's good or I've already done it 10 times.
        int counter = 10;
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
//    funcT get_phi(int indx)
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
//    const std::vector<funcT>& phis()
//    {
//      return _solver->phis();
//    }
//    //*************************************************************************

  };
  //***************************************************************************

}
#define SOLVER_H_


#endif /* SOLVER_H_ */
