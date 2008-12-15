/*
 * solver.h
 *
 *  Created on: Dec 14, 2008
 *      Author: wsttiger
 */

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

    //*************************************************************************
    World& world() {return _world;}
    //*************************************************************************

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
    std::vector<funcT> _eigsa;
    //*************************************************************************

    //*************************************************************************
    std::vector<funcT> _eigsb;
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
       : _world(world), _phisa(phisa), _phisb(phisb), _eigsa(eigsa), _eigsb(eigsb),
       _params(params)
    {
      if (params.periodic)
      {
        Tensor<double> box = FunctionDefaults<NDIM>::get_cell_width();
        _cop = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
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
        _vnuc = apply(_cop, _vnucrhon);
      }
    }
    //*************************************************************************

    //***************************************************************************
    funcT compute_rho(std::vector<phis> phis)
    {
      // Electron density
      funcT rho = FunctionFactory<double,NDIM>(_world);
      // Loop over all wavefunctions to compute density
      for (int j = 0; j < phis.size(); j++)
      {
        // Get phi(j) from iterator
        const funcT& phij = phis[j];
        // Compute the j-th density
        funcT prod = square(phij);
        rho += occs[j]*prod;
      }
      rho.truncate();
      return rho;
    }
    //***************************************************************************

    //***************************************************************************
    std::vector< SeparatedConvolution<T,NDIM> > make_bsh_operators(std::vectors<double> eigs)
    {
      // Make BSH vector
      std::vector< SeparatedConvolution<T,NDIM> > bops;
      // Get defaults
      int k = FunctionDefaults<NDIM>::get_k();
      double tol = FunctionDefaults<NDIM>::get_thresh();
      // Loop through eigenvalues, adding a BSH operator to bops
      // for each eigenvalue
      int sz = _phis.size();
      for (int i = 0; i < sz; i++)
      {
          double eps = _eigs[i];
          if (eps > 0)
          {
              if (_world.rank() == 0)
              {
                  DEBUG_STREAM << "bsh: warning: positive eigenvalue" << i << eps << endl;
              }
              eps = -0.1;
          }
          _bops.push_back(BSHOperator<double,NDIM>(_world, sqrt(-2.0*eps), k, 1e-4, tol));
      }
    }
    //*************************************************************************

    //*************************************************************************
    std::vector<funcT> apply_potential(std::vector<funcT> pfuncsa,
        std::vector<funcT> pfuncsb, std::vector<funcT> phisa,
        std::vector<funcT> phisb, funcT rhoa, funcT rhob, funcT rho)
    {
      // Nuclear and coulomb potentials
      funcT vlocal = _vnuc + apply(_cop, rho);
      // Exchange
      if (_params.functional == 1)
      {
        // LDA, is calculation spin-polarized?
        if (_params.spinpol)
        {
//          funcT vxca = copy(rhoa);
//          funcT vxca = copy(rhob);
//          vxc.unaryop(&::libxc_ldaop);
//          std::vector<funcT> pfuncsa = mul_sparse(_world, vlocal + vxc, phis, _params.thresh * 0.1);
        }
        else
        {
          funcT vxc = copy(rhoa);
          vxc.unaryop(&::libxc_ldaop);
          std::vector<funcT> pfuncsa = mul_sparse(_world, vlocal + vxc, phis, _params.thresh * 0.1);
        }
      }

      return pfuncs;
    }
    //*************************************************************************

    //*************************************************************************
    virtual ~Solver();
    //*************************************************************************

    //*************************************************************************
    void solve()
    {
      for (int it = 0; it < _params.maxits; it++)
      {
        // Compute density
        _rhoa = compute_rho(_phisa);
        _rhob = (params.spinpol) ? compute_rho(_phisb) : rhoa;
        _rho = _rhoa + _rhob;

        // Apply the potentials to the orbitals
        pfa = apply_potential(phisa, phisb, rhoa, rhob, rho);

        // Make BSH Green's function
        std::vector< SeparatedConvolution<T,DIM> > bopsa = make_bsh_operators(_eigsa);
        vector<double> sfactor(pfuncs.size());
        for (unsigned int si = 0; si < sfactor.size(); si++) sfactor[si] = -2.0;
        scale(_world, pfuncsa, sfactor);

        // Apply Green's function to orbitals
        vector<funcT> tmpa = apply(_world, _bops, pfuncs);

        // Gram-Schmidt
        gram_schmidt(tmpa, _phisa);

        // Update orbitals
        truncate(_world, tmpa);
        for (unsigned int ti = 0; ti < tmpa.size(); ti++)
        {
          _phisa[ti] = tmpa[ti].scale(1.0/tmpa[ti].norm2());
        }

        // Update eigenvalues
        update_eigenvalues(tmpa, pfuncsa, phisa, eigsa);

        // Do other spin
        if (_params.spinpol)
        {
          std::vector< SeparatedConvolution<T,DIM> > bopsb = make_bsh_operators(_eigsb);
          scale(_world, pfuncsb, sfactor);
          vector<funcT> tmpa = apply(_world, _bops, pfuncs);
          gram_schmidt(tmpb, _phisb);
          // Update orbitals
          truncate(_world, tmpb);
          for (unsigned int ti = 0; ti < tmpb.size(); ti++)
          {
            _phisb[ti] = tmpb[ti].scale(1.0/tmpb[ti].norm2());
          }
          update_eigenvalues(tmpb, pfuncsb, phisb, eigsb);
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
          a[ti] -= overlap*b[bj];
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    void update_eigenvalues(const std::vector<funcT>& tmp,
        const std::vector<funcT>& pfuncs, const std::vector<funcT>& phis,
        std::vector<double> eigs)
    {
      // Update e
      if (_world.rank() == 0) printf("Updating e ...\n\n");
      for (unsigned int ei = 0; ei < _eigs.size(); ei++)
      {
        funcT r = tmp[ei] - _phis[ei];
        double tnorm = tmp[ei].norm2();
        double rnorm = r.norm2();
        // Compute correction to the eigenvalues
        double ecorrection = -0.5*inner(pfuncs[ei], r) / (tnorm*tnorm);
        double eps_old = _eigs[ei];
        double eps_new = eps_old + ecorrection;
        // Sometimes eps_new can go positive, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
        // I bounce the new eigenvalue back into the negative side of the real axis. I
        // keep doing this until it's good or I've already done it 10 times.
        int counter = 10;
        while (eps_new >= 0.0 && counter < 20)
        {
          // Split the difference between the new and old estimates of the
          // pi-th eigenvalue.
          eps_new = eps_old + 0.33 * (eps_new - eps_old);
          counter++;
        }
        // Still no go, forget about it. (1$ to Donnie Brasco)
        if (eps_new >= 0.0)
        {
          if (_world.rank() == 0) printf("FAILURE OF WST: exiting!!\n\n");
          _exit(0);
        }
        // Set new eigenvalue
        _eigs[ei] = eps_new;
      }
    }
    //*************************************************************************

    //*************************************************************************
    double get_eig(int indx)
    {
      return _solver->get_eig(indx);
    }
    //*************************************************************************

    //*************************************************************************
    funcT get_phi(int indx)
    {
      return _solver->get_phi(indx);
    }
    //*************************************************************************

    //*************************************************************************
    const std::vector<double>& eigs()
    {
      return _solver->eigs();
    }
    //*************************************************************************

    //*************************************************************************
    const std::vector<funcT>& phis()
    {
      return _solver->phis();
    }
    //*************************************************************************

  private:

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    // This variable could either be a nuclear potiential or a nuclear charge
    // density depending on the "ispotential" variable in the
    // ElectronicStructureParams class.
    Function<double,NDIM> _vnucrhon;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

    //*************************************************************************
    World& world() {return _world;}
    //*************************************************************************

  };
  //***************************************************************************

}
#define SOLVER_H_


#endif /* SOLVER_H_ */
