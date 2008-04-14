#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include "hartreefock.h"
#include "util.h"

namespace madness
{

  //*************************************************************************
  HartreeFockNuclearPotentialOp::HartreeFockNuclearPotentialOp(World& world,
    funcT V, double coeff, double thresh) :
  EigSolverOp(world, coeff, thresh)
  {
    _V = V;
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //*************************************************************************

  //*************************************************************************
  HartreeFockCoulombOp::HartreeFockCoulombOp(World& world, double coeff,
      double thresh) : EigSolverOp(world, coeff, thresh)
  {
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //*************************************************************************
  
  //*************************************************************************
  HartreeFockExchangeOp::HartreeFockExchangeOp(World& world, double coeff,
      double thresh) : EigSolverOp(world, coeff, thresh)
  {
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //*************************************************************************
  
  //*************************************************************************
  funcT HartreeFockNuclearPotentialOp::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _V*psi;
    return rfunc;
  }
  //*************************************************************************

  //*************************************************************************
//  funcT HartreeFockCoulombOp::op(const std::vector<funcT>& phis, const funcT& psi)
//  {
//    // Electron density
//    funcT density = FunctionFactory<double,3>(world()).functor(zeros);
//    // Create Coulomb operator
//    SeparatedConvolution<double,3> cop = 
//      CoulombOperator<double,3>(world(), FunctionDefaults<3>::k, 1e-4, thresh());      
//    for (std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
//    {
//      // Get phi(j) from iterator
//      const funcT& phij = (*pj);
//      // Compute the j-th density
//      funcT prod = square(phij);
//      density += prod;
//    }
//    // Transform Coulomb operator into a function (stubbed)
//    if (isPrintingNode()) printf("density.norm2() = %.5f\n\n", density.norm2()); 
//    density.truncate();
//    if (isPrintingNode()) printf("Applying Coulomb operator to density ...\n\n");
//    funcT Vc = apply(cop, density);
//    // Note that we are not using psi
//    // The density is built from all of the wavefunctions. The contribution
//    // psi will be subtracted out later during the exchange.
//    funcT rfunc = Vc*psi;
//    if (isPrintingNode()) printf("Vc.norm2() = %.5f\n\n", Vc.norm2()); 
//    if (isPrintingNode()) printf("pcoulomb.norm2() = %.5f\n\n", rfunc.norm2()); 
//    return  rfunc;
//  }
  //*************************************************************************
  
  //*************************************************************************
  funcT HartreeFockCoulombOp::op_r(const funcT& rho, const funcT& psi)
  {
    // Create Coulomb operator
    SeparatedConvolution<double,3> cop = 
      CoulombOperator<double,3>(world(), FunctionDefaults<3>::k, 1e-4, thresh());      
    // Transform Coulomb operator into a function
    if (isPrintingNode()) printf("rho.norm2() = %.5f\n\n", rho.norm2()); 
    if (isPrintingNode()) printf("Applying Coulomb operator to density ...\n\n");
    // Apply the Coulomb operator
    funcT Vc = apply(cop, rho);
    funcT rfunc = Vc*psi;
    if (isPrintingNode()) printf("Vc.norm2() = %.5f\n\n", Vc.norm2()); 
    if (isPrintingNode()) printf("pcoulomb.norm2() = %.5f\n\n", rfunc.norm2()); 
    return  rfunc;
  }
  //*************************************************************************
  
  //*************************************************************************
  funcT HartreeFockExchangeOp::op_o(const std::vector<funcT>& phis, const funcT& psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,3>(world()).functor(zeros);
    // Create Coulomb operator
    SeparatedConvolution<double,3> cop = CoulombOperator<double, 3>(world(),
        FunctionDefaults<3>::k, 1e-4, thresh());
    // Use the psi and pj wavefunctions to build a product so that the K 
    // operator can be applied to the wavefunction indexed by pj, NOT PSI.
    for (std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
    {
      // Get phi(j) from iterator
      const funcT& phij = (*pj);
      // NOTE that psi is involved in this calculation
      funcT prod = phij*psi;
      if (isPrintingNode()) printf("prod.norm2() = %.5f\n\n", prod.norm2());
      // Transform Coulomb operator into a function (stubbed)
      prod.truncate(thresh());
      funcT Vex = apply(cop, prod);
      if (isPrintingNode()) printf("Vex.norm2() = %.5f\n\n", Vex.norm2());
      // NOTE that the index is j.
      rfunc += Vex*phij;
    }
    return rfunc;
    
  }
  //*************************************************************************
  
  //*************************************************************************
  HartreeFock::HartreeFock(World& world, funcT V, std::vector<funcT> phis,
    std::vector<double> eigs, bool bCoulomb, bool bExchange, double thresh)
   : _world(world), _V(V), _thresh(thresh)
  {
    _bCoulomb = bCoulomb;
    _bExchange = bExchange;
    // Create ops list 
    std::vector<EigSolverOp*> ops;
    // Add nuclear potential to ops list
    ops.push_back(new HartreeFockNuclearPotentialOp(world, V, 1.0, thresh));
    // Check for coulomb and exchange, and add as appropriate
    if (bCoulomb)
    {
      ops.push_back(new HartreeFockCoulombOp(world, 2.0, thresh));
    }
    if (bExchange)
    {
      ops.push_back(new HartreeFockExchangeOp(world, -1.0, thresh));
    }
    // Create solver
    _solver = new EigSolver(world, phis, eigs, ops, thresh);
    _solver->addObserver(this);

    // Misc.
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //***************************************************************************
  
  //***************************************************************************
  HartreeFock::HartreeFock(World& world, funcT V, funcT phi, double eig, 
    bool bCoulomb, bool bExchange, double thresh) : _world(world), _V(V), 
    _thresh(thresh)
  {
    // Create temporary list for eigsolver
    std::vector<funcT> phis;
    std::vector<double> eigs;
    phis.push_back(phi);
    eigs.push_back(eig);
    // Create ops list 
    std::vector<EigSolverOp*> ops;
    // Add nuclear potential to ops list
    ops.push_back(new HartreeFockNuclearPotentialOp(world, V, 1.0, thresh));
    // Check for coulomb and exchange, and add as appropriate
    if (bCoulomb)
    {
      ops.push_back(new HartreeFockCoulombOp(world, 2.0, thresh));
    }
    if (bExchange)
    {
      ops.push_back(new HartreeFockExchangeOp(world, -1.0, thresh));
    }
    // Create solver
    _solver = new EigSolver(world, phis, eigs, ops, thresh);
    _solver->addObserver(this);

    // Misc.
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //***************************************************************************
  
  //***************************************************************************
  HartreeFock::~HartreeFock()
  {
    delete _solver;
  }
  //***************************************************************************
  
  //***************************************************************************
  void HartreeFock::hartree_fock(int maxits)
  {
    _solver->solve(maxits);
  }
  //***************************************************************************

//  //***************************************************************************
//  void HartreeFock::hartree_fock(int maxits)
//  {
//    cout << "bCoulomb is " << _bCoulomb << endl;    
//    cout << "bExchange is " << _bExchange << endl;    
//    for (int it = 0; it < maxits; it++)
//    {
//      printf("//************* iteration #%d *************//\n\n", it);
//      printf("thresh = %.4e\n\n", _tWorld& world() {return _world;}hresh);
//      for (unsigned int pi = 0; pi < _phis.size(); pi++)
//      {
//        // Get psi from collection
//        funcT psi = _phis[pi];
//        // Calculate nuclear contribution
//	//madness::print("THIS IS V");
//        //_V.print_tree();
//	//madness::print("THIS IS PSI");
//        //psi.print_tree();
//        printf("iteration #%d: calc nuclear ...\n\n", it);
//        printf("iteration #%d: psi.norm2() = %.5f\n\n", it, 
//          psi.norm2());
//        funcT pnuclear = _V*psi;
//        printf("iteration #%d: pnuclear.norm2() = %.5f\n\n", it, 
//          pnuclear.norm2());
//        // Calculate the Coulomb contribution to the Fock operator (J)
//        printf("iteration #%d: calc coulomb ...\n\n", it);
//        funcT pcoulomb = calculate_coulomb(psi);
//        // Calculate the Exchange contribution to the Fock operator (K)
//        printf("iteration #%d: calc exchange ...\n\n", it);
//        funcT pexchange = calculate_exchange(psi);
//        // Get new wavefunction
//        funcT pfunc = pnuclear + 2.0 * pcoulomb - pexchange;
//        pfunc.scale(-2.0).truncate(_thresh);
//        printf("iteration #%d: pfunc.norm2() = %.5f\n\n", it, 
//          pfunc.norm2());
//        // Create the free-particle Green's function operator
//        printf("iteration #%d: _eigs[%d] = %.5f\n\n", it, pi, 
//          _eigs[pi]);
//        SeparatedConvolution<double,3> op = 
//          BSHOperator<double,3>(_world, sqrt(-2.0*_eigs[pi]), 
//              FunctionDefaults<3>::k, 1e-3, _thresh);      
//        // Apply the Green's function operator (stubbed)
//        printf("iteration #%d: apply BSH ...\n\n", it);
//        funcT tmp = apply(op, pfunc);
//        printf("iteration #%d (after BSH): tmp.norm2() = %.5f\n\n", it, 
//          tmp.norm2());
//        // (Not sure whether we have to do this mask thing or not!)
//        printf("iteration #%d: doingWorld& world() {return _world;} gram-schmidt ...\n\n", it);
//        for (unsigned int pj = 0; pj < pi; ++pj)
//        {
//          // Project out the lower states
//          // Make sure that pi != pj
//          if (pi != pj)
//          {
//            // Get other wavefunction
//            funcT psij = _phis[pj];
//            double overlap = inner(tmp, psij);
//            tmp -= overlap*psij;
//          }
//        }
//        // Update e
//        funcT r = tmp - psi;
//        double norm = tmp.norm2();
//        double eps_old = _eigs[pi];
//        printf("Updating wavefunction on iteration #%d ...\n\n", it);
//        double ecorrection = -0.5*inner(pfunc, r) / (norm*norm);
//        double eps_new = eps_old + ecorrection;
//        printf("wavefunction #%d: ecorrection = %.5f eps_new = %.5f\n", pi, ecorrection, eps_new);
//        // Sometimes eps_new can go posivite, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
//        // I bounce the new eigenvalue back into the negative side of the real axis. I 
//        // keep doing this until it's good or I've already done it 10 times.
//        int counter = 0;
//        while (eps_new >= 0.0 && counter < 10)
//        {
//          printf("wavefunction #%d: eps_new = %.5f\n", pi, eps_new);
//          // Split the difference between the new and old estimates of the 
//          // pi-th eigenvalue.
//          eps_new = eps_old + 0.5*(eWorld& world() {return _world;}ps_new - eps_old);
//          counter++;
//        }
//        // Still no go, forget about it. (1$ to Donnie Brasco)
//        if (eps_new >= 0.0)
//          {
//            printf("FAILURE OF WST: exiting!!\n\n");
//            _exit(0);
//          }
//        // Update the eigenvalue estimates and wavefunctions.
//        tmp.truncate(_thresh);
//        _eigs[pi] = eps_new;
//        _phis[pi] = tmp.scale(1.0/tmp.norm2());
//        printf("iteration #%d: tmp(/psi).norm2() = %.5f\n\n", it, 
//          tmp.norm2());
//      }
//      // Display energies
//      if (it%3 == 0)
//      {
//        printf("Calculating energies ...\n");
//        printf("Calculating KE ...\n");
//        double ke = 2.0 * calculate_tot_ke_sp();
//        printf("Calculating PE ...\n");
//        double pe = 2.0 * calculate_tot_pe_sp();
//        printf("Calculating CE ...\n");
//        double ce = calculate_tot_coulomb_energy();
//        printf("Calculating EE ...\n");
//        double ee = calculate_tot_exchange_energy();
//        printf("Calculating NE ...\n");
//        double ne = 0.0;
//        printf("Kinetic energy:\t\t\t %.8f\n", ke);
//        printf("Potential energy:\t\t %.8f\n", pe);
//        printf("Two-electron energy:\t\t %.8f\n", 2.0*ce - ee);
//        printf("Total energy:\t\t\t %.8f\n", ke + pe + 2.0*ce - ee + ne);
//      }
//    }World& world() {return _world;}
//  }
//  //***************************************************************************

//  //***************************************************************************
//  funcT HartreeFock::calculate_coulomb(funcT psi)
//  {
//    if (include_coulomb())
//    {
//      // Electron density
//      funcT density = FunctionFactory<double,3>(_world).functor(zeros);
//      // Create Coulomb operator
//      SeparatedConvolution<double,3> op = 
//        CoulombOperator<double,3>(_world, FunctionDefaults<3>::k, 1e-4, _thresh);      
//      for (std::vector<funcT>::iterator pj = _phis.begin(); pj != _phis.end(); ++pj)
//      {
//        // Get phi(j) from iterator
//        funcT& phij = (*pj);
//        // Compute the j-th density
//        funcT prod = square(phij);
//        density += prod;
//      }
//      // Transform Coulomb operator into a function (stubbed)
//      printf("density.norm2() = %.5f\n\n", density.norm2()); 
//      density.truncate();
//      printf("Applying Coulomb operator to density ...\n\n");
//      funcT Vc = apply(op, density);
//      // Note that we are not using psi
//      // The density is built from all of the wavefunctions. The contribution
//      // psi will be subtracted out later during the exchange.
//      funcT rfunc = Vc*psi;
//      printf("Vc.norm2() = %.5f\n\n", Vc.norm2()); 
//      printf("pcoulomb.norm2() = %.5f\n\n", rfunc.norm2()); 
//      return rfunc;
//    }
//    return FunctionFactory<double,3>(_world).functor(zeros);
//  }
//  //***************************************************************************
//
//  //***************************************************************************
//    funcT HartreeFock::calculate_exchange(funcT psi)
//  {
//    // Return value
//    funcT rfunc = FunctionFactory<double,3>(_world).functor(zeros);
//    if (include_exchange())
//    {
//      // Create Coulomb operator
//      SeparatedConvolution<double,3> op = CoulombOperator<double, 3>(_world,
//          FunctionDefaults<3>::k, 1e-4, _thresh);
//      // Use the psi and pj wavefunctions to build a product so that the K 
//      // operator can be applied to the wavefunction indexed by pj, NOT PSI.
//      for (std::vector<funcT>::iterator pj = _phis.begin(); pj != _phis.end(); ++pj)
//      {
//        // Get phi(j) from iterator
//        funcT& phij = (*pj);
//        // NOTE that psi is involved in this calculation
//        funcT prod = phij*psi;
//        printf("prod.norm2() = %.5f\n\n", prod.norm2());
//        // Transform Coulomb operator into a function (stubbed)
//        prod.truncate(_thresh);
//        funcT Vex = apply(op, prod);
//        printf("Vex.norm2() = %.5f\n\n", Vex.norm2());
//        // NOTE that the index is j.
//        rfunc += Vex*phij;
//      }
//    }
//    return rfunc;
//  }
//  //***************************************************************************

  //***************************************************************************
  double HartreeFock::calculate_ke_sp(funcT psi)
  {
    double kenergy = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
      funcT dpsi = diff(psi, axis);
      kenergy += 0.5 * inner(dpsi, dpsi);
    }
    return kenergy;
  }
  //***************************************************************************

  //***************************************************************************
  double HartreeFock::calculate_pe_sp(funcT psi)
  {
    funcT vpsi = _V*psi;
    return vpsi.inner(psi);
  }
  //***************************************************************************

  //***************************************************************************
  double HartreeFock::calculate_coulomb_energy(const std::vector<funcT>& phis, const funcT& psi)
  {
    if (include_coulomb())
    {
      // Electron density
      funcT density = FunctionFactory<double,3>(_world).functor(zeros);
      // Create Coulomb operator
      SeparatedConvolution<double,3> op = 
        CoulombOperator<double,3>(_world, FunctionDefaults<3>::k, 1e-4, _thresh);      
      for (std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
      {
        // Get phi(j) from iterator
        const funcT& phij = (*pj);
        // Compute the j-th density
        funcT prod = phij*phij;
        density += prod;
      }
      // Transform Coulomb operator into a function (stubbed)
      density.truncate(_thresh);
      funcT Vc = apply(op, density);
      // Note that we are not using psi
      // The density is built from all of the wavefunctions. The contribution
      // psi will be subtracted out later during the exchange.
      funcT vpsi = Vc*psi;
      return inner(vpsi, psi);
    }
    return 0.0;
  }
  //***************************************************************************

  //***************************************************************************
  double HartreeFock::calculate_exchange_energy(const std::vector<funcT>& phis, const funcT& psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,3>(world()).functor(zeros);
    if (include_exchange())
    {
      // Create Coulomb operator
      SeparatedConvolution<double,3> op = CoulombOperator<double, 3>(world(),
          FunctionDefaults<3>::k, 1e-4, thresh());
      // Use the psi and pj wavefunctions to build a product so that the K 
      // operator can be applied to the wavefunction indexed by pj, NOT PSI.
      for (std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
      {
        // Get phi(j) from iterator
        const funcT& phij = (*pj);
        // NOTE that psi is involved in this calculation
        funcT prod = phij*psi;
        // Transform Coulomb operator into a function (stubbed)
        funcT Vex = apply(op, prod);
        // NOTE that the index is j.
        rfunc += Vex*phij;
      }
    }
    return inner(rfunc, psi);
  }
  //***************************************************************************

  //***************************************************************************
  double HartreeFock::calculate_tot_ke_sp(const std::vector<funcT>& phis)
  {
    double tot_ke = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate kinetic energy contribution from psi
      tot_ke += calculate_ke_sp(psi);
    }
    return tot_ke;
  }
  //***************************************************************************
  
  //***************************************************************************
  double HartreeFock::calculate_tot_pe_sp(const std::vector<funcT>& phis)
  {
    double tot_pe = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate potential energy contribution from psi
      tot_pe += calculate_pe_sp(psi);
    }
    return tot_pe;
  }
  //***************************************************************************
  
  //***************************************************************************
  double HartreeFock::calculate_tot_coulomb_energy(const std::vector<funcT>& phis)
  {
    double tot_ce = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate coulomb energy contribution from psi
      tot_ce += calculate_coulomb_energy(phis, psi);
    }
    return tot_ce;
  }
  //***************************************************************************
  
  //***************************************************************************
  double HartreeFock::calculate_tot_exchange_energy(const std::vector<funcT>& phis)
  {
    double tot_ee = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate exchange energy contribution from psi
      tot_ee += calculate_exchange_energy(phis, psi);
    }
    return tot_ee;
  }
  //***************************************************************************


  //***************************************************************************
  void HartreeFock::iterateOutput(const std::vector<funcT>& phis,
      const std::vector<double>& eigs, const int& iter)
  {
    if (iter%3 == 0)
    {
      if (world().rank() == 0) printf("Calculating energies ...\n");
      if (world().rank() == 0) printf("Calculating KE ...\n");
      double ke = 2.0 * calculate_tot_ke_sp(phis);
      if (world().rank() == 0) printf("Calculating PE ...\n");
      double pe = 2.0 * calculate_tot_pe_sp(phis);
      if (world().rank() == 0) printf("Calculating CE ...\n");
      double ce = calculate_tot_coulomb_energy(phis);
      if (world().rank() == 0) printf("Calculating EE ...\n");
      double ee = calculate_tot_exchange_energy(phis);
      if (world().rank() == 0) printf("Calculating NE ...\n");
      double ne = 0.0;
      if (world().rank() == 0) printf("Kinetic energy:\t\t\t %.8f\n", ke);
      if (world().rank() == 0) printf("Potential energy:\t\t %.8f\n", pe);
      if (world().rank() == 0) printf("Two-electron energy:\t\t %.8f\n", 2.0*ce - ee);
      if (world().rank() == 0) printf("Total energy:\t\t\t %.8f\n", ke + pe + 2.0*ce - ee + ne);
    }
  }
  //***************************************************************************
}
//*****************************************************************************

