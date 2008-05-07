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
  }
  //*************************************************************************

  //*************************************************************************
  HartreeFockCoulombOp::HartreeFockCoulombOp(World& world, double coeff,
      double thresh) : EigSolverOp(world, coeff, thresh)
  {
  }
  //*************************************************************************
  
  //*************************************************************************
  HartreeFockExchangeOp::HartreeFockExchangeOp(World& world, double coeff,
      double thresh) : EigSolverOp(world, coeff, thresh)
  {
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
  funcT HartreeFockCoulombOp::op_r(const funcT& rho, const funcT& psi)
  {
    // Create Coulomb operator
    SeparatedConvolution<double,3> cop = 
      CoulombOperator<double,3>(world(), FunctionDefaults<3>::get_k(), 1e-4, thresh());      
    // Transform Coulomb operator into a function
//    if (isPrintingNode()) printf("rho.norm2() = %.5f\n\n", rho.norm2()); 
//    if (isPrintingNode()) printf("Applying Coulomb operator to density ...\n\n");
    // Apply the Coulomb operator
    funcT Vc = apply(cop, rho);
    funcT rfunc = Vc*psi;
//    if (isPrintingNode()) printf("Vc.norm2() = %.5f\n\n", Vc.norm2()); 
//    if (isPrintingNode()) printf("pcoulomb.norm2() = %.5f\n\n", rfunc.norm2()); 
    return  rfunc;
  }
  //*************************************************************************
  
  //*************************************************************************
  funcT HartreeFockExchangeOp::op_o(const std::vector<funcT>& phis, const funcT& psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,3>(world());
    // Create Coulomb operator
    SeparatedConvolution<double,3> cop = CoulombOperator<double, 3>(world(),
        FunctionDefaults<3>::get_k(), 1e-4, thresh());
    // Use the psi and pj wavefunctions to build a product so that the K 
    // operator can be applied to the wavefunction indexed by pj, NOT PSI.
    for (std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
    {
      // Get phi(j) from iterator
      const funcT& phij = (*pj);
      // NOTE that psi is involved in this calculation
      funcT prod = phij*psi;
//      if (isPrintingNode()) printf("prod.norm2() = %.5f\n\n", prod.norm2());
      // Transform Coulomb operator into a function (stubbed)
      prod.truncate(thresh());
      funcT Vex = apply(cop, prod);
//      if (isPrintingNode()) printf("Vex.norm2() = %.5f\n\n", Vex.norm2());
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
      funcT density = FunctionFactory<double,3>(_world);
      // Create Coulomb operator
      SeparatedConvolution<double,3> op = 
        CoulombOperator<double,3>(_world, FunctionDefaults<3>::get_k(), 1e-4, _thresh);      
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
    funcT rfunc = FunctionFactory<double,3>(world());
    if (include_exchange())
    {
      // Create Coulomb operator
      SeparatedConvolution<double,3> op = CoulombOperator<double, 3>(world(),
          FunctionDefaults<3>::get_k(), 1e-4, thresh());
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
      const std::vector<double>& eigs,const funcT& rho, const int& iter)
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
      if (world().rank() == 0) printf("gs eigv: = \t\t\t %.4f\n", eigs[0]);
    }
  }
  //***************************************************************************
}
//*****************************************************************************

