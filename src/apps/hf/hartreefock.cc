#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include "hartreefock.h"

namespace madness
{

  //***************************************************************************
  HartreeFock::HartreeFock(World& world, funcT V, std::vector<funcT> phis,
    std::vector<double> eigs, bool bCoulomb, bool bExchange) :
  _world(world), _V(V), _phis(phis), _eigs(eigs), _bCoulomb(bCoulomb),
      _bExchange(bExchange)
  {
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //***************************************************************************
  
  //***************************************************************************
  HartreeFock::HartreeFock(World& world, funcT V, funcT phi, double eig, 
      bool bCoulomb, bool bExchange) : _world(world), _V(V),
      _bCoulomb(bCoulomb), _bExchange(bExchange)
  {
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
    _phis.push_back(phi);
    _eigs.push_back(eig);
  }
  //***************************************************************************
  
  //***************************************************************************
  HartreeFock::~HartreeFock()
  {
  }
  //***************************************************************************
  
  //***************************************************************************
  void HartreeFock::hartree_fock(int maxits)
  {
    for (int it = 0; it < maxits; it++)
    {
      for (int pi = 0; pi < _phis.size(); pi++)
      {
        // Get psi from collection
        funcT psi = _phis[pi];
        // Calculate nuclear contribution
	madness::print("THIS IS V");
        _V.print_tree();
	madness::print("THIS IS PSI");
        psi.print_tree();
        funcT pnuclear = _V*psi;
        // Calculate the Coulomb contribution to the Fock operator (J)
        funcT pcoulomb = calculate_coulomb(psi);
        // Calculate the Exchange contribution to the Fock operator (K)
        funcT pexchange = calculate_exchange(psi);
        // Get new wavefunction
        funcT pfunc = pnuclear + 2.0 * pcoulomb - pexchange;
        // Create the free-particle Green's function operator
        SeparatedConvolution<double,3> op = 
          BSHOperator<double,3>(_world, sqrt(-2.0*_eigs[pi]), 
              FunctionDefaults<3>::k, 1e-4, _thresh);      
        // Apply the Green's function operator (stubbed)
        funcT tmp = apply(op, pfunc);
        // (Not sure whether we have to do this mask thing or not!)
        for (int pj = 0; pj < pi; ++pj)
        {
          // Project out the lower states
          // Make sure that pi != pj
          if (pi != pj)
          {
            // Get other wavefunction
            funcT psij = _phis[pj];
            printf("Calculating overlap on iteration #%d ...\n", it);
            double overlap = inner(tmp, psij);
            tmp -= overlap*psij;
          }
        }
        // Update e
        funcT r = tmp - psi;
        double norm = tmp.norm2();
        double eps_old = _eigs[pi];
        printf("Updating wavefunction on iteration #%d ...\n", it);
        double ecorrection = inner(pfunc, r) / (norm*norm);
        _eigs[pi] += ecorrection;
        _phis[pi] = tmp.scale(1.0/tmp.norm2());
        
      }
    }
  }
  //***************************************************************************

  //***************************************************************************
  funcT HartreeFock::calculate_coulomb(funcT psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,3>(_world).functor(zeros);
    if (include_coulomb())
    {
      // Create Coulomb operator
      SeparatedConvolution<double,3> op = 
        CoulombOperator<double,3>(_world, FunctionDefaults<3>::k, 1e-4, _thresh);      
      for (std::vector<funcT>::iterator pj = _phis.begin(); pj != _phis.end(); ++pj)
      {
        // Get phi(j) from iterator
        funcT& phij = (*pj);
        // Compute the density
        funcT prod = phij*phij;
        // Transform Coulomb operator into a function (stubbed)
        funcT Vc = apply(op, prod);
        // Note that we are not using psi
        // The density is built from all of the wavefunctions. The contribution
        // psi will be subtracted out later during the exchange.
        rfunc += Vc*psi;
      }
    }
    return rfunc;
  }
  //***************************************************************************

  //***************************************************************************
  funcT HartreeFock::calculate_exchange(funcT psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,3>(_world).functor(zeros);
    if (include_exchange())
    {
      // Create Coulomb operator
      SeparatedConvolution<double,3> op = 
        CoulombOperator<double,3>(_world, FunctionDefaults<3>::k, 1e-4, _thresh);      
      // Use the psi and pj wavefunctions to build a product so that the K 
      // operator can be applied to the wavefunction indexed by pj, NOT PSI.
      for (std::vector<funcT>::iterator pi = _phis.begin(); pi != _phis.end(); ++pi)
      {
        for (std::vector<funcT>::iterator pj = _phis.begin(); pj != _phis.end(); ++pj)
        {
          // Get phi(j) from iterator
          funcT& phij = (*pj);
          // NOTE that psi is involved in this calculation
          funcT prod = phij*psi;
          // Transform Coulomb operator into a function (stubbed)
          funcT Vex = apply(op, prod);
          // NOTE that the index is j.
          rfunc += Vex*phij;
        }
      }
    }
    return rfunc;
  }
}
