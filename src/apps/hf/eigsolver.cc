#define WORLD_INSTANTIATE_STATIC_TEMPLATES  

#include "eigsolver.h"
#include "util.h"

using std::cout;
using std::endl;

namespace madness
{

  //***************************************************************************
  EigSolver::EigSolver(World& world, std::vector<funcT> phis, 
      std::vector<double> eigs, std::vector<EigSolverOp*> ops, double thresh, 
      bool periodic)
  : _phis(phis), _eigs(eigs), _ops(ops), _world(world), _thresh(thresh)
  {
    _rho = EigSolver::compute_rho(phis, world);
    _periodic = periodic;
  }
  //***************************************************************************
  
  //***************************************************************************
      EigSolver::~EigSolver()
  {
    // Eigsolver is responsible for deleting the ops
    for (std::vector<EigSolverOp*>::iterator it = _ops.begin(); it != _ops.end(); 
      it++) delete (*it);
    _ops.clear();
    // Clear eigenvectors
    _phis.clear();
    // Clear eigenvalues
    _eigs.clear();
    // Clear observers
    _obs.clear();
  }
  //***************************************************************************
  
  //***************************************************************************
  funcT EigSolver::compute_rho(std::vector<funcT> phis, const World& world)
  {
    // Electron density
    funcT rho = FunctionFactory<double,3>(const_cast<World&>(world));
    // Loop over all wavefunctions to compute density
    for (std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
    {
      // Get phi(j) from iterator
      const funcT& phij = (*pj);
      // Compute the j-th density
      funcT prod = square(phij);
      rho += prod;
    }
    rho.truncate();
    return rho;
  }
  //***************************************************************************

  //***************************************************************************
  double EigSolver::matrix_element(const funcT& phii, const funcT& phij)
  {
    double value = 0.0;
    // Kinetic energy operator
    for (int axis = 0; axis < 3; axis++)
    {
      funcT dpsi_j = diff(phij, axis);
      funcT dpsi_i = diff(phii, axis);
      value += 0.5 * inner(dpsi_i, dpsi_j);
    }
    // Loop through all ops
    for (unsigned int oi = 0; oi < _ops.size(); oi++)
    {
      EigSolverOp* op = _ops[oi];
      // Operate with density-dependent operator
      if (op->is_rd()) value += op->coeff() * phii.inner(op->op_r(_rho, phij));
      // Operate with orbital-dependent operator
      if (op->is_od()) value += op->coeff() * phii.inner(op->op_o(_phis, phij));
    }
    return value;
  }
  //***************************************************************************

  //***************************************************************************
  void EigSolver::print_matrix_elements(const funcT& phii, const funcT& phij)
  {
    double value = 0.0;
    // Kinetic energy operator
    for (int axis = 0; axis < 3; axis++)
    {
      funcT dpsi_j = diff(phij, axis);
      funcT dpsi_i = diff(phii, axis);
      value += 0.5 * inner(dpsi_i, dpsi_j);
    }
    if (_world.rank() == 0)
    {
      cout << "***** Evaluation of matrix elements *****" << endl;
      cout << "KineticEnergyOp:\t\t\t" << value << endl;
    }
    
    // Loop through all ops
    for (unsigned int oi = 0; oi < _ops.size(); oi++)
    {
      value = 0.0;
      EigSolverOp* op = _ops[oi];
      // Operate with density-dependent operator
      if (op->is_rd()) value += op->coeff() * phii.inner(op->op_r(_rho, phij));
      // Operate with orbital-dependent operator
      if (op->is_od()) value += op->coeff() * phii.inner(op->op_o(_phis, phij));
      if (_world.rank() == 0)
      {
        cout << op->messsageME() << ":\t\t\t" << value << endl;
      }
    }
    if (_world.rank() == 0) printf("\n\n");
  }
  //***************************************************************************

  //***************************************************************************
  void EigSolver::solve(int maxits)
  {
    for (int it = 0; it < maxits; it++)
    {
      if (_world.rank() == 0) printf("Iteration #%d\n\n", it);
      for (unsigned int pi = 0; pi < _phis.size(); pi++)
      {
        // Get psi from collection
        funcT psi = _phis[pi];
        funcT pfunc = FunctionFactory<double,3>(_world);
        // Loop through all ops
        if (_world.rank() == 0) printf("Looping through the ops ...\n\n");
        for (unsigned int oi = 0; oi < _ops.size(); oi++)
        {
          EigSolverOp* op = _ops[oi];
          // Operate with density-dependent operator
          if (op->is_rd()) pfunc += op->coeff() * op->op_r(_rho, psi);
          // Operate with orbital-dependent operator
          if (op->is_od()) pfunc += op->coeff() * op->op_o(_phis, psi);
        }
        pfunc.scale(-2.0).truncate(_thresh);
        if (_world.rank() == 0) printf("Creating BSH operator ...\n\n");
        SeparatedConvolution<double,3>* op = 0;
        if (_periodic)
        {
          op = BSHOperatorPtr<double,3>(_world, sqrt(-2.0*_eigs[pi]), 
              FunctionDefaults<3>::get_k(), 1e-3, _thresh);      
        }
        else
        {
          op = BSHOperatorPtr<double,3>(_world, sqrt(-2.0*_eigs[pi]), 
              FunctionDefaults<3>::get_k(), 1e-3, _thresh);      
        }
        // Apply the Green's function operator (stubbed)
          if (_world.rank() == 0) printf("Applying BSH operator ...\n\n");
        funcT tmp = apply(*op, pfunc);
        // (Not sure whether we have to do this mask thing or not!)
        if (_world.rank() == 0) printf("Gram-Schmidt ...\n\n");
        for (unsigned int pj = 0; pj < pi; ++pj)
        {
//          // Make sure that pi != pj
//          MADNESS_ASSERT(pi != pj);
          // Project out the lower states
          // Get other wavefunction
          funcT psij = _phis[pj];
          double overlap = inner(tmp, psij);
          tmp -= overlap*psij;
        }
        // Update e
        if (_world.rank() == 0) printf("Updating e ...\n\n");
        funcT r = tmp - psi;
        double norm = tmp.norm2();
        double eps_old = _eigs[pi];
        double ecorrection = -0.5*inner(pfunc, r) / (norm*norm);
        double eps_new = eps_old + ecorrection;
        // Sometimes eps_new can go positive, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
        // I bounce the new eigenvalue back into the negative side of the real axis. I 
        // keep doing this until it's good or I've already done it 10 times.
        int counter = 0;
        while (eps_new >= 0.0 && counter < 10)
        {
          // Split the difference between the new and old estimates of the 
          // pi-th eigenvalue.
          eps_new = eps_old + 0.5*(eps_new - eps_old);
          counter++;
        }
        // Still no go, forget about it. (1$ to Donnie Brasco)
        if (eps_new >= 0.0)
        {
          printf("FAILURE OF WST: exiting!!\n\n");
          _exit(0);
        }
        // Update the eigenvalue estimates and wavefunctions.
        tmp.truncate(_thresh);
        _eigs[pi] = eps_new;
        _phis[pi] = tmp.scale(1.0/tmp.norm2());
      }
      // Update rho
//      if (_world.rank() == 0) printf("Computing new density for it == #%d\n\n", it);
      _rho = EigSolver::compute_rho(_phis, _world);
      // Output to observables
      for (std::vector<IEigSolverObserver*>::iterator itr = _obs.begin(); itr
        != _obs.end(); ++itr)
      {
        (*itr)->iterateOutput(_phis, _eigs, _rho, it);
      }
    }
  }
  //***************************************************************************
}
