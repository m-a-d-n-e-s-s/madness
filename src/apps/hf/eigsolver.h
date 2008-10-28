#ifndef EIGSOLVER_H_
#define EIGSOLVER_H_
#include <mra/mra.h>
#include <world/world.h>
#include <vector>

namespace madness
{
//***************************************************************************
/// This is the interface the an observer wishing to receive output must
/// implement. This call back gives the current eigenfunctions, eigenvalues,
/// and the density.
/// This is a test LaTeX formula
/// The Pythagorean theorem is
/// \f[
/// c^2 = a^2 + b^2
/// \f]
template <typename T, int NDIM>
class IEigSolverObserver
{
  typedef Function<T,NDIM> funcT;
public:
  virtual void iterateOutput(const std::vector<funcT>& phis,
      const std::vector<double>& eigs, const Function<double, NDIM>& rho, const int& iter, bool periodic) = 0;

  virtual ~IEigSolverObserver() {};
};
//***************************************************************************

class KPoint
{
public:
  KPoint(double kx, double ky, double kz, double weight)
  {
    _kx = kx; _ky = ky; _kz = kz;
    _weight = weight;
  }

  //*************************************************************************
  double kx() {return _kx;}
  double ky() {return _ky;}
  double kz() {return _kz;}
  //*************************************************************************

  //*************************************************************************
  double weight() {return _weight;}
  //*************************************************************************

private:
  //*************************************************************************
  // the actual k-point
  double _kx;
  double _ky;
  double _kz;
  //*************************************************************************

  //*************************************************************************
  // weight
  double _weight;
  //*************************************************************************

};

//***************************************************************************
template <typename T, int NDIM>
class EigSolverOp
{
  // Typedef's
  typedef Function<T,NDIM> funcT;
public:
  //*************************************************************************
  // Constructor
  EigSolverOp(World& world, double coeff, double thresh)
    :  _world(world), _coeff(coeff), _thresh(thresh) {}
  //*************************************************************************

  //*************************************************************************
  // Destructor
  virtual ~EigSolverOp() {}
  //*************************************************************************

  //*************************************************************************
  /// Is there an orbitally-dependent term?
  virtual bool is_od() = 0;
  //*************************************************************************

  //*************************************************************************
  /// Is there a density-dependent term?
  virtual bool is_rd() = 0;
  //*************************************************************************

  //*************************************************************************
  /// Build the potential from a density if a density-dependent operator.
  virtual void prepare_op(funcT rho) {}
  //*************************************************************************

  //*************************************************************************
  /// Orbital-dependent portion of operator
  virtual funcT op_o(const std::vector<funcT>& phis, const funcT& psi)
  {
    funcT func = FunctionFactory<T,NDIM>(_world);
    return func;
  }
  //*************************************************************************

  //*************************************************************************
  /// Density-dependent portion of operator
  virtual funcT op_r(const funcT& rho, const funcT& psi)
  {
    funcT func = FunctionFactory<T,NDIM>(_world);
    return func;
  }
  //*************************************************************************

  //*************************************************************************
  /// Orbital-dependent portion of operator
  virtual std::vector<funcT> multi_op_o(const std::vector<funcT>& phis)
  {
    // Collection of empty functions
    std::vector<funcT> newphis(phis.size(), FunctionFactory<T,NDIM>(_world));
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      newphis[pi] = op_o(phis, phis[pi]);
    }
    _world.gop.fence();
    return newphis;
  }
  //*************************************************************************

  //*************************************************************************
  /// Density-dependent portion of operator
  virtual std::vector<funcT> multi_op_r(const funcT& rho, const std::vector<funcT>& phis)
  {
    std::vector<funcT> newphis(phis.size(), FunctionFactory<T,NDIM>(_world));
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      newphis[pi] = op_r(rho, phis[pi]);
    }
    _world.gop.fence();
    return newphis;
  }
  //*************************************************************************

  //*************************************************************************
  double coeff() {return _coeff;}
  //*************************************************************************

  //*************************************************************************
  std::string messsageME()
  {
    return _messageME;
  }
  //*************************************************************************

  //*************************************************************************
  World& _world;
  //*************************************************************************

protected:
  //*************************************************************************
  double thresh() {return _thresh;}
  //*************************************************************************

  //*************************************************************************
  void messageME(std::string messageME)
  {
    _messageME = messageME;
  }
  //*************************************************************************

private:
  //*************************************************************************
  double _coeff;
  //*************************************************************************

  //*************************************************************************
  double _thresh;
  //*************************************************************************

  //*************************************************************************
  std::string _messageME;
  //*************************************************************************

};
//***************************************************************************

//***************************************************************************
/// The EigSolver class is the class that is the workhorse of both the Hartree
/// Fock and the DFT algorithms. This class relies on the wrapper class to
/// give it a list of operators to implement as its potential. This should
/// allow for much more reuse.
template <typename T, int NDIM>
class EigSolver
{
public:
  //*************************************************************************
  // Typedef's
  typedef Function<T,NDIM> funcT;
//  typedef KPoint<NDIM> kvecT;
  typedef Vector<double,NDIM> kvecT;
  typedef SeparatedConvolution<double,NDIM> operatorT;
  typedef SharedPtr<operatorT> poperatorT;
  //*************************************************************************

  //*************************************************************************
  /// Constructor for periodic system
  EigSolver(World& world, funcT rhon, std::vector<funcT> phis, std::vector<double> eigs,
      std::vector<EigSolverOp<T,NDIM>*> ops, std::vector<kvecT> kpoints,
      double ncharge);
  //*************************************************************************

  //*************************************************************************
  /// Constructor for non-periodic system
  EigSolver(World& world, funcT rhon, std::vector<funcT> phis, std::vector<double> eigs,
      std::vector<EigSolverOp<T,NDIM>*> ops, double ncharge,
      bool periodic = false);
  //*************************************************************************

  //*************************************************************************
  /// Constructor for non-periodic system
  EigSolver(World& world, std::vector<funcT> phis, std::vector<double> eigs,
      std::vector<EigSolverOp<T,NDIM>*> ops, double ncharge);
  //*************************************************************************

  //*************************************************************************
  /// Destructor
  virtual ~EigSolver();
  //*************************************************************************

  //*************************************************************************
  /// This solver has not been optimized for usage in parallel. This solver
  /// processes each eigenfunction in a serial fashion.
  void solve(int maxits);
  //*************************************************************************

  //*************************************************************************
  /// This solver has been optimized for usage in parallel. This solver
  /// processes each eigenfunction in a parallel fashion.
  void multi_solve(int maxits);
  //*************************************************************************

  //*************************************************************************
  double get_eig(int indx)
  {
    return _eigs[indx];
  }
  //*************************************************************************

  //*************************************************************************
  funcT get_phi(int indx)
  {
    return _phis[indx];
  }
  //*************************************************************************

  //*************************************************************************
  const std::vector<funcT>& phis()
  {
    return _phis;
  }
  //*************************************************************************

  //*************************************************************************
  const std::vector<double>& eigs()
  {
    return _eigs;
  }
  //*************************************************************************

  //*************************************************************************
  void addObserver(IEigSolverObserver<T,NDIM>* obs)
  {
    _obs.push_back(obs);
  }
  //*************************************************************************

  //*************************************************************************
  /// Computes a matrix element given the left and right functions.
  T matrix_element(const funcT& phii, const funcT& phij);
  //*************************************************************************

  //*************************************************************************
  /// Prints a matrix element given the left and right functions.
  void print_matrix_elements(const funcT& phii, const funcT& phij);
  //*************************************************************************

  //*************************************************************************
  /// Preprocesses the operators for doing an iteration of "eigensolving".
  void prepare_ops();
  //*************************************************************************

  //*************************************************************************
  /// Makes the BSH Green's functions for the parallel solver (multi_solve()).
  void make_bsh_operators();
  //*************************************************************************

  //*************************************************************************
  void update_occupation();
  //*************************************************************************

  //*************************************************************************
  /// Computes the electronic density
  static funcT compute_rho(std::vector<funcT> phis, std::vector<double> occs,
      const World& world);
  //*************************************************************************

private:
  //*************************************************************************
  /// List of the functions
  std::vector<funcT> _phis;
  //*************************************************************************

  //*************************************************************************
  /// List of the eigenvalues
  std::vector<double> _eigs;
  //*************************************************************************

  //*************************************************************************
  /// List of the ops
  std::vector< EigSolverOp<T,NDIM>* > _ops;
  //*************************************************************************

  //*************************************************************************
  /// List of the ops
  std::vector<kvecT> _kpoints;
  //*************************************************************************

  // Nuclear charge density
  //*************************************************************************
  funcT _rhon;
  //*************************************************************************

  //*************************************************************************
  World& _world;
  //*************************************************************************

  //*************************************************************************
  double _thresh;
  //*************************************************************************

  //*************************************************************************
  // List of the obs
  std::vector<IEigSolverObserver<T,NDIM>*> _obs;
  //*************************************************************************

  // Electronic charge density
  //*************************************************************************
  Function<double,NDIM> _rho;
  //*************************************************************************

  //*************************************************************************
  bool _periodic;
  //*************************************************************************

  //*************************************************************************
  double _ncharge;
  //*************************************************************************

  //*************************************************************************
  // List of the ops
  std::vector<poperatorT> _bops;
  //*************************************************************************

  //*************************************************************************
  // List of the occupation numbers
  std::vector<double> _occs;
  //*************************************************************************

};
//***************************************************************************

}

#endif /*EIGSOLVER_H_*/

