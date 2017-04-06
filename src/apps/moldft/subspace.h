#ifndef SUBSPACE_H
#define SUBSPACE_H

  /*!
   \ingroup periodic_solver

   \brief The Subspace class is a container class holding previous orbitals
   and residuals.
   \par
   The Solver class uses the Krylov Accelerated Inexact Newton Solver (KAIN)
   accelerate the convergence a given calculation. The KAIN solver needs to
   store a subspace of previous orbitals and residuals.
  */

  //***************************************************************************
  class Subspace
  {
    // Typedef's
    typedef std::pair<vector_complex_function_3d,vector_complex_function_3d> pairvecT;
    typedef std::vector<pairvecT> subspaceT;

    //*************************************************************************
    tensor_complex _Q;
    //*************************************************************************

    //*************************************************************************
    subspaceT _subspace;
    //*************************************************************************

    //*************************************************************************
    bool _spinpol;
    //*************************************************************************

    //*************************************************************************
    int _maxsub;
    //*************************************************************************

  public:

    //*************************************************************************
    Subspace(bool spinpol=false, int maxsub=4)
      : _spinpol(spinpol), _maxsub(maxsub)
    {
    }
    //*************************************************************************

    //*************************************************************************
    void update_subspace(World& world,
                         vector_complex_function_3d& awfs_new,
                         vector_complex_function_3d& bwfs_new,
                         const vector_complex_function_3d& awfs_old,
                         const vector_complex_function_3d& bwfs_old,
                         const vector_complex_function_3d& rm)
    {
      // concatenate up and down spins
      vector_complex_function_3d vm = awfs_old;
      if (_spinpol)
      {
        vm.insert(vm.end(), bwfs_old.begin(), bwfs_old.end());
      }

      // Update subspace and matrix Q
      compress(world, vm, false);
      compress(world, rm, false);
      world.gop.fence();
      _subspace.push_back(pairvecT(vm,rm));

      int m = _subspace.size();
      tensor_complex ms(m);
      tensor_complex sm(m);
      for (int s=0; s<m; s++)
      {
          const vector_complex_function_3d& vs = _subspace[s].first;
          const vector_complex_function_3d& rs = _subspace[s].second;
          for (unsigned int i=0; i<vm.size(); i++)
          {
              ms[s] += vm[i].inner_local(rs[i]);
              sm[s] += vs[i].inner_local(rm[i]);
          }
      }
      world.gop.sum(ms.ptr(),m);
      world.gop.sum(sm.ptr(),m);

      tensor_complex newQ(m,m);
      if (m > 1) newQ(Slice(0,-2),Slice(0,-2)) = _Q;
      newQ(m-1,_) = ms;
      newQ(_,m-1) = sm;

      _Q = newQ;
      if (world.rank() == 0) print(_Q);

      // Solve the subspace equations
      tensor_complex c;
      if (world.rank() == 0) {
          double rcond = 1e-12;
          while (1) {
              c = KAIN(_Q,rcond);
              if (abs(c[m-1]) < 3.0) {
                  break;
              }
              else if (rcond < 0.01) {
                  if (world.rank() == 0)
                    print("Increasing subspace singular value threshold ", c[m-1], rcond);
                  rcond *= 100;
              }
              else {
                  if (world.rank() == 0)
                    print("Forcing full step due to subspace malfunction");
                  c = 0.0;
                  c[m-1] = 1.0;
                  break;
              }
          }
      }

      world.gop.broadcast_serializable(c, 0);
      if (world.rank() == 0) {
          print("Subspace solution", c);
      }

      // Form linear combination for new solution
      vector_complex_function_3d phisa_new = zero_functions_compressed<double_complex,3>(world, awfs_old.size());
      vector_complex_function_3d phisb_new = zero_functions_compressed<double_complex,3>(world, bwfs_old.size());
      std::complex<double> one = std::complex<double>(1.0,0.0);
      for (unsigned int m=0; m<_subspace.size(); m++) {
          const vector_complex_function_3d& vm = _subspace[m].first;
          const vector_complex_function_3d& rm = _subspace[m].second;
          const vector_complex_function_3d  vma(vm.begin(),vm.begin()+awfs_old.size());
          const vector_complex_function_3d  rma(rm.begin(),rm.begin()+awfs_old.size());
          const vector_complex_function_3d  vmb(vm.end()-bwfs_old.size(), vm.end());
          const vector_complex_function_3d  rmb(rm.end()-bwfs_old.size(), rm.end());

          gaxpy(world, one, phisa_new, c(m), vma, false);
          gaxpy(world, one, phisa_new,-c(m), rma, false);
          gaxpy(world, one, phisb_new, c(m), vmb, false);
          gaxpy(world, one, phisb_new,-c(m), rmb, false);
      }
      world.gop.fence();

      if (_maxsub <= 1) {
          // Clear subspace if it is not being used
          _subspace.clear();
      }
      else if (_subspace.size() == _maxsub) {
          // Truncate subspace in preparation for next iteration
          _subspace.erase(_subspace.begin());
          _Q = _Q(Slice(1,-1),Slice(1,-1));
      }
      awfs_new = phisa_new;
      bwfs_new = phisb_new;
    }
    //*************************************************************************

    //*************************************************************************
    void update_subspace(World& world,
                         vector_complex_function_3d& awfs_new,
                         const vector_complex_function_3d& awfs_old,
                         const vector_complex_function_3d& rm)
    {
      // concatenate up and down spins
      vector_complex_function_3d vm = awfs_old;

      // Update subspace and matrix Q
      compress(world, vm, false);
      compress(world, rm, false);
      world.gop.fence();
      _subspace.push_back(pairvecT(vm,rm));

      int m = _subspace.size();
      tensor_complex ms(m);
      tensor_complex sm(m);
      for (int s=0; s<m; s++)
      {
          const vector_complex_function_3d& vs = _subspace[s].first;
          const vector_complex_function_3d& rs = _subspace[s].second;
          for (unsigned int i=0; i<vm.size(); i++)
          {
              ms[s] += vm[i].inner_local(rs[i]);
              sm[s] += vs[i].inner_local(rm[i]);
          }
      }
      world.gop.sum(ms.ptr(),m);
      world.gop.sum(sm.ptr(),m);

      tensor_complex newQ(m,m);
      if (m > 1) newQ(Slice(0,-2),Slice(0,-2)) = _Q;
      newQ(m-1,_) = ms;
      newQ(_,m-1) = sm;

      _Q = newQ;
      if (world.rank() == 0) print(_Q);

      // Solve the subspace equations
      tensor_complex c;
      if (world.rank() == 0) {
          double rcond = 1e-12;
          while (1) {
              c = KAIN(_Q,rcond);
              if (abs(c[m-1]) < 3.0) {
                  break;
              }
              else if (rcond < 0.01) {
                  if (world.rank() == 0)
                    print("Increasing subspace singular value threshold ", c[m-1], rcond);
                  rcond *= 100;
              }
              else {
                  if (world.rank() == 0)
                    print("Forcing full step due to subspace malfunction");
                  c = 0.0;
                  c[m-1] = 1.0;
                  break;
              }
          }
      }

      world.gop.broadcast_serializable(c, 0);
      if (world.rank() == 0) {
          print("Subspace solution", c);
      }

      // Form linear combination for new solution
      vector_complex_function_3d phisa_new = zero_functions_compressed<double_complex,3>(world, awfs_old.size());
      std::complex<double> one = std::complex<double>(1.0,0.0);
      for (unsigned int m=0; m<_subspace.size(); m++) {
          const vector_complex_function_3d& vm = _subspace[m].first;
          const vector_complex_function_3d& rm = _subspace[m].second;
          const vector_complex_function_3d  vma(vm.begin(),vm.begin()+awfs_old.size());
          const vector_complex_function_3d  rma(rm.begin(),rm.begin()+awfs_old.size());

          gaxpy(world, one, phisa_new, c(m), vma, false);
          gaxpy(world, one, phisa_new,-c(m), rma, false);
      }
      world.gop.fence();

      if (_maxsub <= 1) {
          // Clear subspace if it is not being used
          _subspace.clear();
      }
      else if (_subspace.size() == _maxsub) {
          // Truncate subspace in preparation for next iteration
          _subspace.erase(_subspace.begin());
          _Q = _Q(Slice(1,-1),Slice(1,-1));
      }
      awfs_new = phisa_new;
    }
    //*************************************************************************

    //*************************************************************************
    void reproject()
    {
      //  //if (world.rank() == 0)
      //    //printf("\n\nreprojecting subspace to wavelet order: %d and thresh: %.5e\n\n",
      //    //FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh());
      //
      //  unsigned int m = _subspace.size();
      //  for (unsigned int s = 0; s < m; s++)
      //  {
      //      vector_complex_function_3d& vs = _subspace[s].first;
      //      vector_complex_function_3d& rs = _subspace[s].second;
      //      reconstruct(world, vs);
      //      reconstruct(world, rs);
      //      unsigned int vm = vs.size();
      //      for (unsigned int i = 0; i < vm; i++)
      //      {
      //        vs[i] = madness::project(vs[i], FunctionDefaults<3>::get_k(),
      //          FunctionDefaults<3>::get_thresh(), false);
      //        rs[i] = madness::project(rs[i], FunctionDefaults<3>::get_k(),
      //          FunctionDefaults<3>::get_thresh(), false);
      //      }
      //      world.gop.fence();
      //      truncate(world, vs);
      //      truncate(world, rs);
      //      normalize(world, vs);
      //  }
      //  world.gop.fence();

    }
    //*************************************************************************

  };

//  template <typename T, int NDIM>
//  struct lbcost {
//      double leaf_value;
//      double parent_value;
//      lbcost(double leaf_value=1.0, double parent_value=0.0) : leaf_value(leaf_value), parent_value(parent_value) {}
//      double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
//          if (key.level() <= 1) {
//              return 100.0*(leaf_value+parent_value);
//          }
//          else if (node.is_leaf()) {
//              return leaf_value;
//          }
//          else {
//              return parent_value;
//          }
//      }
//  };

  //***************************************************************************

  /*! \ingroup periodic_solver
      \brief The main class of the periodic DFT solver
      \f[
      z = frac{x}{1 - y^2}
      \f]
  */

#endif
