/*
 * electronicstructureapp.h
 *
 *  Created on: Nov 5, 2008
 *      Author: eh7
 */

#ifndef ELECTRONICSTRUCTUREAPP_H_
#define ELECTRONICSTRUCTUREAPP_H_

#include "molecularbasis.h"

typedef SharedPtr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef vector<functionT> vecfuncT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef SharedPtr<operatorT> poperatorT;

class MolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
public:
    MolecularPotentialFunctor(const Molecule& molecule)
        : _mentity(molecule)
    {}

    double operator()(const coordT& x) const {
        return _mentity.nuclear_attraction_potential(x[0], x[1], x[2]);
    }
};

class MolecularNuclearChargeDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
public:
    MolecularNuclearChargeDensityFunctor(const Molecule& molecule)
        : _mentity(molecule)
    {}

    double operator()(const coordT& x) const {
        return _mentity.nuclear_charge_density(x[0], x[1], x[2]);
    }
};

class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
    const AtomicBasisSet& _aobasis;
public:
    MolecularGuessDensityFunctor(const MolecularEntity& mentity, const AtomicBasisSet& aobasis)
        : _mentity(mentity), _aobasis(aobasis)
    {}

    double operator()(const coordT& x) const {
        return _aobasis.eval_guess_density(mentity, x[0], x[1], x[2]);
    }
};


class AtomicBasisFunctor : public FunctionFunctorInterface<double,3> {
private:
    const AtomicBasisFunction aofunc;
public:
    AtomicBasisFunctor(const AtomicBasisFunction& aofunc) : aofunc(aofunc)
    {}

    double operator()(const coordT& x) const {
        return aofunc(x[0], x[1], x[2]);
    }
};

/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const coordT& x) const {
        return x[axis];
    }
};

double rsquared(const coordT& r) {
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

/// A MADNESS functor to compute the cartesian moment x^i * y^j * z^k (i, j, k integer and >= 0)
class MomentFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int i, j, k;
public:
    MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
    double operator()(const coordT& r) const {
        double xi=1.0, yj=1.0, zk=1.0;
        for (int p=0; p<i; p++) xi *= r[0];
        for (int p=0; p<j; p++) yj *= r[1];
        for (int p=0; p<k; p++) zk *= r[2];
        return xi*yj*zk;
    }
};

tensorT sqrt(const tensorT& s, double tol=1e-8) {
    int n=s.dim[0], m=s.dim[1];
    MADNESS_ASSERT(n==m);
    tensorT c, e;
    //s.gaxpy(0.5,transpose(s),0.5); // Ensure exact symmetry
    syev(s, &c, &e);
    for (int i=0; i<n; i++) {
        if (e(i) < -tol) {
            MADNESS_EXCEPTION("Matrix square root: negative eigenvalue",i);
        }
        else if (e(i) < tol) { // Ugh ..
            print("Matrix square root: Warning: small eigenvalue ", i, e(i));
            e(i) = tol;
        }
        e(i) = 1.0/sqrt(e(i));
    }
    for (int j=0; j<n; j++) {
        for (int i=0; i<n; i++) {
            c(j,i) *= e(i);
        }
    }
    return c;
}

tensorT energy_weighted_orthog(const tensorT& s, const tensorT eps) {
    int n=s.dim[0], m=s.dim[1];
    MADNESS_ASSERT(n==m);
    tensorT d(n,n);
    for (int i=0; i<n; i++) d(i,i) = eps(i);
    tensorT c, e;
    sygv(d, s, 1, &c, &e);
    return c;
}


template <typename T, int NDIM>
Cost lbcost(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) {
  return 1;
}

class ElectronicStructureApp
{
public:
  ElectronicStructureApp(const std::string& filename)
  {
    init(filename);
  }

  ~ElectronicStructureApp()
  {
    if (_mentity) delete _mentity;
  }

  void init(const std::string& filename)
  {
    read_params(filename);
    read_positions(filename);
    aobasis.read_file("sto-3g");
    mentity->center();
  }

  void read_params(const std::string& filename)
  {

  }

  void read_positions(const std::string& filename)
  {
    _mentity = new MolecularEntity(filename);
  }

  void make_nuclear_potential(World& world)
  {
      //vnuc = factoryT(world).functor(functorT(new MolecularPotentialFunctor(molecule))).thresh(vtol).truncate_on_project();
      //vnuc.truncate();
      //vnuc.reconstruct();
      Function<double, 3> rhon = factoryT(world).functor(
          functorT(new MolecularNuclearChargeDensityFunctor(_mentity))).
          thresh(vtol).truncate_on_project();
  }


private:
  MolecularEntity* _mentity = 0;
  AtomicBasisSet aobasis;
};

#endif /* ELECTRONICSTRUCTUREAPP_H_ */
