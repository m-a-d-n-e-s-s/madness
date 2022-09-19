//
// Created by Florian Bischoff on 5/25/22.
//
#ifndef MADNESS_MOLECULAR_FUNCTORS_H
#define MADNESS_MOLECULAR_FUNCTORS_H

#include<madness/chem/molecule.h>
#include<madness/chem/molecularbasis.h>
#include<madness/mra/mra.h>
#include<madness/mra/functypedefs.h>
namespace madchem {



class MolecularGuessDensityFunctor : public madness::FunctionFunctorInterface<double,3> {
private:
    const madness::Molecule& molecule;
    const madness::AtomicBasisSet& aobasis;
public:
    MolecularGuessDensityFunctor(const madness::Molecule& molecule, const madness::AtomicBasisSet& aobasis)
            : molecule(molecule), aobasis(aobasis) {}

    double operator()(const madness::coord_3d& x) const {
        return aobasis.eval_guess_density(molecule, x[0], x[1], x[2]);
    }

    std::vector<madness::coord_3d> special_points() const {return molecule.get_all_coords_vec();}
};


class AtomicBasisFunctor : public madness::FunctionFunctorInterface<double,3> {
private:
    const madness::AtomicBasisFunction aofunc;

public:
    AtomicBasisFunctor(const madness::AtomicBasisFunction& aofunc)
            : aofunc(aofunc)
    {}

    double operator()(const madness::coord_3d& x) const {
        return aofunc(x[0], x[1], x[2]);
    }

    std::vector<madness::coord_3d> special_points() const {
        return std::vector<madness::coord_3d>(1,aofunc.get_coords_vec());
    }
};


class AtomicAttractionFunctor : public madness::FunctionFunctorInterface<double,3> {
private:
    const madness::Molecule& molecule;
    const int iatom;

public:
    AtomicAttractionFunctor(const madness::Molecule& molecule, int iatom)
            : molecule(molecule), iatom(iatom) {}

    double operator()(const madness::coord_3d& x) const {
        const madness::Atom& atom=molecule.get_atom(iatom);
        const madness::coord_3d coord={atom.x,atom.y,atom.z};
        double r = (x-coord).normf();
        return -atom.q * madness::smoothed_potential(r*molecule.get_rcut()[iatom])
               *molecule.get_rcut()[iatom];
    }

    std::vector<madness::coord_3d> special_points() const {
        return std::vector<madness::coord_3d>(1,molecule.get_atom(iatom).get_coords());
    }
};

class MolecularDerivativeFunctor : public madness::FunctionFunctorInterface<double,3> {
private:
    const madness::Molecule& molecule;
    const int atom;
    const int axis;

public:
    MolecularDerivativeFunctor(const madness::Molecule& molecule, int atom, int axis)
            : molecule(molecule), atom(atom), axis(axis)
    {}

    double operator()(const madness::coord_3d& x) const {
        return molecule.nuclear_attraction_potential_derivative(atom, axis, x[0], x[1], x[2]);
    }

    std::vector<madness::coord_3d> special_points() const {
        return std::vector<madness::coord_3d>(1,molecule.get_atom(atom).get_coords());
    }
};

class MolecularSecondDerivativeFunctor : public madness::FunctionFunctorInterface<double,3> {
private:
    const madness::Molecule& molecule;
    const int atom;
    const int iaxis, jaxis;

public:
    MolecularSecondDerivativeFunctor(const madness::Molecule& molecule, int atom,
                                     int iaxis, int jaxis)
            : molecule(molecule), atom(atom),iaxis(iaxis), jaxis(jaxis)
    {}

    double operator()(const madness::coord_3d& x) const {
        return molecule.nuclear_attraction_potential_second_derivative(atom,
                                                                       iaxis, jaxis, x[0], x[1], x[2]);
    }

    std::vector<madness::coord_3d> special_points() const {
        return std::vector<madness::coord_3d>(1,molecule.get_atom(atom).get_coords());
    }
};


class CorePotentialDerivativeFunctor : public madness::FunctionFunctorInterface<double,3> {
private:
    const madness::Molecule& molecule;
    const int atom;
    const int axis;
    std::vector<madness::coord_3d> specialpt;
public:
    CorePotentialDerivativeFunctor(const madness::Molecule& molecule, int atom, int axis)
            : molecule(molecule), atom(atom), axis(axis) {}

    double operator()(const madness::coord_3d& r) const {
        return molecule.core_potential_derivative(atom, axis, r[0], r[1], r[2]);
    }
};



}
#endif //MADNESS_MOLECULAR_FUNCTORS_H
