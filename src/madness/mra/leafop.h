/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/*
* Leafop.h
 *
 *  Created on: Apr 12, 2016
 *      Author: kottmanj
 *
 *     Definition of Operators which operate on leaf boxes
 *     includes special operators for cuspy funtions (f12 applications) which enforce deeper refinement on the "diagonal" boxes (coordinate r1 == coordinate r2)
 *
 */

#ifndef SRC_MADNESS_MRA_LEAFOP_H_
#define SRC_MADNESS_MRA_LEAFOP_H_

namespace madness {

// forward declarations
template<std::size_t NDIM>
class Key;

template<typename T>
class GenTensor;

template<typename T, std::size_t NDIM>
class SeparatedConvolution;

/// helper structure for the Leaf_op
/// The class has an operator which decides if a given key belongs to a special box (needs further refinement up to a special level)
/// this is the default class which always gives back false
template<typename T, std::size_t NDIM>
struct Specialbox_op {
public:
    Specialbox_op() {}

    virtual ~Specialbox_op() {}

    virtual std::string name() const { return "default special box which only checks for the special points"; }

    /// the operator which deceides if the given box is special
    /// @param[in] the key of the current box
    /// @param[in] the function which will be constructed (if it is needed)
    /// @param[out] bool which states if the box is special or not
    virtual bool operator()(const Key<NDIM>& key, const FunctionImpl<T, NDIM> *const f = NULL) const { return false; }

    // check if on of the special points is in the current box (or at low refinement level in a neighbouring box)
    /// @param[in] the key of the box
    /// @param[in] the function which will be constructed (if it is needed)
    bool
    check_special_points(const Key<NDIM>& key, const FunctionImpl<T, NDIM> *const f) const {
        const std::vector<Vector<double, NDIM> >& special_points = f->get_special_points();
        if (special_points.empty()) return false;
        // level 0 and 1 has only boundary boxes
        if (key.level() > 1 and box_is_at_boundary(key)) return false;
        // check for all special points if they are neighbours of the current box
        BoundaryConditions<NDIM> bc = FunctionDefaults<NDIM>::get_bc();
        const auto bperiodic = bc.is_periodic();
        for (size_t i = 0; i < special_points.size(); ++i) {
            Vector<double, NDIM> simpt;
            user_to_sim(special_points[i], simpt);
            Key<NDIM> specialkey = simpt2key(simpt, key.level());
            // use adaptive scheme: if we are at a low level refine also neighbours
            int ll = get_half_of_special_level(f->get_special_level());
            if (ll < f->get_initial_level()) ll = f->get_initial_level();
            if (key.level() > ll) {
                if (specialkey == key) return true;
                else return false;
            } else {
                if (specialkey.is_neighbor_of(key, bperiodic)) return true;
                else return false;
            }
        }
        return false;
    }


    template<typename Archive>
    void
    serialize(Archive& ar) {}

    /// Determine if a given box is at the boundary of the simulation box
    /// @param[in] the key of the current box
    /// Since boxes of level 0 and 1 are always at the boundary this case should be avoided
    /// if the boundary conditions are periodic then all boxes are boundary boxes
    virtual bool box_is_at_boundary(const Key<NDIM>& key) const {
        // l runs from l=0 to l=2^n-1 for every dimension
        // if one l is = 2^n or 0 we are at the boundary
        // note that boxes of level 0 and 1 are always boundary boxes
        for (size_t i = 0; i < NDIM; i++) {
            if (key.translation()[i] == 0 or key.translation()[i] == pow(2, key.level()) - 1) {
                if (FunctionDefaults<NDIM>::get_bc()(i, 0) != BC_PERIODIC)return true;
            }
        }
        return false;
    }

    // if you want to use milder refinement conditions for special boxes with increasing refinement level
    // Cuspybox_op needs this
    size_t get_half_of_special_level(const size_t& sl = FunctionDefaults<NDIM>::get_special_level()) const {
        size_t ll = sl;
        if (sl % 2 == 0) ll = sl / 2;
        else ll = (sl + 1) / 2;
        return ll;
    }

};

/// the high_dimensional key is broken appart into two lower dimensional keys, if they are neighbours of each other then refinement
/// up to the special_level is enforced
/// For general class description see the base class Specialbox_op
template<typename T, std::size_t NDIM>
struct ElectronCuspyBox_op : public Specialbox_op<T, NDIM> {
public:
    ElectronCuspyBox_op() {}

    ~ElectronCuspyBox_op() {}

    template<typename Archive>
    void serialize(Archive& ar) {}

    std::string name() const { return "Cuspybox_op"; }

    /// Operator which decides if the key belongs to a special box
    /// The key is broken apart in two lower dimensional keys (for electron cusps this is 6D -> 2x3D)
    /// if the keys are neighbours then refinement up to the special level is enforced (means the 6D box is close to the cusp or contains it)
    /// if the refinement level is already beyond half of the special_level then refinement is only enforded if the broken keys are the same (6D box contains cusp)
    bool operator()(const Key<NDIM>& key, const FunctionImpl<T, NDIM> *const f) const {
        // do not treat boundary boxes beyond level 2 as special (level 1 and 0 consist only of boundary boxes)
        if (not(key.level() < 2)) {
            if (this->box_is_at_boundary(key)) return false;
        }

        // Cuspybox_op only valid for even dimensions, if uneven dims are needed just make a new class with NDIM+1/2 and NDIM-LDIM
        if constexpr (NDIM % 2 == 0) {
          BoundaryConditions<NDIM> bc = FunctionDefaults<NDIM>::get_bc();
          const auto bperiodic = bc.template is_periodic<NDIM / 2>();
          Key<NDIM / 2> key1;
          Key<NDIM / 2> key2;
          key.break_apart(key1, key2);
          int ll = this->get_half_of_special_level();
          if (ll < f->get_initial_level())
            ll = f->get_initial_level();
          if (key.level() > ll) {
            if (key1 == key2)
              return true;
            else
              return false;
          } else {
            if (key1.is_neighbor_of(key2, bperiodic))
              return true;
            else
              return false;
          }
        }
        MADNESS_EXCEPTION("We should not end up here (check further of cuspy box)", 1);
        return false;
    }

};

/// This works similar to the Cuspybox_op: The key is broken apart (2N-dimensional -> 2x N-Dimensional)
/// then it is checked if one of the lower dimensional keys contains a lower dimensional special points (which will be the nuclear-coordinates)
template<typename T, std::size_t NDIM>
struct NuclearCuspyBox_op : public Specialbox_op<T, NDIM> {
public:
    NuclearCuspyBox_op() : particle(-1) {
        // failsafe
        if (NDIM % 2 != 0) MADNESS_EXCEPTION("NuclearCuspyBox works only for even dimensions", 1);
    }

    NuclearCuspyBox_op(const size_t p) : particle(p) {
        // failsafe
        if (NDIM % 2 != 0) MADNESS_EXCEPTION("NuclearCuspyBox works only for even dimensions", 1);
        MADNESS_ASSERT(particle == 1 or particle == 2 or particle == 0);
    }

    ~NuclearCuspyBox_op() {}

    // particle==0: check both particles
    // particle==1 or particle==2: check only particle1 or 2
    int particle;

    template<typename Archive>
    void serialize(Archive& ar) { ar & particle; }

    std::string name() const { return "Cuspybox_op for nuclear cusps"; }

    /// Operator which decides if the key belongs to a special box
    /// The key is broken appart in two lower dimensional keys (6D -> 2x3D)
    /// The special points should be given in a format which can also be broken appart
    /// so: 6D special_points = (3D-special-points, 3D-special-points)
    /// if the refinement level is already beyond half of the special_level then refinement is only enforded if the broken keys are the same (6D box contains cusp)
    bool operator()(const Key<NDIM>& key, const FunctionImpl<T, NDIM> *const f) const {
        if (not(key.level() < 2)) {
            if (this->box_is_at_boundary(key)) return false;
        }
        MADNESS_ASSERT(particle == 1 or particle == 2 or particle == 0);
        if (f == NULL) MADNESS_EXCEPTION("NuclearCuspyBox: Pointer to function is NULL", 1);
        const std::vector<Vector<double, NDIM> >& special_points = f->get_special_points();
        if (special_points.empty()) MADNESS_EXCEPTION(
                "Demanded NuclearCuspyBox but the special points of the function are empty", 1);

        // break the special points into 3D points
        std::vector<Vector<double, NDIM / 2> > lowdim_sp;
        for (size_t i = 0; i < special_points.size(); i++) {
            Vector<double, NDIM / 2> lowdim_tmp;
            for (size_t j = 0; j < NDIM / 2; j++) {
                lowdim_tmp[j] = special_points[i][j];
                // check if the formatting is correct
                if (special_points[i][j] != special_points[i][NDIM / 2 + j]) MADNESS_EXCEPTION(
                        "NuclearCuspyBox: Wrong format of special_point: ", 1);
            }
            lowdim_sp.push_back(lowdim_tmp);
        }

        // now break the key appart and check if one if the results is in the neighbourhood of a special point
        BoundaryConditions<NDIM / 2> bc = FunctionDefaults<NDIM / 2>::get_bc();
        const auto bperiodic = bc.is_periodic();
        Key<NDIM / 2> key1;
        Key<NDIM / 2> key2;

        key.break_apart(key1, key2);
        for (size_t i = 0; i < lowdim_sp.size(); ++i) {
            Vector<double, NDIM / 2> simpt;
            user_to_sim(lowdim_sp[i], simpt);
            Key<NDIM / 2> specialkey = simpt2key(simpt, key1.level());
            // use adaptive scheme: if we are at a low level refine also neighbours
            int ll = this->get_half_of_special_level(f->get_special_level());
            if (ll < f->get_initial_level()) ll = f->get_initial_level();
            if (key.level() > ll) {
                if (particle == 1 and specialkey == key1) return true;
                else if (particle == 2 and specialkey == key2) return true;
                else if (particle == 0 and (specialkey == key1 or specialkey == key2)) return true;
                else return false;
            } else {
                if (particle == 1 and specialkey.is_neighbor_of(key1, bperiodic)) return true;
                else if (particle == 2 and specialkey.is_neighbor_of(key2, bperiodic)) return true;
                else if (particle == 0 and
                         (specialkey.is_neighbor_of(key1, bperiodic) or specialkey.is_neighbor_of(key2, bperiodic)))
                    return true;
                else return false;
            }
        }
        return false;
    }


};

template<typename T, std::size_t NDIM, typename opT, typename specialboxT>
class Leaf_op {
public:
    /// the function which the operators use (in most cases this function is also the function that will be constructed)
    const FunctionImpl<T, NDIM> *f;
    /// the operator which is used for screening (null pointer means no screening)
    const opT *op;
    /// special box structure: has an operator which decides if a given key belongs to a special box
    /// the base class Specialbox_op<T,NDIM> just gives back false for every key
    /// the derived class Cuspybox_op<T,NDIM> is used for potentials containing cusps at coalescence points of electrons
    specialboxT specialbox;

    /// pre or post screening (see if this is the general_leaf_op or the leaf_op_other)
    virtual bool
    do_pre_screening() const {
        return false;
    }

    Leaf_op()
            : f(NULL), op(NULL), specialbox(specialboxT()) {
    }

    Leaf_op(const FunctionImpl<T, NDIM> *const tmp)
            : f(tmp), op(NULL), specialbox(specialboxT()) {
    }

    Leaf_op(const FunctionImpl<T, NDIM> *const tmp, specialboxT& sb)
            : f(tmp), op(NULL), specialbox(sb) {
    }

    Leaf_op(const FunctionImpl<T, NDIM> *const tmp, const opT *const ope, specialboxT& sb)
            : f(tmp), op(ope), specialbox(sb) {
    }

    Leaf_op(const Leaf_op& other) : f(other.f), op(other.op), specialbox(other.specialbox) {}

    virtual
    ~Leaf_op() {
    }

    virtual std::string
    name() const {
        return "general_leaf_op";
    }

    /// make sanity check
    virtual void
    sanity() const {
        if (f == NULL) MADNESS_EXCEPTION(("Error in " + name() + " pointer to function is NULL").c_str(), 1);
    }

public:
    /// pre-determination
    /// the decision if the current box  will be a leaf box is made from information of another function
    /// this is needed for on demand function
    /// not that the on-demand function that is being constructed is not the function in this class
    virtual bool
    pre_screening(const Key<NDIM>& key) const {
        MADNESS_EXCEPTION("pre-screening was called for leaf_op != leaf_op_other", 1);
        return false;
    }

    /// post-determination
    /// determination if the box will be a leaf box after the coefficients are made
    /// the function that is being constructed is the function in this class
    /// the function will use the opartor op in order to screen, if op is a NULL pointer the result is always: false
    /// @param[in] key: the key to the current box
    /// @param[in] coeff: Coefficients of the current box
    virtual bool
    post_screening(const Key<NDIM>& key, const GenTensor<T>& coeff) const {
        if (op == NULL) return false;
        if (key.level() < this->f->get_initial_level()) return false;
        sanity();
        const double cnorm = coeff.normf();
        BoundaryConditions<NDIM> bc = FunctionDefaults<NDIM>::get_bc();
        // const auto bperiodic = bc.is_periodic();

        typedef Key<opT::opdim> opkeyT;
        const opkeyT source = op->get_source_key(key);

        const double thresh = (this->f->truncate_tol(this->f->get_thresh(), key));
        const std::vector<opkeyT>& disp = op->get_disp(key.level());
        const opkeyT& d = *disp.begin();         // use the zero-displacement for screening
        const double opnorm = op->norm(key.level(), d, source);
        const double norm = opnorm * cnorm;

        return norm < thresh;
    }

    /// determines if a node is well represented compared to the parents
    /// @param[in]  key the FunctionNode which we want to determine if it's a leaf node
    /// @param[in]  coeff   the coeffs of key
    /// @param[in]  parent  the coeffs of key's parent node
    /// @return is the FunctionNode of key a leaf node?
    bool
    compare_to_parent(const Key<NDIM>& key, const GenTensor<T>& coeff, const GenTensor<T>& parent) const {
        if (key.level() < this->f->get_initial_level()) return false;
        //this->sanity();
        if (parent.has_no_data()) return false;
        if (key.level() < this->f->get_initial_level()) return false;
        GenTensor<T> upsampled = this->f->upsample(key, parent);
        upsampled.scale(-1.0);
        upsampled += coeff;
        const double dnorm = upsampled.normf();
        const bool is_leaf = (dnorm < this->f->truncate_tol(this->f->get_thresh(), key.level()));
        return is_leaf;
    }

    /// check if the box is a special box
    /// first check: Is one of the special points defined in the function f in the current box or a neighbour box
    /// second check: Is the box special according to the specialbox operator (see class description of Specialbox_op)
    /// @param[in] the key of the box
    bool
    special_refinement_needed(const Key<NDIM>& key) const {
        if (key.level() > f->get_special_level()) return false;
        else if (specialbox.check_special_points(key, f)) return true;
        else if (specialbox(key, f)) return true;
        else return false;
    }

    template<typename Archive>
    void
    serialize(Archive& ar) {
        ar & this->f & this->specialbox;
    }

};

/// leaf_op for construction of an on-demand function
/// the function in the class is the function which defines the structure of the on-demand function
/// means: The tree of the function which is to construct will mirror the tree-structure of the function in this class (function variable defined in base-class)
template<typename T, std::size_t NDIM>
class Leaf_op_other : public Leaf_op<T, NDIM, SeparatedConvolution<double, NDIM>, Specialbox_op<T, NDIM> > {
    std::string
    name() const {
        return "leaf_op_other";
    }
    /// we only need the pre-determination based on the already constructed function f
    /// if f at box(key) is a leaf then return true
public:
    bool
    do_pre_screening() const {
        return true;
    }

    Leaf_op_other() : Leaf_op<T, NDIM, SeparatedConvolution<double, NDIM>, Specialbox_op<T, NDIM> >() {
    }

    Leaf_op_other(const FunctionImpl<T, NDIM> *const f_) {
        this->f = f_;
    }

    Leaf_op_other(const Leaf_op_other& other)
            : Leaf_op<T, NDIM, SeparatedConvolution<double, NDIM>, Specialbox_op<T, NDIM> >(other.f) {}

    bool
    pre_screening(const Key<NDIM>& key) const {
        MADNESS_ASSERT(this->f->get_coeffs().is_local(key));
        return (not this->f->get_coeffs().find(key).get()->second.has_children());
    }

    void sanity() const {
        if (this->f == 0) MADNESS_EXCEPTION("Leaf_op_other: f is NULL pointer", 1);
        if (this->op != 0) MADNESS_EXCEPTION("Leaf_op_other: Screening operator was set", 1);
    }

    /// this should never be called for this leaf_op since the function in this class as already constructed
    bool
    post_screening(const Key<NDIM>& key, const GenTensor<T>& G) const {
        MADNESS_EXCEPTION("post-determination was called for leaf_op_other", 1);
        return false;
    }

    /// this should never be called for this leaf_op since the function in this class as already constructed
    bool
    compare_to_parent(const Key<NDIM>& key, const GenTensor<T>& a, const GenTensor<T>& b) const {
        MADNESS_EXCEPTION("compare-to-parent was called for leaf_op_other", 1);
        return false;
    }

    /// this should never be called for this leaf_op since the function in this class as already constructed
    bool
    special_refinement_needed(const Key<NDIM>& key) const {
        MADNESS_EXCEPTION("special_refinement_needed was called for leaf_op_other", 1);
        return false;
    }

    template<typename Archive>
    void
    serialize(Archive& ar) {
        ar & this->f;
    }

};

} /* namespace madness */

#endif /* SRC_MADNESS_MRA_LEAFOP_H_ */
