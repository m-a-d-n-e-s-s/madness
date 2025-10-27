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

#ifndef MADNESS_MRA_FUNCDEFAULTS_H__INCLUDED
#define MADNESS_MRA_FUNCDEFAULTS_H__INCLUDED


/// \file funcdefaults.h
/// \brief Provides FunctionDefaults and utilities for coordinate transformation
/// \ingroup mrabcext

#include <madness/world/MADworld.h>
#include <madness/world/vector.h>
#include <madness/world/worlddc.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/gentensor.h>
#include <madness/mra/bc.h>
#include <madness/mra/key.h>

#include <optional>

namespace madness {
    template <typename T, std::size_t NDIM> class FunctionImpl;

    /// The maximum wavelet order presently supported
    static const int MAXK = 30;

    /// The maximum depth of refinement possible
    static const int MAXLEVEL = 8*sizeof(Translation)-2;

    enum TreeState {
    	reconstructed,				///< s coeffs at the leaves only
		compressed, 				///< d coeffs in internal nodes, s and d coeffs at the root
		nonstandard, 				///< s and d coeffs in internal nodes
    	nonstandard_with_leaves, 	///< like nonstandard, with s coeffs at the leaves
        nonstandard_after_apply, 	///< s and d coeffs, state after operator application
		redundant,					///< s coeffs everywhere
        redundant_after_merge,		///< s coeffs everywhere, must be summed up to yield the result
		on_demand,					///< no coeffs anywhere, but a functor providing if necessary
		unknown
    };

    template<std::size_t NDIM=1>
    std::ostream& operator<<(std::ostream& os, const TreeState treestate) {
        if (treestate==reconstructed) os << "reconstructed";
        else if (treestate==compressed) os << "compressed";
        else if (treestate==nonstandard) os << "nonstandard";
        else if (treestate==nonstandard_with_leaves) os << "nonstandard_with_leaves";
        else if (treestate==nonstandard_after_apply) os << "nonstandard_after_apply";
        else if (treestate==redundant) os << "redundant";
        else if (treestate==redundant_after_merge) os << "redundant_after_merge";
        else if (treestate==on_demand) os << "on_demand";
        else if (treestate==unknown) os << "unknown";
        else {
            MADNESS_EXCEPTION("unknown treestate",1);
        }
        return os;
    }

    /// FunctionDefaults holds default paramaters as static class members

    /// Declared and initialized in mra.cc and/or funcimpl::initialize.
    ///
    /// Currently all functions of the same dimension share the same cell dimensions
    /// since they are stored inside FunctionDefaults ... if you change the
    /// cell dimensions *all* functions of that dimension are affected.
    ///
    /// N.B.  Ultimately, we may need to make these defaults specific to each
    /// world, as should be all global state.
    /// \ingroup mra
    template <std::size_t NDIM>
    class FunctionDefaults {
        MADNESS_PRAGMA_CLANG(diagnostic push)
    	MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

    private:
        static int k;                  ///< Wavelet order
        static double thresh;          ///< Truncation threshold
        static int initial_level;      ///< Initial level for fine scale projection
        static int special_level;      ///< Minimum level for fine scale projection of special boxes
        static int max_refine_level;   ///< Level at which to stop refinement
        static int truncate_mode;    ///< Truncation method
        static bool refine;            ///< Whether to refine new functions
        static bool autorefine;        ///< Whether to autorefine in multiplication, etc.
        static bool debug;             ///< Controls output of debug info
        static bool truncate_on_project; ///< If true initial projection inserts at n-1 not n
        static bool apply_randomize;   ///< If true use randomization for load balancing in apply integral operator
        static bool project_randomize; ///< If true use randomization for load balancing in project/refine
        static std::optional<BoundaryConditions<NDIM>> bc; ///< Default boundary conditions, not initialized by default and must be set explicitly before use
        static Tensor<double> cell ;   ///< cell[NDIM][2] Simulation cell, cell(0,0)=xlo, cell(0,1)=xhi, ...
        static Tensor<double> cell_width;///< Width of simulation cell in each dimension
        static Tensor<double> rcell_width; ///< Reciprocal of width
        static double cell_volume;      ///< Volume of simulation cell
        static double cell_min_width;   ///< Size of smallest dimension
        static TensorType tt;			///< structure of the tensor in FunctionNode
        static std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > pmap; ///< Default mapping of keys to processes
        static int pmap_nproc; ///< Number of processes assumed by pmap, -1 indicates uninitialized pmap

        static Tensor<double> make_default_cell() {
          if (NDIM) {
            Tensor<double> cell(NDIM, 2);
            cell(_, 1) = 1.0;
            return cell;
          } else
            return {};
        }

        static Tensor<double> make_default_cell_width() {
          if (NDIM) {
            Tensor<double> cell(NDIM);
            cell(_) = 1.0;
            return cell;
          } else
            return {};
        }

        static void recompute_cell_info() {
            MADNESS_ASSERT(cell.dim(0)==NDIM && cell.dim(1)==2 && cell.ndim()==2);
            cell_width = cell(_,1)-cell(_,0);
            cell_volume = cell_width.product();
            cell_min_width = cell_width.min();
            rcell_width = copy(cell_width);
            for (std::size_t i=0; i<NDIM; ++i) rcell_width(i) = 1.0/rcell_width(i);
        }

    public:


		/// Used to set defaults to k=7, thresh=1-5, for a unit cube [0,1].
		/// @warning does not reset the boundary conditions if they are already set
		static void set_defaults(World& world);

        static void print();

        /// Returns the default wavelet order
        static int get_k() {
        	return k;
        }

        /// Sets the default wavelet order

        /// Existing functions are unaffacted.
        static void set_k(int value) {
        	k=value;
        	MADNESS_ASSERT(k>0 && k<=MAXK);
        }

        /// Returns the default threshold
        static const double& get_thresh() {
        	return thresh;
        }

        /// Sets the default threshold

        /// Existing functions are unaffected
        static void set_thresh(double value) {
        	thresh=value;
        }

        /// Returns the default initial projection level
        static int get_initial_level() {
        	return initial_level;
        }

        /// Returns the default projection level for special boxes
        static int get_special_level() {
        	return special_level;
        }

        /// Sets the default initial projection level

        /// Existing functions are unaffected
        static void set_initial_level(int value) {
        	initial_level=value;
        	MADNESS_ASSERT(value>0 && value<MAXLEVEL);
        }

        /// Existing functions are unaffected
        static void set_special_level(int value) {
        	special_level=value;
        	MADNESS_ASSERT(value>=0 && value<MAXLEVEL);
        	MADNESS_ASSERT(max_refine_level>=special_level);
        }

        /// Gets the default maximum adaptive refinement level
        static int get_max_refine_level() {
        	return max_refine_level;
        }

        /// Sets the default maximum adaptive refinement level

        /// Existing functions are unaffected
        static void set_max_refine_level(int value) {
        	max_refine_level=value;
        	MADNESS_ASSERT(value>0 && value<MAXLEVEL);
        	MADNESS_ASSERT(max_refine_level>=initial_level);
        	MADNESS_ASSERT(max_refine_level>=special_level);
        }

        /// Gets the default truncation mode
        static int get_truncate_mode() {
        	return truncate_mode;
        }

        /// Sets the default truncation mode

        /// Existing functions are unaffected
        static void set_truncate_mode(int value) {
        	truncate_mode=value;
        	MADNESS_ASSERT(value>=0 && value<4);
        }

        /// Gets the default adaptive refinement flag
        static bool get_refine() {
        	return refine;
        }

        /// Sets the default adaptive refinement flag

        /// Existing functions are unaffected
        static void set_refine(bool value) {
        	refine=value;
        }

        /// Gets the default adaptive autorefinement flag
        static bool get_autorefine() {
        	return autorefine;
        }

        /// Sets the default adaptive autorefinement flag

        /// Existing functions are unaffected
        static void set_autorefine(bool value) {
        	autorefine=value;
        }

        /// Gets the default debug flag (is this used anymore?)
        static bool get_debug() {
        	return debug;
        }

        /// Sets the default debug flag (is this used anymore?)

        /// Not sure if this does anything useful
        static void set_debug(bool value) {
        	debug=value;
        }

        /// Gets the default truncate on project flag
        static bool get_truncate_on_project() {
        	return truncate_on_project;
        }

        /// Sets the default truncate on project flag

        /// Existing functions are unaffected
        static void set_truncate_on_project(bool value) {
        	truncate_on_project=value;
        }

        /// Gets the random load balancing for integral operators flag
        static bool get_apply_randomize() {
        	return apply_randomize;
        }

        /// Sets the random load balancing for integral operators flag
        static void set_apply_randomize(bool value) {
        	apply_randomize=value;
        }


        /// Gets the random load balancing for projection flag
        static bool get_project_randomize() {
        	return project_randomize;
        }

        /// Sets the random load balancing for projection flag
        static void set_project_randomize(bool value) {
        	project_randomize=value;
        }

        /// Returns the default boundary conditions
        static const BoundaryConditions<NDIM>& get_bc() {
          if (!bc.has_value()) {
            const std::string msg = "FunctionDefaults<" + std::to_string(NDIM) + ">::get_bc: must initialize boundary conditions by set_bc or set_defaults or startup";
            MADNESS_EXCEPTION(msg.c_str(), 1);
          }
          return bc.value();
        }

        /// Sets the default boundary conditions
        static void set_bc(const BoundaryConditions<NDIM>& value) {
        	bc=value;
        }

        /// Returns the default tensor type
        static TensorType get_tensor_type() {
        	return tt;
        }

        /// Sets the default tensor type
        static void set_tensor_type(const TensorType& t) {
#if HAVE_GENTENSOR
        	tt=t;
#else
        	tt=TT_FULL;
#endif
        }

        /// adapt the special level to resolve the smallest length scale
        static int set_length_scale(const double lo,const size_t k=get_k()) {
        	const double dk = (double) k;
        	double Lmax=FunctionDefaults<NDIM>::get_cell_width().max();
        	double lo_sim=lo/Lmax;  // lo in simulation coordinates;
        	const int special_level=Level(-log2(lo_sim*dk));
        	return special_level;
        }

        /// Gets the user cell for the simulation
        static const Tensor<double>& get_cell() {
        	return cell;
        }

        /// Gets the user cell for the simulation

        /// Existing functions are probably rendered useless
        static void set_cell(const Tensor<double>& value) {
        	cell=copy(value);
        	recompute_cell_info();
        }

        /// Sets the user cell to be cubic with each dimension having range \c [lo,hi]

        /// Existing functions are probably rendered useless
        static void set_cubic_cell(double lo, double hi) {
        	cell(_,0)=lo;
        	cell(_,1)=hi;
        	recompute_cell_info();
        }

        /// Returns the width of each user cell dimension
        static const Tensor<double>& get_cell_width() {
        	return cell_width;
        }

        /// Returns the reciprocal of the width of each user cell dimension
        static const Tensor<double>& get_rcell_width() {
        	return rcell_width;
        }

        /// Returns the minimum width of any user cell dimension
        static double get_cell_min_width() {
        	return cell_min_width;
        }

        /// Returns the volume of the user cell
        static double get_cell_volume() {
        	return cell_volume;
        }

        /// Returns the default process map that was last initialized via set_default_pmap()
        static std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() {
          if (pmap)
            return pmap;
          else { // try to initialize to default, if not yet
            if (initialized()) {
              set_default_pmap(World::get_default());
              return pmap;
            }
            else
              return pmap; // null ptr if uninitialized
          }
        }

    /// Returns number of the default processes returned by get_pmap()
    static int get_pmap_nproc() {
        return pmap_nproc;
    }

    /// Returns the default process map that can be used with the given world
    static std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > get_pmap(World& world) {
        if (get_pmap_nproc() == world.nproc()) {
            return get_pmap();
        }
        else {
            return make_default_pmap(world);
        }
    }

        /// Sets the default process map (does \em not redistribute existing functions)

        /// Existing functions are probably rendered useless
        static void set_pmap(const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& value) {
        	pmap = value;
        }

        /// Makes a default process map for the given \p world
        static std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > make_default_pmap(World& world);

    	/// Sets the default process map for the use with WorldObjects in \p world
    	/// @note sets default pmap to `make_default_pmap(world)`
        static void set_default_pmap(World& world);

        /// Sets the default process map and redistributes all functions using the old map
        static void redistribute(World& world, const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& newpmap) {
        	pmap->redistribute(world,newpmap);
        	pmap = newpmap;
        }

        MADNESS_PRAGMA_CLANG(diagnostic pop)

    };

    /// Convert user coords (cell[][]) to simulation coords ([0,1]^ndim)
    template <std::size_t NDIM>
    static inline void user_to_sim(const Vector<double,NDIM>& xuser, Vector<double,NDIM>& xsim) {
        for (std::size_t i=0; i<NDIM; ++i)
            xsim[i] = (xuser[i] - FunctionDefaults<NDIM>::get_cell()(i,0)) * FunctionDefaults<NDIM>::get_rcell_width()[i];
    }

    /// Returns the box at level n that contains the given point in simulation coordinates
    /// @param[in] pt point in simulation coordinates
    /// @param[in] n the level of the box
    template <typename T, std::size_t NDIM>
    static inline Key<NDIM> simpt2key(const Vector<T,NDIM>& pt, Level n){
        Vector<Translation,NDIM> l;
        double twon = std::pow(2.0, double(n));
        for (std::size_t i=0; i<NDIM; ++i) {
            l[i] = Translation(twon*pt[i]);
        }
        return Key<NDIM>(n,l);
    }

    /// Convert simulation coords ([0,1]^ndim) to user coords (FunctionDefaults<NDIM>::get_cell())
    template <std::size_t NDIM>
    static void sim_to_user(const Vector<double,NDIM>& xsim, Vector<double,NDIM>& xuser) {
        for (std::size_t i=0; i<NDIM; ++i)
            xuser[i] = xsim[i]*FunctionDefaults<NDIM>::get_cell_width()[i] + FunctionDefaults<NDIM>::get_cell()(i,0);
    }


}
#endif // MADNESS_MRA_FUNCDEFAULTS_H__INCLUDED
