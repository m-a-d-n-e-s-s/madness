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

  $Id$
*/
#ifndef MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
#define MADNESS_MRA_DISPLACEMENTS_H__INCLUDED

#include <madness/mra/indexit.h>
#include <madness/mra/funcdefaults.h>

#include <algorithm>
#include <array>
#include <functional>
#include <iterator>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace madness {

    // How should we treat destinations "extra" to the [0, 2^n) standard domain?
    enum class ExtraDomainPolicy {
        Discard,  // Use case: most computations.
        Keep,     // Use case: PBC w/o lattice sums. Destinations that arise from a source inside the domain and some displacement but are outside [0, 2^n)
                  //           are equivalent to a destination inside [0, 2^n) with the same displacement but a source outside the [0, 2^n).
                  //           That source needs explicit accounting. Keep it. The caller will correct the destination and (if needed) the source.
                  //           We're only responsible for the displacement.
        Translate // Use case: PBC w/ lattice sums. As above, *except* the source outside [0, 2^n) is accounted for by some term source inside [0, 2^n).
                  //           The displacement itself needs changing, so that both source and destination are in the standard domain.
                  //           We're responsible for changing the displacement.
    };

    /// Holds displacements for applying operators to avoid replicating for all operators
    template <std::size_t NDIM>
    class Displacements {

        static std::vector< Key<NDIM> > disp; ///< standard displacements to be used with standard kernels (range-unrestricted, no lattice sum)
        static array_of_bools<NDIM> periodic_axes;  ///< along which axes lattice summation is performed?
        static std::vector< Key<NDIM> > disp_periodic[64];  ///< displacements to be used with lattice-summed kernels

    public:
        static int bmax_default() {
            int bmax;
            if      (NDIM == 1) bmax = 7;
            else if (NDIM == 2) bmax = 5;
            else if (NDIM == 3) bmax = 4;
            else if (NDIM == 4) bmax = 3;
            else if (NDIM == 5) bmax = 3;
            else if (NDIM == 6) bmax = 3;
            else                bmax = 2;
            return bmax;
        }

    private:
        static bool cmp_keys(const Key<NDIM>& a, const Key<NDIM>& b) {
            return a.distsq() < b.distsq();
        }

        static bool cmp_keys_periodic(const Key<NDIM>& a, const Key<NDIM>& b) {
          return a.distsq_bc(periodic_axes) < b.distsq_bc(periodic_axes);
        }

        static void make_disp(int bmax) {
            // Note newer loop structure in make_disp_periodic_sum
            Vector<Translation,NDIM> d(0);

            int num = 1;
            for (std::size_t i=0; i<NDIM; ++i) num *= (2*bmax + 1);
            disp.resize(num,Key<NDIM>(0));

            num = 0;
            if (NDIM == 1) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 2) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 3) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 4) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            for (d[3]=-bmax; d[3]<=bmax; ++d[3])
                                disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 5) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            for (d[3]=-bmax; d[3]<=bmax; ++d[3])
                                for (d[4]=-bmax; d[4]<=bmax; ++d[4])

                                    disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 6) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            for (d[3]=-bmax; d[3]<=bmax; ++d[3])
                                for (d[4]=-bmax; d[4]<=bmax; ++d[4])
                                    for (d[5]=-bmax; d[5]<=bmax; ++d[5])
                                        disp[num++] = Key<NDIM>(0,d);
            }
            else {
                MADNESS_EXCEPTION("make_disp: hard dimension loop",NDIM);
            }

            std::sort(disp.begin(), disp.end(), cmp_keys);
        }

        static void make_disp_periodic(int bmax, Level n) {
            MADNESS_ASSERT(periodic_axes.any());  // else use make_disp
            Translation twon = Translation(1)<<n;

            if (bmax > (twon-1)) bmax=twon-1;

            // Make permissible 1D translations, periodic and nonperiodic (for mixed BC)
            std::vector<Translation> bp(4*bmax+1);
            std::vector<Translation> bnp(2*bmax+1);
            int ip=0;
            int inp=0;
            for (Translation lx=-bmax; lx<=bmax; ++lx) {
                bp[ip++] = lx;
                if ((lx < 0) && (lx+twon > bmax)) bp[ip++] = lx + twon;
                if ((lx > 0) && (lx-twon <-bmax)) bp[ip++] = lx - twon;
                bnp[inp++] = lx;
            }
            MADNESS_ASSERT(ip <= 4*bmax+1);
            MADNESS_ASSERT(inp <= 2*bmax+1);
            const int nbp = ip;
            const int nbnp = inp;

            MADNESS_PRAGMA_CLANG(diagnostic push)
            MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

            disp_periodic[n] = std::vector< Key<NDIM> >();
            Vector<long,NDIM> lim;
            for(size_t i=0; i!=NDIM; ++i) {
              lim[i] = periodic_axes[i] ? nbp : nbnp;
            }
            for (IndexIterator index(lim); index; ++index) {
                Vector<Translation,NDIM> d;
                for (std::size_t i=0; i<NDIM; ++i) {
                  d[i] = periodic_axes[i] ? bp[index[i]] : bnp[index[i]];
                }
                disp_periodic[n].push_back(Key<NDIM>(n,d));
            }

            std::sort(disp_periodic[n].begin(), disp_periodic[n].end(), cmp_keys_periodic);
//             print("KEYS AT LEVEL", n);
//             print(disp_periodic[n]);

            MADNESS_PRAGMA_CLANG(diagnostic pop)

        }


    public:
        /// first time this is called displacements are generated.
        /// if boundary conditions are not periodic, the periodic displacements
        /// are generated for all axes. This allows to support application of
        /// operators with boundary conditions periodic along any axis (including all).
        /// If need to use periodic boundary conditions
        /// for some axes only, make sure to set the boundary conditions appropriately
        /// before the first call to this
        Displacements() {
          MADNESS_PRAGMA_CLANG(diagnostic push)
          MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

          if (disp.empty()) {
                make_disp(bmax_default());
          }

          if constexpr (NDIM <= 3) {
            if (disp_periodic[0].empty()) {  // if not initialized yet
              if (FunctionDefaults<NDIM>::get_bc().is_periodic().any())
                reset_periodic_axes(
                    FunctionDefaults<NDIM>::get_bc().is_periodic());
              else
                reset_periodic_axes(array_of_bools<NDIM>{true});
            }
          }

          MADNESS_PRAGMA_CLANG(diagnostic pop)
        }

        const std::vector< Key<NDIM> >& get_disp(Level n,
                                                 const array_of_bools<NDIM>& kernel_lattice_sum_axes) {
            MADNESS_PRAGMA_CLANG(diagnostic push)
            MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

            if (kernel_lattice_sum_axes.any()) {
                MADNESS_ASSERT(NDIM <= 3);
                MADNESS_ASSERT((std::size_t)n < std::extent_v<decltype(disp_periodic)>);
                if ((kernel_lattice_sum_axes && periodic_axes) != kernel_lattice_sum_axes) {
                  std::string msg =
                      "Displacements<" + std::to_string(NDIM) +
                      ">::get_disp(level, kernel_lattice_sum_axes): kernel_lattice_sum_axes is set for some axes that were not periodic in the FunctionDefault's boundary conditions active at the time when Displacements were initialized; invoke Displacements<NDIM>::reset_periodic_axes(kernel_lattice_sum_axes) to rebuild the periodic displacements";
                  MADNESS_EXCEPTION(msg.c_str(), 1);
                }
                return disp_periodic[n];
            }
            else {
                return disp;
            }

            MADNESS_PRAGMA_CLANG(diagnostic pop)
        }

        /// return the standard displacements appropriate for operators w/o lattice summation
        const std::vector< Key<NDIM> >& get_disp() {
          MADNESS_PRAGMA_CLANG(diagnostic push)
          MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

          return disp;

          MADNESS_PRAGMA_CLANG(diagnostic pop)
        }

        /// rebuilds periodic displacements so that they are optimal for the given set of periodic axes

        /// this must be done while no references to prior periodic displacements are outstanding (i.e. no operator application
        /// tasks in flight)
        /// \param new_periodic_axes the new periodic axes
        static void reset_periodic_axes(const array_of_bools<NDIM>& new_periodic_axes) {
          MADNESS_PRAGMA_CLANG(diagnostic push)
          MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

          MADNESS_ASSERT(new_periodic_axes.any());  // else why call this?
          if (new_periodic_axes != periodic_axes) {

            periodic_axes = new_periodic_axes;
            Level nmax = 8 * sizeof(Translation) - 2;
            for (Level n = 0; n < nmax; ++n)
              make_disp_periodic(bmax_default(), n);
          }
          MADNESS_PRAGMA_CLANG(diagnostic pop)
        }
    };

    template <std::size_t N, std::size_t M>
    constexpr std::enable_if_t<N>=M, std::array<std::size_t, N-M>> iota_array(std::array<std::size_t, M> values_to_skip_sorted) {
      std::array<std::size_t, N - M> result;
      if constexpr (N != M) {
        std::size_t nadded = 0;
        auto value_to_skip_it = values_to_skip_sorted.begin();
        assert(*value_to_skip_it < N);
        auto value_to_skip = *value_to_skip_it;
        for (std::size_t i = 0; i < N; ++i) {
          if (i < value_to_skip) {
            result[nadded++] = i;
          } else if (value_to_skip_it != values_to_skip_sorted.end()) {
            ++value_to_skip_it;
            if (value_to_skip_it != values_to_skip_sorted.end()) {
              value_to_skip = *value_to_skip_it;
            } else
              value_to_skip = N;
          }
        }
      }
      return result;
    }

    /**
     * Generates points at the finite-thickness surface of an N-dimensional box [C1-L1,C1+L1]x...x[CN-LN,CN+LN] centered at point {C1,...CN} in Z^N.
     * For finite thickness T={T1,...,TN} point {x1,...,xN} is at the surface face perpendicular to axis i xi>=Ci-Li-Ti and xi<=Ci-Li+Ti OR xi>=Ci+Li-Ti and xi<=Ci+Li+Ti.
     * For dimensions with unlimited size the point coordinates are limited to [0,2^n], with n being the level of the box.
     * N.B. "points" are really boxes in the standard MADNESS sense, which we'll call "primitive boxes" to disambiguate from box as the product of intervals mentioned above,
     */
    template<std::size_t NDIM>
    class BoxSurfaceDisplacementRange {
    public:
      using Point = Key<NDIM>;
      using PointPattern = Vector<std::optional<Translation>, NDIM>;
      using Displacement = Key<NDIM>;
      /// this callable returns whether a given primitive box (or hyperface if only one coordinate is provided) can be filtered out. if screening a primitive box, the corresponding displacement should be provided both for further screening and for the displacement to be updated, if displacements are translated to connect two cells in the box.
      /// the validator should normally be a BoxSurfaceDisplacementFilter object. anything else is probably a hack.
      using Validator = std::function<bool(Level, const PointPattern&, std::optional<Displacement>&)>;

    private:
      using SurfaceRadius = std::array<std::optional<Translation>, NDIM>;  // null radius = unlimited size
      using SurfaceThickness = std::array<std::optional<Translation>, NDIM>;  // null thickness for dimensions with null radius
      using Box = std::array<std::pair<Translation, Translation>, NDIM>;
      using Hollowness = std::array<bool, NDIM>;  // this can be uninitialized, unlike array_of_bools ... hollow = gap between -radius+thickness and +radius-thickness.
      using Periodicity = array_of_bools<NDIM>;

      Point center_;                          ///< Center point of the surface's box
      SurfaceRadius surface_radius_;          ///< halved size of the surface's box in each dimension, in half-SimulationCells.
      SurfaceThickness
          surface_thickness_;    ///< surface thickness in each dimension, measured in boxes. Real-space surface size is thus n-dependent.
      Box bounds_;               ///< box bounds in each dimension. Does not include thickness.
      Hollowness hollowness_;    ///< does box contain non-surface points along each dimension?
      Periodicity is_lattice_summed_;  ///< which dimensions are lattice summed?
      Validator validator_;      ///< optional validator function
      // TODO: Double-check legitimacy of choosing probing displacement arbitrarily. For isotropic kernels, wouldn't you want it on the face where the boundary is closest (in real-space)?
      //       That depends on both the surface thickness and the side lengths of the simulation cell.
      Displacement probing_displacement_;  ///< displacement to a nearby point on the surface; it may not be able to pass the filter, but is sufficiently representative of the surface displacements to allow screening with isotropic kernels

      /**
     * @brief Iterator class for lazy generation of surface points
     *
     * This iterator generates surface points on-demand by tracking the current fixed
     * dimension and positions in each dimension. It implements the InputIterator concept.
       */
      class Iterator {
      public:
        enum Type {Begin, End};
      private:
        const BoxSurfaceDisplacementRange* parent;  ///< Pointer to parent surface.
        Point point;                                ///< Current point / box. This is always free to leave the simulation cell.
        mutable std::optional<Displacement> disp;   ///< Memoized displacement from parent->center_ to point, computed by displacement(), reset by advance()
        size_t fixed_dim;                           ///< Current fixed dimension (i.e. faces perpendicular to this axis are being iterated over)
        Box unprocessed_bounds;                     ///< The bounds for all *unprocessed* displacements in the finite-thickness surface. Updated as displacements are processed.
                                                    ///  For the dimensions of the parent box, without thickness or regard for displacement processing, use parent->bounds_.
                                                    ///  Tracking `unprocessed_bounds` allows us to avoid double-counting 'edge' boxes that are on multiple hyperfaces.
                                                    ///  e.g., if radius is [5, 5], center is [0, 0] and thickness is [1, 1], the bounds are [-6, 6] x [-6, 6].
                                                    ///  We first evaluate the hyperfaces [-6, -4] x [-5, 5] and then [4, 6] x [-5, 5].
                                                    ///  It remains to evaluate hyperfaces [-5, 5] x [-6, -4] and [-5, 5] x [4, 6], *excluding*
                                                    ///  the edge points shared with the processed hyperfaces. So, we need to evaluate effective hyperfaces
                                                    ///  [-3, 3] x [-6, -4] and [-3, 3] x [4, 6]. The unprocessed_bounds are reset to [-3, 3] x [-6, 6].
        bool done;                                  ///< Flag indicating iteration completion

        // return true if we have another surface layer for the fixed_dim
        // if we do, translate point onto that next surface layer
        bool next_surface_layer() {
          Vector<Translation, NDIM> l = point.translation();
          if (l[fixed_dim] !=
              parent->bounds_[fixed_dim].second +
                  parent->surface_thickness_[fixed_dim].value_or(0)) {
            // if exhausted all layers on the "negative" side of the fixed dimension and there's a gap to the "positive" side,
            // jump to the positive side. otherwise, just take the next layer.
            if (parent->hollowness_[fixed_dim] &&
                l[fixed_dim] ==
                    parent->bounds_[fixed_dim].first +
                        parent->surface_thickness_[fixed_dim].value_or(0)) {
              l[fixed_dim] =
                  parent->bounds_[fixed_dim].second -
                  parent->surface_thickness_[fixed_dim].value_or(0);
            } else
              ++l[fixed_dim];
            point = Point(point.level(), l);
            disp.reset();
            return true;
          } else
            return false;
        };

        /**
         * @brief Advances the iterator to the next surface point
         *
         * This function implements the logic for traversing the box surface by:
         * (1) Incrementing displacement in non-fixed dimensions
         * (2) Switching sides in the fixed dimension when needed
         * (3) Moving to the next fixed dimension when current one is exhausted
         *
         * We filter out layers in (2) but not points within a layer in (1).
         */
        void advance() {
          disp.reset();

          auto increment_along_dim = [this](size_t dim) {
            MADNESS_ASSERT(dim != fixed_dim);
            Vector<Translation, NDIM> unit_displacement(0); unit_displacement[dim] = 1;
            point = point.neighbor(unit_displacement);
          };

          // (1) try all displacements on current NDIM-1 dim layer
          // loop structure is equivalent to NDIM-1 nested, independent for loops
          // over the NDIM-1 dimension of the layer, with last dim as innermost loop
          for (size_t i = NDIM; i > 0; --i) {
            const size_t cur_dim = i - 1;
            if (cur_dim == fixed_dim) continue;

            if (point[cur_dim] < unprocessed_bounds[cur_dim].second) {
              increment_along_dim(cur_dim);
              return;
            }
            reset_along_dim(cur_dim);
          }

          // (2) move to the next surface layer normal to the fixed dimension
          // if we can filter out the entire layer, do so.
          while (next_surface_layer()) {
            const auto filtered_out = [&,this]() {
              bool result = false;
              const auto& validator = this->parent->validator_;
              if (validator) {
                PointPattern point_pattern;
                point_pattern[fixed_dim] = point[fixed_dim];
                std::optional<Displacement> nulldisp;
                result = !validator(point.level(), point_pattern, nulldisp);
              }
              return result;
            };

            if (!filtered_out())
              return;
          }

          // we finished this fixed dimension, so update unprocessed bounds to exclude the layers of the current fixed dimension
          // if box along this dimension is not hollow, the new interval would be [0, 0] - we are done!
          if (parent->hollowness_[fixed_dim]) {
            unprocessed_bounds[fixed_dim] = {
                parent->bounds_[fixed_dim].first +
                    parent->surface_thickness_[fixed_dim].value_or(0) + 1,
                parent->bounds_[fixed_dim].second -
                    parent->surface_thickness_[fixed_dim].value_or(0) - 1};
          }
          else {
            done = true;
            return;
          }
          // (3) switch to next fixed dimension with finite radius
          ++fixed_dim;
          while (!parent->surface_radius_[fixed_dim] && fixed_dim < NDIM) {
            ++fixed_dim;
          }

          // Exit if we've displaced along all dimensions of finite radius
          if (fixed_dim >= NDIM) {
            done = true;
            return;
          }

          // reset our search along all non-fixed dimensions
          // the reset along the fixed_dim returns silently
          for (size_t i = 0; i < NDIM; ++i) {
            reset_along_dim(i);
          }
        }

        /// Perform advance, repeating if you are at a filtered point
        void advance_till_valid() {
          if (parent->validator_) {
            const auto filtered_out = [&]() -> bool {
              this->displacement(); // ensure disp is up to date
              return !parent->validator_(point.level(), point.translation(), disp);
            };

            // if displacement has value, filter has already been applied to it, just advance it
            if (!done && disp) this->advance();

            while (!done && filtered_out()) {
              this->advance();
            }
          }
          else
            this->advance();
        }

        // Recall that the surface is a union of hyperfaces, i.e., direct products of intervals.
        // Reset state on dimension `dim` to initialize for the start of interval `dim` in the the current direct product
        void reset_along_dim(size_t dim) {
          const auto is_fixed_dim = dim == fixed_dim;
          Vector<Translation, NDIM> l = point.translation();
          Translation l_dim_min;
          if (!is_fixed_dim) {
            // This dimension is contiguous boxes on the hyperface.
            // Initialize to the start.
            l_dim_min = unprocessed_bounds[dim].first;
          } else if (!parent->is_lattice_summed_[dim]) {
            // This dimension consists of two finite-thickness hyperfaces, not lattice summed.
            // Initialize to the start of the - hyperface. We trust next_surface_layer()
            // to move to the + hyperface when ready.
            l_dim_min = parent->bounds_[dim].first -
                        parent->surface_thickness_[dim].value_or(0);
          } else {
            // This dimension consists of two finite-thickness hyperfaces, lattice summed.
            // The two hyperfaces are the same interval shifted by parent->surface_radius_[dim]
            // periods. So by lattice summation, the - hyperface is included. Initialize
            // to the start of the + hyperface.
            l_dim_min = parent->bounds_[dim].second -
                        parent->surface_thickness_[dim].value_or(0);
          }
          if (parent->is_lattice_summed_[dim]) {
            // By lattice summation, boxes that differ by a SimulationCell are equivalent.
            // Therefore, we need to sum over equivalence classes and not displacements.
            const auto period = 1 << parent->center_.level();
            const Translation last_equiv_class = is_fixed_dim ? parent->bounds_[dim].second +
              parent->surface_thickness_[dim].value_or(0) : unprocessed_bounds[dim].second;
            const Translation first_equiv_class = last_equiv_class - period + 1;
            l_dim_min = std::max(first_equiv_class, l_dim_min);
          }
          l[dim] = l_dim_min;

          point = Point(point.level(), l);
          disp.reset();

          // if the entire surface layer is filtered out, pick the next one
          if (dim == fixed_dim) {

            const auto filtered_out = [&,this]() {
              bool result = false;
              const auto& validator = this->parent->validator_;
              if (validator) {
                PointPattern point_pattern;
                point_pattern[fixed_dim] = point[fixed_dim];
                std::optional<Displacement> nulldisp;
                result = !validator(point.level(), point_pattern, nulldisp);
              }
              return result;
            };

            if (filtered_out()) {
              bool have_another_surface_layer;
              while ((have_another_surface_layer = next_surface_layer())) {
                if (!filtered_out())
                  break;
              }
              MADNESS_ASSERT(have_another_surface_layer);
            }

          }
        };

        /**
         * @return displacement from the center to the current point
         */
        const std::optional<Displacement>& displacement() const {
          if (!disp) {
            disp = madness::displacement(parent->center_, point);
          }
          return disp;
        }

      public:
        // Iterator type definitions for STL compatibility
        using iterator_category = std::input_iterator_tag;
        using value_type = Point;
        using difference_type = std::ptrdiff_t;
        using pointer = const Point*;
        using reference = const Point&;

        /**
         * @brief Constructs an iterator
         *
         * @param p Pointer to the parent BoxSurfaceDisplacementRange
         * @param type the type of iterator (Begin or End)
         */
        Iterator(const BoxSurfaceDisplacementRange* p, Type type)
            : parent(p), point(parent->center_.level()), fixed_dim(type == End ? NDIM : 0), done(type == End) {
          if (type != End) {
            // skip to first dimensions with limited range
            while (!parent->surface_radius_[fixed_dim] && fixed_dim < NDIM) {
              ++fixed_dim;
            }

            for (size_t d = 0; d != NDIM; ++d) {
              // min/max displacements along this axis ... N.B. take into account surface thickness!
              unprocessed_bounds[d] = parent->surface_radius_[d] ? std::pair{parent->bounds_[d].first -
                                parent->surface_thickness_[d].value_or(0),
                         parent->bounds_[d].second +
                             parent->surface_thickness_[d].value_or(0)} : parent->bounds_[d];
              reset_along_dim(d);
            }
            advance_till_valid();
          }
        }

        /**
         * @brief Dereferences the iterator
         * @return A const reference to the current displacement
         */
        reference operator*() const { return *displacement(); }

        /**
         * @brief Arrow operator for member access
         * @return A const pointer to the current displacement
         */
        pointer operator->() const { return &(*(*this)); }

        /**
         * @brief Pre-increment operator
         * @return Reference to this iterator after advancement
         */
        Iterator& operator++() {
          advance_till_valid();
          return *this;
        }

        /**
         * @brief Post-increment operator
         * @return Copy of the iterator before advancement
         */
        Iterator operator++(int) {
          Iterator tmp = *this;
          ++(*this);
          return tmp;
        }

        /**
         * @brief Equality comparison operator
         * @param a First iterator
         * @param b Second iterator
         * @return true if iterators are equivalent
         */
        friend bool operator==(const Iterator& a, const Iterator& b) {
          if (a.done && b.done) return true;
          if (a.done || b.done) return false;
          return a.fixed_dim == b.fixed_dim &&
                 a.point == b.point;
        }

        /**
         * @brief Inequality comparison operator
         * @param a First iterator
         * @param b Second iterator
         * @return true if iterators are not equivalent
         */
        friend bool operator!=(const Iterator& a, const Iterator& b) {
          return !(a == b);
        }
      };

      friend class Iterator;

    public:
      /**
       * @brief Constructs a box with different radii and thicknesses for each dimension
       *
       * @param center Center primitive box of the box. All displacements will share the `n` of this arg.
       * @param surface_radius Surface radius in each dimension, in half-SimulationCells. Omit for dim `i` to signal that the bound for dim `i` is simply the simulation cell.
       * @param surface_thickness Surface thickness in each dimension, measured in number of addl. boxes *on each half* of the surface box proper. Omit for dim `i` if and only if omitted in `surface_radius`
       * @param is_lattice_summed whether each dimension is periodic; along periodic range-restricted dimensions only one side of the box is iterated over.
       * @param validator Optional filter function (if returns false, displacement is dropped; default: no filter); it may update the displacement to make it valid as needed (e.g. map displacement to the simulation cell)
       * @pre `surface_radius[d]>0 && surface_thickness[d]<=surface_radius[d]`
       *
       */
      explicit BoxSurfaceDisplacementRange(const Key<NDIM>& center,
                                           const std::array<std::optional<std::int64_t>, NDIM>& surface_radius,
                                           const std::array<std::optional<std::int64_t>, NDIM>& surface_thickness,
                                           const array_of_bools<NDIM>& is_lattice_summed,
                                           Validator validator = {})
          : center_(center), surface_radius_(surface_radius),
            surface_thickness_(surface_thickness), is_lattice_summed_(is_lattice_summed), validator_(std::move(validator)) {
        // initialize bounds
        bool has_finite_dimensions = false;
        const auto n = center_.level();
        Vector<Translation, NDIM> probing_displacement_vec(0);
        for (size_t d=0; d!= NDIM; ++d) {
          if (surface_radius_[d]) {
            auto r = *surface_radius_[d];  // in units of 2^{n-1}
            // n = 0 is special b/c << -1 is undefined
            r = (n == 0) ? (r+1)/2 : (r * Translation(1) << (n-1));
            MADNESS_ASSERT(r > 0);
            bounds_[d] = {center_[d] - r, center_[d] + r};
            if (!has_finite_dimensions) // first finite dimension? probing displacement will be nonzero along it, zero along all others
              probing_displacement_vec[d] = r;
            has_finite_dimensions = true;
          } else {
            bounds_[d] = {0, (1 << center_.level()) - 1};
          }
        }
        MADNESS_ASSERT(has_finite_dimensions);
        probing_displacement_ = Displacement(n, probing_displacement_vec);
        for (size_t d=0; d!= NDIM; ++d) {
          // surface thickness should be only given for finite-radius dimensions
          MADNESS_ASSERT(!(surface_radius_[d].has_value() ^ surface_thickness_[d].has_value()));
          MADNESS_ASSERT(surface_thickness_[d].value_or(0) >= 0);
          hollowness_[d] = surface_thickness_[d] ? (bounds_[d].first + surface_thickness_[d].value() < bounds_[d].second - surface_thickness_[d].value()) : false;
        }
      }

      /**
     * @brief Returns an iterator to the beginning of the surface points
     * @return Iterator pointing to the first surface point
       */
      auto begin() const { return Iterator(this, Iterator::Begin); }

      /**
     * @brief Returns an iterator to the end of the surface points
     * @return Iterator indicating the end of iteration
       */
      auto end() const { return Iterator(this, Iterator::End); }

      //      /**
      //     * @brief Returns a view over the surface points
      //     *
      //     * This operator allows the class to be used with C++20 ranges.
      //     *
      //     * @return A view over the surface points
      //       */
      //      auto operator()() const {
      //        return std::ranges::subrange(begin(), end());
      //      }

      /* @return the center of the box
       */
      const Key<NDIM>& center() const { return center_; }

      /**
        * @return the radius of the box in each dimension
       */
      const std::array<std::optional<int64_t>, NDIM>& box_radius() const { return surface_radius_; }

      /**
        * @return the surface thickness in each dimension
       */
      const std::array<std::optional<int64_t>, NDIM>& surface_thickness() const { return surface_thickness_; }

      /**
       * @return flags indicating whether each dimension is lattice summed
       */
      const array_of_bools<NDIM>& is_lattice_summed() const { return is_lattice_summed_; }

      /**
       * @return 'probing" displacement to a nearby point *on* the surface; it may not necessarily be in the range of iteration (e.g., it may not be able to pass the filter) but is representative of the surface displacements for the purposes of screening
       */
      const Displacement& probing_displacement() const {
        return probing_displacement_;
      }
    };  // BoxSurfaceDisplacementRange


    /// This is used to filter out box surface displacements that
    /// - take us outside of the target domain, or
    /// - were already utilized as part of the the standard displacements list.
    /// For dealing with the lattice-summed operators the filter
    /// can adjusts the displacement to make sure that we end up in
    /// the simulation cell.
    template <size_t NDIM>
    class BoxSurfaceDisplacementValidator {
    public:
      using Point = Key<NDIM>;
      using PointPattern = Vector<std::optional<Translation>, NDIM>;
      using Displacement = Key<NDIM>;
      using Periodicity = array_of_bools<NDIM>;
      using DistanceSquaredFunc = std::function<Translation(const Displacement&)>;

      /// \param is_infinite_domain whether the domain along each axis is finite (simulation cell) or infinite (the entire axis); if true for a given axis then any destination coordinate is valid, else only values in [0,2^n) are valid
      /// \param is_lattice_summed if true for a given axis, displacement to x and x+2^n are equivalent, hence will be canonicalized to end up in the simulation cell. Periodic axes imply infinite domain, whatever was passed to `domain_is_infinite`.
      /// \param range the kernel range for each axis
      /// \param default_distance_squared function that converts a displacement to its effective distance squared (effective may be different from the real distance squared due to periodicity)
      /// \param max_distsq_reached max effective distance squared reached by standard displacements
      BoxSurfaceDisplacementValidator(
          const array_of_bools<NDIM>& is_infinite_domain,
          const array_of_bools<NDIM>& is_lattice_summed,
          const std::array<KernelRange, NDIM>& range,
          DistanceSquaredFunc default_distance_squared,
          Translation max_distsq_reached
          ) :
              range_(range),
              default_distance_squared_(default_distance_squared),
              max_distsq_reached_(max_distsq_reached) {
        for (size_t i = 0; i < NDIM; i++) {
          if (is_lattice_summed[i]) {
            domain_policies_[i] = ExtraDomainPolicy::Translate;
          } else if (is_infinite_domain[i]) {
            domain_policies_[i] = ExtraDomainPolicy::Keep;
          } else {
            domain_policies_[i] = ExtraDomainPolicy::Discard;
          }
        }
      }

      /// Apply filter to a displacement ending up at a point or a group of points (point pattern)

      /// @param level the tree level
      /// @param dest the target point (when all elements are nonnull) or point
      /// pattern (when only some are).
      ///        The latter is useful to skip the entire surface layer. The
      ///        point coordinates are only used to determine whether we end up
      ///        in or out of the domain.
      /// @param displacement the optional displacement; if given then will check if it's among
      ///        the standard displacement and whether it was used as part of
      ///        the standard displacement set; if it has not been used and the
      ///        operator is lattice summed, the displacement will be adjusted
      ///        to end up in the simulation cell. Primary use case for omitting `displacement`
      ///        is if `dest` is not equivalent to a point.
      /// @return true if the displacement is to be used
      bool operator()(const Level level, const PointPattern& dest,
          std::optional<Displacement>&  displacement
      ) const {
        // preliminaries
        const auto twon = (static_cast<Translation>(1) << level);  // number of boxes along an axis
        // map_to_range_twon(x) returns for x >= 0 ? x % 2^level : map_to_range_twon(x+2^level)
        // idiv is generally slow, so instead use bit logic that relies on 2's complement representation of integers
        const auto map_to_range_twon = [&, mask = ((~(static_cast<std::uint64_t>(0)) << (64-level)) >> (64-level))](std::int64_t x) -> std::int64_t {
          const std::int64_t x_mapped = x & mask;
          MADNESS_ASSERT(x_mapped >=0 && x_mapped < twon && (std::abs(x_mapped-x)%twon==0));
          return x_mapped;
        };

        const auto out_of_domain = [&](const Translation& t) -> bool {
          return t < 0 || t >= twon;
        };

        // check that dest is in the domain
        const bool dest_is_in_domain = [&]() {
          for(size_t d=0; d!=NDIM; ++d) {
            if (domain_policies_[d] == ExtraDomainPolicy::Discard && dest[d].has_value() && out_of_domain(*dest[d])) return false;
          }
          return true;
        }();

        if (dest_is_in_domain) {
          if (displacement.has_value()) {

            // N.B. avoid duplicates of standard displacements previously included:
            // A displacement has been possibly considered if along EVERY axis the "effective" displacement size
            // fits within the box explored by the standard displacement.
            // If so, skip if <= max magnitude of standard displacements encountered
            // Otherwise this is a new non-standard displacement, consider it
            bool among_standard_displacements = true;
            for(size_t d=0; d!=NDIM; ++d) {
              const auto disp_d = (*displacement)[d];
              auto bmax_standard = Displacements<NDIM>::bmax_default();

              // the effective displacement length depends on whether lattice summation is performed along it
              // compare Displacements::make_disp vs Displacements::make_disp_periodic
              auto disp_d_eff_abs = std::abs(disp_d);
              if (domain_policies_[d] == ExtraDomainPolicy::Translate) {
                MADNESS_ASSERT(range_[d].N() <= 2);  // displacements that exceed 1 whole cell need bit more complex logic

                // for "periodic" displacements the effective disp_d is the shortest of {..., disp_d-twon, disp_d, disp_d+twon, ...} ... see make_disp_periodic
                const std::int64_t disp_d_eff = map_to_range_twon(disp_d);
                disp_d_eff_abs = std::min(disp_d_eff,std::abs(disp_d_eff-twon));

                // IMPORTANT for lattice-summed axes, if the destination is out of the simulation cell map the displacement back to the cell
                // same logic as for disp_d: dest[d] -> dest[d] % twon
                if (dest[d].has_value()) {
                  const Translation dest_d = dest[d].value();
                  const auto dest_d_in_cell = map_to_range_twon(dest_d);
                  MADNESS_ASSERT(!out_of_domain(
                      dest_d_in_cell));
                  // adjust displacement[d] so that it produces dest_d_cell, not dest_d
                  auto t = (*displacement).translation();
                  t[d] += (dest_d_in_cell - dest_d);
                  displacement.emplace(displacement->level(), t);
                }

                // N.B. bmax in make_disp_periodic is clipped in the same way
                if (Displacements<NDIM>::bmax_default() >= twon) bmax_standard = twon-1;
              }

              if (disp_d_eff_abs > bmax_standard) {
                among_standard_displacements = false;
                // Do not break - this loop needs not only to determine among_standard_displacements but to shift the displacement if domain_is_periodic_
                // Therefore, looping over all dim is strictly necessary.
              }
            }
            if (among_standard_displacements) {
              const auto distsq = default_distance_squared_(*displacement);
              if (distsq > max_distsq_reached_) { // among standard displacements => keep if longer than the longest standard displacement considered
                return true;
              } {
                return false;
              }
            }
            else  // not among standard displacements => keep it
              return true;
          }
          else  // skip the displacement-based filter if not given
            return true;
        }
        else
          return false;
      }

    private:
      std::array<ExtraDomainPolicy, 3> domain_policies_;
      std::array<KernelRange, NDIM> range_;
      DistanceSquaredFunc default_distance_squared_;
      Translation max_distsq_reached_;
    };

}  // namespace madness
#endif // MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
