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
            else if (NDIM == 3) bmax = 3;
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
            Translation bp[4*bmax+1];
            Translation bnp[2*bmax+1];
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
            for(int i=0; i!=NDIM; ++i) {
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
                MADNESS_ASSERT(n < std::extent_v<decltype(disp_periodic)>);
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
     */
    template<std::size_t NDIM>
    class BoxSurfaceDisplacementRange {
    private:
      using Point = Key<NDIM>;
      using PointPattern = Vector<std::optional<Translation>, NDIM>;
      using Displacement = Key<NDIM>;
      using BoxRadius = std::array<std::optional<Translation>, NDIM>;  // null radius = unlimited size
      using SurfaceThickness = std::array<std::optional<Translation>, NDIM>;  // null thickness for dimensions with null radius
      using Box = std::array<std::pair<Translation, Translation>, NDIM>;
      using Hollowness = std::array<bool, NDIM>;  // this can be uninitialized, unlike array_of_bools
      using Periodicity = array_of_bools<NDIM>;
      /// this callable filters out points and/or displacements; note that the displacement is optional (this use case supports filtering based on point pattern onlu) and non-const to make it possible for the filter function to update the displacement (e.g. to map it back to the simulation cell)
      using Filter = std::function<bool(Level, const PointPattern&, std::optional<Displacement>&)>;

      Point center_;                          ///< Center point of the box
      BoxRadius box_radius_;                  ///< halved size of the box in each dimension
      SurfaceThickness
          surface_thickness_;    ///< surface thickness in each dimension
      Box box_;                  ///< box bounds in each dimension
      Hollowness hollowness_;    ///< does box contain non-surface points along each dimension?
      Periodicity is_periodic_;  ///< which dimensions are periodic?
      Filter filter_;  ///< optional filter function
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
        const BoxSurfaceDisplacementRange* parent;  ///< Pointer to parent box
        Point point;                                ///< Current point
        mutable std::optional<Displacement> disp;   ///< Memoized displacement from parent->center_ to point, computed by displacement(), reset by advance()
        size_t fixed_dim;                           ///< Current fixed dimension (i.e. faces perpendicular to this axis are being iterated over)
        Box box;                                    ///< updated box bounds in each dimension, used to avoid duplicate displacements by excluding the surface displacements for each processed fixed dim
        bool done;                                  ///< Flag indicating iteration completion

        // return true if we have another surface layer
        bool next_surface_layer() {
          Vector<Translation, NDIM> l = point.translation();
          if (l[fixed_dim] !=
              parent->box_[fixed_dim].second +
                  parent->surface_thickness_[fixed_dim].value_or(0)) {
            // if box is hollow along this dim (has 2 surface layers) and exhausted all layers on the "negative" side of the fixed dimension, move to the first layer on the "positive" side
            if (parent->hollowness_[fixed_dim] &&
                l[fixed_dim] ==
                    parent->box_[fixed_dim].first +
                        parent->surface_thickness_[fixed_dim].value_or(0)) {
              l[fixed_dim] =
                  parent->box_[fixed_dim].second -
                  parent->surface_thickness_[fixed_dim].value_or(0);
            } else
              ++l[fixed_dim];
            point = Point(point.level(), l);
            return true;
          } else
            return false;
        };

        /**
         * @brief Advances the iterator to the next surface point
         *
         * This function implements the logic for traversing the box surface by:
         * 1. Incrementing displacement in non-fixed dimensions
         * 2. Switching sides in the fixed dimension when needed
         * 3. Moving to the next fixed dimension when current one is exhausted
         */
        void advance() {
          disp.reset();

          auto increment_along_dim = [this](size_t dim) {
            MADNESS_ASSERT(dim != fixed_dim);
            Vector<Translation, NDIM> unit_displacement(0); unit_displacement[dim] = 1;
            point = point.neighbor(unit_displacement);
          };

          for (int64_t i = NDIM - 1; i >= 0; --i) {
            if (i == fixed_dim) continue;

            if (point[i] < box[i].second) {
              increment_along_dim(i);
              return;
            }
            reset_along_dim(i);
          }

          // move to the next surface layer normal to the fixed dimension
          while (bool have_another_surface_layer = next_surface_layer()) {
            const auto filtered_out = [&,this]() {
              bool result = false;
              const auto& filter = this->parent->filter_;
              if (filter) {
                PointPattern point_pattern;
                point_pattern[fixed_dim] = point[fixed_dim];
                std::optional<Displacement> nulldisp;
                result = !filter(point.level(), point_pattern, nulldisp);
              }
              return result;
            };

            if (!filtered_out())
              return;
          }

          // ready to switch to next fixed dimension with finite radius
          // but first update box bounds to exclude the surface displacements for the current fixed dimension
          // WARNING if box along this dimension is not hollow we are done!
          if (parent->hollowness_[fixed_dim]) {
            box[fixed_dim] = {
                parent->box_[fixed_dim].first +
                    parent->surface_thickness_[fixed_dim].value_or(0) + 1,
                parent->box_[fixed_dim].second -
                    parent->surface_thickness_[fixed_dim].value_or(0) - 1};
          }
          else {
            done = true;
            return;
          }
          // onto next dimension
          ++fixed_dim;
          while (!parent->box_radius_[fixed_dim] && fixed_dim < NDIM) {
            ++fixed_dim;
          }


          if (fixed_dim >= NDIM) {
            done = true;
            return;
          }

          // reset upon moving to the next fixed dimension
          for (size_t i = 0; i < NDIM; ++i) {
            reset_along_dim(i);
          }
        }

        void advance_till_valid() {
          if (parent->filter_) {
            const auto filtered_out = [&]() -> bool {
              this->displacement(); // ensure disp is up to date
              return !parent->filter_(point.level(), point.translation(), disp);
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

        void reset_along_dim(size_t dim) {
          const auto is_fixed_dim = dim == fixed_dim;
          Vector<Translation, NDIM> l = point.translation();
          // for fixed dimension start with first surface layer, else use box lower bound (N.B. it's updated as fixed dimensions change)
          Translation l_dim_min =
              is_fixed_dim
                  ? parent->box_[dim].first -
                        parent->surface_thickness_[dim].value_or(0)
              : box[dim].first;
          // if dimension is periodic, only include *unique* displacements (modulo period)
          if (parent->is_periodic_[dim]) {
            const auto period = 1 << parent->center_.level();
            const Translation l_dim_max =
                is_fixed_dim ? parent->box_[dim].second +
                          parent->surface_thickness_[dim].value_or(0) :
                 box[dim].second;
            const Translation l_dim_min_unique = l_dim_max - period + 1;
            l_dim_min = std::max(l_dim_min, l_dim_max - period + 1);
            // fixed dim only: l_dim_min may not correspond to a surface layer
            // this can only happen if l_dim_min is in the gap between the surface layers
            if (is_fixed_dim && parent->hollowness_[dim]) {
              if (l_dim_min > parent->box_[dim].first +
                                  parent->surface_thickness_[dim].value_or(0)) {
                l_dim_min = std::max(
                    l_dim_min, parent->box_[dim].second -
                                   parent->surface_thickness_[dim].value_or(0));
              }
            }
          }
          l[dim] = l_dim_min;

          point = Point(point.level(), l);

          // if the entire surface layer is filtered out, pick the next one
          if (dim == fixed_dim) {

            const auto filtered_out = [&,this]() {
              bool result = false;
              const auto& filter = this->parent->filter_;
              if (filter) {
                PointPattern point_pattern;
                point_pattern[fixed_dim] = point[fixed_dim];
                std::optional<Displacement> nulldisp;
                result = !filter(point.level(), point_pattern, nulldisp);
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
            for (size_t d = 0; d != NDIM; ++d) {
              // min/max displacements along this axis ... N.B. take into account surface thickness!
              box[d] = parent->box_radius_[d] ? std::pair{parent->box_[d].first -
                                parent->surface_thickness_[d].value_or(0),
                         parent->box_[d].second +
                             parent->surface_thickness_[d].value_or(0)} : parent->box_[d];
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
       * @brief Constructs a box with different sizes for each dimension
       *
       * @param center Center of the box
       * @param box_radius Box radius in each dimension
       * @param surface_thickness Surface thickness in each dimension
       * @param is_periodic whether each dimension is periodic; along periodic range-restricted dimensions only one side of the box is iterated over.
       * @param filter Optional filter function (if returns false, displacement is dropped; default: no filter); it may update the displacement to make it valid as needed (e.g. map displacement to the simulation cell)
       * @pre `box_radius[d]>0 && surface_thickness[d]<=box_radius[d]`
       *
       */
      explicit BoxSurfaceDisplacementRange(const Key<NDIM>& center,
                                           const std::array<std::optional<std::int64_t>, NDIM>& box_radius,
                                           const std::array<std::optional<std::int64_t>, NDIM>& surface_thickness,
                                           const array_of_bools<NDIM>& is_periodic,
                                           Filter filter = {})
          : center_(center), box_radius_(box_radius),
            surface_thickness_(surface_thickness), is_periodic_(is_periodic), filter_(std::move(filter)) {
        // initialize box bounds
        bool has_finite_dimensions = false;
        const auto n = center_.level();
        Vector<Translation, NDIM> probing_displacement_vec(0);
        for (int d=0; d!= NDIM; ++d) {
          if (box_radius_[d]) {
            auto r = *box_radius_[d];  // in units of 2^{n-1}
            r = (n == 0) ? (r+1)/2 : (r * Translation(1) << (n-1));
            MADNESS_ASSERT(r > 0);
            box_[d] = {center_[d] - r, center_[d] + r};
            if (!has_finite_dimensions) // first finite dimension? probing displacement will be nonzero along it, zero along all others
              probing_displacement_vec[d] = r;
            has_finite_dimensions = true;
          } else {
            box_[d] = {0, (1 << center_.level()) - 1};
          }
        }
        MADNESS_ASSERT(has_finite_dimensions);
        probing_displacement_ = Displacement(n, probing_displacement_vec);
        for (int d=0; d!= NDIM; ++d) {
          // surface thickness should be only given for finite-radius dimensions
          MADNESS_ASSERT(!(box_radius_[d].has_value() ^ surface_thickness_[d].has_value()));
          MADNESS_ASSERT(surface_thickness_[d].value_or(0) >= 0);
          hollowness_[d] = surface_thickness_[d] ? (box_[d].first + surface_thickness_[d].value() < box_[d].second - surface_thickness_[d].value()) : false;
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
      const std::array<std::optional<int64_t>, NDIM>& box_radius() const { return box_radius_; }

      /**
        * @return the surface thickness in each dimension
       */
      const std::array<std::optional<int64_t>, NDIM>& surface_thickness() const { return surface_thickness_; }

      /**
       * @return flags indicating whether each dimension is periodic
       */
      const array_of_bools<NDIM>& is_periodic() const { return is_periodic_; }

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
    class BoxSurfaceDisplacementFilter {
    public:
      using Point = Key<NDIM>;
      using PointPattern = Vector<std::optional<Translation>, NDIM>;
      using Displacement = Key<NDIM>;
      using Periodicity = array_of_bools<NDIM>;
      using DistanceSquaredFunc = std::function<Translation(const Displacement&)>;

      /// \param domain_is_infinite whether the domain along each axis is finite (simulation cell) or infinite (the entire axis); if true for a given axis then any destination coordinate is valid, else only values in [0,2^n) are valid
      /// \param domain_is_periodic if true for a given axis, displacement to x and x+2^n are equivalent, hence will be canonicalized to end up in the simulation cell. Periodic axes imply infinite domain.
      /// \param range the kernel range for each axis
      /// \param default_distance_squared function that converts a displacement to its effective distance squared (effective may be different from the real distance squared due to periodicity)
      /// \param max_distsq_reached max effective distance squared reached by standard displacements
      BoxSurfaceDisplacementFilter(
          const array_of_bools<NDIM>& domain_is_infinite,
          const array_of_bools<NDIM>& domain_is_periodic,
          const std::array<KernelRange, NDIM>& range,
          DistanceSquaredFunc default_distance_squared,
          Translation max_distsq_reached
          ) :
              domain_is_infinite_(domain_is_infinite),
              domain_is_periodic_(domain_is_periodic),
              range_(range),
              default_distance_squared_(default_distance_squared),
              max_distsq_reached_(max_distsq_reached)
      {}

      /// Apply filter to a displacement ending up at a point or a group of points (point pattern)

      /// @param level the tree level
      /// @param dest the target point (when all elements are nonnull) or point pattern (when only some are).
      ///        The latter is useful to skip the entire surface layer. The point coordinates are
      ///        only used to determine whether we end up in or out of the domain.
      /// @param displacement the optional displacement; if given then will check if it's among
      ///        the standard displacement and whether it was used as part of the standard displacement
      ///        set; if it has not been used and the operator is lattice summed, the displacement
      ///        will be adjusted to end up in the simulation cell.
      /// @return true if the displacement is to be used
      bool operator()(
          const Level level,
          const PointPattern& dest,
          std::optional<Displacement>&  displacement
      ) const {
        const auto twon = (1 << level);  // number of boxes along an axis

        const auto out_of_domain = [&](const Translation& t) -> bool {
          return t < 0 || t >= twon;
        };

        // check that desp is in the domain
        const bool dest_is_in_domain = [&]() {
          for(auto d=0; d!=NDIM; ++d) {
            // - if domain is periodic, all displacements will be mapped back to the simulation cell by for_each/neighbor
            // - if kernel is lattice summed mapping back to the simulation cell for standard displacements
            //   is done during their construction, and for boundary displacements manually (see IMPORTANT below in this function)
            //   N.B. due to this mapping back into the cell only half of the surface displacements was generated by BoxSurfaceDisplacementRange
            if (domain_is_infinite_[d] || domain_is_periodic_[d]) continue;
            if (dest[d].has_value() && out_of_domain(*dest[d])) return false;
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
            for(auto d=0; d!=NDIM; ++d) {
              const auto disp_d = (*displacement)[d];
              auto bmax_standard = Displacements<NDIM>::bmax_default();

              // the effective displacement length depends on whether lattice summation is performed along it
              // compare Displacements::make_disp vs Displacements::make_disp_periodic
              auto disp_d_eff_abs = std::abs(disp_d);
              if (domain_is_periodic_[d]) {
                // for "periodic" displacements the effective disp_d is the shortest of {disp_d, disp_d+twon, disp_d-twon} ... see make_disp_periodic
                MADNESS_ASSERT(range_[d].N() <= 2);  // displacements that exceed 1 whole cell need bit more complex logic
                // bmax in make_disp_periodic is wrapped around
                if (Displacements<NDIM>::bmax_default() >= twon) bmax_standard = twon-1;
                disp_d_eff_abs = std::min(std::abs(disp_d < 0 ? disp_d + twon : disp_d - twon), disp_d_eff_abs);

                // IMPORTANT for lattice-summed axes, if the destination is out of the simulation cell map the displacement back to the cell
                if (dest[d].has_value()) {
                  if (dest[d] < 0) {
                    auto t = (*displacement).translation();
                    t[d] += twon;
                    displacement.emplace(displacement->level(), t);
                    MADNESS_ASSERT(!out_of_domain(*dest[d] + twon));
                  }
                  else if (dest[d] >= twon) {
                    auto t = (*displacement).translation();
                    t[d] -= twon;
                    displacement.emplace(displacement->level(), t);
                    MADNESS_ASSERT(!out_of_domain(*dest[d] - twon));
                  }
                }
              }
              if (disp_d_eff_abs > bmax_standard) {
                among_standard_displacements = false;
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
      array_of_bools<NDIM> domain_is_infinite_;
      array_of_bools<NDIM> domain_is_periodic_;
      std::array<KernelRange, NDIM> range_;
      DistanceSquaredFunc default_distance_squared_;
      Translation max_distsq_reached_;
    };

}  // namespace madness
#endif // MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
