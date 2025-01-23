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
            Translation twon = Translation(1)<<n;

            if (bmax > (twon-1)) bmax=twon-1;

            // Make permissible 1D translations
            Translation b[4*bmax+1];
            int i=0;
            for (Translation lx=-bmax; lx<=bmax; ++lx) {
                b[i++] = lx;
                if ((lx < 0) && (lx+twon > bmax)) b[i++] = lx + twon;
                if ((lx > 0) && (lx-twon <-bmax)) b[i++] = lx - twon;
            }
            MADNESS_ASSERT(i <= 4*bmax+1);
            int numb = i;

            MADNESS_PRAGMA_CLANG(diagnostic push)
            MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

            disp_periodic[n] = std::vector< Key<NDIM> >();
            Vector<long,NDIM> lim(numb);
            for (IndexIterator index(lim); index; ++index) {
                Vector<Translation,NDIM> d;
                for (std::size_t i=0; i<NDIM; ++i) {
                    d[i] = b[index[i]];
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
        /// are generated for all axes. If need to use periodic boundary conditions
        /// for some axes only, make sure to set the boundary conditions appropriately
        /// before the first call to this
        Displacements() {
          MADNESS_PRAGMA_CLANG(diagnostic push)
          MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

          if (disp.empty()) {
                make_disp(bmax_default());
          }

          if constexpr (NDIM <= 3) {
            if (disp_periodic[0].empty()) {
              if (FunctionDefaults<NDIM>::get_bc().is_periodic().any())
                periodic_axes = FunctionDefaults<NDIM>::get_bc().is_periodic();
              else
                periodic_axes = decltype(periodic_axes){true};
              Level nmax = 8 * sizeof(Translation) - 2;
              for (Level n = 0; n < nmax; ++n)
                make_disp_periodic(bmax_default(), n);
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
                if (kernel_lattice_sum_axes != periodic_axes) {
                  std::string msg =
                      "Displacements<" + std::to_string(NDIM) +
                      ">::get_disp(level, kernel_lattice_summed): kernel_lattice_summed differs from the boundary conditions FunctionDefault's had when Displacements were initialized; on-demand periodic displacements generation is not supported";
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
      using Displacement = Key<NDIM>;
      using BoxRadius = std::array<std::optional<Translation>, NDIM>;  // null radius = unlimited size
      using SurfaceThickness = std::array<std::optional<Translation>, NDIM>;  // null thickness for dimensions with null radius
      using Box = std::array<std::pair<Translation, Translation>, NDIM>;
      using Hollowness = std::array<bool, NDIM>;
      using Filter = std::function<bool(const Point&, const Displacement&)>;

      Point center_;                          ///< Center point of the box
      BoxRadius box_radius_;                  ///< halved size of the box in each dimension
      SurfaceThickness
          surface_thickness_;    ///< surface thickness in each dimension
      Box box_;                  ///< box bounds in each dimension
      Hollowness hollowness_;    ///< does box contain non-surface points?
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

          // return true if have another surface layer
          auto next_surface_layer = [this]() -> bool {
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

          for (size_t i = NDIM - 1; i > 0; --i) {
            if (i == fixed_dim) continue;

            if (point[i] < box[i].second) {
              increment_along_dim(i);
              return;
            }
            reset_along_dim(i);
          }

          // move to the face on the opposite side of the fixed dimension
          const bool have_another_surface_layer = next_surface_layer();
          if (have_another_surface_layer) return;

          // ready to switch to next fixed dimension with finite radius
          // but first update box bounds to exclude the surface displacements for the current fixed dimension
          // WARNING if box along this dimension is not hollow we are done!
          if (!parent->hollowness_[fixed_dim]) {
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
          while (!parent->box_radius_[fixed_dim] && fixed_dim <= NDIM) {
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
            while (!done && !parent->filter_(point, this->displacement())) {
              ++(*this);
            }
          }
        }

        void reset_along_dim(size_t dim) {
          Vector<Translation, NDIM> l = point.translation();
          if (dim != fixed_dim)
            l[dim] = box[dim].first;
          else
            l[dim] = parent->box_[dim].first - parent->surface_thickness_[dim].value_or(0);
          point = Point(point.level(), l);
        };

        /**
         * @return displacement from the center to the current point
         */
        const Displacement& displacement() const {
          if (!disp) {
            disp = madness::displacement(parent->center_, point);
          }
          return *disp;
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
            : parent(p), point(parent->center_.level()), fixed_dim(type == End ? NDIM : 0), box(parent->box_), done(type == End) {
          if (type != End) {
            for (size_t i = 0; i < NDIM; ++i) {
              reset_along_dim(i);
            }
            advance_till_valid();
          }
        }

        /**
         * @brief Dereferences the iterator
         * @return A const reference to the current displacement
         */
        reference operator*() const { return displacement(); }

        /**
         * @brief Arrow operator for member access
         * @return A const pointer to the current displacement
         */
        pointer operator->() const { return &displacement(); }

        /**
         * @brief Pre-increment operator
         * @return Reference to this iterator after advancement
         */
        Iterator& operator++() {
          advance();
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
       * @param filter Optional filter function (if returns false, displacement is dropped; default: no filter)
       * @throws std::invalid_argument if any size is not positive
       */
      explicit BoxSurfaceDisplacementRange(const Key<NDIM>& center,
                                           const std::array<std::optional<std::int64_t>, NDIM>& box_radius,
                                           const std::array<std::optional<std::int64_t>, NDIM>& surface_thickness,
                                           Filter filter = {})
          : center_(center), box_radius_(box_radius),
            surface_thickness_(surface_thickness), filter_(std::move(filter)) {
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
          MADNESS_ASSERT(!(box_radius[d].has_value() ^ surface_thickness[d].has_value()));
          MADNESS_ASSERT(surface_thickness[d].value_or(0) >= 0);
          hollowness_[d] = surface_thickness[d] ? (box_[d].first + surface_thickness[d].value() < box_[d].second - surface_thickness[d].value()) : false;
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
       * @return 'probing" displacement to a nearby point *on* the surface; it may not necessarily be in the range of iteration (e.g., it may not be able to pass the filter) but is representative of the surface displacements for the purposes of screening
       */
      const Displacement& probing_displacement() const {
        return probing_displacement_;
      }
    };  // BoxSurfaceDisplacementRange

   /**
    * @brief Generates all M-sized combinations of integers from range [0, N)
    *
    * This class generates all possible combinations of M unique integers from the range
    * [0, N) in lexicographic order. It implements the C++ range concept, allowing usage
    * with range-based for loops and standard algorithms.
    *
    * @tparam N The upper bound of the range (exclusive)
    * @tparam M The size of each combination
    *
    * @note N must be greater than or equal to M
    * @internal generated by Claude Sonnet 3.5 model
    */
    template<std::size_t N, std::size_t M>
    class CombinationGenerator {
      static_assert(N >= M, "N must be greater than or equal to M");

    public:
      /**
     * @brief Iterator class for generating combinations
     *
     * Implements input iterator requirements for generating combinations
     * in lexicographic order.
       */
      class iterator {
      public:
        /// Iterator trait definitions
        using iterator_category = std::input_iterator_tag;
        using value_type = std::array<std::size_t, M>;
        using difference_type = std::ptrdiff_t;
        using pointer = const value_type*;
        using reference = const value_type&;

        /**
         * @brief Constructs an end iterator
         */
        constexpr iterator() : done(true) {}

        /**
         * @brief Constructs either a begin or end iterator
         * @param done_ If true, constructs an end iterator; if false, constructs a begin iterator
         */
        constexpr iterator(bool done_) : done(done_) {
          if (!done) {
            // Initialize first combination
            for (std::size_t i = 0; i < M; ++i) {
              current[i] = i;
            }
          }
        }

        /**
         * @brief Dereference operator
         * @return Reference to the current combination
         */
        constexpr reference operator*() const { return current; }

        /**
         * @brief Arrow operator
         * @return Pointer to the current combination
         */
        constexpr pointer operator->() const { return &current; }

        /**
         * @brief Pre-increment operator
         *
         * Advances the iterator to the next combination in lexicographic order.
         *
         * @return Reference to this iterator
         */
        constexpr iterator& operator++() {
          // Find rightmost element that can be incremented
          int i = M - 1;
          while (i >= 0 && current[i] == N - M + i) {
            --i;
          }

          // If no such element exists, we're done
          if (i < 0) {
            done = true;
            return *this;
          }

          // Increment the found element
          ++current[i];

          // Reset all elements to the right
          for (std::size_t j = i + 1; j < M; ++j) {
            current[j] = current[j-1] + 1;
          }

          return *this;
        }

        /**
         * @brief Post-increment operator
         * @return Copy of the iterator before increment
         */
        constexpr iterator operator++(int) {
          iterator tmp = *this;
          ++*this;
          return tmp;
        }

        /**
         * @brief Equality comparison operator
         * @param other Iterator to compare with
         * @return True if both iterators are in the same state
         */
        constexpr bool operator==(const iterator& other) const {
          return done == other.done;
        }

        /**
         * @brief Inequality comparison operator
         * @param other Iterator to compare with
         * @return True if iterators are in different states
         */
        constexpr bool operator!=(const iterator& other) const {
          return !(*this == other);
        }

      private:
        value_type current = {};
        bool done;
      };

      /**
     * @brief Returns an iterator to the beginning of the combination sequence
     * @return Iterator pointing to the first combination
       */
      constexpr iterator begin() const { return iterator(false); }

      /**
     * @brief Returns an iterator to the end of the combination sequence
     * @return Iterator representing the end of combinations
       */
      constexpr iterator end() const { return iterator(true); }
    };

    /**
 * @brief Helper function to create a CombinationGenerator
 *
 * @tparam N The upper bound of the range (exclusive)
 * @tparam M The size of each combination
 * @return CombinationGenerator<N, M> Generator for M-sized combinations from [0, N)
     */
    template<std::size_t N, std::size_t M>
    constexpr auto make_combinations() {
      return CombinationGenerator<N, M>();
    }

    /**
     * Generates points at the finite-thickness surface of an N-dimensional box [C1-L1,C1+L1]x...x[CN-LN,CN+LN] centered at point {C1,...CN} in Z^N.
     * For finite thickness T={T1,...,TN} point {x1,...,xN} is at the surface face perpendicular to axis i xi>=Ci-Li-Ti and xi<=Ci-Li+Ti OR xi>=Ci+Li-Ti and xi<=Ci+Li+Ti.
     * For dimensions with unlimited size the point coordinates are limited to [0,2^n], with n being the level of the box.
     *
     * Surface layer displacements are generated (1) only for NFIXEDDIM>=1 prescribed fixed dimensions, and (2) displacements for the remaining dimensions
     * are produced by Displacements<NDIM-NFIXEDDIM> (sorted in the order of their magnitude);
     * the latter filters out the displacements that touch the surface layers along the finite-range dimensions
     *
     * Examples:
     * - NDIM=3, NFIXEDDIM=1: generates displacements to the surface layer of a 3D box clustered around the centers of the box faces.
     * - NDIM=3, NFIXEDDIM=2: generates displacements to the surface layer of a 3D box clustered around the centers of the box edges.
     * - NDIM=3, NFIXEDDIM=3: generates displacements to the surface layer of a 3D box clustered around the box vertices.
     */
    template<std::size_t NDIM, std::size_t NFIXEDDIM>
    class BoxSurfaceDisplacementRangeV2 {
    private:
      static_assert(NDIM >= NFIXEDDIM, "BoxSurfaceDisplacementRangeV2: # of fixed dimensions cannot exceed the total # of dimensions");
      constexpr static inline std::size_t NFREEDIM = NDIM - NFIXEDDIM;

      using Point = Key<NDIM>;
      using FixedPoint = Key<NFIXEDDIM>;
      using FreePoint = Key<NFREEDIM>;
      using Displacement = Key<NDIM>;
      using BoxRadius = std::array<std::optional<Translation>, NDIM>;  // null radius = unlimited size
      using SurfaceThickness = std::array<std::optional<Translation>, NDIM>;  // null thickness for dimensions with null radius
      using Box = std::array<std::pair<Translation, Translation>, NDIM>;
      using Hollowness = std::array<bool, NDIM>;
      using Filter = std::function<bool(const Point&, const Displacement&)>;

      Point center_;                          ///< Center point of the box
      BoxRadius box_radius_;                  ///< halved size of the box in each dimension
      SurfaceThickness
          surface_thickness_;    ///< surface thickness in each dimension
      Box box_;                  ///< box bounds in each dimension
      Hollowness hollowness_;    ///< does box contain non-surface points?

      std::array<std::size_t, NFIXEDDIM> fixed_dimensions_;  /// indices of the fixed dimensions
      std::array<std::size_t, NFREEDIM> free_dimensions_;  /// indices of the free dimensions
      const std::vector<Key<NFREEDIM>>* free_displacements_ = nullptr;

      Filter filter_;  ///< optional filter function

      bool empty_ = false;  ///< if true, begin() == end(); this is true if for any `d` `box_radius_[fixed_dimensions_[d]]` is null

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
        const BoxSurfaceDisplacementRangeV2* parent;///< Pointer to parent box
        std::size_t free_displacement_ordinal;      ///< index the displacement of free dimensions in Displacements<NFREEDIM>::get_disp()
        FixedPoint point;                           ///< Current point in fixed dimensions
        mutable std::optional<Displacement> disp;   ///< Memoized displacement from parent->center_ to point, computed by displacement(), reset by advance()
        bool done;                                  ///< Flag indicating iteration completion

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

          auto next_free_displacement = [this]() -> bool {
            if constexpr (NFREEDIM == 0) return false;

            MADNESS_ASSERT(parent->free_displacements_);
            ++free_displacement_ordinal;
            // for the dimensions with finite radius make sure we do not touch the box surface layers
            const auto touch_surface_along_finite_radius_dimensions = [&,this]() -> bool {
              for(int d=0; d!=NFREEDIM; ++d) {
                const auto dim = parent->free_dimensions_[d];
                if (parent->box_radius_[dim]) {
                  const auto l = parent->center_[dim] + parent->free_displacements_->operator[](free_displacement_ordinal)[d];
                  if (l <= parent->box_[dim].first +
                               parent->surface_thickness_[dim].value_or(0) ||
                      l >= parent->box_[dim].second -
                               parent->surface_thickness_[dim].value_or(0)) {
                    return true;
                  }
                }
              }
              return false;
            };
            // find the next free displacement
            while (parent->free_displacements_->size() >
                   free_displacement_ordinal) {
              if (touch_surface_along_finite_radius_dimensions()) ++free_displacement_ordinal;
              else return true;
            }
            free_displacement_ordinal = 0;
            return false;
          };

          // return true if have another surface layer along fixed dimension fixed_dim
          auto next_surface_layer = [this](std::size_t fixed_dim) -> bool {
            MADNESS_ASSERT(fixed_dim < NFIXEDDIM);
            const auto dim = parent->fixed_dimensions_.at(fixed_dim);
            Vector<Translation, NFIXEDDIM> l = point.translation();
            if (l[fixed_dim] !=
                parent->box_[dim].second +
                    parent->surface_thickness_[dim].value_or(0)) {
              // if box is hollow along this dim (has 2 surface layers) and exhausted all layers on the "negative" side of the fixed dimension, move to the first layer on the "positive" side
              if (parent->hollowness_[dim] &&
                  l[fixed_dim] ==
                      parent->box_[dim].first +
                          parent->surface_thickness_[dim].value_or(0)) {
                l[fixed_dim] = parent->box_[dim].second -
                               parent->surface_thickness_[dim].value_or(0);
              } else
                ++l[fixed_dim];
              point = FixedPoint(point.level(), l);
              return true;
            } else
              return false;
          };

          // return true if have another displacement along fixed dimensions
          auto next_fixed_displacement = [&]() -> bool {
            for (auto d = 0; d < NFIXEDDIM; ++d) {
              const auto have_next_surface_layer = next_surface_layer(d);
              if (have_next_surface_layer)
                return true;
              else {
                reset_along_dim(d);
              }
            }
            return false;
          };

          const auto have_free_displacement = next_free_displacement();
          if (have_free_displacement)
            return;

          const auto have_fixed_displacement = next_fixed_displacement();
          if (have_fixed_displacement)
            return;

          done = true;
        }

        void advance_till_valid() {
          if (parent->filter_) {
            const auto filtered_out = [&]() -> bool {
              const auto& disp = this->displacement();
              const auto pt = this->parent->center_.neighbor(disp);
              return !parent->filter_(pt, disp);
            };

            // if displacement has value, filter has already been applied to it, just advance it
            if (!done && disp) this->advance();

            while (!done && filtered_out()) {
              this->advance();
            }
          }
        }

        void reset_along_dim(size_t fixed_dim) {
          const auto dim = parent->fixed_dimensions_[fixed_dim];
          Vector<Translation, NFIXEDDIM> l = point.translation();
          l[fixed_dim] = parent->box_[dim].first - parent->surface_thickness_[dim].value_or(0);
          point = FixedPoint(point.level(), l);
        };

        /**
         * @return displacement from the center to the current point
         */
        const Displacement& displacement() const {
          if (!disp) {
            Vector<Translation, NDIM> d;
            for (size_t i = 0; i < NFIXEDDIM; ++i) {
              const auto dim = this->parent->fixed_dimensions_[i];
              d[dim] = point[i] - parent->center_[dim];
            }
            for (size_t i = 0; i < NFREEDIM; ++i) {
              MADNESS_ASSERT(parent->free_displacements_);
              d[this->parent->free_dimensions_[i]] =
                  parent->free_displacements_->operator[](free_displacement_ordinal)[i];
            }
            disp = Displacement(parent->center_.level(), d);
          }
          return *disp;
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
         * @param p Pointer to the parent BoxSurfaceDisplacementRangeV2
         * @param type the type of iterator (Begin or End)
         */
        Iterator(const BoxSurfaceDisplacementRangeV2* p, Type type)
            : parent(p), free_displacement_ordinal(0), point(parent->center_.level()), done(type == End) {
          if (type != End) {
            // since need to avoid displacements along free dimensions with finite radius that touch the surface, make sure the box along such dimensions is hollow
            bool have_hollow_box_along_free_dimensions = true;
            for (int d = 0; d != NFREEDIM; ++d) {
              const auto dim = parent->free_dimensions_[d];
              if (parent->box_radius_[dim] && !parent->hollowness_[dim]) {
                have_hollow_box_along_free_dimensions = false;
                break;
              }
            }
            // if have hollow box along free dimensions, advance to the first valid point
            if (have_hollow_box_along_free_dimensions) {
              for (size_t i = 0; i < NFIXEDDIM; ++i) {
                reset_along_dim(i);
              }
              advance_till_valid();
            } else // else there is nothing to do
              done = true;
          }
        }

        /**
         * @brief Dereferences the iterator
         * @return A const reference to the current displacement
         */
        reference operator*() const { return displacement(); }

        /**
         * @brief Arrow operator for member access
         * @return A const pointer to the current displacement
         */
        pointer operator->() const { return &displacement(); }

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
          return a.free_displacement_ordinal == b.free_displacement_ordinal &&
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
       * @param fixed_dimensions Indices of "fixed" dimensions, i.e. those along which displacements are to the range surface; if for any `d` `box_radius[fixed_dimensions[d]]` is null the range is empty
       * @param filter Optional filter function (if returns false, displacement is dropped; default: no filter)
       * @throws std::invalid_argument if any size is not positive
       */
      explicit BoxSurfaceDisplacementRangeV2(const Key<NDIM>& center,
                                           const std::array<std::optional<std::int64_t>, NDIM>& box_radius,
                                           const std::array<std::optional<std::int64_t>, NDIM>& surface_thickness,
                                           const std::array<std::size_t, NFIXEDDIM>& fixed_dimensions,
                                           Filter filter = {})
          : center_(center), box_radius_(box_radius),
            surface_thickness_(surface_thickness), fixed_dimensions_(fixed_dimensions),
            filter_(std::move(filter)) {

        if constexpr (NFREEDIM > 0) {
          free_displacements_ = &Displacements<NFREEDIM>{}.get_disp();
        }

        std::sort(fixed_dimensions_.begin(), fixed_dimensions_.end());
        // make sure fixed dimension indices are in-range ...
        std::for_each(fixed_dimensions_.begin(), fixed_dimensions_.end(), [](size_t d) { MADNESS_ASSERT(d < NDIM); });
        // ... and unique
        MADNESS_ASSERT(std::unique(fixed_dimensions_.begin(), fixed_dimensions_.end()) == fixed_dimensions_.end());

        // initialize box bounds and free dimensions
        bool has_finite_dimensions = false;
        const auto n = center_.level();
        bool has_fixed_dimensions = false;
        std::size_t free_dimension_count = 0;
        for (int d=0; d!= NDIM; ++d) {
          if (box_radius_[d]) {
            has_finite_dimensions = true;
            auto r = *box_radius_[d];  // in units of 2^{n-1}
            r = (n == 0) ? (r+1)/2 : (r * Translation(1) << (n-1));
            MADNESS_ASSERT(r > 0);
            box_[d] = {center_[d] - r, center_[d] + r};
          } else {
            box_[d] = {0, (1 << center_.level()) - 1};
          }

          // fixed dimension business
          const auto fixed = std::find(fixed_dimensions_.begin(), fixed_dimensions_.end(), d) != fixed_dimensions_.end();
          if (fixed) {
            has_fixed_dimensions = true;
            // fixed dimensions only make sense for finite dimensions, otherwise the range is empty
            if (!box_radius_[d]) {
              empty_ = true;
            }
          } else {
            MADNESS_ASSERT(free_dimension_count < NFREEDIM);
            free_dimensions_[free_dimension_count++] = d;
          }
        }
        MADNESS_ASSERT(has_finite_dimensions);
        MADNESS_ASSERT(has_fixed_dimensions);
        for (int d=0; d!= NDIM; ++d) {
          MADNESS_ASSERT(!(box_radius[d].has_value() ^ surface_thickness[d].has_value()));
          MADNESS_ASSERT(surface_thickness[d].value_or(0) >= 0);
          hollowness_[d] = surface_thickness[d] ? (box_[d].first + surface_thickness[d].value() < box_[d].second - surface_thickness[d].value()) : false;
        }
      }

      /**
       * @brief Returns an iterator to the beginning of the surface points
       * @return Iterator pointing to the first surface point
       */
      auto begin() const { return Iterator(this, empty_ ? Iterator::End : Iterator::Begin); }

      /**
       * @brief Returns an iterator to the end of the surface points
       * @return Iterator indicating the end of iteration
       */
      auto end() const { return Iterator(this, Iterator::End); }

      /// @return true if the range is empty
      bool empty() const { return empty_; }

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
        * @return array listing ordinals of fixed dimensions
       */
      const std::array<std::size_t, NFIXEDDIM>& fixed_dimensions() const { return fixed_dimensions_; }

      /**
        * @return array listing ordinals of free dimensions
       */
      const std::array<std::size_t, NFREEDIM>& free_dimensions() const { return free_dimensions_; }

    };  // BoxSurfaceDisplacementRangeV2

}  // namespace madness
#endif // MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
