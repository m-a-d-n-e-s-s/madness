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
#include <optional>
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
        }
        MADNESS_ASSERT(has_finite_dimensions);
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
    };  // BoxSurfaceDisplacementRange

}  // namespace madness
#endif // MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
