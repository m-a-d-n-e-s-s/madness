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
#include <madness/mra/key.h>
#include <madness/misc/array_of_bools.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace madness {

    template <std::size_t NDIM> class FunctionDefaults;  // funcdefaults.h includes this header at its end

    // How should we treat destinations "extra" to the [0, 2^n) standard domain?
    enum class ExtraDomainPolicy {
        Discard,  // Use case: most computations.
        Keep,     // Use case: PBC w/o lattice sums. Destinations that arise from a source inside the domain and some displacement but are outside [0, 2^n)
                  //           are equivalent to a destination inside [0, 2^n) with the same displacement but a source outside the [0, 2^n).
                  //           That source needs explicit accounting. Keep it. The caller will correct the destination and (if needed) the source.
                  //           We're only responsible for the displacement.
        Translate // Use case: PBC w/ lattice sums. As above, *except* the source outside [0, 2^n) is accounted for by some source inside [0, 2^n).
                  //           The displacement itself needs changing, so that both source and destination are in the standard domain.
                  //           We're responsible for changing the displacement.
    };

    /// Holds displacements for applying operators to avoid replicating for all operators
    template <std::size_t NDIM>
    class Displacements {

        inline static std::vector< Key<NDIM> > disp = {}; ///< standard displacements to be used with standard kernels (range-unrestricted, no lattice sum)
        inline static array_of_bools<NDIM> periodic_axes{false};  ///< along which axes lattice summation is performed?
        inline static std::array<std::vector< Key<NDIM>>, 64 > disp_periodic{};  ///< displacements to be used with lattice-summed kernels
        inline static Tensor<double> widths{NDIM}; ///< cell width, used to order displacements from least to most real space distance

    public:
        static int bmax_default() {
            // Numbers determined by trial and error. The entire idea of bmax is non-adaptive,
            // and the decision to have bmax be isotropic is only valid for hypercubes.
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
            const auto a_width = a.real_distsq(widths);
            const auto b_width = b.real_distsq(widths);
            if (a_width == 0 and a_width == b_width) return a.distsq() < b.distsq();
            else return a_width < b_width;
        }

        static bool cmp_keys_periodic(const Key<NDIM>& a, const Key<NDIM>& b) {
            const auto a_width = a.real_distsq_bc(periodic_axes, widths);
            const auto b_width = b.real_distsq_bc(periodic_axes, widths);
            if (a_width == 0 and a_width == b_width) return a.distsq_bc(periodic_axes) < b.distsq_bc(periodic_axes);
            else return a_width < b_width;
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

          if (widths.normf() < 1e-8) widths = 1;

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
                MADNESS_ASSERT(n < disp_periodic.size());
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

        /// Sets the cell widths used to order the displacements by real-space distance, and reorders them.
        /// Called by FunctionDefaults::recompute_cell_info whenever the cell changes.
        /// @warning not thread safe: must not be called while operators are being applied
        static void set_width(const Tensor<double>& width) {
          widths = copy(width);
          if (!disp.empty()) {
            std::sort(disp.begin(), disp.end(), cmp_keys);
          }
          for (size_t n = 0; n < 64; ++n) {
            if (!disp_periodic[n].empty()) {
              std::sort(disp_periodic[n].begin(), disp_periodic[n].end(), cmp_keys_periodic);
            }
          }
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
    /// Real-space extent of the standard (short-range) displacements, i.e. of what FunctionImpl::do_apply processes
    /// before turning to the surface of the kernel range boundary (see Displacements). BoxSurfaceDisplacementValidator
    /// uses it to skip the surface displacements already processed, and BoxSurfaceDisplacementRange to place its
    /// probing displacements just outside of them; both must therefore see the same reach.
    template <std::size_t NDIM>
    struct StandardDisplacementsReach {
      double max_distsq;                    ///< max real distance squared reached by the standard displacements (see Key::real_distsq_bc)
      std::array<double, NDIM> cell_width;  ///< real-space width of the simulation cell along each axis, as used to compute `max_distsq`
    };

    /// Filters the surface displacements produced by BoxSurfaceDisplacementRange: drops the destinations outside of
    /// the domain, maps the destinations along lattice-summed axes into the simulation cell, and drops the displacements
    /// already processed as standard (short-range) displacements, as described by StandardDisplacementsReach.
    template <size_t NDIM>
    class BoxSurfaceDisplacementValidator {
    public:
      using Point = Key<NDIM>;
      using PointPattern = Vector<std::optional<Translation>, NDIM>;
      using Displacement = Key<NDIM>;
      using Periodicity = array_of_bools<NDIM>;
      using Reach = StandardDisplacementsReach<NDIM>;

      /// \param is_infinite_domain whether the domain along each axis is finite (simulation cell) or infinite (the entire axis); if true for a given axis then any destination coordinate is valid, else only values in [0,2^n) are valid
      /// \param is_lattice_summed if true for a given axis, displacement to x and x+2^n are equivalent, hence will be canonicalized to end up in the simulation cell. Periodic axes imply infinite domain, whatever was passed to `is_infinite_domain`.
      /// \param reach the real-space extent of the standard displacements that have been processed; surface displacements
      ///        within it are filtered out as duplicates. Omit if no standard displacements have been processed (nothing is filtered on that account).
      BoxSurfaceDisplacementValidator(
          const array_of_bools<NDIM>& is_infinite_domain,
          const array_of_bools<NDIM>& is_lattice_summed,
          std::optional<Reach> reach = {}
          ) :
              is_lattice_summed_(is_lattice_summed),
              reach_(std::move(reach)),
              cell_width_(NDIM) {
        for (size_t i = 0; i < NDIM; i++) {
          if (is_lattice_summed[i]) {
            domain_policies_[i] = ExtraDomainPolicy::Translate;
          } else if (is_infinite_domain[i]) {
            domain_policies_[i] = ExtraDomainPolicy::Keep;
          } else {
            domain_policies_[i] = ExtraDomainPolicy::Discard;
          }
          if (reach_) {
            MADNESS_CHECK_THROW(reach_->cell_width[i] > 0, "BoxSurfaceDisplacementValidator: cell widths in StandardDisplacementsReach must be positive");
            cell_width_(i) = reach_->cell_width[i];
          }
        }
      }

      /// @return which axes are lattice summed
      const array_of_bools<NDIM>& is_lattice_summed() const { return is_lattice_summed_; }

      /// @return the real-space extent of the standard displacements this filters out as duplicates; null if none
      const std::optional<Reach>& reach() const { return reach_; }

      /// Apply filter to a displacement ending up at a point or a group of points (point pattern)

      /// @param level the tree level
      /// @param dest the target point (when all elements are nonnull) or point pattern (when only some are).
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
      bool operator()(
          const Level level,
          const PointPattern& dest,
          std::optional<Displacement>&  displacement
      ) const {
        // preliminaries
        const auto twon = (static_cast<Translation>(1) << level);  // number of boxes along an axis
        // map_to_range_twon(x) returns for x >= 0 ? x % 2^level : map_to_range_twon(x+2^level)
        // idiv is generally slow, so instead use bit logic that relies on 2's complement representation of integers
        const auto map_to_range_twon = [&, mask = level == 0 ? std::uint64_t(0) : ((~(static_cast<std::uint64_t>(0)) << (64-level)) >> (64-level))](std::int64_t x) -> std::int64_t {
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
              // N.B. if lattice summation is performed along any axis the standard displacements come from
              // Displacements::make_disp_periodic, which clips bmax to 2^n-1 along *every* axis
              auto bmax_standard = Displacements<NDIM>::bmax_default();
              if (is_lattice_summed_.any() && bmax_standard >= twon) bmax_standard = twon - 1;

              // the effective displacement length depends on whether lattice summation is performed along it
              // compare Displacements::make_disp vs Displacements::make_disp_periodic
              auto disp_d_eff_abs = std::abs(disp_d);
              if (domain_policies_[d] == ExtraDomainPolicy::Translate) {
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
              }

              if (disp_d_eff_abs > bmax_standard) {
                among_standard_displacements = false;
                // Do not break - this loop needs not only to determine among_standard_displacements but to shift the displacement if domain_is_periodic_
                // Therefore, looping over all dim is strictly necessary.
              }
            }
            if (among_standard_displacements) {
              if (!reach_) return true;  // no standard displacements were processed => nothing to duplicate
              // among standard displacements => keep if longer than the longest standard displacement considered
              // N.B. same distance as used to order the standard displacements (see FunctionImpl::do_apply)
              const auto distsq = displacement->real_distsq_bc(is_lattice_summed_, cell_width_);
              return distsq > reach_->max_distsq;
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
      std::array<ExtraDomainPolicy, NDIM> domain_policies_;
      array_of_bools<NDIM> is_lattice_summed_;
      std::optional<Reach> reach_;
      Tensor<double> cell_width_;  ///< reach_->cell_width as a Tensor, for Key::real_distsq_bc
    };


    template<std::size_t NDIM>
    class BoxSurfaceDisplacementRange {
    public:
      using Point = Key<NDIM>;
      using PointPattern = Vector<std::optional<Translation>, NDIM>;
      using Displacement = Key<NDIM>;
      using Validator = BoxSurfaceDisplacementValidator<NDIM>;

    private:
      using BoxRadius = std::array<std::optional<Translation>, NDIM>;  // null radius = unlimited size
      using SurfaceThickness = std::array<std::optional<Translation>, NDIM>;  // null thickness for dimensions with null radius
      using Box = std::array<std::pair<Translation, Translation>, NDIM>;
      using Hollowness = std::array<bool, NDIM>;  // this can be uninitialized, unlike array_of_bools ... hollow = there are boxes between the faces, besides those of the faces themselves
      using Periodicity = array_of_bools<NDIM>;

      Point center_;                          ///< Center point of the box
      BoxRadius box_radius_;          ///< halved size of the box in each dimension, in half-SimulationCells.
      SurfaceThickness
          surface_thickness_;    ///< surface thickness in each dimension, measured in boxes. Real-space surface size is thus n-dependent.
      Box box_;                  ///< box bounds in each dimension.
      Box initial_bounds_;       ///< bounds of the boxes to be iterated over, before any face is processed: the box plus its surface thickness, or, along a lattice-summed dimension, one period ending at the top layer (so that each equivalence class of boxes appears exactly once)
      Hollowness hollowness_;    ///< does box contain non-surface points along each dimension?
      Periodicity is_lattice_summed_;  ///< which dimensions are lattice summed?
      std::optional<Validator> validator_;  ///< optional filter; also the source of the reach of the standard displacements, which the probing displacements are placed outside of
      std::array<std::optional<Displacement>, NDIM> probing_displacements_;  ///< for each finite-radius dimension, a displacement to a nearby point on the faces normal to it (the pair of hyperplanes at -radius and +radius, which lattice summation folds onto each other); it may not be able to pass the filter, but among the displacements those faces contribute it errs toward the largest norm, so that a decaying kernel can be screened with it
      std::array<bool, NDIM> skip_face_{};  ///< faces excluded from iteration (see skip_face())

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
                                                    ///  For the dimensions of the parent box, without thickness or regard for displacement processing, use parent->box_.
                                                    ///  Tracking `unprocessed_bounds` allows us to avoid double-counting 'edge' boxes that are on multiple hyperfaces.
                                                    ///  e.g., if radius is [5, 5], center is [0, 0] and thickness is [1, 1], the bounds are [-6, 6] x [-6, 6].
                                                    ///  We first evaluate the hyperfaces [-6, -4] x [-5, 5] and then [4, 6] x [-5, 5].
                                                    ///  It remains to evaluate hyperfaces [-5, 5] x [-6, -4] and [-5, 5] x [4, 6], *excluding*
                                                    ///  the edge points shared with the processed hyperfaces. So, we need to evaluate effective hyperfaces
                                                    ///  [-3, 3] x [-6, -4] and [-3, 3] x [4, 6]. The unprocessed_bounds are reset to [-3, 3] x [-6, 6].
        bool done;                                  ///< Flag indicating iteration completion
        bool positioned = false;                    ///< whether the iterator is positioned on a point that has been (or is about to be) yielded; false until the first advance_till_valid() completes

        // return true if we have another surface layer for the fixed_dim
        // if we do, translate point onto that next surface layer
        bool next_surface_layer() {
          Vector<Translation, NDIM> l = point.translation();
          if (l[fixed_dim] !=
              parent->box_[fixed_dim].second +
                  parent->surface_thickness_[fixed_dim].value_or(0)) {
            // if exhausted all layers on the "negative" side of the fixed dimension and there's a gap to the "positive" side,
            // jump to the positive side. otherwise, just take the next layer.
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
                result = !(*validator)(point.level(), point_pattern, nulldisp);
              }
              return result;
            };

            if (!filtered_out())
              return;
          }

          // (3) we finished this fixed dimension: move on to the next face
          next_face();
        }

        /// Positions the iterator on the first point of the current face (`fixed_dim`)
        /// @return false if every layer of the face is filtered out, i.e. the face has no point to offer
        bool start_face() {
          bool has_layer = true;
          for (size_t i = 0; i < NDIM; ++i) {
            if (!reset_along_dim(i)) has_layer = false;
          }
          return has_layer;
        }

        /// Leaves the current face (finished, or without any layer to offer) for the next one that has a point
        /// to offer, excluding the layers of the faces left behind from the remaining ones. Sets `done` if none remains.
        void next_face() {
          do {
            if (!exclude_face(fixed_dim)) {
              done = true;
              return;
            }
            select_face(fixed_dim + 1);
            if (done) return;
          } while (!start_face());
        }

        /// Excludes the layers of the faces normal to `dim` from the faces that remain to be processed,
        /// so that the edge boxes shared with them are not visited twice.
        /// @return false if nothing remains, i.e. the box along `dim` is not hollow and every remaining box lies on these faces
        bool exclude_face(size_t dim) {
          if (!parent->hollowness_[dim]) return false;
          // the layers are at both ends of the bounds, or only at the top end if lattice summed (see initial_bounds_)
          const auto nlayers = 2 * parent->surface_thickness_[dim].value_or(0) + 1;
          unprocessed_bounds[dim] = {unprocessed_bounds[dim].first + (parent->is_lattice_summed_[dim] ? 0 : nlayers),
                                     unprocessed_bounds[dim].second - nlayers};
          return true;
        }

        /// Sets `fixed_dim` to the first finite-radius dimension at or after `from` whose faces are not skipped.
        /// Sets `done` if no face remains.
        /// N.B. skipped faces are not excluded from the remaining faces (unlike processed ones), so the
        /// edge boxes they share with them are still visited through them; hence a box is visited iff it
        /// lies on at least one face that is not skipped, regardless of the order in which faces are visited.
        void select_face(size_t from) {
          for (fixed_dim = from; fixed_dim < NDIM; ++fixed_dim) {
            if (parent->box_radius_[fixed_dim] && !parent->skip_face_[fixed_dim]) return;
          }
          done = true;
        }

        /// Leave the current point (if positioned on one) and advance to the next point that passes the filter
        void advance_till_valid() {
          if (positioned && !done) this->advance();
          positioned = true;

          if (parent->validator_) {
            const auto filtered_out = [&]() -> bool {
              this->displacement(); // ensure disp is up to date
              return !(*parent->validator_)(point.level(), point.translation(), disp);
            };

            while (!done && filtered_out()) {
              this->advance();
            }
          }
        }

        // Recall that the surface is a union of hyperfaces, i.e., direct products of intervals.
        // Reset state on dimension `dim` to initialize for the start of interval `dim` in the the current direct product
        // @return false if `dim` is the fixed dimension and every layer of its faces is filtered out
        bool reset_along_dim(size_t dim) {
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
            l_dim_min = parent->box_[dim].first -
                        parent->surface_thickness_[dim].value_or(0);
          } else {
            // This dimension consists of two finite-thickness hyperfaces, lattice summed.
            // The two hyperfaces are the same interval shifted by parent->surface_radius_[dim]
            // periods. So by lattice summation, the - hyperface is included. Initialize
            // to the start of the + hyperface, clipped to one period (= the bounds) in case
            // the layers are thicker than the simulation cell.
            l_dim_min = std::max(parent->box_[dim].second -
                                     parent->surface_thickness_[dim].value_or(0),
                                 unprocessed_bounds[dim].first);
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
                result = !(*validator)(point.level(), point_pattern, nulldisp);
              }
              return result;
            };

            if (filtered_out()) {
              bool have_another_surface_layer;
              while ((have_another_surface_layer = next_surface_layer())) {
                if (!filtered_out())
                  break;
              }
              return have_another_surface_layer;  // false: every layer of this face is filtered out (e.g. lies outside the domain)
            }

          }
          return true;
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
            unprocessed_bounds = parent->initial_bounds_;

            // skip to first dimension with limited range whose faces are not skipped and have a point to offer
            select_face(0);
            if (done) return;
            if (!start_face()) next_face();
            if (done) return;

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
       * @param box_radius Box radius in each dimension, in half-SimulationCells. Omit for dim `i` to signal that the bound for dim `i` is simply the simulation cell.
       * @param surface_thickness Surface thickness in each dimension, measured in number of addl. boxes *on each half* of the surface box proper. Omit for dim `i` if and only if omitted in `box_radius`
       * @param is_lattice_summed whether each dimension is lattice summed; along lattice summed dimensions only one side of the box is iterated over.
       * @param validator Optional filter (if returns false, displacement is dropped; default: no filter); it also maps displacements
       *        along lattice-summed axes into the simulation cell, and carries the real-space reach of the standard displacements
       *        it filters out as duplicates, outside of which the probing displacements are placed. Its lattice-summation flags
       *        must match `is_lattice_summed`. If omitted (or if it carries no reach) nothing is known to be filtered out: the
       *        surface then reaches all the way in to `center`, and the probes fall back to the on-site displacement, which screens nothing.
       * @pre `surface_radius[d]>0 && surface_thickness[d]<=surface_radius[d]`
       */
      explicit BoxSurfaceDisplacementRange(const Key<NDIM>& center,
                                           const std::array<std::optional<std::int64_t>, NDIM>& box_radius,
                                           const std::array<std::optional<std::int64_t>, NDIM>& surface_thickness,
                                           const array_of_bools<NDIM>& is_lattice_summed,
                                           std::optional<Validator> validator = {})
          : center_(center), box_radius_(box_radius),
            surface_thickness_(surface_thickness), is_lattice_summed_(is_lattice_summed), validator_(std::move(validator)) {
        if (validator_) {
          for (size_t d=0; d!= NDIM; ++d)
            MADNESS_CHECK_THROW(validator_->is_lattice_summed()[d] == is_lattice_summed_[d],
                                "BoxSurfaceDisplacementRange: validator and range disagree on which axes are lattice summed");
        }
        // initialize bounds
        bool has_finite_dimensions = false;
        const auto n = center_.level();
        const auto period = Translation(1) << n;
        for (size_t d=0; d!= NDIM; ++d) {
          if (box_radius_[d]) {
            auto r = *box_radius_[d];  // in units of 2^{n-1}
            // n = 0 is special b/c << -1 is undefined
            r = (n == 0) ? (r+1)/2 : (r * Translation(1) << (n-1));
            MADNESS_ASSERT(r > 0);
            box_[d] = {center_[d] - r, center_[d] + r};
            has_finite_dimensions = true;
          } else {
            box_[d] = {0, (1 << center_.level()) - 1};
          }
        }
        MADNESS_ASSERT(has_finite_dimensions);
        for (size_t d=0; d!= NDIM; ++d) {
          if (box_radius_[d]) probing_displacements_[d] = compute_probing_displacement(d);
        }
        for (size_t d=0; d!= NDIM; ++d) {
          // surface thickness should be only given for finite-radius dimensions
          MADNESS_ASSERT(!(box_radius_[d].has_value() ^ surface_thickness_[d].has_value()));
          MADNESS_ASSERT(surface_thickness_[d].value_or(0) >= 0);
          const auto t = surface_thickness_[d].value_or(0);
          if (box_radius_[d]) {
            // the boxes to iterate over: the box plus its surface thickness. Along a lattice-summed dimension the
            // box is at least one simulation cell wide, so instead take one period ending at the top layer: each
            // equivalence class of boxes then appears exactly once, and only the top layers are on the surface.
            initial_bounds_[d] = is_lattice_summed_[d] ? std::pair{box_[d].second + t - period + 1, box_[d].second + t}
                                                       : std::pair{box_[d].first - t, box_[d].second + t};
            // hollow = the bounds hold more boxes than the layers of the faces (both ends, or the top end if lattice summed)
            const auto nlayers = (is_lattice_summed_[d] ? 1 : 2) * (2 * t + 1);
            hollowness_[d] = (initial_bounds_[d].second - initial_bounds_[d].first + 1) > nlayers;
          } else {
            initial_bounds_[d] = box_[d];
            hollowness_[d] = false;
          }
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
       * @return flags indicating whether each dimension is lattice summed
       */
      const array_of_bools<NDIM>& is_lattice_summed() const { return is_lattice_summed_; }

      /**
       * @param face_dimension a dimension with finite radius; its faces are the pair of hyperplanes normal to it at -radius and +radius (which lattice summation folds onto each other)
       * @return "probing" displacement to a nearby point *on* the faces normal to `face_dimension`; it may not necessarily be in the range of iteration (e.g., it may not be able to pass the filter) but, among the displacements those faces contribute, it errs toward the largest norm, so that a decaying kernel can be screened with it. One probe serves both faces since the real-space distance of a displacement depends on its magnitude along each axis only.
       */
      const Displacement& probing_displacement(size_t face_dimension) const {
        MADNESS_ASSERT(face_dimension < NDIM && probing_displacements_[face_dimension].has_value());
        return *probing_displacements_[face_dimension];
      }

      /**
       * @return probing displacements for the faces normal to every dimension; null for dimensions of unlimited size, which have no faces
       * @sa probing_displacement()
       */
      const std::array<std::optional<Displacement>, NDIM>& probing_displacements() const {
        return probing_displacements_;
      }

      /**
       * Excludes the faces normal to `face_dimension` (both, if not lattice summed) from iteration, e.g. because their
       * probing displacement showed their contributions to be negligible. The edge boxes they share with faces that are
       * not skipped are still visited through those faces, i.e. a box is visited iff it lies on at least one face that is not skipped.
       * @param face_dimension a dimension with finite radius
       * @pre no iterator has been obtained from this object yet
       */
      void skip_face(size_t face_dimension) {
        MADNESS_ASSERT(face_dimension < NDIM && box_radius_[face_dimension].has_value());
        skip_face_[face_dimension] = true;
      }

      /**
       * @return whether the faces normal to `face_dimension` are excluded from iteration
       */
      bool face_skipped(size_t face_dimension) const {
        return skip_face_[face_dimension];
      }

    private:
      Displacement compute_probing_displacement(const size_t face_dimension) const {
        // Large boxes we must consider are both those near the center (because 1/r is large
        // for small r), and near the box radius (because going from 1/r to 0 is a sharp change).
        // The probe displacement is a way to screen out cases where the box radius is negligible.
        // Each face (the pair of hyperplanes normal to a dimension with finite box_radius_, or one hyperplane
        // if that dimension is lattice summed) gets its own probe, so that faces at different real-space
        // distances (anisotropic cells, lattice summation along some dimensions only, mixed-parity radii)
        // can be screened independently. Our probe displacement for the face normal to face_dimension must satisfy:
        // (1) It must actually be on that face.
        // To ensure we're probing the box radius effect and not the near-center effect, we require:
        // (2) If at all possible, it must be distinct from the zero displacement and
        //     from the displacements "near" the center, both of which should have already been considered.
        //     Zero displacements are especially pernicious, because self-interaction is always large.
        //  n.b.: Beware that for lattice summed-dimensions, displacements must be distinct in the space
        //     of equivalence classes. For lattice-summed dimensions of an even number of boxes, the origin
        //     of the target face is equivalent the origin.
        //  n.b.: If N even and lattice-summed and 1D, the entire boundary is already equivalent to
        //     the displacements "near" the center.
        // To keep the estimate sharp, we prefer:
        // (3) We want the displacement of minimal real-space r within the above constraints.
        //     Such displacements are more suitable as a heuristic upper bound of the matrix element
        //     controlling the 1/r to 0 change. Not explicitly accounting for this does not seem to
        //     affect whether we're within epsilon, but it's still good practice.
        //     For the same reason, the sigma should matter as well.

        MADNESS_ASSERT(face_dimension < NDIM && box_radius_[face_dimension].has_value());
        const auto face_origin_is_center = [this](size_t d) {
          return is_lattice_summed_[d] && (*box_radius_[d] % 2 == 0);
        };

        // Enforce requirement (1). The faces have finite thickness: their layers span [r-t, r+t], and
        // by (3) the probe goes on the innermost one, which is the nearest to the source.
        Vector<Translation, NDIM> probing_displacement_vec(0);
        const auto n = center_.level();
        auto r = *box_radius_[face_dimension];  // in units of 2^{n-1}
        // n = 0 is special b/c << -1 is undefined
        r = (n == 0) ? (r+1)/2 : (r * Translation(1) << (n-1));
        MADNESS_ASSERT(r > 0);
        probing_displacement_vec[face_dimension] = r - surface_thickness_[face_dimension].value_or(0);
        // Along a lattice-summed dimension fold the probe into the cell, to the representative nearest to the
        // source, as the validator does for the displacements it yields (the operator's norm only sums over a
        // few lattice images of a displacement, so a representative several cells away would be underestimated).
        if (is_lattice_summed_[face_dimension]) {
          const auto period = Translation(1) << n;
          auto& l = probing_displacement_vec[face_dimension];
          l = ((l % period) + period) % period;
          if (l > period / 2) l -= period;
        }

        // In these cases, requirement (2) is already satisfied or unsatisfiable.
        // Choosing 0 for all other dimensions satisfies requirement (3).
        if (!face_origin_is_center(face_dimension) || n == 0 || NDIM == 1)
          return Displacement(n, probing_displacement_vec);

        // If nothing is known to be filtered out, none of the surface points have been processed,
        // so the surface reaches all the way in to center_ and the on-site probe is the only safe choice.
        if (!validator_ || !validator_->reach())
          return Displacement(n, probing_displacement_vec);
        const auto& reach = *validator_->reach();

        // Else, we still need to satisfy requirement (2) while trying to obey (3). We need to displace along
        // a different dimension.

        // The offset along axis d is the least number of boxes that takes us out of the region covered by
        // the standard displacements, which the validator filters out (see BoxSurfaceDisplacementValidator):
        // either beyond bmax boxes (see Displacements::make_disp), or, within bmax, beyond sqrt(max_distsq) in real space.
        // Key::real_distsq_bc measures cell_width*(|l|-1) along an axis, so invert that.
        // Cap the offset at half a cell. If our dimension is lattice-summed, it's even, and half a cell
        // is where it's furthest from the origin. Else, half a cell is the furthest away we can
        // guarantee we can displace to, in the case of an open dimension and the center_ is the origin.
        const Translation half_cell = Translation(1) << (n-1);
        const Translation bmax = Displacements<NDIM>::bmax_default();
        const auto offset_along = [&](size_t d) -> Translation {
          const double width = reach.cell_width[d];  // positive, checked by the validator
          const Translation nboxes = 1 + static_cast<Translation>(std::sqrt(reach.max_distsq) / width);
          return std::min(std::min(nboxes, bmax) + 1, half_cell);
        };
        // real-space distance of the offset; since the face axis folds to zero this is the probe's real distance
        const auto offset_distance = [&](size_t d) -> double {
          return reach.cell_width[d] * (offset_along(d) - 1);
        };

        // choose the dimension to displace along: the one with the least real-space offset (requirement (3)),
        // which for anisotropic cells need not be the narrowest one in boxes. Break ties in favor of
        // finite dimensions (offset is guaranteed to stay on the face) with the smallest radius, then by index.
        const auto offset_sort_key = [&](size_t d) {
          return std::make_tuple(offset_distance(d), !box_radius_[d].has_value(), box_radius_[d].value_or(0));
        };
        size_t offset_dimension = NDIM;
        for (size_t d=0; d != NDIM; ++d) {
          if (d == face_dimension) continue;
          if (offset_dimension == NDIM || offset_sort_key(d) < offset_sort_key(offset_dimension))
            offset_dimension = d;
        }
        MADNESS_ASSERT(offset_dimension != NDIM);  // NDIM > 1, so some dimension was found

        const auto d = offset_dimension;
        const Translation offset = offset_along(d);
        if (box_radius_[d]) {
          // the offset stays on the face: box_radius_ >= 1 means the box spans at least a half
          // simulation cell along this dimension, and offset <= half_cell
          probing_displacement_vec[d] = offset;
        } else {
          // we're bounded by the simulation cell; displace toward whichever side of center_ has more room
          const auto left_distance = center_[d] - box_[d].first;
          const auto right_distance = box_[d].second - center_[d];
          const auto sign = right_distance >= left_distance ? +1 : -1;
          probing_displacement_vec[d] = sign * offset;
        }
        return Displacement(n, probing_displacement_vec);
      }   // compute_probing_displacement
    };  // BoxSurfaceDisplacementRange


    /// This is used to filter out box surface displacements that
    /// - take us outside of the target domain, or
    /// - were already utilized as part of the the standard displacements list.
    /// For dealing with the lattice-summed operators the filter
    /// can adjusts the displacement to make sure that we end up in
    /// the simulation cell.
}  // namespace madness
#endif // MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
