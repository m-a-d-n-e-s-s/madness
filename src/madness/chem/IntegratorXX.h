//
// Created by Florian Bischoff on 4/11/24.
//

#ifndef INTEGRATORXX_H
#define INTEGRATORXX_H

#include<madness/world/world.h> // for madness::print
#include<madness/world/vector.h>

#ifdef MADNESS_HAS_INTEGRATORXX
#ifndef INTEGRATORXX_INCLUDED
#define INTEGRATORXX_INCLUDED
#include <integratorxx/generators/impl/impl.hpp>
#endif
#endif


/// Builder for a molecular integration grid, using IntegratorXX
///
/// usage:
///     GridBuilder builder;
///     builder.set_nradial(nrad);
///     builder.set_origin(gridcenter);
///     builder.set_angular_order(order);
///     builder.make_grid();
///     std::vector<Vector<double,3>> points=builder.get_points();
/// see test_IntegratorXX.cc for an example
class GridBuilder {
private:
    typedef std::array<double,3> pointT;
    bool debug=false;
    std::vector<pointT> points;         // final grid points (after make_grid())
    madness::Vector<double,3> origin={0,0,0}; // origin/center of the grid
    std::size_t nradial=25; // the number of radial quadrature points
    std::size_t nangular=0; // to be set from angular order, depending on angular_scheme
    std::size_t angular_order=7; // a quadrature of a give order will integrate a spherical harmonic of that order exactly
    std::string radial_scheme="TA";
    std::string angular_scheme="LebedevLaikov";
#ifdef MADNESS_HAS_INTEGRATORXX
    IntegratorXX::SphericalGridFactory::spherical_grid_ptr grid;
#else
    struct dummygrid {
        std::vector<pointT> points() const {return std::vector<pointT>();};
        std::vector<double> weights() const {return std::vector<double>();}
    };
    std::shared_ptr<dummygrid> grid;
#endif

public:
    GridBuilder() {
#ifndef MADNESS_HAS_INTEGRATORXX
        madness::print("\nno GridBuilder without IntegratorXX\n");
        MADNESS_EXCEPTION("no GridBuilder without IntegratorXX",1);
#endif
    }

    /// a quadrature of a give order will integrate a spherical harmonic of that order exactly
    void set_angular_order(const std::size_t order) {
        angular_order=order;
    }

    /// a quadrature of a give order will integrate a spherical harmonic of that order exactly
    std::size_t get_angular_order() const {
        return angular_order;
    }

    /// number of radial gridpoints on the interval [0,\inf)
    void set_nradial(const std::size_t nradial1) {
        nradial=nradial1;
    }

    /// number of radial gridpoints on the interval [0,\inf)
    std::size_t get_nradial() const {
        return nradial;
    }

    /// the origin/center of the grid
    void set_origin(const madness::Vector<double,3>& o) {
        origin=o;
    }

    void make_grid() {
        set_angular_points_from_order(angular_order);
        if (debug) {
            madness::print("generating grid with nradial=",nradial,"nangular=",nangular,
                "angular_order=",angular_order,"npoints=",nradial*nangular);
        }

#ifdef MADNESS_HAS_INTEGRATORXX
        auto radial_quadrature=IntegratorXX::radial_from_string("MuraKnowles");
        auto radial_traits=IntegratorXX::make_radial_traits(radial_quadrature,nradial,1.0);
        auto angular_quadrature=IntegratorXX::angular_from_string("LebedevLaikov");
        IntegratorXX::UnprunedSphericalGridSpecification unp{
            radial_quadrature, *radial_traits,
            angular_quadrature, nangular
          };

        auto pruning_spec = IntegratorXX::robust_psi4_pruning_scheme(unp);
        grid=IntegratorXX::SphericalGridFactory::generate_grid(pruning_spec);
#endif
    }

    /// get the grid points
    std::vector<madness::Vector<double,3>> get_points() const {
        MADNESS_CHECK_THROW(grid.get(),"grid not initialized");
        auto points0=grid->points();
        std::vector<madness::Vector<double,3>> points;
        for (auto p : points0) points.push_back(madness::Vector<double,3>(p)+origin);
        return points;
    }

    std::vector<double> get_weights() const {
        MADNESS_CHECK_THROW(grid.get(),"grid not initialized");
        return grid->weights();
    }

private:
    /// number of angular gridpoints, derived from the angular order
    void set_angular_points_from_order(const std::size_t order) {
#ifdef MADNESS_HAS_INTEGRATORXX
        // this is stupid: why have a runtime switch for the angular scheme if the traits are templated???
        if (angular_scheme=="LebedevLaikov") {
            using traits=IntegratorXX::quadrature_traits< IntegratorXX::LebedevLaikov<double>>;
            int nextorder=traits::next_algebraic_order(order);
            nangular = traits::npts_by_algebraic_order(nextorder);
        } else if (angular_scheme=="AhrensBeylkin") {
            using traits=IntegratorXX::quadrature_traits< IntegratorXX::AhrensBeylkin<double>>;
            int nextorder=traits::next_algebraic_order(order);
            nangular = traits::npts_by_algebraic_order(nextorder);
        } else {
            MADNESS_EXCEPTION("unknown angular quadrature scheme",1);
        }
#endif
    }
};



#endif //INTEGRATORXX_H
