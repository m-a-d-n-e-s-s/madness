//
// Created by Florian Bischoff on 6/27/22.
//

#ifndef MADNESS_CCPAIRFUNCTION_H
#define MADNESS_CCPAIRFUNCTION_H



#include <madness/mra/mra.h>
#include <chem/commandlineparser.h>
#include <chem/QCCalculationParametersBase.h>
#include <algorithm>
#include <iomanip>
#include <madness/mra/macrotaskq.h>

namespace madness {

class CCConvolutionOperator;

/// FuncTypes used by the CC_function_6d structure
/// Types of Functions used by CC_function class
enum FuncType { UNDEFINED, HOLE, PARTICLE, MIXED, RESPONSE };

/// FuncTypes used by the CC_function_6d structure
enum PairFormat { PT_UNDEFINED, PT_FULL, PT_DECOMPOSED, PT_OP_DECOMPOSED };



/// structure for a CC Function 3D which holds an index and a type
// the type is defined by the enum FuncType (definition at the start of this file)
struct CCFunction {
    CCFunction() : current_error(99), i(99), type(UNDEFINED) {};

    CCFunction(const real_function_3d& f) : current_error(99), function(f), i(99), type(UNDEFINED) {};

    CCFunction(const real_function_3d& f, const size_t& ii) : current_error(99), function(f), i(ii), type(UNDEFINED) {};

    CCFunction(const real_function_3d& f, const size_t& ii, const FuncType& type_) : current_error(99), function(f),
                                                                                     i(ii), type(type_) {};

    CCFunction(const CCFunction& other) : current_error(other.current_error), function(other.function), i(other.i),
                                          type(other.type) {};
    double current_error;
    real_function_3d function;

    real_function_3d get() const { return function; }

    real_function_3d f() const { return function; }

    void set(const real_function_3d& other) { function = other; }

    size_t i;
    FuncType type;

    void info(World& world, const std::string& msg = " ") const;

    std::string name() const;

    double inner(const CCFunction& f) const {
        return inner(f.function);
    }

    double inner(const real_function_3d& f) const {
        return function.inner(f);
    }

    /// scalar multiplication
    CCFunction operator*(const double& fac) const {
        real_function_3d fnew = fac * function;
        return CCFunction(fnew, i, type);
    }

    // for convenience
    bool operator==(const CCFunction& other) const {
        if (i == other.i and type == other.type) return true;
        else return false;
    }

    /// plotting
    void plot(const std::string& msg = "") const {
        plot_plane(function.world(), function, msg + name());
    }
};




/// Helper structure for the coupling potential of CC Singles and Doubles
/// because of the regularization of the CC-Wavefunction (for CC2: |tauij> = |uij> + Qt12*f12*|titj>)
/// we have 6D-functions in std format |u> : type==pure_
/// we have 6D-functions in sepparated format: type==decomposed_ (e.g O1*f12*|titj> = |xy> with x=|k> and y=<k|f12|ti>*|tj>)
/// we have 6D-function like f12|xy> which are not needed to be represented on the 6D MRA-Grid, type==op_decomposed_
    struct CCPairFunction {

    public:
        CCPairFunction(World& world, const real_function_6d& ket) : world(world), type(PT_FULL), a(), b(), op(0),
                                                                    u(ket) {}

        CCPairFunction(World& world, const vector_real_function_3d& f1, const vector_real_function_3d& f2) : world(
                world),
                                                                                                             type(PT_DECOMPOSED),
                                                                                                             a(f1),
                                                                                                             b(f2),
                                                                                                             op(0),
                                                                                                             u() {}

        CCPairFunction(World& world, const std::pair<vector_real_function_3d, vector_real_function_3d>& f) : world(
                world),
                                                                                                             type(PT_DECOMPOSED),
                                                                                                             a(f.first),
                                                                                                             b(f.second),
                                                                                                             op(0),
                                                                                                             u() {}

        CCPairFunction(World& world, const CCConvolutionOperator *op_, const CCFunction& f1, const CCFunction& f2)
                : world(
                world), type(PT_OP_DECOMPOSED), a(), b(), op(op_), x(f1), y(f2), u() {}

        CCPairFunction(const CCPairFunction& other) : world(other.world), type(other.type), a(other.a), b(other.b),
                                                      op(other.op), x(other.x), y(other.y), u(other.u) {}

        CCPairFunction
        operator=(const CCPairFunction& other);

        void info() const {
            if (world.rank() == 0) std::cout << "Information about Pair " << name() << "\n";
            print_size();
        }

        /// deep copy
        CCPairFunction
        copy() const;


        /// make a deep copy and invert the sign
        /// deep copy necessary otherwise: shallow copy errors
        CCPairFunction
        invert_sign();

        CCPairFunction operator*(const double fac) const {
            if (type == PT_FULL) return CCPairFunction(world, fac * u);
            else if (type == PT_DECOMPOSED) return CCPairFunction(world, fac * a, b);
            else if (type == PT_OP_DECOMPOSED) return CCPairFunction(world, op, x * fac, y);
            else MADNESS_EXCEPTION("wrong type in CCPairFunction scale", 1);
        }

        /// print the size of the functions
        void print_size() const;

        std::string name() const;

        /// @param[in] f: a 3D-CC_function
        /// @param[in] particle: the particle on which the operation acts
        /// @param[out] <f|u>_particle (projection from 6D to 3D)
        real_function_3d project_out(const CCFunction& f, const size_t particle) const;

        // result is: <x|op12|f>_particle
        /// @param[in] x: a 3D-CC_function
        /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
        /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
        /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
        real_function_3d
        dirac_convolution(const CCFunction& x, const CCConvolutionOperator& op, const size_t particle) const;

        /// @param[out] particles are interchanged, if the function was u(1,2) the result is u(2,1)
        CCPairFunction swap_particles() const;

        double
        make_xy_u(const CCFunction& xx, const CCFunction& yy) const;

    public:
        /// the 3 types of 6D-function that occur in the CC potential which coupled doubles to singles

        World& world;
        /// the type of the given 6D-function
        const PairFormat type;
        /// if type==decomposed this is the first particle
        vector_real_function_3d a;
        /// if type==decomposed this is the second particle
        vector_real_function_3d b;
        /// if type==op_decomposed_ this is the symmetric 6D-operator (g12 or f12) in u=op12|xy>
        const CCConvolutionOperator *op;
        /// if type==op_decomposed_ this is the first particle in u=op12|xy>
        CCFunction x;
        /// if type==op_decomposed_ this is the second particle in u=op12|xy>
        CCFunction y;
        /// if type=pure_ this is just the MRA 6D-function
        real_function_6d u;

        /// @param[in] f: a 3D-CC_function
        /// @param[in] particle: the particle on which the operation acts
        /// @param[out] <f|u>_particle (projection from 6D to 3D) for the case that u=|ab> so <f|u>_particle = <f|a>*|b> if particle==1
        real_function_3d project_out_decomposed(const real_function_3d& f, const size_t particle) const;

        /// @param[in] f: a 3D-CC_function
        /// @param[in] particle: the particle on which the operation acts
        /// @param[out] <f|u>_particle (projection from 6D to 3D) for the case that u=op|xy> so <f|u>_particle = <f|op|x>*|y> if particle==1
        real_function_3d project_out_op_decomposed(const CCFunction& f, const size_t particle) const;

        /// @param[in] x: a 3D-CC_function
        /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
        /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
        /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
        /// in this case u=|ab> and the result is <x|op|u>_1 = <x|op|a>*|b> for particle==1
        real_function_3d
        dirac_convolution_decomposed(const CCFunction& x, const CCConvolutionOperator& op, const size_t particle) const;

        /// small helper function that gives back (a,b) or (b,a) depending on the value of particle
        const std::pair<vector_real_function_3d, vector_real_function_3d> assign_particles(const size_t particle) const;

        /// swap particle function if type==pure_
        CCPairFunction swap_particles_pure() const;

        /// swap particle function if type==decomposed_
        CCPairFunction swap_particles_decomposed() const;

        /// swap particle function if type==op_decomposed_ (all ops are assumed to be symmetric)
        CCPairFunction swap_particles_op_decomposed() const;
    };
} // namespace madness

#endif //MADNESS_CCPAIRFUNCTION_H
