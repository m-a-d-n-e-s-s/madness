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


class TwoBodyFunctionComponentBase {
public:
    virtual void swap_particles_inplace() = 0;
    virtual bool is_pure() const {return false;}
    virtual bool is_decomposed() const {return false;}
    virtual bool has_operator() const {return false;}
    virtual void print_size() const = 0;
    virtual std::string name() const = 0;
};

/// a two-body, explicitly 6-dimensional function
template<typename T>
class TwoBodyFunctionPureComponent : public TwoBodyFunctionComponentBase {

public:

    TwoBodyFunctionPureComponent(const Function<T,6>& f) : u(f) {}

    template<typename Q>
    TwoBodyFunctionPureComponent& operator*=(const Q fac) {
        u.scale(fac);
        return *this;
    }

    bool is_pure() const override {return true;}

    void serialize() {}

    void print_size() const override {
        u.print_size(name());
    }

    std::string name() const override {
        return "|u>";
    }


//    template<typename Q, std::size_t MDIM>
//    TwoBodyFunctionPureComponent operator*()(const Function<Q,MDIM>& g, const int particle=0) {}
//
    template<typename Q, std::size_t MDIM>
    TwoBodyFunctionPureComponent apply(const SeparatedConvolution<Q,MDIM>* op, const int particle=0) {}

    /// return f(2,1)
    void swap_particles_inplace() override {
        u=swap_particles(u);
    }

    real_function_6d& get_function() {
        return u;
    }

private:
    /// pure 6D function
    real_function_6d u;

};

template<typename T>
class TwoBodyFunctionSeparatedComponent : public TwoBodyFunctionComponentBase {

public:
    TwoBodyFunctionSeparatedComponent(const std::vector<Function<T,3>>& a,
                                      const std::vector<Function<T,3>>& b) : a(a), b(b), op(nullptr) {};

    TwoBodyFunctionSeparatedComponent(const std::vector<Function<T,3>>& a,
                                      const std::vector<Function<T,3>>& b,
                                      const CCConvolutionOperator* op) : a(a), b(b), op(op) {};

    template<typename Q>
    TwoBodyFunctionSeparatedComponent& operator*=(const Q fac) {
        if (a.size()>0 and a.front().is_initialized()) scale(a.front().world(),a,fac);
        return *this;
    }

    bool is_decomposed() const override {return true;}
    bool has_operator() const override {return op!=nullptr;}

    void print_size() const override {
        if (a.size() > 0) {
            World& world = a.front().world();
            madness::print_size(world, a, "a from " + name());
            madness::print_size(world, b, "b from " + name());
        }
    }


    std::string name() const override;

    void serialize() {}

    template<typename Q, std::size_t MDIM>
    TwoBodyFunctionPureComponent<T> apply(const SeparatedConvolution<Q,MDIM>* op, const int particle=0) {}

    /// return f(2,1)
    void swap_particles_inplace() override {
        std::swap(a,b);
    }

    vector_real_function_3d get_a() const {return a;}
    vector_real_function_3d get_b() const {return b;}
    CCConvolutionOperator* get_operator_ptr() const {return op;};

private:

    vector_real_function_3d a;
    vector_real_function_3d b;
    const CCConvolutionOperator* op;

};




/// Helper structure for the coupling potential of CC Singles and Doubles
/// because of the regularization of the CC-Wavefunction (for CC2: |tauij> = |uij> + Qt12*f12*|titj>)
/// we have 6D-functions in std format |u> : type==pure_
/// we have 6D-functions in sepparated format: type==decomposed_ (e.g O1*f12*|titj> = |xy> with x=|k> and y=<k|f12|ti>*|tj>)
/// we have 6D-function like f12|xy> which are not needed to be represented on the 6D MRA-Grid, type==op_decomposed_
    struct CCPairFunction {

        using T=double;
    public:
        CCPairFunction(World& world, const real_function_6d& ket) {
            component.reset(new TwoBodyFunctionPureComponent<T>(ket));
        }

        CCPairFunction(World& world, const vector_real_function_3d& f1, const vector_real_function_3d& f2) {
            component.reset(new TwoBodyFunctionSeparatedComponent<T>(f1,f2));
        }

        CCPairFunction(World& world, const std::pair<vector_real_function_3d, vector_real_function_3d>& f) {
            component.reset(new TwoBodyFunctionSeparatedComponent<T>(f.first,f.second));
        }

        CCPairFunction(World& world, const CCConvolutionOperator *op_, const CCFunction& f1, const CCFunction& f2) {
            component.reset(new TwoBodyFunctionSeparatedComponent<T>({f1.function},{f2.function},op_));
        }

        CCPairFunction(const CCPairFunction& other) = default;

        CCPairFunction
        operator=(const CCPairFunction& other);

        void info() const { print_size(); }

        Function<T,6>& get_function() {
            MADNESS_CHECK(component and (component->is_pure()));
            return pure().get_function();
        }

        Function<T,6>& get_function() const {
            MADNESS_CHECK(component and (component->is_pure()));
            return pure().get_function();
        }

        std::vector<Function<T,3>> get_function_vector_a() const {
            MADNESS_CHECK(component and (component->is_decomposed()));
            return decomposed().get_a();
        }

        std::vector<Function<T,3>> get_function_vector_b() const {
            MADNESS_CHECK(component and (component->is_decomposed()));
            return decomposed().get_b();
        }

        /// make a deep copy and invert the sign
        /// deep copy necessary otherwise: shallow copy errors
        CCPairFunction invert_sign();

        CCPairFunction operator*(const double fac) const {
            CCPairFunction result(*this);
            result*=fac;
            return result;
        }

        TwoBodyFunctionPureComponent<T>& pure() const {
            if (auto ptr=dynamic_cast<TwoBodyFunctionPureComponent<T>*>(component.get())) return *ptr;
            MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
        }

        TwoBodyFunctionSeparatedComponent<T>& decomposed() const {
            if (auto ptr=dynamic_cast<TwoBodyFunctionSeparatedComponent<T>*>(component.get())) return *ptr;
            MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
        }

        vector_real_function_3d get_a() const {
            MADNESS_CHECK(component->is_decomposed());
            return decomposed().get_a();
        }
        vector_real_function_3d get_b() const {
            MADNESS_CHECK(component->is_decomposed());
            return decomposed().get_b();
        }

        CCPairFunction& operator*=(const double fac) {
            if (component->is_pure()) pure()*=fac;
            if (component->is_decomposed()) decomposed()*=fac;
            return *this;
        }

        /// print the size of the functions
        void print_size() const {
            if (component) component->print_size();
        };

        std::string name() const {
            if (not component) return "empty";
            return component->name();
        };

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
        CCPairFunction swap_particles() const {
            CCPairFunction result(*this);
            result.component->swap_particles_inplace();
            return result;
        };

        double
        make_xy_u(const CCFunction& xx, const CCFunction& yy) const;

        double inner_internal(const CCPairFunction& other, const real_function_3d& R2) const;

        friend double inner(const CCPairFunction& a, const CCPairFunction& b, const real_function_3d& R2) {
            return a.inner_internal(b,R2);
        }

    public:
        /// the 3 types of 6D-function that occur in the CC potential which coupled doubles to singles
        std::shared_ptr<TwoBodyFunctionComponentBase> component;

//         World& world;
//         /// the type of the given 6D-function
//         const PairFormat type;
//         /// if type==decomposed this is the first particle
//         vector_real_function_3d a;
//         /// if type==decomposed this is the second particle
//         vector_real_function_3d b;
//         /// if type==op_decomposed_ this is the symmetric 6D-operator (g12 or f12) in u=op12|xy>
//         const CCConvolutionOperator *op;
//         /// if type==op_decomposed_ this is the first particle in u=op12|xy>
//         CCFunction x;
//         /// if type==op_decomposed_ this is the second particle in u=op12|xy>
//         CCFunction y;
//         /// if type=pure_ this is just the MRA 6D-function
//         real_function_6d u;

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
    };
} // namespace madness

#endif //MADNESS_CCPAIRFUNCTION_H
