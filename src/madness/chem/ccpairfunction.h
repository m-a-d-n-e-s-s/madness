//
// Created by Florian Bischoff on 6/27/22.
//

#ifndef MADNESS_CCPAIRFUNCTION_H
#define MADNESS_CCPAIRFUNCTION_H



#include <madness/mra/mra.h>
#include<madness/chem/commandlineparser.h>
#include<madness/chem/QCCalculationParametersBase.h>
#include <algorithm>
#include <iomanip>
#include <madness/mra/macrotaskq.h>

namespace madness {

class CCConvolutionOperator;
class ProjectorBase;

/// FuncTypes used by the CC_function_6d structure
/// Types of Functions used by CC_function class
enum FuncType { UNDEFINED, HOLE, PARTICLE, MIXED, RESPONSE };

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
    virtual bool has_operator() const = 0;
    virtual void set_operator(const std::shared_ptr<CCConvolutionOperator> op) = 0;
    virtual const std::shared_ptr<CCConvolutionOperator> get_operator_ptr() const = 0;
    virtual void print_size() const = 0;
    virtual std::string name(const bool transpose=false) const = 0;
    virtual World& world() const =0;
    virtual std::shared_ptr<TwoBodyFunctionComponentBase> clone() = 0;
};

/// a two-body, explicitly 6-dimensional function
template<typename T>
class TwoBodyFunctionPureComponent : public TwoBodyFunctionComponentBase {

public:
    TwoBodyFunctionPureComponent() = default;

    explicit TwoBodyFunctionPureComponent(const Function<T,6>& f) : u(f) {}
    explicit TwoBodyFunctionPureComponent(const std::shared_ptr<CCConvolutionOperator> op, const Function<T,6>& f)
            : u(f), op(op) {}

    /// deep copy
    std::shared_ptr<TwoBodyFunctionComponentBase> clone() override {
        TwoBodyFunctionPureComponent<T> result(op,madness::copy(u));
        return std::make_shared<TwoBodyFunctionPureComponent<T>>(result);
    }

    template<typename Q>
    TwoBodyFunctionPureComponent& operator*=(const Q fac) {
        u.scale(fac);
        return *this;
    }

    bool is_pure() const override {return true;}

    bool has_operator() const override {return op!=nullptr;}

    World& world() const override {return u.world();};

    void serialize() {}

    void print_size() const override {
        u.print_size(name(false));
    }

    std::string name(const bool transpose) const override {
        if (transpose) {
            if (has_operator()) return "< u |"+get_operator_ptr()->name();
            return "< u |";
        }
        if (has_operator()) return get_operator_ptr()->name() + "| u >";
        return "| u >";
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

    const std::shared_ptr<CCConvolutionOperator> get_operator_ptr() const override {return op;};

    void set_operator(const std::shared_ptr<CCConvolutionOperator> op1) override {op=op1;}

    real_function_6d& get_function() {
        return u;
    }

private:
    /// pure 6D function
    real_function_6d u;
    std::shared_ptr<CCConvolutionOperator> op;

};

/// holds two vectors a and b of low-dimensional functions forming a high-dim function by a sum of outer products
/// f(1,2) = sum_i |a_i b_i >
template<typename T>
class TwoBodyFunctionSeparatedComponent : public TwoBodyFunctionComponentBase {

public:
    TwoBodyFunctionSeparatedComponent() = default;

    TwoBodyFunctionSeparatedComponent(const std::vector<Function<T,3>>& a,
                                      const std::vector<Function<T,3>>& b) : a(a), b(b), op(nullptr) {};

    TwoBodyFunctionSeparatedComponent(const std::vector<Function<T,3>>& a,
                                      const std::vector<Function<T,3>>& b,
                                      const std::shared_ptr<CCConvolutionOperator> op) : a(a), b(b), op(op) {};

    TwoBodyFunctionSeparatedComponent(const TwoBodyFunctionSeparatedComponent& other) = default;

    /// deep copy
    std::shared_ptr<TwoBodyFunctionComponentBase> clone() override {
        TwoBodyFunctionSeparatedComponent<T> result(madness::copy(world(),a),madness::copy(world(),b),op);
        return std::make_shared<TwoBodyFunctionSeparatedComponent<T>>(result);
    }

    template<typename Q>
    TwoBodyFunctionSeparatedComponent& operator*=(const Q fac) {
        if (a.size()>0 and a.front().is_initialized()) scale(a.front().world(),a,fac);
        return *this;
    }

    bool is_decomposed() const override {return true;}
    bool has_operator() const override {return op!=nullptr;}

    World& world() const override {
        MADNESS_ASSERT(a.size()>0 and a.front().is_initialized());
        return a.front().world();
    };

    void print_size() const override {
        if (a.size() > 0) {
            World& world = a.front().world();
            madness::print_size(world, a, "a from " + name(false));
            madness::print_size(world, b, "b from " + name(false));
        }
    }

    std::string name(const bool transpose) const override {
        if (transpose) {
            if (has_operator()) return "<ab |"+get_operator_ptr()->name();
            return "<ab |";
        }
        if (has_operator()) return get_operator_ptr()->name() + "| ab>";
        return "| ab>";
    };

    void serialize() {}

    template<typename Q, std::size_t MDIM>
    TwoBodyFunctionPureComponent<T> apply(const SeparatedConvolution<Q,MDIM>* op, const int particle=0) {
        MADNESS_EXCEPTION("TwoBodyFunctionPureComponent<T> apply not yet implemented",1);
    }

    /// return f(2,1)
    void swap_particles_inplace() override {
        std::swap(a,b);
    }

    long rank() const {
        MADNESS_CHECK(a.size()==b.size());
        return a.size();
    }

    std::vector<Function<T,3>> get_a() const {return a;}
    std::vector<Function<T,3>> get_b() const {return b;}
    std::vector<Function<T,3>> get_vector(const int i) const {
        MADNESS_CHECK(i==0 or i==1);
        if (i==0) return a;
        else if (i==1) return b;
        else {
            MADNESS_EXCEPTION("confused index in TwoBodyFunctionSeparatedComponent",1);
        }
    }

    const std::shared_ptr<CCConvolutionOperator> get_operator_ptr() const override {return op;};

    void set_operator(const std::shared_ptr<CCConvolutionOperator> op1) override {op=op1;}

private:

    std::vector<Function<T,3>> a;
    std::vector<Function<T,3>> b;
    std::shared_ptr<CCConvolutionOperator> op;

};




/// Helper structure for the coupling potential of CC Singles and Doubles
/// because of the regularization of the CC-Wavefunction (for CC2: |tauij> = |uij> + Qt12*f12*|titj>)
/// we have 6D-functions in std format |u> : type==pure_
/// we have 6D-functions in separated format: type==decomposed_ (e.g O1*f12*|titj> = |xy> with x=|k> and y=<k|f12|ti>*|tj>)
/// we have 6D-function like f12|xy> which are not needed to be represented on the 6D MRA-Grid, type==op_decomposed_


/** functionality
 *
 *  - ctor
 *  - assignment
 *  - add
 *  - scalar multiplication
 *  - inner
 *  - inner_partial
 *  - swap_particles
 *  - apply
 *  - apply_partial (i.e. exchange)
 *  - serialize
 *  - callapse_to_pure (excl g!)
 *  - mul_partial
 */

/// a 6D function, either in full or low rank form, possibly including an 2-particle function

/**
 * the function is stored as
 *  - pure: full rank form, 6D
 *  - op_pure: full rank form, 6D with an 2-particle function f(1,2) |u>
 *  - decomposed: sum of two vectors of 3D functions \sum_i |a_i(1) b_i(2)>
 *  - op_decomposed: as above, with an 2-particle function: f(1,2) \sum_i |a_i b_i>
 *
**/
struct CCPairFunction {

using T=double;
using pureT=Function<T,6>;

public:

    /// empty ctor
    CCPairFunction() = default;

    /// takes a deep copy of the argument function
    explicit CCPairFunction(const real_function_6d& ket) {
        component.reset(new TwoBodyFunctionPureComponent<T>(copy(ket)));
    }

    /// takes a deep copy of the argument function
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator> op_, const real_function_6d& ket) {
        component.reset(new TwoBodyFunctionPureComponent<T>(op_,copy(ket)));
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const vector_real_function_3d& f1, const vector_real_function_3d& f2) {
        World& world=f1.front().world();
        component.reset(new TwoBodyFunctionSeparatedComponent<T>(copy(world,f1),copy(world,f2)));
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const real_function_3d& f1, const real_function_3d& f2) :
            CCPairFunction(std::vector<real_function_3d>({f1}),std::vector<real_function_3d>({f2})) {
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::pair<vector_real_function_3d, vector_real_function_3d>& f) :
            CCPairFunction(f.first,f.second) {
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator> op_, const CCFunction& f1, const CCFunction& f2) :
            CCPairFunction(op_,std::vector<real_function_3d>({f1.function}),std::vector<real_function_3d>({f2.function})) {
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator> op_, const std::vector<real_function_3d>& f1,
                            const std::vector<real_function_3d>& f2) {
        World& world=f1.front().world();
        component.reset(new TwoBodyFunctionSeparatedComponent<T>(copy(world,f1),copy(world,f2),op_));
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator> op_, const real_function_3d& f1,
                            const real_function_3d& f2) : CCPairFunction(op_,std::vector<real_function_3d>({f1}),
                                                                         std::vector<real_function_3d>({f2})) {
    };

    /// shallow assignment operator
    CCPairFunction& operator()(const CCPairFunction& other) {
        component=other.component;
        return *this;
    }

    /// copy ctor -- shallow
    CCPairFunction(const CCPairFunction& other) = default;

    /// deep copy
    friend CCPairFunction copy(const CCPairFunction& other) {
        CCPairFunction result;
        result.component=other.component->clone();
        return result;
    }

    void info() const { print_size(); }

    World& world() const {
        MADNESS_CHECK(component);
        return component->world();
    }

    Function<T,6>& get_function() {
        MADNESS_CHECK(component and (component->is_pure()));
        return pure().get_function();
    }

    Function<T,6>& get_function() const {
        MADNESS_CHECK(component and (component->is_pure()));
        return pure().get_function();
    }

    /// make a deep copy and invert the sign
    /// deep copy necessary otherwise: shallow copy errors
    CCPairFunction invert_sign();

    CCPairFunction operator*(const double fac) const {
        CCPairFunction result=copy(*this);
        result*=fac;
        return result;
    }

    friend CCPairFunction operator*(const double fac, const CCPairFunction& f) {
        return fac*f;
    }

    bool has_operator() const {return component->has_operator();}
    bool is_pure() const {return component->is_pure();}
    bool is_op_pure() const {return is_pure() and has_operator();}
    bool is_pure_no_op() const {return is_pure() and (not has_operator());}
    bool is_decomposed() const {return component->is_decomposed();}
    bool is_op_decomposed() const {return component->is_decomposed() and component->has_operator();}
    bool is_decomposed_no_op() const {return component->is_decomposed() and (not component->has_operator());}

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

    std::vector<Function<T,3>> get_vector(const int i) const {
        MADNESS_CHECK(component->is_decomposed());
        return decomposed().get_vector(i);
    }

    const CCConvolutionOperator& get_operator() const {
        MADNESS_CHECK(component and component->has_operator());
        return *(component->get_operator_ptr());
    }

    const std::shared_ptr<CCConvolutionOperator> get_operator_ptr() const {
        MADNESS_CHECK(component);
        return component->get_operator_ptr();
    }

    /// can this be converted to a pure representation (depends on the operator, if present)
    bool is_convertible_to_pure_no_op() const;

    /// convert this into a pure hi-dim function
    void convert_to_pure_no_op_inplace();

    CCPairFunction& operator*=(const double fac) {
        if (component->is_pure()) pure()*=fac;
        if (component->is_decomposed()) decomposed()*=fac;
        return *this;
    }

    /// print the size of the functions
    void print_size() const {
        if (component) component->print_size();
    };

    std::string name(const bool transpose=false) const {
        if (not component) return "empty";
        return component->name(transpose);
    }

    /// multiply CCPairFunction with a 3D function of one of the two particles
    friend CCPairFunction multiply(const CCPairFunction& other, const real_function_3d& f, const std::array<int, 3>& v1);

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D)
    real_function_3d project_out(const CCFunction& f, const size_t particle) const;

    /// result is: <x|op12|f>_particle

    /// @param[in] x: a 3D-CC_function
    /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
    /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
    /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
    real_function_3d
    dirac_convolution(const CCFunction& x, const CCConvolutionOperator& op, const size_t particle) const;

    /// @param[out] particles are interchanged, if the function was u(1,2) the result is u(2,1)
    CCPairFunction swap_particles() const {
        CCPairFunction result=copy(*this);
        result.component->swap_particles_inplace();
        return result;
    };

    double
    make_xy_u(const CCFunction& xx, const CCFunction& yy) const;

    /// compute the inner product of this and other
    double inner_internal(const CCPairFunction& other, const real_function_3d& R2) const;

    friend double inner(const CCPairFunction& a, const CCPairFunction& b, const real_function_3d& R2) {
        return a.inner_internal(b,R2);
    }

    friend double inner(const CCPairFunction& a, const CCPairFunction& b) {
        real_function_3d R2;
        return a.inner_internal(b,R2);
    }

    friend double inner(const std::vector<CCPairFunction>& va, const std::vector<CCPairFunction>& vb) {
        double wall0=cpu_time();
        real_function_3d R2;
        double result=0.0;
        for (auto& a : va) {
            for (auto& b : vb) {
                double tmp=a.inner_internal(b,R2);
                double wall1=cpu_time();
//                print("result from inner",a.name(true),b.name(),tmp,wall1-wall0,"s");
                wall0=wall1;
                result+=tmp;
            }
        }
        return result;
    }


public:
    /// the 3 types of 6D-function that occur in the CC potential which coupled doubles to singles
    std::shared_ptr<TwoBodyFunctionComponentBase> component;

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

    template<typename Q, std::size_t MDIM>
    friend CCPairFunction apply(const SeparatedConvolution<Q,MDIM>& G, const CCPairFunction& argument) {
        if (argument.is_pure()) {
            return CCPairFunction(G(argument.get_function()));
        } else if (argument.is_decomposed_no_op()) {
            real_function_6d result=real_factory_6d(argument.world()).compressed();

            MADNESS_ASSERT(argument.get_a().size() == argument.get_b().size());

            for (size_t k = 0; k < argument.get_a().size(); k++) {
                const real_function_6d tmp = G(argument.get_a()[k], argument.get_b()[k]);
                result += tmp;
            }
            return CCPairFunction(result);
        } else {
            MADNESS_EXCEPTION("unknown type in CCPairFunction::apply",1);
        }
        return CCPairFunction();
    };

    real_function_3d partial_inner(const real_function_3d& f,
                                   const std::array<int, 3>& v1,
                                   const std::array<int, 3>& v2) const;

    CCPairFunction partial_inner(const CCPairFunction& other,
                                   const std::array<int, 3>& v1,
                                   const std::array<int, 3>& v2) const;

};

/// apply the projector on the argument function, potentially yielding a vector of CCPairfunctions as result

/// result can be
///  Q12 f12 |ij> = (1 - O1) (1 - O2) f12 i(1) j(2)
///              = f12 ij - \sum_k k(1) f_ik(2) j(2) - \sum_k k(2) f_ij(1)j(1)
/// which is a pure function and a decomposed function
std::vector<CCPairFunction> apply(const ProjectorBase& P, const std::vector<CCPairFunction>& argument);

/// convenience function
CCPairFunction apply(const ProjectorBase& P, const CCPairFunction& argument);

real_function_3d inner(const CCPairFunction& c, const real_function_3d& f,
                       const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2={0,1,2});

CCPairFunction inner(const CCPairFunction& c1, const CCPairFunction& c2,
                       const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2);

} // namespace madness

#endif //MADNESS_CCPAIRFUNCTION_H
