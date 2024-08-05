//
// Created by Florian Bischoff on 6/27/22.
//

#ifndef MADNESS_CCPAIRFUNCTION_H
#define MADNESS_CCPAIRFUNCTION_H



#include <madness/mra/mra.h>
#include<madness/world/timing_utilities.h>
#include<madness/mra/commandlineparser.h>
#include<madness/mra/QCCalculationParametersBase.h>
#include <algorithm>
#include <iomanip>
#include <madness/mra/macrotaskq.h>

namespace madness {

template<typename T, std::size_t NDIM>
class CCConvolutionOperator;
class ProjectorBase;

/// FuncTypes used by the CC_function_6d structure
/// Types of Functions used by CC_function class
enum FuncType { UNDEFINED, HOLE, PARTICLE, MIXED, RESPONSE };

/// structure for a CC Function 3D which holds an index and a type
// the type is defined by the enum FuncType (definition at the start of this file)
template<typename T=double, std::size_t NDIM=3>
class CCFunction : public archive::ParallelSerializableObject {
public:
    CCFunction() : current_error(99), i(99), type(UNDEFINED) {};

    CCFunction(const Function<T,NDIM>& f) : current_error(99), function(f), i(99), type(UNDEFINED) {};

//    CCFunction(const Function<T,NDIM>& f, const size_t& ii) : current_error(99), function(f), i(ii), type(UNDEFINED) {};
//
    CCFunction(const Function<T,NDIM>& f, const size_t& ii, const FuncType& type_) : current_error(99), function(f),
                                                                                     i(ii), type(type_) {};

    CCFunction(const CCFunction& other) : current_error(other.current_error), function(other.function), i(other.i),
                                          type(other.type) {};
    double current_error;
    Function<T,NDIM> function;

    Function<T,NDIM> get() const { return function; }

    Function<T,NDIM> f() const { return function; }

    void set(const Function<T,NDIM>& other) { function = other; }

    size_t i;
    FuncType type;

    void info(World& world, const std::string& msg = " ") const {
        if (world.rank() == 0) {
            std::cout << "Information about 3D function: " << name() << " " << msg << std::endl;
            std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) << " |f|    : " << function.norm2()
                      << std::endl;
            std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) << " |error|: " << current_error << std::endl;
        }
    };

    std::string name() const {
        if (type == HOLE) {
            return "phi" + stringify(i);
        } else if (type == PARTICLE) {
            return "tau" + stringify(i);
        } else if (type == MIXED) {
            return "t" + stringify(i);
        } else if (type == RESPONSE) {
            return "x" + stringify(i);
        } else {
            return "function" + stringify(i);
        }
    };

    double inner(const CCFunction& f) const {
        return inner(f.function);
    }

    double inner(const Function<T,NDIM>& f) const {
        return function.inner(f);
    }

    /// scalar multiplication
    CCFunction operator*(const double& fac) const {
        Function<T,NDIM> fnew = fac * function;
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

    template<typename Archive>
    void serialize(const Archive& ar) {
        ar & type & i & current_error & function;
    }
};


class TwoBodyFunctionComponentBase {
public:
    virtual void swap_particles_inplace() = 0;
    virtual bool is_pure() const {return false;}
    virtual bool is_decomposed() const {return false;}
    virtual bool has_operator() const = 0;
//    virtual void set_operator(const std::shared_ptr<CCConvolutionOperator> op) = 0;
//    virtual const std::shared_ptr<CCConvolutionOperator> get_operator_ptr() const = 0;
    virtual void print_size() const = 0;
    virtual std::string name(const bool transpose=false) const = 0;
    virtual World& world() const =0;
    virtual std::shared_ptr<TwoBodyFunctionComponentBase> clone() = 0;
    virtual ~TwoBodyFunctionComponentBase() {}
    virtual hashT hash() const = 0;
};

/// a two-body, explicitly 6-dimensional function
template<typename T, std::size_t NDIM>
class TwoBodyFunctionPureComponent : public TwoBodyFunctionComponentBase {
    static constexpr std::size_t LDIM=NDIM/2;
    static_assert(NDIM%2==0,"NDIM must be even");

public:
    TwoBodyFunctionPureComponent() = default;

    explicit TwoBodyFunctionPureComponent(const Function<T,NDIM>& f) : u(f) {}
    explicit TwoBodyFunctionPureComponent(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op, const Function<T,NDIM>& f)
            : u(f), op(op) {}

    /// deep copy
    std::shared_ptr<TwoBodyFunctionComponentBase> clone() override {
        TwoBodyFunctionPureComponent<T,NDIM> result(op,madness::copy(u));
        return std::make_shared<TwoBodyFunctionPureComponent<T,NDIM>>(result);
    }

    template<typename Q>
    TwoBodyFunctionPureComponent& operator*=(const Q fac) {
        u.scale(fac);
        return *this;
    }

    bool is_pure() const override {return true;}

    bool has_operator() const override {return op!=nullptr;}

    World& world() const override {return u.world();};

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

    const std::shared_ptr<CCConvolutionOperator<T,LDIM>> get_operator_ptr() const {return op;};

    void set_operator(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op1) {op=op1;}

    Function<T,NDIM>& get_function() {
        return u;
    }

    hashT hash() const override {
        hashT h1=hash_value(u.get_impl());
        if (op) hash_combine(h1,hash_value(*op));
        return h1;
    }

private:
    /// pure 6D function
    Function<T,NDIM> u;
    std::shared_ptr<CCConvolutionOperator<T,LDIM>> op;

};

/// holds two vectors a and b of low-dimensional functions forming a high-dim function by a sum of outer products
/// f(1,2) = sum_i |a_i b_i >
template<typename T, std::size_t NDIM>
class TwoBodyFunctionSeparatedComponent : public TwoBodyFunctionComponentBase {
    static constexpr std::size_t LDIM=NDIM/2;
    static_assert(NDIM%2==0,"NDIM must be even");

public:
    TwoBodyFunctionSeparatedComponent() = default;

    TwoBodyFunctionSeparatedComponent(const std::vector<Function<T,LDIM>>& a,
                                      const std::vector<Function<T,LDIM>>& b) : a(a), b(b), op(nullptr) {};

    TwoBodyFunctionSeparatedComponent(const std::vector<Function<T,LDIM>>& a,
                                      const std::vector<Function<T,LDIM>>& b,
                                      const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op) : a(a), b(b), op(op) {};

    TwoBodyFunctionSeparatedComponent(const TwoBodyFunctionSeparatedComponent& other) = default;

    /// deep copy
    std::shared_ptr<TwoBodyFunctionComponentBase> clone() override {
        TwoBodyFunctionSeparatedComponent<T,NDIM> result(madness::copy(world(),a),madness::copy(world(),b),op);
        return std::make_shared<TwoBodyFunctionSeparatedComponent<T,NDIM>>(result);
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

    hashT hash() const override {
        hashT h1=0;
        for (const auto& aa : a) hash_combine(h1,hash_value(aa.get_impl()));
        for (const auto& bb : b) hash_combine(h1,hash_value(bb.get_impl()));
        // print("hashvalue of TwoBodyFunctionSeparatedComponent: ",h1);

        if (op) hash_combine(h1,hash_value(*op));
        return h1;
    }


        template<typename Q, std::size_t MDIM>
    TwoBodyFunctionPureComponent<T,NDIM> apply(const SeparatedConvolution<Q,MDIM>* op, const int particle=0) {
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

    std::vector<Function<T,LDIM>> get_a() const {return a;}
    std::vector<Function<T,LDIM>> get_b() const {return b;}
    std::vector<Function<T,LDIM>> get_vector(const int i) const {
        MADNESS_CHECK(i==0 or i==1);
        if (i==0) return a;
        else if (i==1) return b;
        else {
            MADNESS_EXCEPTION("confused index in TwoBodyFunctionSeparatedComponent",1);
        }
    }

    const std::shared_ptr<CCConvolutionOperator<T,LDIM>> get_operator_ptr() const {return op;};

    void set_operator(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op1) {op=op1;}

private:

    std::vector<Function<T,LDIM>> a;
    std::vector<Function<T,LDIM>> b;
    std::shared_ptr<CCConvolutionOperator<T,LDIM>> op;

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
template<typename T=double, std::size_t NDIM=6>
class CCPairFunction : public archive::ParallelSerializableObject {
public:
    static constexpr std::size_t LDIM=NDIM/2;
    static_assert(NDIM%2==0,"NDIM must be even");

using pureT=Function<T,NDIM>;

public:

    /// empty ctor
    CCPairFunction() = default;

    /// takes a shallow copy of the argument function
    explicit CCPairFunction(const Function<T,NDIM>& ket) {
        component.reset(new TwoBodyFunctionPureComponent<T,NDIM>(ket));
    }

    /// takes a shallow copy of the argument function
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op_, const Function<T,NDIM>& ket) {
        component.reset(new TwoBodyFunctionPureComponent<T,NDIM>(op_,ket));
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::vector<Function<T,LDIM>>& f1, const std::vector<Function<T,LDIM>>& f2) {
        World& world=f1.front().world();
        component.reset(new TwoBodyFunctionSeparatedComponent<T,NDIM>(copy(world,f1),copy(world,f2)));
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const Function<T,LDIM>& f1, const Function<T,LDIM>& f2) :
            CCPairFunction(std::vector<Function<T,LDIM>>({f1}),std::vector<Function<T,LDIM>>({f2})) {
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::pair<std::vector<Function<T,LDIM>>, std::vector<Function<T,LDIM>>>& f) :
            CCPairFunction(f.first,f.second) {
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op_, const CCFunction<T,LDIM>& f1, const CCFunction<T,LDIM>& f2) :
            CCPairFunction(op_,std::vector<Function<T,LDIM>>({f1.function}),std::vector<Function<T,LDIM>>({f2.function})) {
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op_, const std::vector<Function<T,LDIM>>& f1,
                            const std::vector<Function<T,LDIM>>& f2) {
        World& world=f1.front().world();
        component.reset(new TwoBodyFunctionSeparatedComponent<T,NDIM>(copy(world,f1),copy(world,f2),op_));
    }

    /// takes a deep copy of the argument functions
    explicit CCPairFunction(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op_, const Function<T,LDIM>& f1,
                            const Function<T,LDIM>& f2) : CCPairFunction(op_,std::vector<Function<T,LDIM>>({f1}),
                                                                         std::vector<Function<T,LDIM>>({f2})) {
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

    bool is_assigned() const {
        return component.get();
    }

    friend hashT hash_value(const CCPairFunction& f) {
        if (not f.is_assigned()) { return hashT(); }
        return f.component->hash();
    }

private:
    std::vector<CCPairFunction> consolidate(const std::vector<CCPairFunction>& other,
                                                   const std::vector<std::string>& options,
                                                   const std::vector<Vector<double,LDIM>>& centers) const;

    /// turn decomposed functions with operator into decomposed functions using LowRankFunction
    static std::vector<CCPairFunction> op_dec_to_dec(const std::vector<CCPairFunction>& other,
                                                     const std::vector<Vector<double,LDIM>>& centers);

    /// turn pure functions with operator into pure functions without operators
    static std::vector<CCPairFunction> op_pure_to_pure(const std::vector<CCPairFunction>& other);

    /// turn decomposed functions with operator into pure functions without operators
    static std::vector<CCPairFunction> op_dec_to_pure(const std::vector<CCPairFunction>& other);

    /// remove linear dependent terms in the low-rank parts
    static std::vector<CCPairFunction> remove_linearly_dependent_terms(const std::vector<CCPairFunction>& other,
        double thresh=-1.0);

    static std::vector<CCPairFunction> collect_same_types(const std::vector<CCPairFunction>& other);

public:
    // check if all types (pure. op_pure, decomposed, op_decomposed, with various ops) occur only once
    static bool is_collected(const std::vector<CCPairFunction<T,NDIM>>& other);

    /// collect the terms into a compact format

    /// @param[in] other: a vector of CCPairFunctions
    /// @param[in] options: a vector of strings which can be "one_term", "op_pure_to_pure", "svd"
    /// @param[in] centers: a vector of 3D-vectors which are the centers of the grid for low-rank functions
    /// TODO: implement a function for removing linearly dependent terms without orthonormalization
    friend std::vector<CCPairFunction> consolidate(const std::vector<CCPairFunction>& other,
                                                   const std::vector<std::string> options,
                                                   const std::vector<Vector<double,LDIM>> centers=std::vector<Vector<double,LDIM>>()) {

        if (other.size()>0) return other.front().consolidate(other,options,centers); // workaround
        return other;
    };


    void info() const { print_size(); }

    World& world() const {
        MADNESS_CHECK(component);
        return component->world();
    }

    Function<T,NDIM>& get_function() {
        MADNESS_CHECK(component and (component->is_pure()));
        return pure().get_function();
    }

    Function<T,NDIM>& get_function() const {
        MADNESS_CHECK(component and (component->is_pure()));
        return pure().get_function();
    }

    /// make a deep copy and invert the sign
    /// deep copy necessary otherwise: shallow copy errors
    CCPairFunction invert_sign();

    /// scalar multiplication: f*fac
    CCPairFunction operator*(const double fac) const {
        CCPairFunction result=copy(*this);
        result*=fac;
        return result;
    }

    /// scalar multiplication: fac*f
    friend CCPairFunction operator*(const double fac, const CCPairFunction& f) {
        return f*fac;
    }

    /// multiplication with a 2-particle function
    friend CCPairFunction operator*(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op, const CCPairFunction& f) {
        CCPairFunction result=copy(f);
        return result.multiply_with_op_inplace(op);
    }

    /// multiplication with a 2-particle function
    friend std::vector<CCPairFunction> operator*(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op,
            const std::vector<CCPairFunction>& f) {
        std::vector<CCPairFunction> result;
        for (auto& ff : f) {
            result.push_back(copy(ff));
            result.back().multiply_with_op_inplace(op);
        }
        return result;
    }

    friend std::vector<CCPairFunction> multiply(const std::vector<CCPairFunction>& other, const Function<T,LDIM> f,
                                                const std::array<int, LDIM>& v1) {
        std::vector<CCPairFunction> result;
        for (auto& o : other) {
//            double cpu0=cpu_time();
//            std::cout << "multiply " << o.name();
            result.push_back(multiply(o,f,v1));
//            double cpu1=cpu_time();
//            std::cout << " done after " << cpu1-cpu0 << std::endl;
        }
        return result;
    }

    /// multiplication with a 2-particle function
    CCPairFunction operator*(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op) {
        CCPairFunction result=copy(*this);
        return result.multiply_with_op_inplace(op);
    }

    CCPairFunction<T,NDIM>& multiply_with_op_inplace(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op);


    bool has_operator() const {return component->has_operator();}
    bool is_pure() const {return component->is_pure();}
    bool is_op_pure() const {return is_pure() and has_operator();}
    bool is_pure_no_op() const {return is_pure() and (not has_operator());}
    bool is_decomposed() const {return component->is_decomposed();}
    bool is_op_decomposed() const {return component->is_decomposed() and component->has_operator();}
    bool is_decomposed_no_op() const {return component->is_decomposed() and (not component->has_operator());}

    TwoBodyFunctionPureComponent<T,NDIM>& pure() const {
        if (auto ptr=dynamic_cast<TwoBodyFunctionPureComponent<T,NDIM>*>(component.get())) return *ptr;
        MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
    }

    TwoBodyFunctionSeparatedComponent<T,NDIM>& decomposed() const {
        if (auto ptr=dynamic_cast<TwoBodyFunctionSeparatedComponent<T,NDIM>*>(component.get())) return *ptr;
        MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
    }

    std::vector<Function<T,LDIM>> get_a() const {
        MADNESS_CHECK(component->is_decomposed());
        return decomposed().get_a();
    }

    std::vector<Function<T,LDIM>> get_b() const {
        MADNESS_CHECK(component->is_decomposed());
        return decomposed().get_b();
    }

    std::vector<Function<T,LDIM>> get_vector(const int i) const {
        MADNESS_CHECK(component->is_decomposed());
        return decomposed().get_vector(i);
    }

    const CCConvolutionOperator<T,LDIM>& get_operator() const {
        MADNESS_CHECK(component and component->has_operator());
        if (is_pure()) return *(pure().get_operator_ptr());
        if (is_decomposed()) return *(decomposed().get_operator_ptr());
        MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
        return *(decomposed().get_operator_ptr());
    }

    const std::shared_ptr<CCConvolutionOperator<T,LDIM>> get_operator_ptr() const {
        MADNESS_CHECK(component);
        if (is_pure()) return (pure().get_operator_ptr());
        if (is_decomposed()) return (decomposed().get_operator_ptr());
        MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
//        return component->get_operator_ptr();
    }

    void reset_operator(const std::shared_ptr<CCConvolutionOperator<T,LDIM>> op) {
        MADNESS_CHECK(component);
        if (is_pure()) pure().set_operator(op);
        else if (is_decomposed()) decomposed().set_operator(op);
        else {
            MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
        }
    }

    /// can this be converted to a pure representation (depends on the operator, if present)
    bool is_convertible_to_pure_no_op() const;

    /// out-of-place conversion to pure function
    CCPairFunction to_pure() const {
        auto tmp=copy(*this);
        MADNESS_CHECK(tmp.is_convertible_to_pure_no_op());
        tmp.convert_to_pure_no_op_inplace();
        return tmp;
    }

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

    typename Tensor<T>::scalar_type norm2() const {
        if (component->is_pure()) return pure().get_function().norm2();
        if (component->is_decomposed()) {
            Function<T,LDIM> R2;
            auto tmp= inner_internal(*this,R2);
            typename Tensor<T>::scalar_type result=std::real(tmp);
            typename Tensor<T>::scalar_type imag=std::imag(tmp);
            if ((imag>1.e-14) or (result<-1.e-14)) {
                MADNESS_EXCEPTION("bad norm in TwoBodyFunction",1);
            }
            return sqrt(std::abs(result));
        }
        MADNESS_EXCEPTION("bad cast in TwoBodyFunction",1);
        return 0.0;
    }

    /// multiply CCPairFunction with a 3D function of one of the two particles
    friend CCPairFunction<T,NDIM> multiply(const CCPairFunction<T,NDIM>& other, const Function<T,LDIM>& f,
                                           const std::array<int, LDIM>& v1) {
        auto a012=std::array<int,LDIM>();
        auto a345=std::array<int,LDIM>();
        for (int i=0; i<LDIM; ++i) {
            a012[i]=i;
            a345[i]=i+LDIM;
        }
        int particle=-1;
        if (v1== a012) particle=0;
        if (v1== a345) particle=1;
        MADNESS_CHECK(particle==0 or particle==1);
        World& world=other.world();

        if (other.is_decomposed()) {
            if (particle == 0) {
                return CCPairFunction<T,NDIM>(other.get_operator_ptr(), f * other.get_a(), copy(world, other.get_b()));
            } else {
                return CCPairFunction<T,NDIM>(other.get_operator_ptr(), copy(world, other.get_a()), f * other.get_b());
            }
        } else if (other.is_pure()) {
            auto tmp=multiply(other.get_function(),f,particle+1);
            return CCPairFunction<T,NDIM>(other.get_operator_ptr(),tmp);
        } else  {
            MADNESS_EXCEPTION("confused CCPairFunction<T,NDIM> in multiply",1);
        }
    };

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D)
    Function<T,LDIM> project_out(const CCFunction<T,LDIM>& f, const size_t particle) const;

    /// result is: <x|op12|f>_particle

    /// @param[in] x: a 3D-CC_function
    /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
    /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
    /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
    Function<T,LDIM>
    dirac_convolution(const CCFunction<T,LDIM>& x, const CCConvolutionOperator<T,LDIM>& op, const size_t particle) const;

    /// @param[out] particles are interchanged, if the function was u(1,2) the result is u(2,1)
    CCPairFunction swap_particles() const {
        CCPairFunction result=copy(*this);
        result.component->swap_particles_inplace();
        return result;
    };

    double
    make_xy_u(const CCFunction<T,LDIM>& xx, const CCFunction<T,LDIM>& yy) const;

    /// compute the inner product of this and other
    double inner_internal(const CCPairFunction& other, const Function<T,LDIM>& R2) const;

    friend double inner(const CCPairFunction& a, const CCPairFunction& b, const Function<T,LDIM>& R2) {
        return a.inner_internal(b,R2);
    }

    friend double inner(const CCPairFunction& a, const CCPairFunction& b) {
        Function<T,LDIM> R2;
        return a.inner_internal(b,R2);
    }

    friend double inner(const std::vector<CCPairFunction>& va, const std::vector<CCPairFunction>& vb,
                        const Function<T,LDIM> R2=Function<T,LDIM>()) {
        double wall0=cpu_time();
//        Function<T,LDIM> R2;
        double result=0.0;
        for (auto& a : va) {
            for (auto& b : vb) {
                double tmp=a.inner_internal(b,R2);
                double wall1=cpu_time();
                constexpr std::size_t bufsize=256;
                char buf[bufsize];
                snprintf(buf,bufsize,"result from inner %10s %10s %12.8f %4.1fs",a.name(true).c_str(),b.name().c_str(),tmp,wall1-wall0);
                print(std::string(buf));
                wall0=wall1;
                result+=tmp;
            }
        }
        return result;
    }


    friend std::vector<CCPairFunction> swap_particles(const std::vector<CCPairFunction>& argument) {
        std::vector<CCPairFunction> result;
        for (auto& a : argument) result.push_back(a.swap_particles());
        return result;
    };

public:
    /// the 3 types of 6D-function that occur in the CC potential which coupled doubles to singles
    std::shared_ptr<TwoBodyFunctionComponentBase> component;

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D) for the case that u=|ab> so <f|u>_particle = <f|a>*|b> if particle==1
    Function<T,LDIM> project_out_decomposed(const Function<T,LDIM>& f, const size_t particle) const;

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D) for the case that u=op|xy> so <f|u>_particle = <f|op|x>*|y> if particle==1
    Function<T,LDIM> project_out_op_decomposed(const CCFunction<T,LDIM>& f, const size_t particle) const;

    /// @param[in] x: a 3D-CC_function
    /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
    /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
    /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
    /// in this case u=|ab> and the result is <x|op|u>_1 = <x|op|a>*|b> for particle==1
    Function<T,LDIM>
    dirac_convolution_decomposed(const CCFunction<T,LDIM>& x, const CCConvolutionOperator<T,LDIM>& op, const size_t particle) const;

    /// small helper function that gives back (a,b) or (b,a) depending on the value of particle
    const std::pair<std::vector<Function<T,LDIM>>, std::vector<Function<T,LDIM>>> assign_particles(const size_t particle) const;

    static std::vector<CCPairFunction<T,NDIM>> apply(const ProjectorBase& P, const std::vector<CCPairFunction<T,NDIM>>& argument);

    /// apply the operator on a CCPairfunction, both with the same dimension

    /// note there is another function, where the operator works only on some dimensions of the CCPairFunction!
    /// @return result(x) = \int op(x,x') arg(x') dx': a CCPairfunction with the same dimension as the argument
    friend CCPairFunction<T,NDIM> apply(const SeparatedConvolution<T,NDIM>& G, const CCPairFunction<T,NDIM>& argument) {
        CCPairFunction result;
        timer t1(argument.world());
        if (argument.is_pure()) {
            result=CCPairFunction(G(argument.get_function()));
        } else if (argument.is_decomposed_no_op()) {
            Function<T,NDIM> result1=real_factory_6d(argument.world()).compressed();

            MADNESS_ASSERT(argument.get_a().size() == argument.get_b().size());
            MADNESS_CHECK_THROW(G.particle()==-1,"G must be a two-particle operator in apply(CCPairFunction)");

            for (size_t k = 0; k < argument.get_a().size(); k++) {
                const Function<T,NDIM> tmp = G(argument.get_a()[k], argument.get_b()[k]);
                result1 += tmp;
            }
            result=CCPairFunction(result1);
        } else {
            MADNESS_EXCEPTION("unknown type in CCPairFunction::apply",1);
        }
        t1.end("applying G to " + argument.name());
        return result;
    };


    Function<T,LDIM> partial_inner(const Function<T,LDIM>& f,
                                   const std::array<int, LDIM>& v1,
                                   const std::array<int, LDIM>& v2) const;

    CCPairFunction partial_inner(const CCPairFunction& other,
                                   const std::array<int, LDIM>& v1,
                                   const std::array<int, LDIM>& v2) const;

};

namespace archive {
template <class archiveT, class T, std::size_t NDIM>
struct ArchiveLoadImpl< ParallelInputArchive<archiveT>, CCPairFunction<T,NDIM> > {
    static inline void load(const ParallelInputArchive<archiveT>& ar, CCPairFunction<T,NDIM>& p) {
        constexpr std::size_t LDIM=CCPairFunction<T,NDIM>::LDIM;
        bool exists=false;
        bool is_pure=false;
        bool has_operator=false;
        ar & exists;
        if (exists) {
            ar & is_pure & has_operator;
            if (is_pure) {
                Function<T,NDIM> f;
                ar & f;
                p=CCPairFunction<T,NDIM>(f);
            } else {
                std::size_t sz=0;
                ar & sz;
                std::vector<Function<T,LDIM>> a(sz),b(sz);
                for (auto& aa : a) ar & aa;
                for (auto& bb : b) ar & bb;
                p=CCPairFunction<T,NDIM>(a,b);
            }

            // store construction parameters of the operator, not the operator itself
            if (has_operator) {
                auto param=typename CCConvolutionOperator<T,LDIM>::Parameters();
                OpType type=OT_UNDEFINED;
                ar & param & type;
                auto op=std::make_shared<CCConvolutionOperator<T,LDIM>>(*ar.get_world(),type,param);
                p.reset_operator(op);
            }
        }
    }
};

template <class archiveT, class T, std::size_t NDIM>
struct ArchiveStoreImpl< ParallelOutputArchive<archiveT>, CCPairFunction<T,NDIM> > {
    static inline void store(const ParallelOutputArchive<archiveT>& ar, const CCPairFunction<T,NDIM>& f) {
        bool exists=f.is_assigned();
        ar & exists;
        if (exists) {
            ar & f.is_pure() & f.has_operator();
            if (f.is_pure()) ar & f.get_function();
            if (f.is_decomposed()) {
                auto avec=f.get_a();
                auto bvec=f.get_b();
                ar & avec.size();
                for (const auto& a : avec) ar & a;
                for (const auto& b : bvec) ar & b;
            }
            // store construction parameters of the operator, not the operator itself
            if (f.has_operator()) {
                ar & f.get_operator().parameters & f.get_operator().type();
            }
        }
    }
};
}

/// apply the operator to the argument

/// the operator is applied to one particle only, the other one is left untouched
/// note the ordering of the particles, cf the corresponding comment in mra.h
/// op.particle==1 :  op(f(x,y)) = op(x,x') f(x',y) = result(x,y);
/// op.particle==2 :  op(f(x,y)) = op(y,y') f(x,y') = result(y,x);
template<typename T, std::size_t NDIM>
CCPairFunction<T,NDIM> apply(const SeparatedConvolution<T,NDIM/2>& op, const CCPairFunction<T,NDIM>& arg) {
    bool convert_to_pure=(arg.has_operator() or arg.is_pure());
    CCPairFunction<T, NDIM> result;
    World& world = arg.world();

    if (convert_to_pure) {
        auto tmp=arg.to_pure().get_function();
        tmp=op(tmp);

        result=(CCPairFunction<T,NDIM>(tmp));
        // !! confusing ordering of the result variables!!
        if (op.particle()==2) result=result.swap_particles();

    } else if (arg.is_decomposed_no_op()) {
        MADNESS_CHECK(op.particle()==1 or op.particle()==2);
        if (op.particle()==1) {
            auto tmp= madness::apply(world,op,arg.get_a());
            result=(CCPairFunction<T,NDIM>(tmp,arg.get_b()));
        } else if (op.particle()==2) {
            auto tmp= madness::apply(world,op,arg.get_b());
            result=(CCPairFunction<T,NDIM>(tmp,arg.get_a()));
        }

    } else {
        MADNESS_CHECK_THROW(false,"confused type in apply(CCPairFunction<T,NDIM>)");
    }

    return result;
}

/// apply the operator to the argument

/// the operator is applied to one particle only, the other one is left untouched
/// note the ordering of the particles, cf the corresponding comment in mra.h
/// op.particle==1 :  op(f(x,y)) = op(x,x') f(x',y) = result(x,y);
/// op.particle==2 :  op(f(x,y)) = op(y,y') f(x,y') = result(y,x);
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> apply(const SeparatedConvolution<T,NDIM/2>& op, const std::vector<CCPairFunction<T,NDIM>>& argument) {
    std::vector<CCPairFunction<T, NDIM>> result;
    for (const auto& arg : argument) result.push_back(madness::apply(op, arg));
    return result;
}

template<typename T, std::size_t NDIM>
CCPairFunction<T,NDIM> apply(const ProjectorBase& projector, const CCPairFunction<T,NDIM>& argument) {
    auto result=madness::apply(projector,std::vector<CCPairFunction<T,NDIM>> (1,argument));
    MADNESS_CHECK(result.size()==1);
    return result[0];
}

/// apply the projector on the argument function, potentially yielding a vector of CCPairfunctions as result

/// result can be
///  Q12 f12 |ij> = (1 - O1) (1 - O2) f12 i(1) j(2)
///              = f12 ij - \sum_k k(1) f_ik(2) j(2) - \sum_k k(2) f_ij(1)j(1)
/// which is a pure function and a decomposed function
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> apply(const ProjectorBase& projector, const std::vector<CCPairFunction<T,NDIM>>& argument) {
    return CCPairFunction<T,NDIM>::apply(projector,argument);
}


template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM>inner(const CCPairFunction<T,NDIM>& c, const Function<T,CCPairFunction<T,NDIM>::LDIM>& f,
                                              const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2) {
    constexpr std::size_t LDIM=CCPairFunction<T,NDIM>::LDIM;
    auto v11=std::array<int,LDIM>({std::get<0>(v1),std::get<1>(v1),std::get<2>(v1)});
    auto v22=std::array<int,LDIM>({std::get<0>(v2),std::get<1>(v2),std::get<2>(v2)});

    return c.partial_inner(f,v11,v22);
}

template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM>inner(const CCPairFunction<T,NDIM>& c, const Function<T,CCPairFunction<T,NDIM>::LDIM>& f,
                                              const std::array<int,CCPairFunction<T,NDIM>::LDIM>& v1,
                                              const std::array<int,CCPairFunction<T,NDIM>::LDIM>& v2) {
    return c.partial_inner(f,v1,v2);
}

template<typename T, std::size_t NDIM>
CCPairFunction<T,NDIM> inner(const CCPairFunction<T,NDIM>& c1, const CCPairFunction<T,NDIM>& c2,
                             const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2) {
    constexpr std::size_t LDIM=CCPairFunction<T,NDIM>::LDIM;
    auto v11=std::array<int,LDIM>({std::get<0>(v1),std::get<1>(v1),std::get<2>(v1)});
    auto v22=std::array<int,LDIM>({std::get<0>(v2),std::get<1>(v2),std::get<2>(v2)});

    return c1.partial_inner(c2,v11,v22);
}

template<typename T, std::size_t NDIM>
CCPairFunction<T,NDIM> inner(const CCPairFunction<T,NDIM>& c1, const CCPairFunction<T,NDIM>& c2,
                             const std::array<int,CCPairFunction<T,NDIM>::LDIM>& v1,
                             const std::array<int,CCPairFunction<T,NDIM>::LDIM>& v2) {
    return c1.partial_inner(c2,v1,v2);
}

template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> inner(const std::vector<CCPairFunction<T,NDIM>>& c1,
                                          const std::vector<CCPairFunction<T,NDIM>>& c2,
                                          const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2) {
    constexpr std::size_t LDIM=CCPairFunction<T,NDIM>::LDIM;
    auto v11=std::array<int,LDIM>({std::get<0>(v1),std::get<1>(v1),std::get<2>(v1)});
    auto v22=std::array<int,LDIM>({std::get<0>(v2),std::get<1>(v2),std::get<2>(v2)});
    return inner(c1,c2,v11,v22);
}

template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> inner(const std::vector<CCPairFunction<T,NDIM>>& c1,
                                          const std::vector<CCPairFunction<T,NDIM>>& c2,
                                          const std::array<int,CCPairFunction<T,NDIM>::LDIM>& v1,
                                          const std::array<int,CCPairFunction<T,NDIM>::LDIM>& v2) {
    std::vector<CCPairFunction<T,NDIM>> result;
    for (const auto& cc1 : c1) {
        for (const auto& cc2 : c2) {
            print("inner of ",cc1.name(), cc2.name());
            result.push_back(inner(cc1,cc2,v1,v2));
        }
    }
    return result;
}


template <typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM> >& operator+=(std::vector<CCPairFunction<T,NDIM> >& rhs,
        const std::vector<CCPairFunction<T,NDIM> >& lhs) {
    for (const auto& l : lhs) rhs.push_back(l);
    return rhs;
}

template <typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM> >& operator-=(std::vector<CCPairFunction<T,NDIM> >& rhs,
        const std::vector<CCPairFunction<T,NDIM> >& lhs) {
    for (const auto& l : lhs) rhs.push_back(-1.0*l);
    return rhs;
}

template <typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM> > operator*(const double fac, const std::vector<CCPairFunction<T,NDIM> >& arg) {
    std::vector<CCPairFunction<T,NDIM>> result;
    for (const auto& l : arg) result.push_back(fac*l);
    return result;
}


template<typename T, std::size_t NDIM>
bool is_collected(const std::vector<CCPairFunction<T,NDIM>>& other) {
    return CCPairFunction<T,NDIM>::is_collected(other);

}
} // namespace madness

#endif //MADNESS_CCPAIRFUNCTION_H
