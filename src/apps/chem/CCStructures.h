/*
 * CCStructures.h
 *
 *  Created on: Sep 3, 2015
 *      Author: kottmanj
 */


/// File holds all helper structures necessary for the CC_Operator and CC2 class
#ifndef CCSTRUCTURES_H_
#define CCSTRUCTURES_H_

#include <madness/mra/mra.h>
#include <algorithm>
#include <iomanip>

namespace madness {
/// FuncTypes used by the CC_function_6d structure
enum PairFormat {
    PT_UNDEFINED, PT_FULL, PT_DECOMPOSED, PT_OP_DECOMPOSED
};
/// Operatortypes used by the CCConvolutionOperator Class
enum OpType {
    OT_UNDEFINED, OT_G12, OT_F12
};
/// Calculation Types used by CC2
enum CalcType {
    CT_UNDEFINED, CT_MP2, CT_CC2, CT_LRCCS, CT_LRCC2, CT_CISPD, CT_ADC2, CT_TDHF, CT_TEST
};
/// Types of Functions used by CC_function class
enum FuncType {
    UNDEFINED, HOLE, PARTICLE, MIXED, RESPONSE
};
/// Type of Pairs used by CC_Pair2 class
enum CCState {
    CCSTATE_UNDEFINED, GROUND_STATE, EXCITED_STATE
};
/// CC2 Singles Potentials
enum PotentialType {
    POT_UNDEFINED,
    POT_F3D_,
    POT_s3a_,
    POT_s3b_,
    POT_s3c_,
    POT_s5a_,
    POT_s5b_,
    POT_s5c_,
    POT_s2b_,
    POT_s2c_,
    POT_s4a_,
    POT_s4b_,
    POT_s4c_,
    POT_s6_,
    POT_ccs_,
    POT_cis_,
    POT_singles_
};

/// Assigns strings to enums for formated output
std::string
assign_name(const PairFormat &input);

/// Assigns strings to enums for formated output
std::string
assign_name(const CCState &input);

/// Assigns strings to enums for formated output
std::string
assign_name(const OpType &input);

/// Assigns enum to string
CalcType
assign_calctype(const std::string name);

/// Assigns strings to enums for formated output
std::string
assign_name(const CalcType &inp);

/// Assigns strings to enums for formated output
std::string
assign_name(const PotentialType &inp);

/// Assigns strings to enums for formated output
std::string
assign_name(const FuncType &inp);

// Little structure for formated output and to collect warnings
// much room to improve
struct CCMessenger {
    CCMessenger(World &world) : world(world), output_prec(10), scientific(true), debug(false), os(std::cout) {}

    World &world;
    size_t output_prec;
    bool scientific;
    bool debug;

    void operator()(const std::string &msg) const { output(msg); }

    void debug_output(const std::string &msg) const {
        if (debug) output(msg);
    }

    void
    output(const std::string &msg) const;

    void
    section(const std::string &msg) const;

    void
    subsection(const std::string &msg) const;

    void
    warning(const std::string &msg) const;

    void print_warnings() const {
        for (const auto &x:warnings) if (world.rank() == 0) std::cout << x << "\n";
    }

    template<class T>
    CCMessenger operator<<(const T &t) const {
        using madness::operators::operator<<;
        if (world.rank() == 0) os << t;
        return *this;
    }

    /// collect all warnings that occur to print out at the end of the job
    mutable std::vector<std::string> warnings;
    /// output stream
    std::ostream &os;
};

/// Timer Structure
struct CCTimer {
    /// TDA_TIMER contructor
    /// @param[in] world the world
    /// @param[in] msg	a string that contains the desired printout when info function is called
    CCTimer(World &world, std::string msg) : world(world), start_wall(wall_time()), start_cpu(cpu_time()),
                                             operation(msg), end_wall(0.0), end_cpu(0.0), time_wall(-1.0),
                                             time_cpu(-1.0) {}

    World &world;
    double start_wall;
    double start_cpu;
    std::string operation;
    double end_wall;
    double end_cpu;
    double time_wall;
    double time_cpu;

    void update_time() {
        time_wall = wall_time() - start_wall;
        time_cpu = cpu_time() - start_cpu;
    }

public:
    /// print out information about the passed time since the CC_TIMER object was created
    void
    info(const bool debug = true, const double norm = 12345.6789);

    CCTimer start() {
        start_wall = wall_time();
        start_cpu = cpu_time();
        return *this;
    }

    CCTimer stop() {
        end_wall = wall_time();
        end_cpu = cpu_time();
        time_wall = end_wall - start_wall;
        time_cpu = end_cpu - start_cpu;
        return *this;
    }


    double get_wall_time_diff() const { return end_wall; }

    double get_cpu_time_diff() const { return end_cpu; }

    std::pair<double, double> current_time(bool printout = false) {
        if (time_wall < 0.0 or time_cpu < 0.0) stop();
        return std::make_pair(time_wall, time_cpu);
    }

    double current_wall() { return current_time().first; }

    double current_cpu() { return current_time().second; }

    void print() {
        print(current_time());
    }

    void print() const {
        print(std::make_pair(time_wall, time_cpu));
    }

    void print(const std::pair<double, double> &times) const {
        if (world.rank() == 0) {
            std::cout << std::setfill(' ') << std::scientific << std::setprecision(2)
                      << "Timer: " << times.first << " (Wall), " << times.second << " (CPU)" << ", (" + operation + ")"
                      << "\n";
        }
    }


};

/// Calculation TDHFParameters for CC2 and TDA calculations
/// Maybe merge this with calculation_parameters of SCF at some point, or split into TDA and CC
struct CCParameters {

    CCParameters() {};

    const double uninitialized = 123.456;

    /// copy constructor
    CCParameters(const CCParameters &other);

    /// ctor reading out the input file
    CCParameters(const std::string &input, const double &low);

    // the demanded calculation: possibilities are MP2_, CC2_, CIS_, CCS_ (same as CIS), CISpD_
    CalcType calculation = CT_LRCC2;
    double lo = 1.e-7;
    // the finest length to be resolved by 6D operators which needs special refinement
    // this will define the depth of the special level (default is 1.0 bohr)
    double dmin = 1.0;
    // function thresh 3D
    double thresh_3D = FunctionDefaults<3>::get_thresh();
    double tight_thresh_3D = FunctionDefaults<3>::get_thresh() * 0.1;
    // function thresh 6D
    double thresh_6D = FunctionDefaults<6>::get_thresh();
    double tight_thresh_6D = FunctionDefaults<3>::get_thresh() * 0.1;
    // BSH thresh
    double thresh_bsh_3D = std::min(1.e-4, FunctionDefaults<3>::get_thresh());
    double thresh_bsh_6D = std::min(1.e-4, FunctionDefaults<3>::get_thresh());
    // Poisson thresh
    double thresh_poisson = std::min(1.e-4, FunctionDefaults<3>::get_thresh());
    // f12 thresh
    double thresh_f12 = std::min(1.e-4, FunctionDefaults<3>::get_thresh());
    // Ue thresh
    double thresh_Ue = std::min(1.e-4, FunctionDefaults<3>::get_thresh());
    // Convergence for Correlation Energy (overall and pairs)
    double econv = FunctionDefaults<6>::get_thresh();
    double econv_pairs = FunctionDefaults<6>::get_thresh();
    // Convergence for CC-singles
    double dconv_3D = FunctionDefaults<6>::get_thresh();
    // Convergence for CC-Doubles
    double dconv_6D = FunctionDefaults<6>::get_thresh();;
    // iterations
    size_t iter_max = 10;
    size_t iter_max_3D = 10;
    size_t iter_max_6D = 10;
    // restart
    bool restart = false;
    bool no_compute = false;
    std::pair<int, int> only_pair = std::make_pair(-1, -1);
    bool no_compute_gs = false;
    bool no_compute_response = false;
    bool no_compute_mp2 = false;
    bool no_compute_cc2 = false;
    bool no_compute_cispd = false;
    bool no_compute_lrcc2 = false;
    // Exponent for the correlation factor
    double corrfac_gamma = 1.0;
    // for formated output
    double output_prec = 1.e-8;
    // debug mode
    bool debug = false;
    // make additional plots
    bool plot = false;
    // use kain
    bool kain = true;
    size_t kain_subspace = 5;
    // freeze MOs
    size_t freeze = 0;

    // Gamma of the correlation factor
    double gamma() const {
        if (corrfac_gamma < 0) MADNESS_EXCEPTION("ERROR in CC_PARAMETERS: CORRFAC_GAMMA WAS NOT INITIALIZED", 1);
        return corrfac_gamma;
    }

    bool test = false;
    // choose if Q for the constant part of MP2 and related calculations should be decomposed: GQV or GV - GO12V
    bool decompose_Q = false;
    // if true the ansatz for the CC2 ground state pairs is |tau_ij> = |u_ij> + Qtf12|titj>, with Qt = Q - |tau><phi|
    // if false the ansatz is the same with normal Q projector
    // the response ansatz is the corresponding response of the gs ansatz
    bool QtAnsatz = true;

    /// a vector containing the excitations which shall be optizmized later (with CIS(D) or CC2)
    std::vector<size_t> excitations_;

    // TDHFParameters for the TDA Algorithm

    /// print out the parameters
    void information(World &world) const;

    /// check if parameters are set correct
    void sanity_check(World &world) const;

    void error(World &world, const std::string &msg) const {
        if (world.rank() == 0)
            std::cout << "\n\n\n\n\n!!!!!!!!!\n\nERROR IN CC_PARAMETERS:\n    ERROR MESSAGE IS: " << msg
                      << "\n\n\n!!!!!!!!" << std::endl;
        MADNESS_EXCEPTION("ERROR IN CC_PARAMETERS", 1);
    }

    size_t warning(World &world, const std::string &msg) const {
        if (world.rank() == 0) std::cout << "WARNING IN CC_PARAMETERS!: " << msg << std::endl;
        return 1;
    }
};


/// POD holding all electron pairs with easy access
/// Similar strucutre than the Pair structure from MP2 but with some additional features (merge at some point)
/// This structure will also be used for intermediates
template<typename T>
struct Pairs {


    typedef std::map<std::pair<int, int>, T> pairmapT;
    pairmapT allpairs;


    /// getter
    const T &operator()(int i, int j) const {
        return allpairs.at(std::make_pair(i, j));
    }

    /// getter
    // at instead of [] operator bc [] inserts new element if nothing is found while at throws out of range error
    T &operator()(int i, int j) {
        return allpairs.at(std::make_pair(i, j));
    }

    /// setter
    /// can NOT replace elements (for this construct new pair map and swap the content)
    void insert(int i, int j, const T &pair) {
        std::pair<int, int> key = std::make_pair(i, j);
        allpairs.insert(std::make_pair(key, pair));
    }

    /// swap the contant of the pairmap
    void swap(Pairs<T> &other) {
        allpairs.swap(other.allpairs);
    }

    bool empty() const {
        if (allpairs.size() == 0) return true;
        else return false;
    }
};

/// f12 and g12 intermediates of the form <f1|op|f2> (with op=f12 or op=g12) will be saved using the pair structure
typedef Pairs<real_function_3d> intermediateT;

/// Returns the size of an intermediate
double
size_of(const intermediateT &im);


/// structure for a CC Function 3D which holds an index and a type
// the type is defined by the enum FuncType (definition at the start of this file)
struct CCFunction {
    CCFunction() : current_error(99), i(99), type(UNDEFINED) {};

    CCFunction(const real_function_3d &f) : current_error(99), function(f), i(99), type(UNDEFINED) {};

    CCFunction(const real_function_3d &f, const size_t &ii) : current_error(99), function(f), i(ii), type(UNDEFINED) {};

    CCFunction(const real_function_3d &f, const size_t &ii, const FuncType &type_) : current_error(99), function(f),
                                                                                     i(ii), type(type_) {};

    CCFunction(const CCFunction &other) : current_error(other.current_error), function(other.function), i(other.i),
                                          type(other.type) {};
    double current_error;
    real_function_3d function;

    real_function_3d get() const { return function; }

    real_function_3d f() const { return function; }

    void set(const real_function_3d &other) { function = other; }

    size_t i;
    FuncType type;

    void info(World &world, const std::string &msg = " ") const;

    std::string name() const;

    double inner(const CCFunction &f) const {
        return inner(f.function);
    }

    double inner(const real_function_3d &f) const {
        return function.inner(f);
    }

    /// scalar multiplication
    CCFunction operator*(const double &fac) const {
        real_function_3d fnew = fac * function;
        return CCFunction(fnew, i, type);
    }

    // for convenience
    bool operator==(const CCFunction &other) const {
        if (i == other.i and type == other.type) return true;
        else return false;
    }

    /// plotting
    void plot(const std::string &msg = "") const {
        plot_plane(function.world(), function, msg + name());
    }
};


// structure for CC Vectorfunction
/// A helper structure which holds a map of functions
struct CC_vecfunction {

    CC_vecfunction() : type(UNDEFINED), omega(0.0), excitation(-1), current_error(99.9), delta(0.0) {}

    CC_vecfunction(const FuncType type_) : type(type_), omega(0.0), excitation(-1), current_error(99.9), delta(0.0) {}

    CC_vecfunction(const vector_real_function_3d &v) : type(UNDEFINED), omega(0.0), excitation(-1), current_error(99.9),
                                                       delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction tmp(v[i], i, type);
            functions.insert(std::make_pair(i, tmp));
        }
    }

    CC_vecfunction(const std::vector<CCFunction> &v) : type(UNDEFINED), omega(0.0), excitation(-1), current_error(99.9),
                                                       delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            functions.insert(std::make_pair(v[i].i, v[i]));
        }
    }

    CC_vecfunction(const vector_real_function_3d &v, const FuncType &type) : type(type), omega(0.0), excitation(-1),
                                                                             current_error(99.9), delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction tmp(v[i], i, type);
            functions.insert(std::make_pair(i, tmp));
        }
    }

    CC_vecfunction(const vector_real_function_3d &v, const FuncType &type, const size_t &freeze) : type(type),
                                                                                                   omega(0.0),
                                                                                                   excitation(-1),
                                                                                                   current_error(99.9),
                                                                                                   delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction tmp(v[i], freeze + i, type);
            functions.insert(std::make_pair(freeze + i, tmp));
        }
    }

    CC_vecfunction(const std::vector<CCFunction> &v, const FuncType type_)
            : type(type_), omega(0.0), excitation(-1), current_error(99.9), delta(0.0) {
        for (auto x:v) functions.insert(std::make_pair(x.i, x));
    }

    /// copy ctor (shallow)
    CC_vecfunction(const CC_vecfunction &other)
            : functions(other.functions), type(other.type), omega(other.omega),
              excitation(other.excitation), current_error(other.current_error),
              delta(other.delta), irrep(other.irrep) {
    }

    /// assignment operator
//    CC_vecfunction& operator=(const CC_vecfunction& other) = default;
    CC_vecfunction &operator=(const CC_vecfunction &other) {
        if (this == &other) return *this;
        functions = other.functions;
        type = other.type;
        omega = other.omega;
        excitation = other.excitation;
        current_error = other.current_error;
        delta = other.delta;
        irrep = other.irrep;
        return *this;
    }

    typedef std::map<std::size_t, CCFunction> CC_functionmap;
    CC_functionmap functions;

    /// returns a deep copy (void shallow copy errors)
    CC_vecfunction
    copy() const;


    FuncType type;
    double omega; /// excitation energy
    int excitation; /// the excitation number
    double current_error;
    double delta; // Last difference in Energy
    std::string irrep = "null";    /// excitation irrep (direct product of x function and corresponding orbital)

    std::string
    name() const;

    /// getter
    const CCFunction &operator()(const CCFunction &i) const {
        return functions.find(i.i)->second;
    }

    /// getter
    const CCFunction &operator()(const size_t &i) const {
        return functions.find(i)->second;
    }

    /// getter
    CCFunction &operator()(const CCFunction &i) {
        return functions[i.i];
    }

    /// getter
    CCFunction &operator()(const size_t &i) {
        return functions[i];
    }

    /// setter
    void insert(const size_t &i, const CCFunction &f) {
        functions.insert(std::make_pair(i, f));
    }

    /// setter
    void set_functions(const vector_real_function_3d &v, const FuncType &type, const size_t &freeze) {
        functions.clear();
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction tmp(v[i], freeze + i, type);
            functions.insert(std::make_pair(freeze + i, tmp));
        }
    }

    /// Returns all the functions of the map as vector
    vector_real_function_3d get_vecfunction() const {
        vector_real_function_3d tmp;
        for (auto x:functions) tmp.push_back(x.second.function);
        return tmp;
    }

    /// Get the size vector (number of functions in the map)
    size_t size() const {
        return functions.size();
    }

    /// Print the memory of which is used by all the functions in the map
    void
    print_size(const std::string &msg = "!?not assigned!?") const;

    /// scalar multiplication
    CC_vecfunction operator*(const double &fac) const {
        vector_real_function_3d vnew = fac * get_vecfunction();
        const size_t freeze = functions.cbegin()->first;
        return CC_vecfunction(vnew, type, freeze);
    }

    /// scaling (inplace)
    void scale(const double &factor) {
        for (auto &ktmp:functions) {
            ktmp.second.function.scale(factor);
        }
    }

    /// operator needed for sort operation (sorted by omega values)
    bool operator<=(const CC_vecfunction &b) const { return omega <= b.omega; }

    /// operator needed for sort operation (sorted by omega values)
    bool operator<(const CC_vecfunction &b) const { return omega < b.omega; }

    /// store functions on disc
    void save_functions(const std::string msg = "") const {
        std::string pre_name = "";
        if (msg != "") pre_name = msg + "_";
        for (const auto &tmp:functions) save<double, 3>(tmp.second.function, pre_name + tmp.second.name());
    }

    // plotting
    void plot(const std::string &msg = "") const {
        for (auto &ktmp:functions) {
            ktmp.second.plot(msg);
        }
    }

};

/// Helper Structure that carries out operations on CC_functions
/// The structure can hold intermediates for g12 and f12 of type : <mo_bra_k|op|type> with type=HOLE,PARTICLE or RESPONSE
/// some 6D operations are also included
/// The structure does not know if nuclear correlation facors are used, so the corresponding bra states have to be prepared beforehand
struct CCConvolutionOperator {

    /// parameter class
    struct Parameters {
        Parameters() {};
        Parameters(const Parameters& other) :
            thresh_op(other.thresh_op),
            lo(other.lo),
            freeze(other.freeze),
            gamma(other.gamma) {
        }

        Parameters(const CCParameters &param) : thresh_op(param.thresh_poisson), lo(param.lo), freeze(param.freeze),
                                                gamma(param.gamma()) {};
        double thresh_op=FunctionDefaults<3>::get_thresh();
        double lo = 1.e-6;
        int freeze = 0;
        double gamma = 1.0; /// f12 exponent
    };


    /// @param[in] world
    /// @param[in] optype: the operatortype (can be g12_ or f12_)
    /// @param[in] param: the parameters of the current CC-Calculation (including function and operator thresholds and the exponent for f12)
    CCConvolutionOperator(World &world, const OpType type, Parameters param) : parameters(param), world(world),
                                                                                      operator_type(type),
                                                                                      op() {
    }

    CCConvolutionOperator(const CCConvolutionOperator& other) = default;

    /// @param[in] f: a 3D function
    /// @param[out] the convolution op(f), no intermediates are used
    real_function_3d operator()(const real_function_3d &f) const { return ((*op)(f)).truncate(); }

    /// @param[in] bra a CC_vecfunction
    /// @param[in] ket a CC_function
    /// @param[out] vector[i] = <bra[i]|op|ket>
    vector_real_function_3d operator()(const CC_vecfunction &bra, const CCFunction &ket) const {
        vector_real_function_3d result;
        if (bra.type == HOLE) {
            for (const auto &ktmp:bra.functions) {
                const CCFunction &brai = ktmp.second;
                const real_function_3d tmpi = this->operator()(brai, ket);
                result.push_back(tmpi);
            }
        } else {
            vector_real_function_3d tmp = mul(world, ket.function, bra.get_vecfunction());
            result = apply(world, (*op), tmp);
            truncate(world, result);
        }
        return result;
    }

    // @param[in] f: a vector of 3D functions
    // @param[out] the convolution of op with each function, no intermeditates are used
    vector_real_function_3d operator()(const vector_real_function_3d &f) const {
        return apply<double, double, 3>(world, (*op), f);
    }

    // @param[in] bra: a 3D CC_function, if nuclear-correlation factors are used they have to be applied before
    // @param[in] ket: a 3D CC_function,
    // @param[in] use_im: default is true, if false then no intermediates are used
    // @param[out] the convolution <bra|op|ket> = op(bra*ket), if intermediates were calculated before the operator uses them
    real_function_3d operator()(const CCFunction &bra, const CCFunction &ket, const bool use_im = true) const;

    // @param[in] u: a 6D-function
    // @param[out] the convolution \int g(r,r') u(r,r') dr' (if particle==2) and g(r,r') u(r',r) dr' (if particle==1)
    // @param[in] particle: specifies on which particle of u the operator will act (particle ==1 or particle==2)
    real_function_6d operator()(const real_function_6d &u, const size_t particle) const;

    // @param[in] bra: a 3D-CC_function, if nuclear-correlation factors are used they have to be applied before
    // @param[in] u: a 6D-function
    // @param[in] particle: specifies on which particle of u the operator will act (particle ==1 or particle==2)
    // @param[out] the convolution <bra|g12|u>_particle
    real_function_3d operator()(const CCFunction &bra, const real_function_6d &u, const size_t particle) const;

    /// @param[in] bra: a vector of CC_functions, the type has to be HOLE
    /// @param[in] ket: a vector of CC_functions, the type can be HOLE,PARTICLE,RESPONSE
    /// updates intermediates of the type <bra|op|ket>
    void update_elements(const CC_vecfunction &bra, const CC_vecfunction &ket);

    /// @param[out] prints the name of the operator (convenience) which is g12 or f12 or maybe other things like gf in the future
    std::string name() const { return assign_name(operator_type); }

    /// @param[in] the type of which intermediates will be deleted
    /// e.g if(type==HOLE) then all intermediates of type <mo_bra_k|op|HOLE> will be deleted
    void clear_intermediates(const FuncType &type);

    /// name speaks for itself
    void clear_all_intermediates() {
        clear_intermediates(HOLE);
        clear_intermediates(PARTICLE);
        clear_intermediates(RESPONSE);
    }

    /// prints out information (operatorname, number of stored intermediates ...)
    size_t info() const;

    /// sanity check .. doens not do so much
    void sanity() const { print_intermediate(HOLE); }

    /// @param[in] type: the type of intermediates which will be printed, can be HOLE,PARTICLE or RESPONSE
    void print_intermediate(const FuncType type) const {
        if (type == HOLE)
            for (const auto &tmp:imH.allpairs)
                tmp.second.print_size("<H" + std::to_string(tmp.first.first) + "|" + assign_name(operator_type) + "|H" +
                                      std::to_string(tmp.first.second) + "> intermediate");
        else if (type == PARTICLE)
            for (const auto &tmp:imP.allpairs)
                tmp.second.print_size("<H" + std::to_string(tmp.first.first) + "|" + assign_name(operator_type) + "|P" +
                                      std::to_string(tmp.first.second) + "> intermediate");
        else if (type == RESPONSE)
            for (const auto &tmp:imR.allpairs)
                tmp.second.print_size("<H" + std::to_string(tmp.first.first) + "|" + assign_name(operator_type) + "|R" +
                                      std::to_string(tmp.first.second) + "> intermediate");
    }

    /// create a TwoElectronFactory with the operatorkernel
    TwoElectronFactory get_kernel() const {
        if (type() == OT_G12) return TwoElectronFactory(world).dcut(1.e-7);
        else if (type() == OT_F12) return TwoElectronFactory(world).dcut(1.e-7).f12().gamma(parameters.gamma);
        else error("no kernel of type " + name() + " implemented");
        return TwoElectronFactory(world);
    }

    OpType type() const { return operator_type; }

    const Parameters parameters;
private:
    /// the world
    World &world;
    /// the operatortype, currently this can be g12_ or f12_
    const OpType operator_type = OT_UNDEFINED;

    /// @param[in] optype: can be f12_ or g12_ depending on which operator shall be intitialzied
    /// @param[in] parameters: parameters (thresholds etc)
    /// initializes the operators
    SeparatedConvolution<double, 3> *init_op(const OpType &type, const Parameters &parameters) const;

    std::shared_ptr<real_convolution_3d> op;
    intermediateT imH;
    intermediateT imP;
    intermediateT imR;

    /// @param[in] msg: output message
    /// the function will throw an MADNESS_EXCEPTION
    void error(const std::string &msg) const {
        if (world.rank() == 0)
            std::cout << "\n\n!!!!ERROR in CCConvolutionOperator " << assign_name(operator_type) << ": " << msg
                      << "!!!!!\n\n" << std::endl;
        MADNESS_EXCEPTION(msg.c_str(), 1);
    }
};

/// Helper structure for the coupling potential of CC Singles and Doubles
/// because of the regularization of the CC-Wavefunction (for CC2: |tauij> = |uij> + Qt12*f12*|titj>)
/// we have 6D-functions in std format |u> : type==pure_
/// we have 6D-functions in sepparated format: type==decomposed_ (e.g O1*f12*|titj> = |xy> with x=|k> and y=<k|f12|ti>*|tj>)
/// we have 6D-function like f12|xy> which are not needed to be represented on the 6D MRA-Grid, type==op_decomposed_
struct CCPairFunction {

public:
    CCPairFunction(World &world, const real_function_6d &ket) : world(world), type(PT_FULL), a(), b(), op(0), u(ket) {}

    CCPairFunction(World &world, const vector_real_function_3d &f1, const vector_real_function_3d &f2) : world(world),
                                                                                                         type(PT_DECOMPOSED),
                                                                                                         a(f1), b(f2),
                                                                                                         op(0), u() {}

    CCPairFunction(World &world, const std::pair<vector_real_function_3d, vector_real_function_3d> &f) : world(world),
                                                                                                         type(PT_DECOMPOSED),
                                                                                                         a(f.first),
                                                                                                         b(f.second),
                                                                                                         op(0), u() {}

    CCPairFunction(World &world, const CCConvolutionOperator *op_, const CCFunction &f1, const CCFunction &f2) : world(
            world), type(PT_OP_DECOMPOSED), a(), b(), op(op_), x(f1), y(f2), u() {}

    CCPairFunction(const CCPairFunction &other) : world(other.world), type(other.type), a(other.a), b(other.b),
                                                  op(other.op), x(other.x), y(other.y), u(other.u) {}

    CCPairFunction
    operator=(CCPairFunction &other);

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
    void
    print_size() const;

    std::string name() const {
        if (type == PT_FULL) return "|u>";
        else if (type == PT_DECOMPOSED) return "|ab>";
        else if (type == PT_OP_DECOMPOSED) return op->name() + "|xy>";
        return "???";
    }


    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D)
    real_function_3d project_out(const CCFunction &f, const size_t particle) const;

    // result is: <x|op12|f>_particle
    /// @param[in] x: a 3D-CC_function
    /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
    /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
    /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
    real_function_3d
    dirac_convolution(const CCFunction &x, const CCConvolutionOperator &op, const size_t particle) const;

    /// @param[out] particles are interchanged, if the function was u(1,2) the result is u(2,1)
    CCPairFunction swap_particles() const;

    /// @param[in] the Greens operator
    /// @param[out] the Greens operator is applied to the function: G(u)
    real_function_6d apply_G(const real_convolution_6d &G) const;


    double
    make_xy_u(const CCFunction &xx, const CCFunction &yy) const;

public:
    /// the 3 types of 6D-function that occur in the CC potential which coupled doubles to singles

    World &world;
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
    real_function_3d project_out_decomposed(const real_function_3d &f, const size_t particle) const;

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D) for the case that u=op|xy> so <f|u>_particle = <f|op|x>*|y> if particle==1
    real_function_3d project_out_op_decomposed(const CCFunction &f, const size_t particle) const;

    /// @param[in] x: a 3D-CC_function
    /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
    /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
    /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
    /// in this case u=|ab> and the result is <x|op|u>_1 = <x|op|a>*|b> for particle==1
    real_function_3d
    dirac_convolution_decomposed(const CCFunction &x, const CCConvolutionOperator &op, const size_t particle) const;

    /// small helper function that gives back (a,b) or (b,a) depending on the value of particle
    const std::pair<vector_real_function_3d, vector_real_function_3d> assign_particles(const size_t particle) const;

    /// swap particle function if type==pure_
    CCPairFunction swap_particles_pure() const;

    /// swap particle function if type==decomposed_
    CCPairFunction swap_particles_decomposed() const;

    /// swap particle function if type==op_decomposed_ (all ops are assumed to be symmetric)
    CCPairFunction swap_particles_op_decomposed() const;
};

class CCPair : public archive::ParallelSerializableObject {
public:
    CCPair(const size_t ii, const size_t jj, const CCState t, const CalcType c) : type(t), ctype(c), i(ii), j(jj),
                                                                                  bsh_eps(12345.6789) {};

    CCPair(const size_t ii, const size_t jj, const CCState t, const CalcType c, const std::vector<CCPairFunction> &f)
            : type(t), ctype(c), i(ii), j(jj), functions(f), bsh_eps(12345.6789) {};

    CCPair(const CCPair &other) : type(other.type), ctype(other.ctype), i(other.i), j(other.j),
                                  functions(other.functions), constant_part(other.constant_part),
                                  bsh_eps(other.bsh_eps) {};

    const CCState type;
    const CalcType ctype;
    const size_t i;
    const size_t j;
    int excitation = -1;

    /// gives back the pure 6D part of the pair function
    real_function_6d function() const {
        MADNESS_ASSERT(not functions.empty());
        MADNESS_ASSERT(functions[0].type == PT_FULL);
        return functions[0].u;
    }

    /// updates the pure 6D part of the pair function
    void update_u(const real_function_6d &u) {
        MADNESS_ASSERT(not functions.empty());
        MADNESS_ASSERT(functions[0].type == PT_FULL);
        CCPairFunction tmp(u.world(), u);
        functions[0] = tmp;
    }

    /// the functions which belong to the pair
    std::vector<CCPairFunction> functions;

    /// the constant part
    real_function_6d constant_part;

    /// Energy for the BSH Operator
    /// Ground State: e_i + e_j
    /// Excited State: e_i + e_j + omega
    double bsh_eps;

    std::string name() const {
        std::string name = "???";
        if (type == GROUND_STATE) name = assign_name(ctype) + "_pair_u_";
        if (type == EXCITED_STATE) name = assign_name(ctype) + "_pair_x_";
        return name + stringify(i) + stringify(j);
    }

    void
    info() const;

};

/// little helper structure which manages the stored singles potentials
struct CCIntermediatePotentials {
    CCIntermediatePotentials(World &world, const CCParameters &p) : world(world), parameters(p) {};

    /// fetches the correct stored potential or throws an exception
    vector_real_function_3d
    operator()(const CC_vecfunction &f, const PotentialType &type) const;

    /// fetch the potential for a single function
    real_function_3d
    operator()(const CCFunction &f, const PotentialType &type) const;

    vector_real_function_3d
    get_unprojected_cc2_projector_response() const { return unprojected_cc2_projector_response_; }

    void add_unprojected_cc2_projector_response(
            const vector_real_function_3d &tmp) { unprojected_cc2_projector_response_ = copy(world, tmp); }

    /// deltes all stored potentials
    void clear_all() {
        current_singles_potential_gs_.clear();
        current_singles_potential_ex_.clear();
        current_s2b_potential_gs_.clear();
        current_s2b_potential_ex_.clear();
        current_s2c_potential_gs_.clear();
        current_s2c_potential_ex_.clear();
    }

    /// clears only potentials of the response
    void clear_response() {
        current_singles_potential_ex_.clear();
        current_s2b_potential_ex_.clear();
        current_s2c_potential_ex_.clear();
    }

    /// insert potential
    void
    insert(const vector_real_function_3d &potential, const CC_vecfunction &f, const PotentialType &type);

private:
    World &world;
    const CCParameters &parameters;
    /// whole ground state singles potential without fock-residue
    vector_real_function_3d current_singles_potential_gs_;
    /// whole excited state singles potential without fock-residue
    vector_real_function_3d current_singles_potential_ex_;
    /// s2b_potential for the pure 6D-part of the ground-state (expensive and constant during singles iterations)
    vector_real_function_3d current_s2b_potential_gs_;
    /// s2b_potential for the pure 6D-part of the excited-state (expensive and constant during singles iterations)
    vector_real_function_3d current_s2b_potential_ex_;
    /// s2c_potential for the pure 6D-part of the ground-state (expensive and constant during singles iterations)
    vector_real_function_3d current_s2c_potential_gs_;
    /// s2c_potential for the pure 6D-part of the excited_state (expensive and constant during singles iterations)
    vector_real_function_3d current_s2c_potential_ex_;
    /// unprojected S3c + S5c + S2b + S2c potential of CC2 singles
    /// for the projector response of the CC2 singles potential
    vector_real_function_3d unprojected_cc2_projector_response_;

    /// structured output
    void output(const std::string &msg) const {
        if (world.rank() == 0 and parameters.debug)
            std::cout << "Intermediate Potential Manager: " << msg << "\n";
    }
};


}//namespace madness

#endif /* CCSTRUCTURES_H_ */
