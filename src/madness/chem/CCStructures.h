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
#include<madness/mra/commandlineparser.h>
#include<madness/chem/ccpairfunction.h>
#include<madness/mra/QCCalculationParametersBase.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <madness/mra/macrotaskq.h>

#include "lowrankfunction.h"

namespace madness {

/// Calculation Types used by CC2
enum CalcType {
    CT_UNDEFINED, CT_MP2, CT_MP3, CT_CC2, CT_LRCCS, CT_LRCC2, CT_CISPD, CT_ADC2, CT_TDHF, CT_TEST
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
assign_name(const CCState& input);

/// Assigns enum to string
CalcType
assign_calctype(const std::string name);

/// Assigns strings to enums for formated output
std::string
assign_name(const CalcType& inp);

/// Assigns strings to enums for formated output
std::string
assign_name(const PotentialType& inp);

/// Assigns strings to enums for formated output
std::string
assign_name(const FuncType& inp);

/// check memory usage using getrusage
inline void print_memory_usage(const World& world) {
    long mem=get_memory_usage();
    std::string hostname=get_hostname();
    std::stringstream ss;
    ss << "memory usage of process "<< world.rank()<< " on "<< hostname<< ": "<< mem/1024/1024<<"MB";
    std::string msg=ss.str();
    auto memusage=world.gop.concat0(std::vector<std::string>(1,msg));
    std::sort(memusage.begin(),memusage.end());
    if (world.rank()==0) for (const auto& msg : memusage) print(msg);
}

// Little structure for formated output and to collect warnings
// much room to improve
struct CCMessenger {
    CCMessenger(World& world) : world(world), output_prec(10), scientific(true), debug(false), os(std::cout) {}

    World& world;
    size_t output_prec;
    bool scientific;
    bool debug;

    void operator()(const std::string& msg) const { output(msg); }

    void debug_output(const std::string& msg) const {
        if (debug) output(msg);
    }

    void
    output(const std::string& msg) const;

    void
    section(const std::string& msg) const;

    void
    subsection(const std::string& msg) const;

    void
    warning(const std::string& msg) const;

    void print_warnings() const {
        for (const auto& x:warnings) if (world.rank() == 0) std::cout << x << "\n";
    }

    template<class T>
    CCMessenger operator<<(const T& t) const {
        using madness::operators::operator<<;
        if (world.rank() == 0) os << t;
        return *this;
    }

    /// collect all warnings that occur to print out at the end of the job
    mutable std::vector<std::string> warnings;
    /// output stream
    std::ostream& os;
};


/// Timer Structure
struct CCTimer {
    /// TDA_TIMER constructor
    /// @param[in] world the world
    /// @param[in] msg	a string that contains the desired printout when info function is called
    CCTimer(World& world, std::string msg) : world(world), start_wall(wall_time()), start_cpu(cpu_time()),
                                             operation(msg), end_wall(0.0), end_cpu(0.0), time_wall(-1.0),
                                             time_cpu(-1.0) {}

    World& world;
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

    double reset() {
        stop();
        double wtime=time_wall;
        start();
        return wtime;
    }


    double get_wall_time_diff() const { return end_wall; }

    double get_cpu_time_diff() const { return end_cpu; }

    std::pair<double, double> current_time(bool printout = false) {
        if (time_wall < 0.0 or time_cpu < 0.0) stop();
        return std::make_pair(time_wall, time_cpu);
    }

    void print() {
        print(current_time());
    }

    void print() const {
        print(std::make_pair(time_wall, time_cpu));
    }

    void print(const std::pair<double, double>& times) const {
        if (world.rank() == 0) {
            std::cout << std::setfill(' ') << std::scientific << std::setprecision(2)
                      << "Timer: " << times.first << " (Wall), " << times.second << " (CPU)" << ", (" + operation + ")"
                      << "\n";
        }
    }
};

/// Calculation TDHFParameters for CC2 and TDA calculations
/// Maybe merge this with calculation_parameters of SCF at some point, or split into TDA and CC
struct CCParameters : public QCCalculationParametersBase {

    CCParameters() {
        initialize_parameters();
    };

    /// copy constructor
    CCParameters(const CCParameters& other) =default;

    /// ctor reading out the input file
    CCParameters(World& world, const commandlineparser& parser) {
        initialize_parameters();
        read_input_and_commandline_options(world,parser,"cc2");
        set_derived_values();
    };


    void initialize_parameters() {
        double thresh=1.e-3;
        double thresh_operators=1.e-6;
        initialize < std::string > ("calc_type", "mp2", "the calculation type", {"mp2", "mp3", "cc2", "cis", "lrcc2", "cispd", "adc2", "test"});
        initialize < double > ("lo", 1.e-7, "the finest length scale to be resolved by 6D operators");
        initialize < double > ("dmin", 1.0, "defines the depth of the special level");
        initialize < double > ("thresh_6d", thresh, "threshold for the 6D wave function");
        initialize < double > ("tight_thresh_6d", 0.1*thresh, "tight threshold for the 6D wave function");
        initialize < double > ("thresh_3d", 0.01*thresh, "threshold for the 3D reference wave function");
        initialize < double > ("tight_thresh_3d", 0.001*thresh, "tight threshold for the 3D reference wave function");
        initialize < double > ("thresh_bsh_3d", thresh_operators, "threshold for BSH operators");
        initialize < double > ("thresh_bsh_6d", thresh_operators, "threshold for BSH operators");
        initialize < double > ("thresh_poisson", thresh_operators, "threshold for Poisson operators");
        initialize < double > ("thresh_f12", thresh_operators, "threshold for Poisson operators");
        initialize < double > ("thresh_Ue", thresh_operators, "ue threshold");
        initialize < double > ("econv", thresh, "overal convergence threshold ");
        initialize < double > ("econv_pairs", 0.1*thresh, "convergence threshold for pairs");
        initialize < double > ("dconv_3d", 0.3*thresh, "convergence for cc singles");
        initialize < double > ("dconv_6d", 3.0*thresh, "convergence for cc doubles");
        initialize < std::size_t > ("iter_max", 10, "max iterations");
        initialize < std::size_t > ("iter_max_3d", 10, "max iterations for singles");
        initialize < std::size_t > ("iter_max_6d", 10, "max iterations for doubles");
        initialize < std::pair<int, int>> ("only_pair", {-1, -1}, "compute only a single pair");
        initialize < bool > ("restart", false, "restart");
        initialize < bool > ("no_compute", false, "no compute");
        initialize < bool > ("no_compute_gs", false, "no compute");
        initialize < bool > ("no_compute_mp2_constantpart", false, "no compute");
        initialize < bool > ("no_compute_response", false, "no compute");
        initialize < bool > ("no_compute_mp2", false, "no compute");
        initialize < bool > ("no_compute_cc2", false, "no compute");
        initialize < bool > ("no_compute_cispd", false, "no compute");
        initialize < bool > ("no_compute_lrcc2", false, "no compute");
        initialize < double > ("corrfac_gamma", 1.0, "exponent for the correlation factor");
        initialize < std::size_t > ("output_prec", 8, "for formatted output");
        initialize < bool > ("debug", false, "");
        initialize < bool > ("plot", false, "");
        initialize < bool > ("kain", true, "");
        initialize < std::size_t > ("kain_subspace", 3, "");
        initialize < long > ("freeze", -1, "number of frozen orbitals: -1: automatic");
        initialize < bool > ("test", false, "");
        // choose if Q for the constant part of MP2 and related calculations should be decomposed: GQV or GV - GO12V
        initialize < bool > ("decompose_Q", true, "always true",{true});
        // if true the ansatz for the CC2 ground state pairs is |tau_ij> = |u_ij> + Qtf12|titj>, with Qt = Q - |tau><phi|
        // if false the ansatz is the same with normal Q projector
        // the response ansatz is the corresponding response of the gs ansatz
        initialize < bool > ("QtAnsatz", true, "always true",{true});
        // a vector containing the excitations which shall be optizmized later (with CIS(D) or CC2)
        initialize < std::vector<size_t>>
        ("excitations", {}, "vector containing the excitations");
    }

    void set_derived_values();

    CalcType calc_type() const {
        std::string value = get<std::string>("calc_type");
        if (value == "mp2") return CT_MP2;
        if (value == "mp3") return CT_MP3;
        if (value == "cc2") return CT_CC2;
        if (value == "cis") return CT_LRCCS;
        if (value == "lrcc2") return CT_LRCC2;
        if (value == "cispd") return CT_CISPD;
        if (value == "adc2") return CT_ADC2;
        if (value == "test") return CT_TEST;
        MADNESS_EXCEPTION("faulty CalcType", 1);
    }

    bool response() const {return calc_type()==CT_ADC2 or calc_type()==CT_CISPD or calc_type()==CT_LRCC2 or calc_type()==CT_LRCCS;}
    double lo() const { return get<double>("lo"); }

    double dmin() const { return get<double>("dmin"); }

    double thresh_3D() const { return get<double>("thresh_3d"); }

    double tight_thresh_3D() const { return get<double>("tight_thresh_3d"); }

    double thresh_6D() const { return get<double>("thresh_6d"); }

    double tight_thresh_6D() const { return get<double>("tight_thresh_6d"); }

    double thresh_bsh_3D() const { return get<double>("thresh_bsh_3d"); }

    double thresh_bsh_6D() const { return get<double>("thresh_bsh_6d"); }

    double thresh_poisson() const { return get<double>("thresh_poisson"); }

    double thresh_f12() const { return get<double>("thresh_f12"); }

    double thresh_Ue() const { return get<double>("thresh_ue"); }

    double econv() const { return get<double>("econv"); }

    double econv_pairs() const { return get<double>("econv_pairs"); }

    double dconv_3D() const { return get<double>("dconv_3d"); }

    double dconv_6D() const { return get<double>("dconv_6d"); }

    std::size_t iter_max() const { return get<std::size_t>("iter_max"); }

    std::size_t iter_max_3D() const { return get<std::size_t>("iter_max_3d"); }

    std::size_t iter_max_6D() const { return get<std::size_t>("iter_max_6d"); }

    std::pair<int, int> only_pair() const { return get<std::pair<int, int>>("only_pair"); }

    bool restart() const { return get<bool>("restart"); }

    bool no_compute() const { return get<bool>("no_compute"); }

    bool no_compute_gs() const { return get<bool>("no_compute_gs"); }

    bool no_compute_mp2_constantpart() const { return get<bool>("no_compute_mp2_constantpart"); }

    bool no_compute_response() const { return get<bool>("no_compute_response"); }

    bool no_compute_mp2() const { return get<bool>("no_compute_mp2"); }

    bool no_compute_cc2() const { return get<bool>("no_compute_cc2"); }

    bool no_compute_cispd() const { return get<bool>("no_compute_cispd"); }

    bool no_compute_lrcc2() const { return get<bool>("no_compute_lrcc2"); }

    bool debug() const { return get<bool>("debug"); }

    bool plot() const { return get<bool>("plot"); }

    bool kain() const { return get<bool>("kain"); }

    bool test() const { return get<bool>("test"); }

    bool decompose_Q() const { return get<bool>("decompose_q"); }

    bool QtAnsatz() const { return get<bool>("qtansatz"); }

    std::size_t output_prec() const { return get<std::size_t>("output_prec"); }

    std::size_t kain_subspace() const { return get<std::size_t>("kain_subspace"); }

    long freeze() const { return get<long>("freeze"); }

    std::vector<std::size_t> excitations() const { return get<std::vector<std::size_t>>("excitations"); }

    double gamma() const {return get<double>("corrfac_gamma");}

    /// print out the parameters
    void information(World& world) const;

    /// check if parameters are set correct
    void sanity_check(World& world) const;

    void error(World& world, const std::string& msg) const {
        if (world.rank() == 0)
            std::cout << "\n\n\n\n\n!!!!!!!!!\n\nERROR IN CC_PARAMETERS:\n    ERROR MESSAGE IS: " << msg
                      << "\n\n\n!!!!!!!!" << std::endl;
        MADNESS_EXCEPTION("ERROR IN CC_PARAMETERS", 1);
    }

    size_t warning(World& world, const std::string& msg) const {
        if (world.rank() == 0) std::cout << "WARNING IN CC_PARAMETERS!: " << msg << std::endl;
        return 1;
    }

};

struct PairVectorMap {

    std::vector<std::pair<int, int>> map; ///< maps pair index (i,j) to vector index k
    PairVectorMap() = default;
    PairVectorMap(const std::vector<std::pair<int, int>> map1) : map(map1) {}

    static PairVectorMap triangular_map(const int nfreeze, const int nocc) {
        std::vector<std::pair<int, int>> map; ///< maps pair index (i,j) to vector index k
        for (int i=nfreeze; i<nocc; ++i) {
            for (int j=i; j<nocc; ++j) {
                map.push_back(std::make_pair(i,j));
            }
        }
        return PairVectorMap(map);
    }

    static PairVectorMap quadratic_map(const int nfreeze, const int nocc) {
        std::vector<std::pair<int, int>> map; ///< maps pair index (i,j) to vector index k
        for (int i=nfreeze; i<nocc; ++i) {
            for (int j=nfreeze; j<nocc; ++j) {
                map.push_back(std::make_pair(i,j));
            }
        }
        return PairVectorMap(map);
    }

    void print(const std::string msg="PairVectorMap") const {
        madness::print(msg);
        madness::print("vector element <-> pair index");
        for (size_t i=0; i<map.size(); ++i) {
            madness::print(i, " <-> ",map[i]);
        }
    }

};

/// POD holding all electron pairs with easy access
/// Similar strucutre than the Pair structure from MP2 but with some additional features (merge at some point)
/// This structure will also be used for intermediates
template<typename T>
struct Pairs {

    typedef std::map<std::pair<int, int>, T> pairmapT;
    pairmapT allpairs;

    /// convert Pairs<T> to another type

    /// opT op takes an object of T and returns the result type
    template<typename R, typename opT>
    Pairs<R> convert(const Pairs<T> arg, const opT op) const {
        Pairs<R> result;
        for (auto& p : arg.allpairs) {
            int i=p.first.first;
            int j=p.first.second;
            result.insert(i,j,op(p.second));
        }
        return result;
    }

    static Pairs vector2pairs(const std::vector<T>& argument, const PairVectorMap map) {
        Pairs<T> pairs;
        for (int i=0; i<argument.size(); ++i) {
            pairs.insert(map.map[i].first,map.map[i].second,argument[i]);
        }
        return pairs;
    }

    static std::vector<T> pairs2vector(const Pairs<T>& argument, const PairVectorMap map) {
        std::vector<T> vector;
        for (size_t i=0; i<argument.allpairs.size(); ++i) {
            vector.push_back(argument(map.map[i].first,map.map[i].second));
        }
        return vector;
    }

    /// getter
    const T& operator()(int i, int j) const {
        return allpairs.at(std::make_pair(i, j));
    }

    /// getter
    // at instead of [] operator bc [] inserts new element if nothing is found while at throws out of range error
    // back to before
    T& operator()(int i, int j) {
        // return allpairs.at(std::make_pair(i, j));
        return allpairs[std::make_pair(i, j)];
    }

    /// setter
    /// can NOT replace elements (for this construct new pair map and swap the content)
    void insert(int i, int j, const T& pair) {
        std::pair<int, int> key = std::make_pair(i, j);
        allpairs.insert(std::make_pair(key, pair));
    }

    /// swap the contant of the pairmap
    void swap(Pairs<T>& other) {
        allpairs.swap(other.allpairs);
    }

    bool empty() const {
        if (allpairs.size() == 0) return true;
        else return false;
    }
};

/// f12 and g12 intermediates of the form <f1|op|f2> (with op=f12 or op=g12) will be saved using the pair structure
template <typename T, std::size_t NDIM>
using intermediateT = Pairs<Function<T,NDIM>>;

/// Returns the size of an intermediate
//double
//size_of(const intermediateT& im);
/// Returns the size of an intermediate
template<typename T, std::size_t NDIM>
double
size_of(const intermediateT<T,NDIM>& im) {
    double size = 0.0;
    for (const auto& tmp : im.allpairs) {
        size += get_size<T, NDIM>(tmp.second);
    }
    return size;
}



// structure for CC Vectorfunction
/// A helper structure which holds a map of functions
struct CC_vecfunction : public archive::ParallelSerializableObject {

    CC_vecfunction() : type(UNDEFINED), omega(0.0), current_error(99.9), delta(0.0) {}

    CC_vecfunction(const FuncType type_) : type(type_), omega(0.0), current_error(99.9), delta(0.0) {}

    CC_vecfunction(const vector_real_function_3d& v) : type(UNDEFINED), omega(0.0), current_error(99.9), delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction<double,3> tmp(v[i], i, type);
            functions.insert(std::make_pair(i, tmp));
        }
    }

    CC_vecfunction(const std::vector<CCFunction<double,3>>& v) : type(UNDEFINED), omega(0.0), current_error(99.9), delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            functions.insert(std::make_pair(v[i].i, v[i]));
        }
    }

    CC_vecfunction(const vector_real_function_3d& v, const FuncType& type) : type(type), omega(0.0),
                                                                             current_error(99.9), delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction<double,3> tmp(v[i], i, type);
            functions.insert(std::make_pair(i, tmp));
        }
    }

    CC_vecfunction(const vector_real_function_3d& v, const FuncType& type, const size_t& freeze) : type(type),
                                                                                                   omega(0.0),
                                                                                                   current_error(99.9),
                                                                                                   delta(0.0) {
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction<double,3> tmp(v[i], freeze + i, type);
            functions.insert(std::make_pair(freeze + i, tmp));
        }
    }

    CC_vecfunction(const std::vector<CCFunction<double,3>>& v, const FuncType type_)
            : type(type_), omega(0.0), current_error(99.9), delta(0.0) {
        for (auto x:v) functions.insert(std::make_pair(x.i, x));
    }

    /// copy ctor (shallow)
    CC_vecfunction(const CC_vecfunction& other)
            : functions(other.functions), type(other.type), omega(other.omega),
              current_error(other.current_error),
              delta(other.delta), irrep(other.irrep) {
    }

    /// assignment operator, shallow wrt the functions
//    CC_vecfunction& operator=(const CC_vecfunction& other) = default;
    CC_vecfunction& operator=(const CC_vecfunction& other) {
        if (this == &other) return *this;
        functions = other.functions;
        type = other.type;
        omega = other.omega;
        current_error = other.current_error;
        delta = other.delta;
        irrep = other.irrep;
        return *this;
    }


    /// returns a deep copy (void shallow copy errors)
    friend CC_vecfunction
    copy(const CC_vecfunction& other) {
        CC_vecfunction tmp=other;
        tmp.functions.clear();
        for (const auto& x : other.functions) {
            tmp.functions.insert(std::make_pair(x.first, copy(x.second)));
        }
        return tmp;
    }

    void reconstruct() const {
        for (auto& x : functions) x.second.function.reconstruct();
    }

//madness::CC_vecfunction
//CC_vecfunction::copy() const {
//    std::vector<CCFunction<double,3>> vn;
//    for (auto x : functions) {
//        const CCFunction<double,3> fn(madness::copy(x.second.function), x.second.i, x.second.type);
//        vn.push_back(fn);
//    }
//    CC_vecfunction result(vn, type);
//    result.irrep = irrep;
//    return result;
//}
//

    static CC_vecfunction load_restartdata(World& world, std::string filename) {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, filename.c_str());
        CC_vecfunction tmp;
        ar & tmp;
        return tmp;
    }

    void save_restartdata(World& world, std::string filename) const {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, filename.c_str());
        ar & *this;
    }

    template<typename Archive>
    void serialize(const Archive& ar) {
        typedef std::vector<std::pair<std::size_t, CCFunction<double,3>>> CC_functionvec;

        auto map2vector = [] (const CC_functionmap& map) {
            return CC_functionvec(map.begin(), map.end());
        };
        auto vector2map = [] (const CC_functionvec& vec) {
            return CC_functionmap(vec.begin(), vec.end());
        };

        ar & type & omega & current_error & delta & irrep ;
        if (ar.is_input_archive) {
	    std::size_t size=0; // set to zero to silence compiler warning
            ar & size;
            CC_functionvec tmp(size);

            for (auto& t : tmp) ar & t.first & t.second;
            functions=vector2map(tmp);
        } else if (ar.is_output_archive) {
            auto tmp=map2vector(functions);
            ar & tmp.size();
            for (auto& t : tmp) ar & t.first & t.second;
        }
    }

    hashT hash() const {
        hashT hashval = std::hash<FuncType>{}(type);
        for (const auto& f : functions) hash_combine(hashval, hash_value(f.second.f().get_impl()->id()));

        return hashval;
    }

    typedef std::map<std::size_t, CCFunction<double,3>> CC_functionmap;
    CC_functionmap functions;

    FuncType type;
    double omega; /// excitation energy
    double current_error;
    double delta; // Last difference in Energy
    std::string irrep = "null";    /// excitation irrep (direct product of x function and corresponding orbital)

    std::string
    name(const int ex) const {
        return madness::name(type,ex);
    };

    bool is_converged(const double econv, const double dconv) const {
        return (current_error<dconv) and (std::fabs(delta)<econv);
    }

    /// getter
    const CCFunction<double,3>& operator()(const CCFunction<double,3>& i) const {
        return functions.find(i.i)->second;
    }

    /// getter
    const CCFunction<double,3>& operator()(const size_t& i) const {
        return functions.find(i)->second;
    }

    /// getter
    CCFunction<double,3>& operator()(const CCFunction<double,3>& i) {
        return functions[i.i];
    }

    /// getter
    CCFunction<double,3>& operator()(const size_t& i) {
        return functions[i];
    }

    /// setter
    void insert(const size_t& i, const CCFunction<double,3>& f) {
        functions.insert(std::make_pair(i, f));
    }

    /// setter
    void set_functions(const vector_real_function_3d& v, const FuncType& type, const size_t& freeze) {
        functions.clear();
        for (size_t i = 0; i < v.size(); i++) {
            CCFunction<double,3> tmp(v[i], freeze + i, type);
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
    print_size(const std::string& msg = "!?not assigned!?") const;

    /// scalar multiplication
    CC_vecfunction operator*(const double& fac) const {
        vector_real_function_3d vnew = fac * get_vecfunction();
        const size_t freeze = functions.cbegin()->first;
        return CC_vecfunction(vnew, type, freeze);
    }

    /// scaling (inplace)
    void scale(const double& factor) {
        for (auto& ktmp:functions) {
            ktmp.second.function.scale(factor);
        }
    }

    /// operator needed for sort operation (sorted by omega values)
    bool operator<(const CC_vecfunction& b) const { return omega < b.omega; }

    // plotting
    void plot(const std::string& msg = "") const {
        for (auto& ktmp:functions) {
            ktmp.second.plot(msg);
        }
    }
public:
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(CC_vecfunction, omega, irrep, current_error)

};

/// Helper Structure that carries out operations on CC_functions
/// The structure can hold intermediates for g12 and f12 of type : <mo_bra_k|op|type> with type=HOLE,PARTICLE or RESPONSE
/// some 6D operations are also included
/// The structure does not know if nuclear correlation facors are used, so the corresponding bra states have to be prepared beforehand
template<typename T=double, std::size_t NDIM=3>
class CCConvolutionOperator {
public:

    /// parameter class
    struct Parameters {
        Parameters() {};

        Parameters(const Parameters& other) :
                thresh_op(other.thresh_op),
                lo(other.lo),
                freeze(other.freeze),
                gamma(other.gamma) {
        }

        Parameters(const CCParameters& param) : thresh_op(param.thresh_poisson()), lo(param.lo()),
                                                freeze(param.freeze()),
                                                gamma(param.gamma()) {};
        double thresh_op = FunctionDefaults<3>::get_thresh();
        double lo = 1.e-6;
        int freeze = 0;
        double gamma = 1.0; /// f12 exponent

        template<typename archiveT>
        void serialize(archiveT& ar) {
            ar & thresh_op & lo & freeze & gamma;
        }
    };


    /// @param[in] world
    /// @param[in] optype: the operatortype (can be g12_ or f12_)
    /// @param[in] param: the parameters of the current CC-Calculation (including function and operator thresholds and the exponent for f12)
    CCConvolutionOperator(World& world, const OpType type, Parameters param) : parameters(param), world(world),
                                                                               op(init_op(type, parameters)) {
    }

    CCConvolutionOperator(const CCConvolutionOperator& other) = default;

    static inline
    std::shared_ptr<CCConvolutionOperator> CCConvolutionOperatorPtr(World& world, const OpType type, Parameters param) {
        return std::shared_ptr<CCConvolutionOperator>(new CCConvolutionOperator(world, type, param));
    }

protected:

    friend CCConvolutionOperator combine(const CCConvolutionOperator& a, const CCConvolutionOperator& b) {
        auto info= SeparatedConvolution<T,NDIM>::combine_OT((*a.get_op()),(*b.get_op()));
        Parameters param;
        param.gamma=info.mu;
        param.thresh_op=info.thresh;
        param.lo=info.lo;
        param.freeze=a.parameters.freeze;
        return CCConvolutionOperator(a.world, info.type, param);
    }

    friend std::shared_ptr<CCConvolutionOperator> combine(const std::shared_ptr<CCConvolutionOperator>& a,
                                         const std::shared_ptr<CCConvolutionOperator>& b) {
        if (a and (not b)) return a;
        if ((not a) and b) return b;
        if ((not a) and (not b)) return nullptr;
        return std::shared_ptr<CCConvolutionOperator>(new CCConvolutionOperator(combine(*a,*b)));
    }

public:
    /// @param[in] f: a 3D function
    /// @param[out] the convolution op(f), no intermediates are used
    Function<T,NDIM> operator()(const Function<T,NDIM>& f) const {
        if (op) return ((*op)(f)).truncate();
        return f;
    }

    /// @param[in] bra a CC_vecfunction
    /// @param[in] ket a CC_function
    /// @param[out] vector[i] = <bra[i]|op|ket>
    std::vector<Function<T,NDIM>> operator()(const CC_vecfunction& bra, const CCFunction<T,NDIM>& ket) const {
        MADNESS_CHECK(op);
        std::vector<Function<T, NDIM>> result;
        if constexpr (NDIM == 3) {
            if (bra.type == HOLE) {
                for (const auto& ktmp: bra.functions) {
                    const CCFunction<T, NDIM>& brai = ktmp.second;
                    const Function<T, NDIM> tmpi = this->operator()(brai, ket);
                    result.push_back(tmpi);
                }
            } else {
                std::vector<Function<T, NDIM>> tmp = mul(world, ket.function, bra.get_vecfunction());
                result = apply(world, (*op), tmp);
                truncate(world, result);
            }
        } else {
            MADNESS_EXCEPTION("not implemented", 1);
        }

        return result;
    }

    // @param[in] f: a vector of 3D functions
    // @param[out] the convolution of op with each function, no intermeditates are used
    std::vector<Function<T,NDIM>> operator()(const std::vector<Function<T,NDIM>>& f) const {
        MADNESS_CHECK(op);
        return apply<T,T,NDIM,NDIM>(world, (*op), f);
    }

    // @param[in] bra: a 3D CC_function, if nuclear-correlation factors are used they have to be applied before
    // @param[in] ket: a 3D CC_function,
    // @param[in] use_im: default is true, if false then no intermediates are used
    // @param[out] the convolution <bra|op|ket> = op(bra*ket), if intermediates were calculated before the operator uses them
    Function<T,NDIM> operator()(const CCFunction<T,NDIM>& bra, const CCFunction<T,NDIM>& ket, const bool use_im = true) const;

    // @param[in] u: a 6D-function
    // @param[out] the convolution \int g(r,r') u(r,r') dr' (if particle==2) and g(r,r') u(r',r) dr' (if particle==1)
    // @param[in] particle: specifies on which particle of u the operator will act (particle ==1 or particle==2)
    Function<T,2*NDIM> operator()(const Function<T,2*NDIM>& u, const size_t particle) const;

    // @param[in] bra: a 3D-CC_function, if nuclear-correlation factors are used they have to be applied before
    // @param[in] u: a 6D-function
    // @param[in] particle: specifies on which particle of u the operator will act (particle ==1 or particle==2)
    // @param[out] the convolution <bra|g12|u>_particle
    Function<T,NDIM> operator()(const CCFunction<T,NDIM>& bra, const Function<T,2*NDIM>& u, const size_t particle) const;

    /// @param[in] bra: a vector of CC_functions, the type has to be HOLE
    /// @param[in] ket: a vector of CC_functions, the type can be HOLE,PARTICLE,RESPONSE
    /// updates intermediates of the type <bra|op|ket>
    void update_elements(const CC_vecfunction& bra, const CC_vecfunction& ket);

    /// @param[out] prints the name of the operator (convenience) which is g12 or f12 or maybe other things like gf in the future
    std::string name() const {
        std::stringstream ss;
        ss << type();
        return ss.str();
    }

    /// @param[in] the type of which intermediates will be deleted
    /// e.g if(type==HOLE) then all intermediates of type <mo_bra_k|op|HOLE> will be deleted
    void clear_intermediates(const FuncType& type);

    /// prints out information (operatorname, number of stored intermediates ...)
    size_t info() const;

    friend hashT hash_value(CCConvolutionOperator<T,NDIM>& op) {
        hashT h;
        hash_combine(h, op.parameters.thresh_op);
        hash_combine(h, op.parameters.lo);
        hash_combine(h, op.parameters.freeze);
        hash_combine(h, op.parameters.gamma);
        hash_combine(h, int(op.type()));
        return h;
    }

    /// sanity check .. doens not do so much
    void sanity() const { print_intermediate(HOLE); }

    /// @param[in] type: the type of intermediates which will be printed, can be HOLE,PARTICLE or RESPONSE
    void print_intermediate(const FuncType type) const {
        if (type == HOLE)
            for (const auto& tmp:imH.allpairs)
                tmp.second.print_size("<H" + std::to_string(tmp.first.first) + "|" + name() + "|H" +
                                      std::to_string(tmp.first.second) + "> intermediate");
        else if (type == PARTICLE)
            for (const auto& tmp:imP.allpairs)
                tmp.second.print_size("<H" + std::to_string(tmp.first.first) + "|" + name() + "|P" +
                                      std::to_string(tmp.first.second) + "> intermediate");
        else if (type == RESPONSE)
            for (const auto& tmp:imR.allpairs)
                tmp.second.print_size("<H" + std::to_string(tmp.first.first) + "|" + name() + "|R" +
                                      std::to_string(tmp.first.second) + "> intermediate");
    }

    /// create a TwoElectronFactory with the operatorkernel
    TwoElectronFactory<T,2*NDIM> get_kernel() const {
        auto factory=TwoElectronFactory<T,2*NDIM>(world);
        factory.set_info(op->info);
        return factory;
    }

    OpType type() const { return get_op()->info.type; }

    const Parameters parameters;

    std::shared_ptr<SeparatedConvolution<T,NDIM>> get_op() const {return op;};

private:
    /// the world
    World& world;

    /// @param[in] optype: can be f12_ or g12_ depending on which operator shall be intitialzied
    /// @param[in] parameters: parameters (thresholds etc)
    /// initializes the operators
    SeparatedConvolution<T,NDIM> *init_op(const OpType& type, const Parameters& parameters) const;

    std::shared_ptr<SeparatedConvolution<T,NDIM>> op;
    intermediateT<T,NDIM> imH;
    intermediateT<T,NDIM> imP;
    intermediateT<T,NDIM> imR;

    /// @param[in] msg: output message
    /// the function will throw an MADNESS_EXCEPTION
    void error(const std::string& msg) const {
        if (world.rank() == 0)
            std::cout << "\n\n!!!!ERROR in CCConvolutionOperator " << name() << ": " << msg
                      << "!!!!!\n\n" << std::endl;
        MADNESS_EXCEPTION(msg.c_str(), 1);
    }
public:
};

template<typename T, std::size_t NDIM>
std::shared_ptr<CCConvolutionOperator<T,NDIM>> CCConvolutionOperatorPtr(World& world, const OpType type,
                                                                       typename CCConvolutionOperator<T,NDIM>::Parameters param) {
    return std::shared_ptr<CCConvolutionOperator<T,NDIM>>(new CCConvolutionOperator<T,NDIM>(world,type,param));
}

/// little helper structure which manages the stored singles potentials
struct CCIntermediatePotentials {
    CCIntermediatePotentials() = default;
    CCIntermediatePotentials(const CCParameters& p) : parameters(p) {};
    CCIntermediatePotentials(const CCIntermediatePotentials& other) = default;
    CCIntermediatePotentials& operator=(const CCIntermediatePotentials& other) = default;

    /// check if the intermediate potential exists
    bool potential_exists(const CC_vecfunction& f, const PotentialType& type) const {
        return potential_exists(type,f.type);
    }

    /// check if the intermediate potential exists
    bool potential_exists(const PotentialType& type,const FuncType& ftype) const {
        bool exists=get_potential(type,ftype,false).size()>0;
        return exists;
    }

    /// return a vector of the intermediate potentials

    /// @param[in] ptype: the potential type (POT_SINGLES, POT_S2B, ..)
    /// @param[in] ftype: the function type (HOLE, PARTICLE, RESPONSE)
    vector_real_function_3d
    get_potential(const PotentialType& ptype, const FuncType& ftype, const bool throw_if_empty) const;

    /// fetches the correct stored potential or throws an exception
    vector_real_function_3d
    operator()(const CC_vecfunction& f, const PotentialType& type, const bool throw_if_empty) const;

    /// fetch the potential for a single function
    Function<double,3>
    operator()(const CCFunction<double,3>& f, const PotentialType& type, const bool throw_if_empty) const;

    void reconstruct() const {
        madness::reconstruct(current_s2b_potential_ex_);
        madness::reconstruct(current_s2b_potential_gs_);
        madness::reconstruct(current_s2c_potential_ex_);
        madness::reconstruct(current_s2c_potential_gs_);
        madness::reconstruct(current_singles_potential_ex_);
        madness::reconstruct(current_singles_potential_gs_);
        madness::reconstruct(unprojected_cc2_projector_response_);
    }

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
    insert(const vector_real_function_3d& potential, const CC_vecfunction& f, const PotentialType& type);

    Recordlist<Cloud::keyT> cloud_store(World& world, Cloud& cloud) const {
        Recordlist<Cloud::keyT> records;
        records+=cloud.store(world,parameters);
        records+=cloud.store(world,current_s2b_potential_ex_);
        records+=cloud.store(world,current_s2b_potential_gs_);
        records+=cloud.store(world,current_s2c_potential_ex_);
        records+=cloud.store(world,current_s2c_potential_gs_);
        records+=cloud.store(world,current_singles_potential_ex_);
        records+=cloud.store(world,current_singles_potential_gs_);
        records+=cloud.store(world,unprojected_cc2_projector_response_);
        return records;
    }

    void cloud_load(World& world, const Cloud& cloud, Recordlist<Cloud::keyT>& recordlist) {
        parameters=cloud.forward_load<CCParameters>(world,recordlist);
        current_s2b_potential_ex_=cloud.forward_load<vector_real_function_3d>(world,recordlist);
        current_s2b_potential_gs_=cloud.forward_load<vector_real_function_3d>(world,recordlist);
        current_s2c_potential_ex_=cloud.forward_load<vector_real_function_3d>(world,recordlist);
        current_s2c_potential_gs_=cloud.forward_load<vector_real_function_3d>(world,recordlist);
        current_singles_potential_ex_=cloud.forward_load<vector_real_function_3d>(world,recordlist);
        current_singles_potential_gs_=cloud.forward_load<vector_real_function_3d>(world,recordlist);
        unprojected_cc2_projector_response_=cloud.forward_load<vector_real_function_3d>(world,recordlist);
    }

    friend hashT hash_value(const CCIntermediatePotentials& ip) {
        auto hash_vector_of_functions =[](const vector_real_function_3d& v) {
            hashT h;
            for (const auto& f : v) {
                hash_combine(h, hash_value(f.get_impl()->id()));
            }
            return h;
        };
        hashT h;
        hash_combine(h, hash_vector_of_functions(ip.current_s2b_potential_ex_));
        hash_combine(h, hash_vector_of_functions(ip.current_s2b_potential_gs_));
        hash_combine(h, hash_vector_of_functions(ip.current_s2c_potential_ex_));
        hash_combine(h, hash_vector_of_functions(ip.current_s2c_potential_gs_));
        hash_combine(h, hash_vector_of_functions(ip.current_singles_potential_ex_));
        hash_combine(h, hash_vector_of_functions(ip.current_singles_potential_gs_));
        hash_combine(h, hash_vector_of_functions(ip.unprojected_cc2_projector_response_));
        return h;
    }

    CCParameters parameters;
private:
    // World& world;
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
    void output(const std::string& msg) const {
        if (parameters.debug())
            std::cout << "Intermediate Potential Manager: " << msg << "\n";
    }
};

/// POD holding some basic functions and some intermediates for the CC2 calculation

/// the class is cloud-serializable and can be used in MacroTasks
struct Info {
    std::vector<Function<double,3>> mo_ket;
    std::vector<Function<double,3>> mo_bra;
    std::vector<madness::Vector<double,3>> molecular_coordinates;
    CCParameters parameters;
    std::vector<double> orbital_energies;
    Tensor<double> fock;
    CCIntermediatePotentials intermediate_potentials;
    Function<double,3> R_square, U2, R;;
    std::vector<Function<double,3>> U1;

    vector_real_function_3d get_active_mo_ket() const {
        vector_real_function_3d result;
        for (size_t i = parameters.freeze(); i < mo_ket.size(); i++) result.push_back(mo_ket[i]);
        return result;
    }

    vector_real_function_3d get_active_mo_bra() const {
        vector_real_function_3d result;
        for (size_t i = parameters.freeze(); i < mo_bra.size(); i++) result.push_back(mo_bra[i]);
        return result;
    }

    void reconstruct() const {
        madness::reconstruct(mo_bra);
        madness::reconstruct(mo_ket);
        R_square.reconstruct();
        madness::reconstruct(U1);
        U2.reconstruct();
        intermediate_potentials.reconstruct();

    }

    /// customized function to store this to the cloud

    /// functions and constant_part can be very large and we want to split them and store them in different records
    Recordlist<Cloud::keyT> cloud_store(World& world, Cloud& cloud) const {
        Recordlist<Cloud::keyT> records;
        records+=cloud.store(world,mo_bra);
        records+=cloud.store(world,mo_ket);
        records+=cloud.store(world,parameters);
        records+=cloud.store(world,orbital_energies);
        records+=cloud.store(world,fock);
        records+=cloud.store(world,intermediate_potentials);
        records+=cloud.store(world,R_square);
        records+=cloud.store(world,molecular_coordinates);
        records+=cloud.store(world,U2);
        records+=cloud.store(world,U1);
        return records;
    }

    /// customized function to load this from the cloud

    /// functions and constant_part can be very large and we want to split them and store them in different records
    /// @param[inout] recordlist: containing the keys of the member variables -> will be reduced by the keys which are used
    void cloud_load(World& world, const Cloud& cloud, Recordlist<Cloud::keyT>& recordlist) {
        // load bookkeeping stuff in a vector
        mo_bra=cloud.forward_load<std::vector<Function<double,3>>>(world,recordlist);
        mo_ket=cloud.forward_load<std::vector<Function<double,3>>>(world,recordlist);
        parameters=cloud.forward_load<CCParameters>(world,recordlist);
        orbital_energies=cloud.forward_load<std::vector<double>>(world,recordlist);
        fock=cloud.forward_load<Tensor<double>>(world,recordlist);
        intermediate_potentials=cloud.forward_load<CCIntermediatePotentials>(world,recordlist);
        R_square=cloud.forward_load<Function<double,3>>(world,recordlist);
        molecular_coordinates=cloud.forward_load<std::vector<madness::Vector<double,3>>>(world,recordlist);
        U2=cloud.forward_load<Function<double,3>>(world,recordlist);
        U1=cloud.forward_load<std::vector<Function<double,3>>>(world,recordlist);
    }

};


class CCPair : public archive::ParallelSerializableObject {
public:
    CCPair() = default;

    CCPair(const size_t ii, const size_t jj, const CCState t, const CalcType c)
    : type(t), ctype(c), i(ii), j(jj), bsh_eps(12345.6789) {};

    CCPair(const size_t ii, const size_t jj, const CCState t, const CalcType c,
        const std::vector<CCPairFunction<double,6>>& f)
            : type(t), ctype(c), i(ii), j(jj), functions(f), bsh_eps(12345.6789) {};

    CCPair(const CCPair& other) : type(other.type), ctype(other.ctype), i(other.i), j(other.j),
                                  functions(other.functions), constant_part(other.constant_part),
                                  bsh_eps(other.bsh_eps) {};

    CCState type;
    CalcType ctype;
    size_t i;
    size_t j;

    /// customized function to store this to the cloud

    /// functions and constant_part can be very large and we want to split them and store them in different records
    /// *NOTE* only the 6d function and the constant part are stored in the cloud, not the 3d functions *NOTE*
    Recordlist<Cloud::keyT> cloud_store(World& world, Cloud& cloud) const {
        // save bookkeeping stuff in a vector
        std::vector<unsigned char> v;
        archive::VectorOutputArchive arout(v);
        bool function_is_assigned=(functions.size()>0 && functions[0].is_assigned());
        arout & type & ctype & i & j & bsh_eps & function_is_assigned & constant_part.is_initialized();

        Recordlist<Cloud::keyT> records;
        records+=cloud.store(world,v);
        if (function_is_assigned) records+=cloud.store(world,functions[0]);
        if (constant_part.is_initialized()) records+=cloud.store(world,constant_part);
        return records;
   }

    /// customized function to load this from the cloud

    /// functions and constant_part can be very large and we want to split them and store them in different records
    /// @param[inout] recordlist: containing the keys of the member variables -> will be reduced by the keys which are used
    void cloud_load(World& world, const Cloud& cloud, Recordlist<Cloud::keyT>& recordlist) {
        // load bookkeeping stuff in a vector
        std::vector<unsigned char> v=cloud.forward_load<std::vector<unsigned char>>(world,recordlist);
        archive::VectorInputArchive arin(v);
        bool function_is_assigned = false, constant_part_is_initialized=false;
        arin & type & ctype & i & j & bsh_eps & function_is_assigned & constant_part_is_initialized;
        functions.clear();
        constant_part.clear();

        if (function_is_assigned) functions.emplace_back(cloud.forward_load<CCPairFunction<double,6>>(world,recordlist));
        if (constant_part_is_initialized) constant_part=cloud.forward_load<real_function_6d>(world,recordlist);
   }

    bool function_exists() const {
        return (functions.size()>0 and functions[0].is_assigned() and functions[0].is_pure());
    }

    /// gives back the pure 6D part of the pair function
    real_function_6d function() const {
        MADNESS_CHECK_THROW(not functions.empty(), "no function assigned in CCPair::function()");
        MADNESS_CHECK_THROW(functions[0].is_pure(),"function is not pure in CCPair::function()");
        return functions[0].get_function();
    }

    /// updates the pure 6D part of the pair function
    void update_u(const real_function_6d& u) {
        // print("updating u(",i,j,")");
        CCPairFunction tmp(u);
        if (functions.size() == 0) functions.push_back(tmp);
        else { //(functions.size() > 1) {
            MADNESS_CHECK_THROW(functions[0].is_pure(),"function is not pure in CCPair::update_u()");
            functions[0]=tmp;
        }
    }

    template<typename Archive>
    void serialize(const Archive& ar) {
        size_t f_size = functions.size();
        bool fexist = (f_size > 0) && (functions[0].get_function().is_initialized());
        bool cexist = constant_part.is_initialized();
        ar & type & ctype & i & j & bsh_eps & fexist & cexist & f_size;
        if constexpr (Archive::is_input_archive) {
            if (fexist) {
                real_function_6d func;
                ar & func;
                CCPairFunction f1(func);
                functions.push_back(f1);
            }
        } else {
            if (fexist) ar & functions[0].get_function();
        }
        if (cexist) ar & constant_part;
    }

    /// reconstruct constant part and all functions
    void reconstruct() const {
        constant_part.reconstruct();
        for (auto& f : functions) {
            if (f.is_assigned() and f.is_pure()) f.get_function().reconstruct();
        }
    }

    bool load_pair(World& world, const bool verbose=false) {
        std::string fname=this->name();
        if (verbose and world.rank()==0) print("loading pair from file", fname);
        bool exists = archive::ParallelInputArchive<archive::BinaryFstreamInputArchive>::exists(world, fname.c_str());
        if (exists) {
            archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, fname.c_str(), 1);
            ar & *this;
            if (functions[0].get_function().is_initialized()) functions[0].get_function().set_thresh(FunctionDefaults<6>::get_thresh());
            if (constant_part.is_initialized()) constant_part.set_thresh(FunctionDefaults<6>::get_thresh());
        }
        return exists;
    }

    void store_pair(World& world, const bool verbose=false) {
        std::string fname =this->name();
        if (verbose and world.rank()==0) print("loading pair from file", fname);
        this->reconstruct();
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, fname.c_str(), 1);
        ar & *this;
    }

    hashT hash() const {
        hashT hash_i = std::hash<std::size_t>{}(i);
        hash_combine(hash_i, std::hash<std::size_t>{}(j));
        if (constant_part.is_initialized()) {
            hash_combine(hash_i, hash_value(constant_part.get_impl()->id()));
        }
        return hash_i;
    }

    /// the functions which belong to the pair
    std::vector<CCPairFunction<double,6>> functions;

    /// the constant part
    real_function_6d constant_part;

    /// Energy for the BSH Operator
    /// Ground State: e_i + e_j
    /// Excited State: e_i + e_j + omega
    /// default to positive value to make sure this is set somewhere
    double bsh_eps=1.0;

    /// return the base name like "MP2_pair_u" or "CC2_pair_x"
    std::string basename() const {
        std::string name = "???";
        if (type == GROUND_STATE) name = assign_name(ctype) + "_pair_u";
        if (type == EXCITED_STATE) name = assign_name(ctype) + "_pair_x";
        return name;

    }
    std::string name() const {
        return basename() +"_" + stringify(i) + stringify(j);
    }

    void
    info() const;

};

/// build an MP2 or CC2 or LRCC2 etc pair, possibly including the lo-rank parts
class CCPairBuilder {
public:
    CCPairBuilder(World& world, const Info& info) : world(world), info(info) {};

    /// provide ground-state singles, needed for CC2 and LRCC2 wave function ansatz
    CCPairBuilder& set_gs_singles(const CC_vecfunction& gs) {
        gs_singles = gs;
        return *this;
    }

    /// provide excited state singles, needed for CC2 and LRCC2 wave function ansatz
    CCPairBuilder& set_ex_singles(const CC_vecfunction& ex) {
        ex_singles = ex;
        return *this;
    }

    CCPairBuilder& set_ctype(const CalcType& type) {
        ctype = type;
        return *this;
    }

    /// make a CCPair without the 6d function and some bookkeeping information
    CCPair make_bare_pair(const int i, const int j) const;

    /// make a CCPair with the 6d function only and some bookkeeping information
    CCPair make_bare_pair_from_file(const int i, const int j) const;

    /// make a CCPair
    CCPair make_pair(const int i, const int j, const std::vector<CCPairFunction<double,6>>& u) const;

    inline static CCState cc_state(const CalcType& type) {
        if (type==CT_MP2 or type==CT_CC2 or type==CT_MP3) return GROUND_STATE;
        else if (type==CT_LRCC2 or type==CT_ADC2 or type==CT_CISPD) return EXCITED_STATE;
        else {
            MADNESS_EXCEPTION("unknown cc-state",1);
        }
    }

    /// make a CCPair with the 6d function only and some bookkeeping information
    Pairs<CCPair> make_all_bare_pairs() const {
        Pairs<CCPair> pairs;
        for (size_t i = info.parameters.freeze(); i < info.mo_bra.size(); i++) {
            for (size_t j = info.parameters.freeze(); j < info.mo_ket.size(); j++) {
                pairs.insert(i,j,make_bare_pair(i, j));
            }
        }
        return pairs;
    }

    /// complete the given pair with the low-rank parts

    /// will use pair's ctype, while builder's ctype is ignored
    CCPair complete_pair_with_low_rank_parts(const CCPair& pair) const;


    /// Function to load a function from disk

    /// @param[in] name of the file in which the function was stored
    /// @param do_print
    /// @return the function, possibly not initialized if not found on disk
    template <typename T, size_t NDIM>
    Function<T,NDIM> load_function(const std::string name, bool do_print) const {
        Function<T,NDIM> f;
        bool exists = archive::ParallelInputArchive<
            archive::BinaryFstreamInputArchive>::exists(world, name.c_str());
        if (exists) {
            if ((world.rank() == 0) and do_print) print("loading function", name);
            archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str());
            ar & f;
            if (do_print) f.print_size(name);
            if (f.is_compressed()) {
                if (world.rank()==0 and do_print) print("function is compressed -- reconstructing");
                f.change_tree_state(reconstructed);
                if (do_print) f.print_size(name+" reconstructed");
                save(f, name);
            }
            f.set_thresh(FunctionDefaults<NDIM>::get_thresh());
            f.truncate();
            f.print_size(name);
        } else {
            if ((world.rank()==0) and do_print) print("could not find function",name);
        }
        return f;
    }

    World& world;
    const Info& info;
    CC_vecfunction gs_singles, ex_singles;
    CalcType ctype=CT_UNDEFINED;

};


/// print accumulated size of all functions
struct CCSize {
    double size_local=0;

    CCSize() = default;

    template<typename T, std::size_t NDIM>
    void add_helper(const std::vector<Function<T,NDIM>>& v) {
        if (v.size()>0) size_local+=get_size_local(v.front().world(),v);
    }

    void add_helper(const std::vector<CCPair>& vp) {
        if (vp.empty()) return;
        for (const auto& p : vp) {
            size_local+=get_size(p.constant_part);
            if (p.function_exists()) size_local+=get_size_local(p.function());
        }
    }

    /// variadic template parameters to add the size of all functions and pairs
    template<typename... Args>
    void add(const Args&... args) {
        (add_helper(args), ...);
    }

    void print(World& world, const std::string msg="") const {
        double size_global=size_local;
        world.gop.sum(size_global);
        if (msg.size()>0 and world.rank()==0) madness::print(msg);
        world.gop.fence();
        madness::print("size of all functions on rank",world.rank(),size_local);
        world.gop.fence();
        if (world.rank()==0) madness::print("total size of all functions",size_global);

    }
};


class MacroTaskMp2ConstantPart : public MacroTaskOperationBase {

    class ConstantPartPartitioner : public MacroTaskPartitioner {
    public:
        ConstantPartPartitioner() {};

        partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                   const std::string policy) const override {
            partitionT p;
            for (size_t i = 0; i < vsize1; i++) {
                Batch batch(Batch_1D(i,i+1), Batch_1D(i,i+1));
                p.push_back(std::make_pair(batch,1.0));
            }
            return p;
        }
    };

public:
    MacroTaskMp2ConstantPart(){partitioner.reset(new ConstantPartPartitioner());}

    // typedef std::tuple<const std::vector<CCPair>&, const std::vector<Function<double,3>>&,
            // const std::vector<Function<double,3>>&, const CCParameters&, const Function<double,3>&,
            // const std::vector<Function<double,3>>&, const std::vector<std::string>& > argtupleT;
    typedef std::tuple<const std::vector<CCPair>&, const madness::Info&, const std::vector<std::string>& > argtupleT;

    using resultT = std::vector<real_function_6d>;

    resultT allocator(World& world, const argtupleT& argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        resultT result = zero_functions_auto_tree_state<double, 6>(world, n);
        return result;
    }

//    resultT operator() (const std::vector<CCPair>& pair, const std::vector<Function<double,3>>& mo_ket,
//                        const std::vector<Function<double,3>>& mo_bra, const CCParameters& parameters,
//                        const Function<double,3>& Rsquare, const std::vector<Function<double,3>>& U1,
//                        const std::vector<std::string>& argument) const;
    resultT operator() (const std::vector<CCPair>& pair, const Info& info, const std::vector<std::string>& argument) const;
};

/// compute the "constant" part of MP2, CC2, or LR-CC2
///
/// the constant part is
/// result = G [F,f] |ij>  for MP2
/// result = G [F,f] |t_i t_j>  for CC2
/// result = G [F,f] |t_i x_j> + |x_i t_j>  for LR-CC2
class MacroTaskConstantPart : public MacroTaskOperationBase {

    class ConstantPartPartitioner : public MacroTaskPartitioner {
    public:
        ConstantPartPartitioner() {};

        partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                   const std::string policy) const override {
            partitionT p;
            for (size_t i = 0; i < vsize1; i++) {
                Batch batch(Batch_1D(i,i+1), Batch_1D(i,i+1));
                p.push_back(std::make_pair(batch,1.0));
            }
            return p;
        }
    };

public:
    MacroTaskConstantPart()  {
        partitioner.reset(new ConstantPartPartitioner());
        name="ConstantPart";
    }

    // typedef std::tuple<const std::vector<CCPair>&, const std::vector<Function<double,3>>&,
    // const std::vector<Function<double,3>>&, const CCParameters&, const Function<double,3>&,
    // const std::vector<Function<double,3>>&, const std::vector<std::string>& > argtupleT;
    typedef std::tuple<const std::vector<CCPair>&,
                       const std::vector<Function<double,3>>&, const std::vector<Function<double,3>>&,
                       const madness::Info&> argtupleT;

    using resultT = std::vector<real_function_6d>;

    resultT allocator(World& world, const argtupleT& argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        resultT result = zero_functions_auto_tree_state<double, 6>(world, n);
        return result;
    }
    resultT operator() (const std::vector<CCPair>& pair,
        const std::vector<Function<double,3>>& gs_singles,
        const std::vector<Function<double,3>>& ex_singles,
        const Info& info) const;
};

class MacroTaskMp2UpdatePair : public MacroTaskOperationBase {

    class UpdatePairPartitioner : public MacroTaskPartitioner {
    public :
        UpdatePairPartitioner() {
            set_dimension(2);
        }

        partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                   const std::string policy) const override {
            partitionT p;
            for (size_t i = 0; i < vsize1; i++) {
                Batch batch(Batch_1D(i, i+1), Batch_1D(i, i+1), Batch_1D(i,i+1));
                p.push_back(std::make_pair(batch, 1.0));
            }
            return p;
        }
    };
public:
    MacroTaskMp2UpdatePair() {
        partitioner.reset(new UpdatePairPartitioner());
        name="MP2UpdatePair";
    }

    // typedef std::tuple<const std::vector<CCPair>&, const std::vector<real_function_6d>&, const CCParameters&,
                        // const std::vector< madness::Vector<double,3> >&,
                       // const std::vector<Function<double,3>>&, const std::vector<Function<double,3>>&,
                       // const std::vector<Function<double,3>>&, const Function<double,3>&> argtupleT;
    typedef std::tuple<const std::vector<CCPair>&, const std::vector<real_function_6d>&,
                    const std::vector<madness::Vector<double,3>>&, const Info& > argtupleT;

    using resultT = std::vector<real_function_6d>;

    resultT allocator(World& world, const argtupleT& argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        resultT result = zero_functions_auto_tree_state<double, 6>(world, n);
        return result;
    }

//    resultT operator() (const std::vector<CCPair>& pair, const std::vector<real_function_6d>& mp2_coupling, const CCParameters& parameters,
//                        const std::vector< madness::Vector<double,3> >& all_coords_vec,
//                        const std::vector<Function<double,3>>& mo_ket, const std::vector<Function<double,3>>& mo_bra,
//                        const std::vector<Function<double,3>>& U1, const Function<double,3>& U2) const;
    resultT operator() (const std::vector<CCPair>& pair, const std::vector<real_function_6d>& mp2_coupling,
                        const std::vector< madness::Vector<double,3> >& all_coords_vec, const Info& info) const;
};


class MacroTaskIteratePair : public MacroTaskOperationBase {

    class IteratePairPartitioner : public MacroTaskPartitioner {
    public :
        IteratePairPartitioner() = default;

        partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                   const std::string policy) const override {
            partitionT p;
            for (size_t i = 0; i < vsize1; i++) {
                Batch batch(Batch_1D(i, i+1), Batch_1D(i, i+1), Batch_1D(i,i+1));
                p.push_back(std::make_pair(batch,1.0));
            }
            return p;
        }
    };
public:
    MacroTaskIteratePair() {
        partitioner.reset(new IteratePairPartitioner());
        name="IteratePair";
    }

    typedef std::tuple<
        const std::vector<CCPair>&,      // pair
        const std::vector<real_function_6d>&,   // local coupling
        const CC_vecfunction&,          // gs singles
        const CC_vecfunction&,          // ex singles
        const Info&,
        const std::size_t&
        > argtupleT;

    using resultT = std::vector<real_function_6d>;

    resultT allocator(World& world, const argtupleT& argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        resultT result = zero_functions_auto_tree_state<double, 6>(world, n);
        return result;
    }

    /// iterate a given pair of the MP2, CC2 or LRCC2 calculation

    /// will *NOT* compute the local coupling,
    /// will apply the Fock operators (J-K+V)|pair> and use
    /// the (excited) singles vectors to update the pair
    /// @param[in] pair: the pair which will be updated
    /// @param[in] gs_singles: the ground state singles, may be dummy for MP2
    /// @param[in] ex_singles: the excited state singles, may be dummy for MP2, CC2
    /// @param[in] all_coords_vec: the coordinates of the atoms
    /// @param[in] info: the info structure
    /// @param[in] maxiter: the maximal number of iterations
    resultT operator() (const std::vector<CCPair>& pair,
        const std::vector<real_function_6d>& local_coupling,
        const CC_vecfunction& gs_singles,
        const CC_vecfunction& ex_singles,
        const Info& info,
        const std::size_t& maxiter) const;
};


class MacroTaskSinglesPotentialEx : public MacroTaskOperationBase {
public:
    std::string basename="SinglesPotentialEx";
    MacroTaskSinglesPotentialEx() {
        name="SinglesPotentialEx";
        partitioner->max_batch_size=2;
        partitioner->min_batch_size=2;
    }

    typedef std::tuple<
        const std::vector<int>&,    // result_index,
        const CC_vecfunction&,      // singles_gs,
        const std::vector<CCPair>&,       // doubles_gs,
        const CC_vecfunction&,      // singles_ex,
        const std::vector<CCPair>&,       // doubles_ex,
        const int&,       // name,
        const Info&                 // info
    > argtupleT;

    using resultT = std::tuple<std::vector<real_function_3d>,std::vector<real_function_3d>>;

    resultT allocator(World& world, const argtupleT& argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        std::vector<real_function_3d> result = zero_functions_auto_tree_state<double, 3>(world, n);
        std::vector<real_function_3d> intermediate = zero_functions_auto_tree_state<double, 3>(world, n);
        const_cast<std::string&>(name) =basename+"_"+assign_name(PotentialType(std::get<5>(argtuple)));
        return std::make_tuple(result,intermediate);
    }

    resultT operator() (const std::vector<int>& result_index,
                        const CC_vecfunction& singles_gs,
                        const std::vector<CCPair>& doubles_gs,
                        const CC_vecfunction& singles_ex,
                        const std::vector<CCPair>& doubles_ex,
                        const int& name,
                        const Info& info);
};

class MacroTaskSinglesPotentialGs : public MacroTaskOperationBase {
public:
    std::string basename="SinglesPotentialGs";
    MacroTaskSinglesPotentialGs() {
        name="SinglesPotentialGs";
    }

    typedef std::tuple<
        const std::vector<int>&,    // result_index,
        const CC_vecfunction&,      // singles_gs,
        const std::vector<CCPair>&,       // doubles_gs,
        const int&,       // name,
        const Info&                 // info
    > argtupleT;

    /// first vector is the potential, second is an intermediate (if applicable, e.g. for s2b and s2c potentials)
    using resultT = std::tuple<std::vector<real_function_3d>,std::vector<real_function_3d>>;

    /// allocate the result and set the name of this task
    resultT allocator(World& world, const argtupleT& argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        std::vector<real_function_3d> result = zero_functions_auto_tree_state<double, 3>(world, n);
        std::vector<real_function_3d> intermediate = zero_functions_auto_tree_state<double, 3>(world, n);
        const_cast<std::string&>(name) =basename+"_"+assign_name(PotentialType(std::get<3>(argtuple)));
        return std::make_tuple(result,intermediate);
    }

    resultT operator() (const std::vector<int>& result_index,
                        const CC_vecfunction& singles_gs,
                        const std::vector<CCPair>& doubles_gs,
                        const int& name,
                        const Info& info);
};


class MacroTaskComputeCorrelationEnergy : public MacroTaskOperationBase {
public:
    std::string basename="CorrelationEnergy";
    MacroTaskComputeCorrelationEnergy() {
        name="CorrelationEnergy";
    }

    typedef std::tuple<
        const std::vector<CCPair>&,
        const CC_vecfunction&,
        const Info&
    > argtupleT;

    /// first vector is the potential, second is an intermediate (if applicable, e.g. for s2b and s2c potentials)
    typedef std::vector<ScalarResult<double>> resultT;


    /// allocate the result and set the name of this task
    resultT allocator(World &world, const argtupleT &argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        return scalar_result_vector<double>(world,n);
    }

    resultT operator() (const std::vector<CCPair>& pairs,
                        const CC_vecfunction& singles_gs,
                        const Info& info) const;
};

}//namespace madness

#endif /* CCSTRUCTURES_H_ */
