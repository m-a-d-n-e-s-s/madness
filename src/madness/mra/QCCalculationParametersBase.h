/*
 * QCCalculationParametersBase.h
 *
 *  Created on: 27 Jun 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_QCCALCULATIONPARAMETERSBASE_H_
#define SRC_APPS_CHEM_QCCALCULATIONPARAMETERSBASE_H_

#include "madness/external/nlohmann_json/json.hpp"
#include "madness/misc/misc.h"
#include "madness/mra/commandlineparser.h"
#include "madness/world/archive.h"
#include "madness/world/world.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <typeindex>
#include <typeinfo>


namespace madness {

    using json = nlohmann::json;

    template<typename T>
    static typename std::enable_if<std::is_floating_point<T>::value, void>::type check_for_inf(const std::string &str, T &arg) {
        std::string sinf;
        std::stringstream ss;
        ss << std::numeric_limits<T>::infinity();
        sinf = ss.str();

        if (sinf == str) arg = std::numeric_limits<T>::infinity();
    }

    template<typename T>
    static typename std::enable_if<!std::is_floating_point<T>::value, void>::type check_for_inf(const std::string &str, T &arg) {
        return;
    }

    /// inverting the print method from print.h for std::vector
    /// TODO: move this where it belongs (into print.h ??)
    template<typename T, typename A = std::allocator<T>>
    std::istream &operator>>(std::istream &is, std::vector<T, A> &v) {

        // get the full line from opening to closing brackets [ .. ]
        std::string word, line = "";
        while (is >> word) {
            line += word;
            if (word.find(']') != std::string::npos) break;
        }
        if (line.size() != 0) is.clear();

        // remove enclosing brackets and commas
        auto find_c = [](char &c) { return ((c == ',') or (c == '[') or (c == ']')); };
        std::replace_if(line.begin(), line.end(), find_c, ' ');// 0 2 0 4 0 6 0 8 0

        // stream the values into the container
        std::stringstream sline(line);
        T tmp;
        while (sline >> word) {
            std::stringstream sword(word);
            sword >> tmp;
            check_for_inf(word, tmp);
            v.push_back(tmp);
        }
        if (sline.bad()) {
            madness::print("error while reading vector from istream: ");
            madness::print(line, "\n");
            throw std::runtime_error("IO error");
        }

        return is;
    }


    /// inverting the print method from print.h for std::vector
    /// TODO: move this where it belongs (into print.h ??)
    template<typename Q, typename T>
    std::istream &operator>>(std::istream &is, std::pair<T, Q> &p) {

        // get all words from the line
        std::string word, line = "";
        while (is >> word) { line += word + " "; };
        // have to set this here to account for the error handling later..
        is.clear();

        // remove enclosing brackets and commas
        auto find_c = [](char &c) { return ((c == ',') or (c == '(') or (c == ')')); };
        std::replace_if(line.begin(), line.end(), find_c, ' ');// 0 2 0 4 0 6 0 8 0

        // stream the values into the container
        std::stringstream sline(line);
        T tmp1;
        Q tmp2;
        sline >> tmp1 >> tmp2;
        if (sline.bad() or sline.fail()) {
            madness::print("error while reading vector from istream: ");
            madness::print(line, "\n");
            throw std::runtime_error("IO error");
        }
        p = std::pair<T, Q>(tmp1, tmp2);

        return is;
    }


    /// structure holding the value for a given parameter

    /// keeps logic about default, derived and user-defined settings (with increasing priority),
    /// as well as comments for the user (e.g. a recommended range for a given parameter).
    ///
    /// might be extended to hold allowed values for certain parameters, e.g. localization procedures.
    ///
    /// all values are stored as strings and must be converted to their respective types by the
    /// QCCalculationParametersBase class (see below)
    struct QCParameter {
    public:
        QCParameter(){};

        QCParameter(const std::string v, const std::string t, const std::string comment = "", const std::vector<std::string> allowed_values1 = {})
            : default_value(v), type(t), comment(comment), allowed_values(allowed_values1) {
            static int i = 0;
            print_order = i++;
            set_all();
        }

        void set_derived_value(const std::string val) {
            precedence = std::max(derived, precedence);
            derived_value = val;
            set_all();
        }
        void set_user_defined_value(const std::string val) {
            precedence = std::max(defined, precedence);
            user_defined_value = val;
            set_all();
        }

        bool is_user_defined() const { return (precedence == defined); }

        std::string get_value() const { return value; }
        std::string get_type() const { return type; }
        std::string get_comment() const { return comment; }
        std::string print_precedence() const {
            if (precedence == def) return "default";
            if (precedence == derived) return "derived";
            if (precedence == defined) return "defined";
            std::stringstream sprecedence;
            sprecedence << precedence;
            throw std::runtime_error("unknown precedence in QCParameter" + sprecedence.str());
            return "darn";
        }

        int get_print_order() const { return print_order; }

        std::string print_line(const std::string &key) const {

            auto fill_left = [](const int size, const std::string word) {
                int nspaces = std::max(int(0), size - int(word.length()));
                return std::string(nspaces, ' ') + word;
            };
            auto fill_right = [](const int size, const std::string word) {
                int nspaces = std::max(int(0), size - int(word.length()));
                return word + std::string(nspaces, ' ');
            };

            // key-value block
            std::string keyval = fill_left(20, key) + "  " + fill_right(10, get_value()) + " # " + fill_right(10, print_precedence());
            std::string empty_keyval(keyval.size(), ' ');
            empty_keyval[33] = '#';

            std::string allowed_val;
            if (allowed_values.size() > 0) {
                using madness::operators::operator<<;
                std::stringstream ss;
                ss << allowed_values;
                allowed_val += ss.str();
            }

            // split comment into several lines: split onto words and add linebreak
            bool leave_space_for_allowed_values = (allowed_val.size() > 0);

            // first line breaks after 80 characters, all other lines after 120 (leave space f
            long keyvalsize = keyval.size();// length of key, value, precedence
            std::string comment1 = get_comment();
            auto commentwords = commandlineparser::split(comment1, " ");
            std::vector<std::string> commentlines(1);
            long nchar = 0;
            for (auto word: commentwords) {

                bool is_first_line = commentlines.size() == 1;
                long thislinebreak = 120;
                if (is_first_line and leave_space_for_allowed_values) thislinebreak = 80;
                long commentsize = thislinebreak - keyvalsize;

                nchar += word.size() + 1;
                if (nchar > commentsize) {// start newline
                    commentlines.push_back("");
                    nchar = word.size() + 1;
                }
                commentlines.back() += word + " ";
            }

            std::string result;
            for (size_t i = 0; i < commentlines.size(); ++i) {
                if (i == 0) result = keyval + fill_right(40, commentlines[i]) + allowed_val;
                else
                    result += "\n" + empty_keyval + commentlines[i];
            }

            // trim result
            std::size_t last = result.find_last_not_of(' ');
            return result.substr(0, last + 1);
        }

        template<typename Archive>
        void serialize(Archive &ar) {
            ar & value & default_value & derived_value & user_defined_value & type & null & comment & allowed_values & print_order & precedence;
        }

        enum { def, derived, defined } precedence = def;

        hashT hash() const { return hash_value(value); }

    private:
        void set_all() {
            value = default_value;
            if (derived_value != null) value = derived_value;
            if (user_defined_value != null) value = user_defined_value;
            if (not check_allowed()) throw std::invalid_argument(not_allowed_errmsg());
        }

        bool check_allowed() {
            if (allowed_values.size() == 0) return true;
            auto it = std::find(allowed_values.begin(), allowed_values.end(), value);
            return (it != allowed_values.end());
        }

        std::string not_allowed_errmsg() const {
            using madness::operators::operator<<;
            std::stringstream ss;
            ss << allowed_values;
            std::string errmsg = "\ntrying to assign a value that's not allowed\n\n";
            errmsg += "\tuser-defined value: " + value + "\n";
            errmsg += "\tallowed values:     " + ss.str() + "\n\n";
            return errmsg;
        }


        std::string value;
        std::string default_value = "";
        std::string derived_value = "";
        std::string user_defined_value = "";
        std::string type = "";
        std::string null = "";
        std::string comment = "";
        std::vector<std::string> allowed_values = std::vector<std::string>();
        int print_order = 0;// use this for printing the parameters in the same order as they are defined
    };

    /// class for holding the parameters for calculation

    /// Actual parameter classes will be derived from this class with a simple constructor
    /// (see test_QCCalculationParametersBase.cc for an example) and convenience
    /// getters for the parameters of each parameter class.
    /// Having the base class will allow consistent parameter input/output handling
    /// for all madness programs and reuse of the parsing methods.
    /// Even if the same parameter is used in different programs, default might differ
    /// (i.e. econv for 3D/6D calculations). The parameter class will effectively serve as
    /// a factory for the calculation classes (SCF, nemo, mp2, etc)
    ///
    /// parameters are kept in a map with key (string) and value (QCParameter),
    /// types are converted whenever a parameter is accessed (i.e. should not be
    /// done too frequently in an inner loop)
    class QCCalculationParametersBase {

    public:
        /// print all parameters
        void print(const std::string header = "", const std::string footer = "") const;

        std::string print_to_string(bool non_defaults_only = false) const;

        template<typename T>
        T get(const std::string key) const {
            const QCParameter &parameter = get_parameter(key);
            MADNESS_ASSERT(check_type<T>(key, parameter));
            if (std::is_same<T, std::string>::value) { return fromstring<T>(add_quotes(parameter.get_value())); }
            return fromstring<T>(parameter.get_value());
        }

        bool is_user_defined(std::string key) const { return get_parameter(key).is_user_defined(); }

        template<typename Archive>
        void serialize(Archive &ar) {
            ar & parameters & print_debug;
        }

        hashT hash() const { return hash_range(parameters.begin(), parameters.end()); }

    protected:
        typedef std::map<std::string, QCParameter> ParameterContainerT;
        ParameterContainerT parameters;

        virtual void read_input_and_commandline_options(World &world, const commandlineparser &parser, const std::string tag) {
            try {
                // check that user-defined input files actually exist
                bool file_ok = true;
                if (parser.key_exists("user_defined_input_file")) file_ok = file_exists(world, parser.value("input"));
                if (file_ok) read_input(world, parser.value("input"), tag);
                else {
                    std::string msg = "could not find user-defined input file: " + parser.value("input") + "\n";
                    throw std::invalid_argument(msg);
                }
            } catch (std::invalid_argument &e) { throw; } catch (std::exception &e) {
                madness::print(e.what());
            }
            read_commandline_options(world, parser, tag);
        }

    public:
        bool file_exists(World &world, std::string filename) const;

    private:
        /// read the parameters from file

        /// only world.rank()==0 reads the input file and broadcasts to all other nodes,
        /// so we don't need to serialize the ParameterMap
        void read_input(World &world, const std::string filename, const std::string tag);

        void read_commandline_options(World &world, const commandlineparser &parser, const std::string tag);

    protected:
        bool print_debug = false;
        bool ignore_unknown_keys = true;
        bool ignore_unknown_keys_silently = false;
        bool throw_if_datagroup_not_found = true;

        /// ctor for testing
        QCCalculationParametersBase() {}

        /// copy ctor
        QCCalculationParametersBase(const QCCalculationParametersBase &other) : parameters(other.parameters), print_debug(other.print_debug) {}

        /// destructor
        virtual ~QCCalculationParametersBase() {}

        template<typename T>
        void initialize(const std::string &key, const T &value, const std::string comment = "", const std::vector<T> allowed_values = {}) {

            if (parameters.find(key) != parameters.end()) {
                madness::print("you cannot initialize a parameter twice: ", key);
                throw std::runtime_error("initialization error");
            }

            std::string svalue = tostring(value);
            std::string type = std::type_index(typeid(T)).name();

            // transform everything to lower case
            std::string key_lower = key;
            std::transform(key_lower.begin(), key_lower.end(), key_lower.begin(), ::tolower);
            std::transform(svalue.begin(), svalue.end(), svalue.begin(), ::tolower);
            std::vector<std::string> av_lower_vec;
            for (auto av: allowed_values) {
                std::string av_lower = tostring(av);
                std::transform(av_lower.begin(), av_lower.end(), av_lower.begin(), ::tolower);
                av_lower_vec.push_back(av_lower);
            }

            parameters.insert(std::make_pair<std::string, QCParameter>(std::string(key_lower), QCParameter(svalue, type, comment, av_lower_vec)));
        }

    public:
        template<typename T>
        void set_derived_value(const std::string &key, const T &value) {

            QCParameter &parameter = get_parameter(key);
            if (not check_type_silent<T>(parameter)) { throw std::runtime_error("type error in set_derived_value for key " + key); }
            parameter.set_derived_value(tostring(value));
        }

        ParameterContainerT get_all_parameters() const { return parameters; }


    protected:
        template<typename T>
        bool try_setting_user_defined_value(const std::string &key, const std::string &val) {

            if (not check_type_silent<T>(get_parameter(key))) return false;

            if (print_debug) ::madness::print("key:", key, "will set type", std::type_index(typeid(T)).name());
            T value = fromstring<T>(val);
            set_user_defined_value<T>(key, value);
            return true;
        }

    public:
        void from_json(const json &j) {
            for (const auto & [key, value]: j.items()) {
                QCParameter &parameter = get_parameter(key);
                // ::print("key: ", key, " value: ", value);
                parameter.set_user_defined_value(tostring(value));
            }
        }
        [[nodiscard]] json to_json() const {
            json j_params = {};
            // TODO Is there a way to the get member for every parameter even though get is a template function?
            for (auto &p: parameters) {
                auto param_type = p.second.get_type();
                if (param_type == "i") {
                    j_params[p.first] = get<int>(p.first);
                    // if vector of double
                } else if (param_type == "d") {
                    j_params[p.first] = get<double>(p.first);

                    // if vector of bool
                } else if (param_type == "b") {
                    j_params[p.first] = get<bool>(p.first);

                    // if vector of doubles?
                } else if (param_type == "St6vectorIdSaIdEE") {
                    j_params[p.first] = get<std::vector<double>>(p.first);
                } else if (param_type == "NSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE") {
                    auto sval = get<std::string>(p.first);
                    if (!sval.empty()) continue;
                    j_params[p.first] = sval;
                    // size t
                } else if (p.second.get_type() == "m") {
                    j_params[p.first] = get<size_t>(p.first);
                }
            }
            return j_params;
        }

        [[nodiscard]] json to_json_if_precedence(const std::string &precedence) const {
            json j_params = {};
            for (auto &p: parameters) {
                auto param_type = p.second.get_type();
                if (p.second.print_precedence() == precedence) {

                    if (param_type == "i") {
                        j_params[p.first] = get<int>(p.first);
                        // if vector of double
                    } else if (param_type == "d") {
                        j_params[p.first] = get<double>(p.first);

                        // if vector of bool
                    } else if (param_type == "b") {
                        j_params[p.first] = get<bool>(p.first);

                        // if vector of doubles?
                    } else if (param_type == "St6vectorIdSaIdEE") {
                        j_params[p.first] = get<std::vector<double>>(p.first);
                    } else if (param_type == "NSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE") {
                        auto sval = get<std::string>(p.first);
                        if (!sval.empty()) continue;
                        j_params[p.first] = sval;
                        // size t
                    } else if (p.second.get_type() == "m") {
                        j_params[p.first] = get<size_t>(p.first);
                    }
                }
            }
            return j_params;
        }
        /**
         * Adds the response parameters to an existing json under key "parameters"
         * @param j
         */
        void to_json(json &j) const {
            json j_params = {};
            // TODO Is there a way to the get member for every parameter even though get is a template function?
            for (auto &p: parameters) {
                auto param_type = p.second.get_type();
                if (param_type == "i") {
                    j_params[p.first] = get<int>(p.first);
                    // if vector of double
                } else if (param_type == "d") {
                    j_params[p.first] = get<double>(p.first);

                    // if vector of bool
                } else if (param_type == "b") {
                    j_params[p.first] = get<bool>(p.first);

                    // if vector of doubles?
                } else if (param_type == "St6vectorIdSaIdEE") {
                    j_params[p.first] = get<std::vector<double>>(p.first);
                } else if (param_type == "NSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE") {
                    auto sval = get<std::string>(p.first);
                    if (!sval.empty()) continue;
                    j_params[p.first] = sval;
                    // size t
                } else if (p.second.get_type() == "m") {
                    j_params[p.first] = get<size_t>(p.first);
                }
            }
            j["parameters"] = j_params;
        }

        bool operator==(const QCCalculationParametersBase &other) const {

            for (auto &p: parameters) {
                auto param_type = p.second.get_type();

                if (param_type == "i") {
                    if (get<int>(p.first) != other.get<int>(p.first)) { return false; };
                    // if vector of double
                } else if (param_type == "d") {
                    if (get<double>(p.first) != other.get<double>(p.first)) { return false; };
                    // if vector of bool
                } else if (param_type == "b") {
                    if (get<bool>(p.first) != other.get<bool>(p.first)) { return false; };
                    // if vector of doubles?
                } else if (param_type == "St6vectorIdSaIdEE") {
                    if (get<std::vector<double>>(p.first) != other.get<std::vector<double>>(p.first)) { return false; };
                } else if (param_type == "NSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE") {
                    if (get<std::string>(p.first) != other.get<std::string>(p.first)) { return false; };
                } else if (p.second.get_type() == "m") {
                    if (get<size_t>(p.first) != other.get<size_t>(p.first)) { return false; };
                }
            }
            return true;
        }

        template<typename T>
        void set_user_defined_value(const std::string &key, const T &value) {

            QCParameter &parameter = get_parameter(key);
            if (not check_type_silent<T>(parameter)) { throw std::runtime_error("type error in set_user_defined_value"); }

            parameter.set_user_defined_value(tostring(value));
        }

    protected:
        const QCParameter &get_parameter(const std::string & key) const {
            if (not parameter_exists(key)) { throw std::runtime_error("could not find parameter for key " + key); }
            const QCParameter &parameter = parameters.find(key)->second;
            return parameter;
        }

    public:
        QCParameter &get_parameter(const std::string key) {
            if (not parameter_exists(key)) {
                madness::print("\ncould not find parameter for key", key, "\n");
                throw std::runtime_error("could not find parameter for key " + key);
            }
            QCParameter &parameter = parameters.find(key)->second;
            return parameter;
        }

        bool parameter_exists(const std::string &key) const { return (parameters.find(key) != parameters.end()); }


        template<typename T>
        static bool check_type(const std::string key, const QCParameter &parameter) {
            if (check_type_silent<T>(parameter)) return true;

            madness::print("trying to get the wrong type in QCCalculationParametersBase");
            madness::print("key             ", key);
            madness::print("parameter type  ", parameter.get_type());
            madness::print("setting type    ", std::type_index(typeid(T)).name());
            madness::print("value           ", parameter.get_value());
            return false;
        }

        template<typename T>
        static bool check_type_silent(const QCParameter &parameter) {
            return (parameter.get_type() == std::type_index(typeid(T)).name());
        }

        /// read the stream, starting from tag

        /// only parameters that are defined in the constructor will be processed,
        /// all others will be discarded.
        virtual void read_internal(World &world, std::string &filecontents, std::string tag);


        static std::string tostring(const bool &arg) {
            std::ostringstream ss;
            ss << std::boolalpha << arg;
            return ss.str();
        }

        template<typename T>
        static std::string tostring(const T &arg) {
            using madness::operators::operator<<;
            std::ostringstream ss;

            ss << std::scientific << std::setprecision(4) << arg;
            std::string str = ss.str();
            std::transform(str.begin(), str.end(), str.begin(), ::tolower);

            overwrite_if_inf(str, arg);
            return str;
        }

        template<typename T>
        static typename std::enable_if<std::is_floating_point<T>::value, void>::type overwrite_if_inf(std::string &str, const T &arg) {
            if (std::isinf(arg)) {
                std::stringstream ss;
                ss << std::numeric_limits<T>::infinity();
                str = ss.str();
            }
        }

        template<typename T>
        static typename std::enable_if<!std::is_floating_point<T>::value, void>::type overwrite_if_inf(std::string &str, const T &arg) {
            return;
        }

        template<typename T>
        static typename std::enable_if<!std::is_same<T, bool>::value, T>::type fromstring(const std::string &arg) {

            std::stringstream ssvalue(arg);

            // if argument is std::string read the everything between possible double quotes
            T result = read_quotes<T>(ssvalue);

            bool type_conversion_failed = ssvalue.fail();

            // check for infinity in floating point conversions
            if (type_conversion_failed and (std::is_floating_point<T>::value)) {

                const static T inf = std::numeric_limits<T>::infinity();
                std::string sinf = tostring(inf);// repeat type conversion from above
                if (sinf == arg) result = inf;
                type_conversion_failed = false;
            }

            if (type_conversion_failed) {

                std::string errmsg = "error in type conversion for argument >> " + arg + " << to type " + std::type_index(typeid(T)).name();
                throw std::runtime_error(errmsg);
            }

            // check for trailing characters
            std::string word;
            while (ssvalue >> word) {
                std::string errmsg = "trailing characters in arguement >> " + arg + " <<";
                throw std::runtime_error(errmsg);
            }
            return result;
        }

        template<typename T>
        static typename std::enable_if<std::is_same<T, std::string>::value, T>::type read_quotes(std::stringstream &ssvalue) {
            T arg = ssvalue.str();
            T result;

            if (arg.find("\"") == std::string::npos) {// no double quotes found
                ssvalue >> result;

            } else {// found double quotes
                int counter = 0;
                while (counter < 2) {
                    T tmp;
                    ssvalue >> tmp;
                    if (ssvalue.fail()) {
                        std::string errmsg = "missing closing double quote in line >> " + arg;
                        throw std::runtime_error(errmsg);
                    }
                    result += " " + tmp;
                    counter = std::count(result.begin(), result.end(), '"');
                }

                // use only the text between the double quotes
                result = trim_blanks(trim_quotes(result));
            }
            return result;
        }

        static std::string trim_blanks(const std::string arg) {
            std::size_t first = arg.find_first_not_of(' ');
            std::size_t last = arg.find_last_not_of(' ');
            return arg.substr(first, last - first + 1);
        }

        static std::string trim_quotes(const std::string arg) {
            std::size_t first = arg.find_first_of('"');
            std::size_t last = arg.find_last_of('"');
            return arg.substr(first + 1, last - first - 1);
        }

        static std::string add_quotes(const std::string arg) { return "\"" + arg + "\""; }


        template<typename T>
        static typename std::enable_if<!std::is_same<T, std::string>::value, T>::type read_quotes(std::stringstream &ssvalue) {
            T result;
            ssvalue >> result;
            return result;
        }

        template<typename T>
        static typename std::enable_if<std::is_same<T, bool>::value, T>::type fromstring(const std::string &arg) {
            std::string str = arg;
            std::transform(str.begin(), str.end(), str.begin(), ::tolower);
            if (str == "true" or str == "1" or str == "yes") return true;
            if (str == "false" or str == "0" or str == "no") return false;
            std::string errmsg = "error in type conversion for argument >> " + arg + " << to type " + std::type_index(typeid(T)).name();
            throw std::runtime_error(errmsg);
            return 0;
        }
    };


    bool operator!=(const QCCalculationParametersBase &p1, const QCCalculationParametersBase &p2);


} /* namespace madness */

#endif /* SRC_APPS_CHEM_QCCALCULATIONPARAMETERSBASE_H_ */
