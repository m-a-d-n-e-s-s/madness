//
// Created by Florian Bischoff on 2/10/21.
//

#ifndef MADNESS_COMMANDLINEPARSER_H
#define MADNESS_COMMANDLINEPARSER_H

#include<map>
#include<algorithm>
namespace madness {
/// very simple command line parser

/// parser reads out key/value pairs from the command line of the form
/// `--key=val` or `--key`. By default values are case-folded so that
/// chemistry-style enum knobs (e.g. `--xc=HF`, `--localize=BOYS`) match
/// the lowercased `allowed_values` of `QCCalculationParametersBase`-derived
/// classes.
///
/// Two accessors are exposed for retrieving values:
///   * `value(key)`     — case-folded value (back-compat). Use for enum-like
///                        chemistry parameters.
///   * `value_raw(key)` — original-case value. Use for **paths and free-form
///                        identifiers** (archive filenames, JSON paths, …).
///
/// The bare-positional input-file form (e.g. `madness_app /Some/Path/in.in`)
/// is also case-preserving — both `value("input")` and `value_raw("input")`
/// return the path with original case.
struct commandlineparser {

    std::map<std::string, std::string> keyval;        // case-folded values
    std::map<std::string, std::string> keyval_raw;    // original-case values

    commandlineparser() {
        set_defaults();
    }

    // parse command line arguments
    // mp2 --mp2='maxiter 10; freeze 1' --dft:maxiter=20 --Xmpi:debug=true
    commandlineparser(int argc, char **argv) {
        set_defaults();
        std::vector<std::string> allArgs_raw(argv, argv + argc);
        allArgs_raw.erase(allArgs_raw.begin());     // first argument is the name of the binary
        for (auto &a : allArgs_raw) {
            // special treatment for the input file: no hyphens
            a=check_for_input_file(a);
            a= remove_first_equal(remove_front_hyphens(a));
            std::replace_copy(a.begin(), a.end(), a.begin(), '=', ' ');
            std::string key, val;
            std::stringstream sa(a);
            sa >> key;
            val=a.substr(key.size());
            if (key=="input") set_keyval("user_defined_input_file","1");
            if (key=="prefix") set_keyval("user_defined_prefix","1");

            // Path-like keys: keep case in BOTH maps so the back-compat
            // `value("input"|"file")` form also returns the original case
            // (already true today for "file"; extending to "input" so the
            // bare-positional form `app /Some/Path/in.in` survives).
            if (key=="file" || key=="input") {
                set_keyval_keep_case(key,val);
            } else {
                set_keyval(key,val);
            }
        }
    }

    /// set default values from the command line
    void set_defaults() {
        keyval["input"]="input";
        keyval["prefix"]="mad";
//        keyval["geometry"]="input_file";
    }

    void print_map() const {
        for (auto&[key, val] : keyval) {
            printf("%20s %20s \n", key.c_str(), val.c_str());
        }
    }

    bool key_exists(std::string key) const {
        return (keyval.count(tolower(key))==1);
    }

    std::string value(const std::string key) const {
        std::string msg= "key not found: " + key;
        MADNESS_CHECK_THROW(key_exists(key), msg.c_str());
        return keyval.find(tolower(key))->second;
    }

    /// Original-case value lookup — use for paths and free-form identifiers
    /// (archive filenames, JSON paths, etc.). For enum-like chemistry knobs
    /// matched against `allowed_values`, prefer `value(key)` instead.
    std::string value_raw(const std::string key) const {
        std::string msg = "key not found: " + key;
        MADNESS_CHECK_THROW(key_exists(key), msg.c_str());
        auto it = keyval_raw.find(tolower(key));
        if (it == keyval_raw.end()) {
            // Defensive fallback for keys set via the default-constructor
            // path (defaults), which only seed `keyval`. Return the
            // case-folded version rather than throwing.
            return keyval.find(tolower(key))->second;
        }
        return it->second;
    }

    void set_keyval(const std::string key, const std::string value) {
        const std::string trimmed = trim_blanks(value);
        keyval[tolower(key)]    = tolower(trimmed);
        keyval_raw[tolower(key)] = trimmed;
    }

    void set_keyval_keep_case(const std::string key, const std::string value) {
        const std::string trimmed = trim_blanks(value);
        keyval[tolower(key)]    = trimmed;
        keyval_raw[tolower(key)] = trimmed;
    }

public:

    /// special option: the input file has no hyphens in front and is just a value
    std::string check_for_input_file(std::string line) {
        if (line[0]=='-') return line;
        auto words=split(line,"=");
        if (words.size()==1) line="input="+line;
        return line;
    }
    /// make lower case
    static std::string tolower(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        return s;
    }

    /// split a string s into a vector of strings, using delimiter

    /// @param[in]  s   the string (pass by value!)
    static std::vector<std::string> split(std::string s, const std::string delimiter) {
        std::size_t pos = 0;
        std::string token;
        std::vector<std::string> result;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            result.push_back(token);
            s.erase(0, pos + delimiter.length());
        }
        result.push_back(s);
        return result;
    }

    static std::string remove_front_hyphens(const std::string arg) {
        std::size_t first=arg.find_first_not_of('-');
        return arg.substr(first);
    }

    static std::string remove_first_equal(const std::string arg) {
        std::string result=arg;
        const std::string item="=";
        const std::string blank=" ";
        auto it=std::find_first_of(result.begin(),result.end(),item.begin(),item.end());
        std::replace(it,it+1,item.front(),blank.front());
        return result;
    }

    /// remove all blanks
    static std::string remove_blanks(const std::string arg) {
        std::string str2 = arg;
        str2.erase(std::remove_if(str2.begin(), str2.end(),
                                  [](unsigned char x){return std::isspace(x);}),str2.end());
        return str2;
    }

    /// remove blanks at the beginning and the end only
    static std::string trim_blanks(const std::string arg) {
        if (arg.size()==0) return arg;
        std::size_t first=arg.find_first_not_of(' ');
        std::size_t last=arg.find_last_not_of(' ');
        return arg.substr(first,last-first+1);
    }

    static std::string base_name(std::string const & path, std::string const & delims = "/")
    {
        return path.substr(path.find_last_of(delims) + 1);
    }

    static std::string remove_extension(std::string const & filename)
    {
        std::size_t p=filename.find_last_of('.');
        return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
    }

};
}
#endif //MADNESS_COMMANDLINEPARSER_H
