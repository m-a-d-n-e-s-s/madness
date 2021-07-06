//
// Created by Florian Bischoff on 2/10/21.
//

#ifndef MADNESS_COMMANDLINEPARSER_H
#define MADNESS_COMMANDLINEPARSER_H

#include<map>
namespace madness {
/// very simple command line parser

/// parser reads out key/value pairs from the command line of the from  --key=val or --key
/// currently no error handling, feel free to add.
struct commandlineparser {

    std::map<std::string, std::string> keyval;

    commandlineparser() {
        set_defaults();
    }

    // parse command line arguments
    commandlineparser(int argc, char **argv) {
        set_defaults();
        std::vector<std::string> allArgs(argv, argv + argc);
        for (auto &a : allArgs) {
            std::replace_copy(a.begin(), a.end(), a.begin(), '=', ' ');
            std::replace_copy(a.begin(), a.end(), a.begin(), '-', ' ');
            std::string key, val;
            std::stringstream sa(a);
            sa >> key >> val;

            keyval[tolower(key)] = tolower(val);
        }
    }

    /// set default values from the command line
    void set_defaults() {
        keyval["input"]="input";
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
        return keyval.find(tolower(key))->second;
    }

    void set_keyval(const std::string key, const std::string value) {
        keyval[tolower(key)]=tolower(value);
    }

private:
    /// make lower case
    std::string tolower(std::string s) const {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        return s;
    }

};
}
#endif //MADNESS_COMMANDLINEPARSER_H
