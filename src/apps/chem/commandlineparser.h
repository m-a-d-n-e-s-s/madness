//
// Created by Florian Bischoff on 2/10/21.
//

#ifndef MADNESS_COMMANDLINEPARSER_H
#define MADNESS_COMMANDLINEPARSER_H

#include<map>
namespace madness {
struct commandlineparser {

    std::map<std::string, std::string> keyval;

    // parse command line arguments
    commandlineparser(int argc, char **argv) {
        std::vector<std::string> allArgs(argv, argv + argc);
        for (auto &a : allArgs) {
            std::replace_copy(a.begin(), a.end(), a.begin(), '=', ' ');
            std::replace_copy(a.begin(), a.end(), a.begin(), '-', ' ');
            std::string key, val;
            std::stringstream sa(a);
            sa >> key >> val;
            keyval[key] = val;
        }
    }

    void print_map() const {
        for (auto&[key, val] : keyval) {
            printf("%20s %20s \n", key.c_str(), val.c_str());
        }
    }

    bool key_exists(std::string key) const {
        return (keyval.count(key)==1);
    }

};
}
#endif //MADNESS_COMMANDLINEPARSER_H
