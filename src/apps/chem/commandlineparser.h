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
    // mp2 --mp2='maxiter 10; freeze 1' --dft:maxiter=20 --Xmpi:debug=true
    commandlineparser(int argc, char **argv) {
        set_defaults();
        std::vector<std::string> allArgs_raw(argv, argv + argc);
        for (auto &a : allArgs_raw) {
            a= remove_first_equal(remove_front_hyphens(a));
            std::string key, val;
            std::stringstream sa(a);
            sa >> key >> val;
            val=a.substr(key.size());
            std::replace_copy(val.begin(), val.end(), val.begin(), '=', ' ');

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
        MADNESS_CHECK(key_exists(key));
        return keyval.find(tolower(key))->second;
    }

    void set_keyval(const std::string key, const std::string value) {
        keyval[tolower(key)]=tolower(value);
    }

private:
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
        auto it=find_first_of(result.begin(),result.end(),item.begin(),item.end());
        replace(it,it+1,item.front(),blank.front());
        return result;
    }
    static std::string remove_blanks(const std::string arg) {
        std::string str2 = arg;
        str2.erase(std::remove_if(str2.begin(), str2.end(),
                                  [](unsigned char x){return std::isspace(x);}),str2.end());
        return str2;
    }

};
}
#endif //MADNESS_COMMANDLINEPARSER_H
