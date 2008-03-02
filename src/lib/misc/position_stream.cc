#include <misc/misc.h>
#include <cstring>

namespace madness {
    std::istream& position_stream(std::istream& f, const std::string& tag) {
        std::string s;
        while (std::getline(f,s)) {
            std::string::size_type loc = s.find(tag, 0);
            if(loc != std::string::npos) return f;
        }
        std::string errmsg = std::string("position_stream: failed to locate ") + tag;
        MADNESS_EXCEPTION(errmsg.c_str(),0);
    }

    std::string lowercase(const std::string& s) {
        std::string r(s);
        for (unsigned int i=0; i<r.size(); i++) r[i] = tolower(r[i]);
        return r;
    }

}

