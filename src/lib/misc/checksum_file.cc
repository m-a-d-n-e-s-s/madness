#include <cstdio>
using std::fopen;
using std::fgetc;
using std::fclose;

/// \file checksum_file.cc
/// \brief Miscellaneous useful stuff

namespace madness {
    
    /// Simple checksum for ASCII characters in file
    unsigned long checksum_file(const char* filename) {
        FILE *file = fopen(filename,"r");
        if (!file) return 0;
        
        unsigned long sum = 0;
        int c;
        while ((c = fgetc(file)) != EOF) {
            sum = (sum*31u + ((unsigned) c)) & 0xffffff;
        }
        fclose(file);
        return sum;
    }
}
