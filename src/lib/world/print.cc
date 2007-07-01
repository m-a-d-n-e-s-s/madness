#include <world/world.h>

namespace madness {

    void print_justified(const char* s, int column, bool underline) {
        for (int i=0; i<column; i++) std::cout << " ";
        std::cout << s << ENDL;
        if (underline) {
            for (int i=0; i<column; i++) std::cout << " ";
            for (unsigned int i=0; i<std::strlen(s); i++) std::cout << "-";
            std::cout << s << ENDL;
        }
        FLUSH();
    }
    
    void print_centered(const char* s, int column, bool underline) {
        print_justified(s, column-std::strlen(s)/2, underline);
    }

}
