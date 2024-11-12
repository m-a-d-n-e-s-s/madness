#include <random>
#include <iostream>

int main() {
    std::mt19937 gen;    
    std::uniform_int_distribution r(0,3);
    for (int i=0; i<20; i++) std::cout << r(gen) << std::endl;
    
    return 0;
}
