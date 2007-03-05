#include <iostream>

using namespace std;

template<typename T, size_t N>
void tag(T (&x)[N]) {
    std::cout << "array object of size " << sizeof(x)
    << " and length " << N << std::endl;
}

template <size_t n>
template <class T[n]>
class is_an_array {
public:
    static const bool value = false;
};

template <>
class is_an_array<int[3]> {
public:
    static const bool value = true;
};


int main() {



    cout << is_an_array<int>::value << " " << is_an_array<int*> ::value<< " " << is_an_array<int [3]>::value << endl;

    return 0;
}
