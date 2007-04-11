#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#ifndef POWER_H
#define POWER_H

using namespace std;

namespace madness {

template <int D>
inline int power(int base = 2) {
    return (int) pow((double) base, (int) D);
};

template <>
inline int power<0>(int base) {
    return 1;
};

template <>
inline int power<1>(int base) {
    return base;
};

template <>
inline int power<2>(int base) {
    return (int) (base*base);
};

template <>
inline int power<3>(int base) {
    return (int) (base*base*base);
};

template <>
inline int power<4>(int base) {
    return power<2>(power<2>(base));
};

template <>
inline int power<5>(int base) {
    return (power<2>(base)*power<3>(base));
};

template <>
inline int power<6>(int base) {
    return power<2>(power<3>(base));
};

template <>
inline int power<7>(int base) {
    return (power<4>(base)*power<3>(base));
};

template <>
inline int power<8>(int base) {
    return (power<2>(power<4>(base)));
};

template <>
inline int power<9>(int base) {
    return (power<3>(power<3>(base)));
};

template <>
inline int power<10>(int base) {
    return (power<2>(power<5>(base)));
};

template <>
inline int power<12>(int base) {
    return (power<3>(power<4>(base)));
};
}

#endif
