#ifndef __MY_QD_H
#define __MY_QD_H

#include <qd/qd_real.h>
#include <type_traits>
#include <string>

// Code using qd/dd CANNOT be compiled with options that fully relax floating point standards.
// -O3 works 
// -Ofast breaks
// -O3 -ffast-math breaks
// -Ofast -fno-associative-math works but should be more extensively tested

// Inject math routines into std name space to enable templated types
namespace std {
    qd_real abs(const qd_real& x) {return ::abs(x);}
    qd_real cos(const qd_real& x) {return ::cos(x);}
    qd_real sin(const qd_real& x) {return ::sin(x);}
    qd_real pow(const qd_real& x, int n) {return ::pow(x,n);}
    qd_real pow(const qd_real& x, const qd_real& n) {return ::pow(x,n);}
    qd_real acos(const qd_real& x) {return ::acos(x);}
    qd_real sqrt(const qd_real& x) {return ::sqrt(x);}
    
    dd_real abs(const dd_real& x) {return ::abs(x);}
    dd_real cos(const dd_real& x) {return ::cos(x);}
    dd_real sin(const dd_real& x) {return ::sin(x);}
    dd_real pow(const dd_real& x, int n) {return ::pow(x,n);}
    dd_real pow(const dd_real& x, const dd_real& n) {return ::pow(x,n);}
    dd_real acos(const dd_real& x) {return ::acos(x);}
    dd_real sqrt(const dd_real& x) {return ::sqrt(x);}
}

template <typename T> T from_str(const char* s);

template <>
qd_real from_str<qd_real>(const char* s) {
    return qd_real(s);
}

template <>
dd_real from_str<dd_real>(const char* s) {
    return dd_real(s);
}

template <>
double from_str<double>(const char* s) {
    double d;
    auto status = sscanf(s, "%lf", &d);
    if (status == EOF) throw "EOF reading double from string";
    if (status != 1) throw "failed to read double from string";
    return d;
}

template <>
float from_str<float>(const char* s) {
    float d;
    auto status = sscanf(s, "%f", &d);
    if (status == EOF) throw "EOF reading float from string";
    if (status != 1) throw "failed to read float from string";
    return d;
}

std::string to_str(const qd_real& t) {return t.to_string();}
std::string to_str(const dd_real& t) {return t.to_string();}
std::string to_str(double t) {char buf[256]; sprintf(buf,"%.17f",t); return buf;}

double to_double(double d) {return d;}

double to_double(float f) {return f;}

// Needed since double(qd_real) not implemented and I fear what would break if I added it.
// Also, qd_real(unsigned or size_t) does not compile due to ambiguity ... just force to double before final conversion
// Also dd_real(qd_real) is missing ... have to go via string

template <typename T, typename R>
T convert(const R& r) {
    if constexpr(std::is_same<T,double>::value) {
        return to_double(r);
    }
    else if constexpr((std::is_same<T,qd_real>::value || std::is_same<T,dd_real>::value) && std::is_integral<R>::value) {
        return T(double(r)); // int to dd or qd is ambiguous conversion
    }
    else if constexpr(std::is_same<T,dd_real>::value && std::is_same<R,qd_real>::value) {
        return T(r.to_string().c_str());
    }
    else {
        return T(r);
    }
};

#endif
