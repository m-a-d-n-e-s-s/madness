#ifndef __MY_QD_H
#define __MY_QD_H

#include <qd/qd_real.h>
#include <type_traits>
#include <string>
#include <cstdint>
#include <algorithm>
#include <cmath>

// Code using qd/dd CANNOT be compiled with options that fully relax floating point standards.
// -O3 works 
// -Ofast breaks
// -O3 -ffast-math breaks
// -Ofast -fno-associative-math works but should be more extensively tested

template <typename T, typename R> T convert(const R& r);

// Inject math routines into std name space to enable templated types
namespace std {
    qd_real abs(const qd_real& x) {return ::abs(x);}
    qd_real exp(const qd_real& x) {return ::exp(x);}
    qd_real ldexp(const qd_real &a, int exp) {return ::ldexp(a,exp);}
    qd_real cos(const qd_real& x) {return ::cos(x);}
    qd_real sin(const qd_real& x) {return ::sin(x);}
    qd_real pow(const qd_real& x, int n) {return ::pow(x,n);}
    qd_real pow(const qd_real& x, const qd_real& n) {return ::pow(x,n);}
    qd_real acos(const qd_real& x) {return ::acos(x);}
    qd_real sqrt(const qd_real& x) {return ::sqrt(x);}
    
    dd_real abs(const dd_real& x) {return ::abs(x);}
    dd_real exp(const dd_real& x) {return ::exp(x);}
    dd_real ldexp(const dd_real &a, int exp) {return ::ldexp(a,exp);}
    dd_real cos(const dd_real& x) {return ::cos(x);}
    dd_real sin(const dd_real& x) {return ::sin(x);}
    dd_real pow(const dd_real& x, int n) {return ::pow(x,n);}
    dd_real pow(const dd_real& x, const dd_real& n) {return ::pow(x,n);}
    dd_real acos(const dd_real& x) {return ::acos(x);}
    dd_real sqrt(const dd_real& x) {return ::sqrt(x);}
}

// Option 1: 15-digit block parsing into exact integers / doubles
template <typename T>
T myread(const std::string& s) {
    int sx = 1, sexpnt = 1, expnt = 0;
    int decimal_point = -1;
    bool doingexpt = false;
    std::string digits;
    digits.reserve(80);

    for (char c : s) {
        if (doingexpt) {
            if (c == '-') sexpnt = -1;
            else if (c == '+') continue;
            else if (c >= '0' && c <= '9') expnt = expnt * 10 + (c - '0');
        } else {
            if (c == '-') sx = -1;
            else if (c == '.') decimal_point = digits.size();
            else if (c == 'e' || c == 'E') doingexpt = true;
            else if (c >= '0' && c <= '9') digits.push_back(c);
        }
    }

    size_t first_nz = digits.find_first_not_of('0');
    if (first_nz == std::string::npos) return T(0.0);

    if (decimal_point == -1) decimal_point = digits.size();
    int net_exp = sexpnt * expnt - (digits.size() - decimal_point);

    // Limit digits to maximum needed for T
    size_t max_digits = 32;
    if constexpr (std::is_same_v<T, qd_real>) max_digits = 64;
    else if constexpr (std::is_same_v<T, double>) max_digits = 17;
    else if constexpr (std::is_same_v<T, float>) max_digits = 9;

    if (digits.size() > max_digits) {
        net_exp += static_cast<int>(digits.size() - max_digits);
        digits.resize(max_digits);
    }

    // Process in 15-digit blocks (exact in uint64_t / double)
    T result = 0.0;
    int cur_exp = net_exp;
    int len = static_cast<int>(digits.size());

    for (int i = len; i > 0; i -= 15) {
        int start = std::max(0, i - 15);
        int block_len = i - start;
        uint64_t val = 0;
        for (int j = start; j < i; ++j) {
            val = val * 10 + (digits[j] - '0');
        }
        
        T block_val = T(static_cast<double>(val));
        
        T scale;
        if constexpr (std::is_same_v<T, dd_real> || std::is_same_v<T, qd_real>) {
            qd_real p10 = std::pow(qd_real(10.0), cur_exp);
            scale = convert<T>(p10);
        } else {
            scale = std::pow(T(10.0), cur_exp);
        }

        result += block_val * scale;
        cur_exp += block_len;
    }

    return (sx < 0) ? -result : result;
}

template <typename T>
T myread(const char* s) {
    return myread<T>(std::string(s));
}

template <typename T> T from_str(const char* s);

template <>
qd_real from_str<qd_real>(const char* s) {
    return myread<qd_real>(s);
}

template <>
dd_real from_str<dd_real>(const char* s) {
    return myread<dd_real>(s);
}

template <>
double from_str<double>(const char* s) {
    return atof(s);
}

template <>
float from_str<float>(const char* s) {
    float d;
    auto status = sscanf(s, "%f", &d);
    if (status == EOF) throw "EOF reading float from string";
    if (status != 1) throw "failed to read float from string";
    return d;
}

template <typename T> T from_str(const std::string& s) {
    return from_str<T>(s.c_str());
}

std::string to_str(const qd_real& t) {return t.to_string();}
std::string to_str(const dd_real& t) {return t.to_string();}
std::string to_str(double t) {char buf[256]; sprintf(buf,"%.19e",t); return buf;}
std::string to_str(float t) {char buf[256]; sprintf(buf,"%.9e",t); return buf;}

double to_double(double d) {return d;}
double to_double(float f) {return f;}

// Needed since double(qd_real) not implemented and I fear what would break if I added it.
// Also, qd_real(unsigned or size_t) does not compile due to ambiguity ... just force to double before final conversion
// Also dd_real(qd_real) is missing ... use direct dd_real(r[0], r[1]) conversion

template <typename T, typename R>
T convert(const R& r) {
    if constexpr(std::is_same<T,double>::value) {
        return to_double(r);
    }
    else if constexpr((std::is_same<T,qd_real>::value || std::is_same<T,dd_real>::value) && std::is_integral<R>::value) {
        return T(double(r)); // int to dd or qd is ambiguous conversion
    }
    else if constexpr(std::is_same<T,dd_real>::value && std::is_same<R,qd_real>::value) {
        return dd_real(r[0], r[1]);
    }
    else {
        return T(r);
    }
};

#endif
