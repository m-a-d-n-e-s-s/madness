#ifndef PRINT_H
#define PRINT_H

#include <iostream>

#include <typestuff.h>

// Braindead Jaguar
#define ENDL "\n"
static inline void FLUSH() {};

/// \file print.h
/// \brief Defines simple templates for printing to cout a la Python

namespace madness {
    
    /// Print a fixed dimension array to cout terminating with a new line
    template <typename T, std::size_t N>
    void print_array(const T (&v)[N]) {
        for (int i=0; i<(int)N; i++) std::cout << v[i] << " ";
        std::cout<<ENDL;
        FLUSH();
    };

    /// Print a single item to std::cout terminating with new line
    template <class A>
    void print(const A& a) {
        std::cout << a << ENDL;
        FLUSH();
    }

    /// Print two items separated by spaces to std::cout terminating with new line
    template <class A, class B>
    void print(const A& a, const B& b) {
        std::cout << a << " " << b << ENDL;
        FLUSH();
    }

    /// Print three items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C>
    void print(const A& a, const B& b, const C& c) {
        std::cout << a << " " << b << " " << c << ENDL;
        FLUSH();
    }

    /// Print four items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D>
    void print(const A& a, const B& b, const C& c, const D& d) {
        std::cout << a << " " << b << " " << c << " " << d << ENDL;
        FLUSH();
    }

    /// Print five items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << ENDL;
        FLUSH();
    }

    /// Print six items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << ENDL;
        FLUSH();
    }

    /// Print seven items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << ENDL;
        FLUSH();
    }

    /// Print eight items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << ENDL;
        FLUSH();
    }
    
    /// Print nine items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << ENDL;
        FLUSH();
    }
    
    /// Print ten items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I, class J>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << " " << j << ENDL;
        FLUSH();
    }
}
#endif
