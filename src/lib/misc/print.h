#ifndef PRINT_H
#define PRINT_H

#include <iostream>

/// \file print.h
/// \brief Defines simple templates for printing to cout a la Python

namespace madness {

    /// Print a single item to cout terminating with new line
    template <class A>
    void print(const A& a) {
        std::cout << a << std::endl;
        std::cout.flush();
    }

    /// Print two items separated by spaces to cout terminating with new line
    template <class A, class B>
    void print(const A& a, const B& b) {
        std::cout << a << " " << b << std::endl;
        std::cout.flush();
    }

    /// Print three items separated by spaces to cout terminating with new line
    template <class A, class B, class C>
    void print(const A& a, const B& b, const C& c) {
        std::cout << a << " " << b << " " << c << std::endl;
        std::cout.flush();
    }

    /// Print four items separated by spaces to cout terminating with new line
    template <class A, class B, class C, class D>
    void print(const A& a, const B& b, const C& c, const D& d) {
        std::cout << a << " " << b << " " << c << " " << d << std::endl;
        std::cout.flush();
    }

    /// Print five items separated by spaces to cout terminating with new line
    template <class A, class B, class C, class D, class E>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << std::endl;
        std::cout.flush();
    }

    /// Print six items separated by spaces to cout terminating with new line
    template <class A, class B, class C, class D, class E, class F>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
        std::cout.flush();
    }

    /// Print seven items separated by spaces to cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << std::endl;
        std::cout.flush();
    }

    /// Print eight items separated by spaces to cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << std::endl;
        std::cout.flush();
    }
    
    /// Print nine items separated by spaces to cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << std::endl;
        std::cout.flush();
    }
    
    /// Print ten items separated by spaces to cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I, class J>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << " " << j << std::endl;
        std::cout.flush();
    }
}
#endif
