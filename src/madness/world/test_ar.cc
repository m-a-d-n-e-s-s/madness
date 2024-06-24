/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#include <iostream>
using std::cout;
using std::endl;

#include <cmath>

// make more fundamental types ostreammable
#include <madness/world/array_addons.h>
#include <madness/world/print.h>

template <typename T>
struct type_printer;

/// \file test.cc
/// \brief Tests serialization by some of the archives

#include <madness/world/text_fstream_archive.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;

#include <madness/world/binary_fstream_archive.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;

#include <madness/world/vector_archive.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

#include <madness/world/buffer_archive.h>
using madness::archive::BufferInputArchive;
using madness::archive::BufferOutputArchive;

#include <madness/world/cereal_archive.h>
#ifdef MADNESS_HAS_CEREAL
#include <cereal/archives/binary.hpp>
using CerealBinaryInputArchive = madness::archive::CerealInputArchive<cereal::BinaryInputArchive>;
using CerealBinaryOutputArchive = madness::archive::CerealOutputArchive<cereal::BinaryOutputArchive>;
static_assert(madness::is_archive_v<CerealBinaryInputArchive>);
static_assert(madness::is_input_archive_v<CerealBinaryInputArchive>);
static_assert(!madness::is_output_archive_v<CerealBinaryInputArchive>);
static_assert(!madness::is_text_archive_v<CerealBinaryInputArchive>, "ouch");
static_assert(!madness::is_text_archive_v<CerealBinaryInputArchive>, "ouch");
#include <cereal/archives/portable_binary.hpp>
using CerealPortableBinaryInputArchive = madness::archive::CerealInputArchive<cereal::PortableBinaryInputArchive>;
using CerealPortableBinaryOutputArchive = madness::archive::CerealOutputArchive<cereal::PortableBinaryOutputArchive>;
static_assert(!madness::is_text_archive_v<CerealPortableBinaryInputArchive>, "ouch");
static_assert(!madness::is_text_archive_v<CerealPortableBinaryInputArchive>, "ouch");
#include <cereal/archives/json.hpp>
using CerealJSONInputArchive = madness::archive::CerealInputArchive<cereal::JSONInputArchive>;
using CerealJSONOutputArchive = madness::archive::CerealOutputArchive<cereal::JSONOutputArchive>;
static_assert(madness::is_text_archive_v<CerealJSONInputArchive>, "ouch");
static_assert(madness::is_text_archive_v<CerealJSONOutputArchive>, "ouch");
#include <cereal/archives/xml.hpp>
using CerealXMLInputArchive = madness::archive::CerealInputArchive<cereal::XMLInputArchive>;
using CerealXMLOutputArchive = madness::archive::CerealOutputArchive<cereal::XMLOutputArchive>;
static_assert(madness::is_text_archive_v<CerealXMLInputArchive>, "ouch");
static_assert(madness::is_text_archive_v<CerealXMLOutputArchive>, "ouch");
#endif // MADNESS_HAS_CEREAL

#include <madness/world/world.h>
#include <madness/world/worldgop.h>

template <typename T>
struct type_printer;

///// basic traits tests

static_assert(!madness::has_member_serialize_v<int, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_member_serialize_v<int, madness::archive::BufferInputArchive>);
static_assert(!madness::has_nonmember_serialize_v<int, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_nonmember_serialize_v<int, madness::archive::BufferInputArchive>);
static_assert(!madness::has_freestanding_serialize_v<int, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_freestanding_serialize_v<int, madness::archive::BufferInputArchive>);
static_assert(!madness::has_freestanding_serialize_with_size_v<int*, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_freestanding_serialize_with_size_v<int*, madness::archive::BufferInputArchive>);
static_assert(madness::has_freestanding_default_serialize_v<int, madness::archive::BufferOutputArchive>);
static_assert(madness::has_freestanding_default_serialize_v<int, madness::archive::BufferInputArchive>);
static_assert(madness::has_freestanding_default_serialize_with_size_v<int*, madness::archive::BufferOutputArchive>);
static_assert(madness::has_freestanding_default_serialize_with_size_v<int*, madness::archive::BufferInputArchive>);
static_assert(madness::has_freestanding_default_serialize_with_size_v<int*, madness::archive::BinaryFstreamOutputArchive>);
static_assert(madness::has_nonmember_store_v<int, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_load_v<int, madness::archive::BufferInputArchive>);
static_assert(!madness::is_user_serializable_v<int, madness::archive::BufferOutputArchive>);
static_assert(!madness::is_user_serializable_v<int, madness::archive::BufferInputArchive>);

// A is a class that provides a symmetric serialize method
class A {
public:
    float a;
    A(float a = 0.0) : a(a) {};

    template <class Archive>
    inline void serialize(Archive& ar) {
        ar & a;
    }
};

static_assert(madness::has_member_serialize_v<A, madness::archive::BufferOutputArchive>);
static_assert(madness::has_member_serialize_v<A, madness::archive::BufferInputArchive>);
// N.B. nonmember serialize is provided for types with serialize member
static_assert(madness::has_nonmember_serialize_v<A, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_serialize_v<A, madness::archive::BufferInputArchive>);
// N.B. nonmember load/store is provided for types with serialize member
static_assert(madness::has_nonmember_store_v<A, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_load_v<A, madness::archive::BufferInputArchive>);
static_assert(!madness::has_nonmember_load_and_store_v<A, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_nonmember_load_and_store_v<A, madness::archive::BufferInputArchive>);

// B is a class without a serialize method but with symmetric serialization.
class B {
public:
    bool b;
    B(bool b = false) : b(b) {};
};

namespace madness {
    namespace archive {

        template <class Archive>
        struct ArchiveSerializeImpl<Archive,B> {
            static inline void serialize(const Archive& ar, B& b) {
                ar & b.b;
            };
        };
    }
}

static_assert(!madness::has_member_serialize_v<B, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_member_serialize_v<B, madness::archive::BufferInputArchive>);
static_assert(madness::has_nonmember_serialize_v<B, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_serialize_v<B, madness::archive::BufferInputArchive>);
// N.B. nonmember load/store is provided for types with nonmember serialize
static_assert(madness::has_nonmember_store_v<B, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_load_v<B, madness::archive::BufferInputArchive>);
static_assert(!madness::has_nonmember_load_and_store_v<B, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_nonmember_load_and_store_v<B, madness::archive::BufferInputArchive>);

// C is a class with asymmetric load/store.
class C {
public:
    long c;
    C(long c = 0) : c(c) {};
};

namespace madness {
    namespace archive {
        template <class Archive>
        struct ArchiveLoadImpl<Archive,C> {
            static inline void load(const Archive& ar, C& c) {
                ar >> c.c;
            };
        };

        template <class Archive>
        struct ArchiveStoreImpl<Archive,C> {
            static inline void store(const Archive& ar, const C& c) {
                ar << c.c;
            };
        };
    }
}

static_assert(!madness::has_member_serialize_v<C, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_member_serialize_v<C, madness::archive::BufferInputArchive>);
static_assert(!madness::has_nonmember_serialize_v<C, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_nonmember_serialize_v<C, madness::archive::BufferInputArchive>);
static_assert(madness::has_nonmember_store_v<C, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_load_v<C, madness::archive::BufferInputArchive>);
static_assert(madness::has_nonmember_load_and_store_v<C, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_load_and_store_v<C, madness::archive::BufferInputArchive>);

// POD can be serialized to most (except text stream) archives without any effort
struct D {
  int32_t i;
  int64_t l;
};

static_assert(!madness::has_member_serialize_v<D, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_member_serialize_v<D, madness::archive::BufferInputArchive>);
static_assert(!madness::has_nonmember_serialize_v<D, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_nonmember_serialize_v<D, madness::archive::BufferInputArchive>);
// N.B. nonmember load/store is provided for default-serializable types
static_assert(madness::has_nonmember_store_v<D, madness::archive::BufferOutputArchive>);
static_assert(madness::has_nonmember_load_v<D, madness::archive::BufferInputArchive>);
static_assert(!madness::has_nonmember_load_and_store_v<D, madness::archive::BufferOutputArchive>);
static_assert(!madness::has_nonmember_load_and_store_v<D, madness::archive::BufferInputArchive>);

// to serialize a POD to a text stream just overload stream redirection operators
struct F {
  int32_t i;
  int64_t l;
};
std::ostream& operator<<(std::ostream& os, const F& data) {
  os << data.i << " " << data.l;
  return os;
}
std::istream& operator>>(std::istream& os, F& data) {
  os >> data.i >> data.l;
  return os;
}

// make sure streaming ops have correct return types
struct NotStreammable {
};
int operator<<(std::ostream& os, const NotStreammable& data) {
  return 0;
}
void operator>>(std::istream& os, NotStreammable& data) {
}

static_assert(!madness::is_ostreammable_v<A>, "A is not ostreammable");
static_assert(madness::is_ostreammable_v<F>, "F is ostreammable");
static_assert(madness::is_ostreammable_v<bool>, "bool is ostreammable");
static_assert(madness::is_ostreammable_v<int>, "int is ostreammable");
static_assert(madness::is_ostreammable_v<std::array<int,3>>, "std::array<int,3> is ostreammable");
static_assert(madness::is_ostreammable_v<std::vector<int>>, "std::vector<int> is ostreammable");
static_assert(!madness::is_ostreammable_v<NotStreammable>, "NotStreammable is ostreammable");
static_assert(!madness::is_istreammable_v<A>, "A is not istreammable");
static_assert(madness::is_istreammable_v<F>, "F is istreammable");
static_assert(madness::is_istreammable_v<bool>, "bool is istreammable");
static_assert(madness::is_istreammable_v<int>, "int is istreammable");
static_assert(!madness::is_istreammable_v<std::array<int,3>>, "std::array<int,3> is not istreammable");
static_assert(!madness::is_istreammable_v<std::vector<int>>, "std::vector<int> is not istreammable");
static_assert(!madness::is_istreammable_v<NotStreammable>, "NotStreammable is not istreammable");

static_assert(!madness::is_default_serializable_v<madness::archive::TextFstreamOutputArchive,std::vector<std::vector<int>>>);
static_assert(!madness::is_default_serializable_v<madness::archive::TextFstreamInputArchive,std::vector<std::vector<int>>>);

// serialization of trivially-serializable types can be overloaded
struct G1 {
  int32_t i;
  int64_t l;
  template <typename Archive>
  void serialize(Archive& ar) {
    ar & i & l;
    if constexpr (madness::is_input_archive_v<Archive>) {
      int junk=0;
      ar >> junk;
      MADNESS_ASSERT(junk == 17);
    }
    else {
      static_assert(madness::is_output_archive_v<Archive>);
      ar << 17;
    }
  }
};

// A better example of a class with asym load/store
class linked_list {
    int value;
    linked_list *next;
public:
    linked_list(int value = 0) : value(value), next(0) {};

    void append(int value) {
        if (next) next->append(value);
        else next = new linked_list(value);
    };

    void set_value(int val) {
        value = val;
    };

    int get_value() const {
        return value;
    };

    linked_list* get_next() const {
        return next;
    };

    void print() {
        if (next) {
            cout << value << " --> ";
            next->print();
        }
        else cout << value << endl;
    };
};

namespace {
void free_fn() {}
struct Member {
  static void static_fn(){}
  void fn(){}
  virtual void virtual_fn() {}
};
}

namespace madness {
    namespace archive {
        template <class Archive>
        struct ArchiveStoreImpl<Archive,linked_list> {
            static void store(const Archive& ar, const linked_list& c) {
                ar & c.get_value() & bool(c.get_next());
                if (c.get_next()) ar & *c.get_next();
            };
        };

        template <class Archive>
        struct ArchiveLoadImpl<Archive,linked_list> {
            static void load(const Archive& ar, linked_list& c) {
                int value = 0;
                bool flag = false;
                ar & value & flag;
                c.set_value(value);
                if (flag) {
                    c.append(0);
                    ar & *c.get_next();
                }
            };
        };
    }
}

namespace madness {
    namespace archive {
        typedef std::map< short,std::complex<double> > map_short_complex_double;
        typedef std::pair< short,std::complex<double> > pair_short_complex_double;
        typedef std::pair<int,double> pair_int_double;
        ARCHIVE_REGISTER_TYPE_AND_PTR(A,128);
        ARCHIVE_REGISTER_TYPE_AND_PTR(B,129);
        ARCHIVE_REGISTER_TYPE_AND_PTR(C,130);
        ARCHIVE_REGISTER_TYPE_AND_PTR(linked_list,131);
        ARCHIVE_REGISTER_TYPE_AND_PTR(pair_int_double,132);
        ARCHIVE_REGISTER_TYPE_AND_PTR(map_short_complex_double,133);
        ARCHIVE_REGISTER_TYPE_AND_PTR(pair_short_complex_double, 134);
    }
}

using namespace std;

using madness::archive::wrap;

typedef std::complex<double> double_complex;
typedef std::tuple<int,double,std::complex<float>> tuple_int_double_complexfloat;

template <bool DoSerialize>
struct pod_serialize_dispatch;
template <>
struct pod_serialize_dispatch<false> {
  template <typename Archive, typename POD>
  void operator()(Archive&& ar, const POD& pod) {
  }
};
template <>
struct pod_serialize_dispatch<true> {
  template <typename Archive, typename POD>
  void operator()(Archive&& ar, const POD& pod) {
    ar &pod;
    ar << pod;
  }
};

template <bool DoSerialize>
struct pod_deserialize_dispatch;
template <>
struct pod_deserialize_dispatch<false> {
  template <typename Archive, typename POD>
  void operator()(Archive&& ar, const POD& pod) {
  }
};
template <>
struct pod_deserialize_dispatch<true> {
  template <typename Archive, typename POD>
  void operator()(Archive&& ar, const POD& pod) {
    ar &pod;
    ar >> pod;
  }
};

template <class OutputArchive>
void test_out(const OutputArchive& oar) {
    constexpr const int n = 3;
    A a, an[n];
    B b, bn[n];
    C c, cn[n];
    D d, dn[n];
    F f, fn[n];
    G1 g1, g1n[n];
    int i, in[n];
    double *p = new double[n];
    A *q = new A[n];
    vector<int> v(n);
    vector<vector<int>> vv(n);
    std::array<int64_t,n> arr;
    pair<int,double> pp(33,99.0);
    map<short,double_complex> m;
    const char* teststr = "hello \n dude !";
    string str(teststr);
    linked_list list(0);
    double pi = atan(1.0)*4.0;
    double e = exp(1.0);
    tuple_int_double_complexfloat t = std::make_tuple(1,2.0,std::complex<float>(3.0f,4.0f));
    std::set<int> s{1, 33, 6, 352};

    // Initialize data
    a.a = 1; b.b = 1; c.c = 1; i = 1;
    d.i = 1;  d.l = 2;
    f.i = 1;  f.l = 2;
    g1.i = 1; g1.l = 2;
    for (int k=0; k<n; ++k) {
        p[k] = q[k].a = an[k].a = v[k] = arr[k] = cn[k].c = in[k] = k;
        dn[k].i = k+1;
        dn[k].l = k+2;
        fn[k].i = k+3;
        fn[k].l = k+4;
        g1n[k].i = k+5;
        g1n[k].l = k+6;
        vv[k] = {k+1, k+2, k+3, k+4};
        bn[k].b = k&1;
        m[k] = double_complex(k,k);
        list.append(k+1);
    }

    // test example list code
    list.print();
    oar & list;

    oar & pi & e;

    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " constant double value" << std::endl);
    oar & 1.0;
    oar << 1.0;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int" << std::endl);
    oar & i;
    oar << i;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A" << std::endl);
    oar & a;
    oar << a;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B" << std::endl);
    oar & b;
    oar << b;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C" << std::endl);
    oar & c;
    oar << c;
    // only binary archives can serialize an un-annotated POD
    constexpr const bool D_is_serializable = !madness::is_text_archive_v<OutputArchive>;
    if (D_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " D" << std::endl);
      pod_serialize_dispatch<D_is_serializable>{}(oar, d);
    }
    // only MADNESS text archive can use operator<<(std::ostream) to serialize un-annotated POD
    constexpr const bool F_is_serializable = !madness::is_text_archive_v<OutputArchive> || std::is_same<OutputArchive, TextFstreamOutputArchive>::value;
    if (F_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " F" << std::endl);
      pod_serialize_dispatch<F_is_serializable>{}(oar, f);
    }
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " G1" << std::endl);
    oar & g1;
    oar << g1;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int[]" << std::endl);
    oar & in;
    oar << in;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A[]" << std::endl);
    oar & an;
    oar << an;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B[]" << std::endl);
    oar << bn;
    oar & bn;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C[]" << std::endl);
    oar << cn;
    oar & cn;
    // only binary archives can serialize an un-annotated POD
    if (D_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " D[]" << std::endl);
      pod_serialize_dispatch<D_is_serializable>{}(oar, dn);
    }
    // only MADNESS text archive can use operator<<(std::ostream) to serialize un-annotated POD
    if (F_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " F[]" << std::endl);
      pod_serialize_dispatch<F_is_serializable>{}(oar, fn);
    }
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " G1[]" << std::endl);
    oar << g1n;
    oar & g1n;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " double *p wrapped" << std::endl);
    oar << wrap(p,n);
    oar & wrap(p,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A *q wrapped" << std::endl);
    oar << wrap(q,n);
    oar & wrap(q,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<int>" << std::endl);
    oar << v;
    oar & v;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<vector<int>>" << std::endl);
    oar << vv;
    oar & vv;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " array<int64_t," << n << ">" << std::endl);
    oar << arr;
    oar & arr;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " pair<int,double>" << std::endl);
    oar << pp;
    oar & pp;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " map<short,complex<double>>" << std::endl);
    oar << m;
    oar & m;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " tuple<int,double,complex<float>>" << std::endl);
    oar << t;
    oar & t;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl<< " set<int>" << std::endl);
    oar << s;
    oar & s;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " string" << std::endl);
    oar << str;
    oar & str;

    constexpr const bool ptr_is_serializable = !madness::is_cereal_archive<OutputArchive>::value;
    if constexpr(ptr_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " function pointer" << std::endl);
      oar & free_fn;  // no difference between storing function ref and function pointer
      oar << (&free_fn);
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " static member function pointer"
                                  << std::endl);
      oar & Member::static_fn;  // no difference between storing function ref and function pointer
      oar << (&Member::static_fn);
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " non-static member function pointer"
                                  << std::endl);
      oar & (&Member::fn);
      oar << (&Member::fn);
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " non-static virtual member function pointer"
                                  << std::endl);
      oar & (&Member::virtual_fn);
      oar << (&Member::virtual_fn);
    }

    oar & 1.0 & i & a & b & c & in & an & bn & cn & wrap(p,n) & wrap(q,n) & pp & m & t & s & str;
    if constexpr(ptr_is_serializable) {
      oar & free_fn &(&free_fn) & Member::static_fn & (&Member::static_fn) & (&Member::fn) & (&Member::virtual_fn);
    }

    if (D_is_serializable) {
      pod_serialize_dispatch<D_is_serializable>{}(oar, d);
      pod_serialize_dispatch<D_is_serializable>{}(oar, dn);
    }
    if (F_is_serializable) {
      pod_serialize_dispatch<F_is_serializable>{}(oar, f);
      pod_serialize_dispatch<F_is_serializable>{}(oar, fn);
    }
}

template <class InputArchive>
void test_in(const InputArchive& iar) {
    constexpr const int n = 3;
    A a, an[n];
    B b, bn[n];
    C c, cn[n];
    D d, dn[n];
    F f, fn[n];
    G1 g1, g1n[n];
    int i, in[n];
    double *p = new double[n];
    A *q = new A[n];
    vector<int> v(n);
    vector<vector<int>> vv(n);
    std::array<int64_t,n> arr;
    pair<int,double> pp(33,99.0);
    map<short,double_complex> m;
    std::set<int> s;
    std::set<int> s2 = {1, 33, 6, 352};
    const char *teststr = "hello \n dude !";
    string str(teststr);
    linked_list list;
    double pi = 0.0, e = 0.0;
    tuple_int_double_complexfloat t;
    decltype(&free_fn) free_fn_ptr1=nullptr, free_fn_ptr2=nullptr;
    decltype(&Member::static_fn) static_member_fn_ptr1=nullptr, static_member_fn_ptr2=nullptr;
    decltype(&Member::fn) member_fn_ptr=nullptr, virtual_member_fn_ptr=nullptr;

    // Destroy in-core data
    a.a = 0; b.b = 0; c.c = 0; i = 0;
    d.i = -1;  d.l = -1;
    f.i = -1;  f.l = -1;
    g1.i = -1; g1.l = -1;
    for (int k=0; k<n; ++k) {
        p[k] = q[k].a = an[k].a = v[k] = arr[k] = cn[k].c = in[k] = -1;
        dn[k].i = -1;
        dn[k].l = -1;
        fn[k].i = -1;
        fn[k].l = -1;
        g1n[k].i = -1;
        g1n[k].l = -1;
        vv[k] = {};
        bn[k].b = (k+1)&1;
        m[k] = double_complex(0,0);
    }
    pp = pair<int,double>(0,0.0);
    str = string("");

    iar & list;
    list.print();

    iar & pi & e;
    cout.setf(std::ios::scientific);
    cout << "error in pi " << (pi - 4.0*atan(1.0)) << endl;
    cout << "error in e " << (e - exp(1.0)) << endl;

    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " constant double value" << std::endl);
    double val;
    iar & val;
    iar >> val;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int" << std::endl);
    iar & i;
    iar >> i;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A" << std::endl);
    iar & a;
    iar >> a;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B" << std::endl);
    iar & b;
    iar >> b;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C" << std::endl);
    iar & c;
    iar >> c;
    // only binary archives can serialize an un-annotated POD
    constexpr const bool D_is_serializable = !madness::is_text_archive_v<InputArchive>;
    if (D_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " D" << std::endl);
      pod_deserialize_dispatch<D_is_serializable>{}(iar, d);
    }
    // only MADNESS text archive can use operator<<(std::ostream) to serialize un-annotated POD
    constexpr const bool F_is_serializable = !madness::is_text_archive_v<InputArchive> || std::is_same<InputArchive, TextFstreamInputArchive>::value;
    if (F_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " F" << std::endl);
      pod_deserialize_dispatch<F_is_serializable>{}(iar, f);
    }
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " G1" << std::endl);
    iar & g1;
    iar >> g1;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int[]" << std::endl);
    iar & in;
    iar >> in;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A[]" << std::endl);
    iar & an;
    iar >> an;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B[]" << std::endl);
    iar & bn;
    iar >> bn;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C[]" << std::endl);
    iar & cn;
    iar >> cn;
    if (D_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " D[]" << std::endl);
      pod_deserialize_dispatch<D_is_serializable>{}(iar, dn);
    }
    if (F_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl << " F[]" << std::endl);
      pod_deserialize_dispatch<F_is_serializable>{}(iar, fn);
    }
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " G1[]" << std::endl);
    iar & g1n;
    iar >> g1n;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " double *p wrapped" << std::endl);
    iar & wrap(p,n);
    iar >> wrap(p,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A *q wrapped" << std::endl);
    iar & wrap(q,n);
    iar >> wrap(q,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<int>" << std::endl);
    iar & v;
    iar >> v;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<vector<int>>" << std::endl);
    iar & vv;
    iar >> vv;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " array<int64_t," << n << ">" << std::endl);
    iar & arr;
    iar >> arr;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " pair<int,double>" << std::endl);
    iar & pp;
    iar >> pp;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " map<short,complex<double>>" << std::endl);
    iar & m;
    iar >> m;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " tuple<int,double,complex<float>>" << std::endl);
    iar & t;
    iar >> t;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl<< " set<int>" << std::endl);
    iar >> s;
    iar & s;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " string" << std::endl);
    iar & str;
    iar >> str;

    constexpr const bool ptr_is_serializable = true && !madness::is_cereal_archive<InputArchive>::value;
    if constexpr(ptr_is_serializable) {
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " function pointer" << std::endl);
      static_assert(!madness::is_istreammable_v<void(*)()>);
      static_assert(madness::is_ostreammable_v<void(*)()>);
      static_assert(madness::is_default_serializable_v<madness::archive::TextFstreamInputArchive, void (*)()>);
      static_assert(madness::is_default_serializable_v<madness::archive::TextFstreamInputArchive, void (*)()>);
      iar & free_fn_ptr1;
      iar >> free_fn_ptr2;
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " static member function pointer"
                                  << std::endl);
      iar & static_member_fn_ptr1;
      iar >> static_member_fn_ptr2;
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " non-static member function pointer"
                                  << std::endl);
      decltype(&Member::fn) nonstatic_member_fn_ptr;
      iar & nonstatic_member_fn_ptr;
      iar >> nonstatic_member_fn_ptr;
      MAD_ARCHIVE_DEBUG(std::cout << std::endl
                                  << " non-static virtual member function pointer"
                                  << std::endl);
      decltype(&Member::virtual_fn) nonstatic_virtual_member_fn_ptr;
      iar & nonstatic_virtual_member_fn_ptr;
      iar >> nonstatic_virtual_member_fn_ptr;
    }

    iar & 1.0 & i & a & b & c & in & an & bn & cn & wrap(p,n) & wrap(q,n) & pp & m & t & s & str;
    if constexpr(ptr_is_serializable) {
      iar &free_fn_ptr1 &free_fn_ptr2 &static_member_fn_ptr1
          &static_member_fn_ptr2 & member_fn_ptr & virtual_member_fn_ptr;
    }
    if (D_is_serializable) {
      pod_deserialize_dispatch<D_is_serializable>{}(iar, d);
      pod_deserialize_dispatch<D_is_serializable>{}(iar, dn);
    }
    if (F_is_serializable) {
      pod_deserialize_dispatch<F_is_serializable>{}(iar, f);
      pod_deserialize_dispatch<F_is_serializable>{}(iar, fn);
    }
    // Test data
    bool status = true;

#define TEST(cond) status &= cond;  \
                   if (!(cond)) std::cout << #cond << " failed" << std::endl
    TEST(a.a == 1);
    TEST(b.b == 1);
    TEST(c.c == 1);
    if (D_is_serializable) {
      TEST(d.i == 1);
      TEST(d.l == 2);
    }
    if (F_is_serializable) {
      TEST(f.i == 1);
      TEST(f.l == 2);
    }
    TEST(g1.i == 1);
    TEST(g1.l == 2);
    TEST(i == 1);
    for (int k=0; k<n; ++k) {
        TEST(an[k].a == k);
        TEST(bn[k].b == (k&1));
        TEST(cn[k].c == k);
        if (D_is_serializable) {
          TEST(dn[k].i == k+1);
          TEST(dn[k].l == k+2);
        }
        if (F_is_serializable) {
          TEST(fn[k].i == k + 3);
          TEST(fn[k].l == k + 4);
        }
        TEST(g1n[k].i == k + 5);
        TEST(g1n[k].l == k + 6);
        TEST(in[k] == k);
        TEST(p[k] == k);
        TEST(q[k].a == k);
        TEST(v[k] == k);
        TEST(vv[k].size() == 4);
        TEST(arr[k] == k);
        TEST(vv[k][0] == k+1);
        TEST(vv[k][1] == k+2);
        TEST(vv[k][2] == k+3);
        TEST(vv[k][3] == k+4);
        TEST(m[k] == double_complex(k,k));
    }
    TEST(pp.first==33 && pp.second==99.0);
    TEST(str == string(teststr));
    TEST(t == std::make_tuple(1,2.0,std::complex<float>(3.0f,4.0f)));
    TEST(s == s2);
    if constexpr (ptr_is_serializable)
    {
      TEST(free_fn_ptr1 == &free_fn);
      TEST(free_fn_ptr2 == &free_fn);
      TEST(static_member_fn_ptr1 == &Member::static_fn);
      TEST(static_member_fn_ptr2 == &Member::static_fn);
      TEST(member_fn_ptr == &Member::fn);
      TEST(virtual_member_fn_ptr == &Member::virtual_fn);
    }

#undef TEST

    if (status)
        std::cout << "Serialization appears to work." << std::endl;
    else
        std::cout << "Sorry, back to the drawing board.\n";
}

int main(int argc, char* argv[]) {

  auto &world = madness::initialize(argc, argv);

  const auto me = world.rank();
  const auto nproc = world.size();
  const auto is_write_proc = (me == 0);
  const auto is_read_proc = (nproc == 1) ? (me == 0) : (me == 1);
  const int nworkers = (nproc == 1) ? 1 : 2;
  for(int rank=0; rank!=nworkers; ++rank) {
    if (rank == me) {
      cout << "rank=" << rank << " fn_ptr_origin=" << std::hex
           << madness::archive::fn_ptr_origin() << endl;
    }
    world.gop.barrier();
  }

  madness::archive::archive_initialize_type_names();
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(A);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(B);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(C);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(linked_list);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(madness::archive::pair_int_double);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(
      madness::archive::pair_short_complex_double);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(
      madness::archive::map_short_complex_double);

  {
    const char *f = "test.dat";
    if (is_write_proc) {
      cout << endl << "testing binary fstream archive" << endl;
      BinaryFstreamOutputArchive oar(f);
      test_out(oar);
      oar.close();
    }
    world.gop.barrier();

    if (is_read_proc) {
      BinaryFstreamInputArchive iar(f);
      test_in(iar);
      iar.close();
    }
    world.gop.barrier();
  }

  if (me == 0) {
    cout << endl << "testing vector archive" << endl;
    std::vector<unsigned char> f;
    VectorOutputArchive oar(f);
    test_out(oar);
    oar.close();

    VectorInputArchive iar(f);
    test_in(iar);
    iar.close();
  }
  world.gop.barrier();

  {
    const char *f = "test.dat";
    if (is_write_proc) {
      cout << endl << "testing text fstream archive" << endl;
      TextFstreamOutputArchive oar(f);
      test_out(oar);
      oar.close();
    }
    world.gop.barrier();

    if (is_read_proc) {
      TextFstreamInputArchive iar(f);
      test_in(iar);
      iar.close();
    }
    world.gop.barrier();
  }

  // buffer archive only tested on rank 0
  if (me == 0) {
    cout << endl << "testing buffer archive" << endl;
    unsigned char buf[32768];
    BufferOutputArchive oar(buf, sizeof(buf));
    test_out(oar);
    std::size_t nbyte = oar.size();
    oar.close();

    BufferInputArchive iar(buf, nbyte);
    test_in(iar);
    iar.close();
  }
  world.gop.barrier();

#ifdef MADNESS_HAS_CEREAL
  {
    const char *f = "test.dat";
    if (is_write_proc) {
      cout << endl << "testing binary Cereal archive" << endl;
      std::ofstream fout(f, std::ios_base::binary | std::ios_base::out |
                                std::ios_base::trunc);
      CerealBinaryOutputArchive oar(fout);
      test_out(oar);
      oar.close();
    }
    world.gop.barrier();

    if (is_read_proc) {
      std::ifstream fin(f, std::ios_base::binary | std::ios_base::in);
      CerealBinaryInputArchive iar(fin);
      test_in(iar);
      iar.close();
    }
    world.gop.barrier();

    if (is_write_proc) {
      cout << endl << "testing portable binary Cereal archive" << endl;
      std::ofstream fout(f, std::ios_base::binary | std::ios_base::out |
                                std::ios_base::trunc);
      CerealPortableBinaryOutputArchive oar(fout);
      test_out(oar);
      oar.close();
    }
    world.gop.barrier();

    if (is_read_proc) {
      std::ifstream fin(f, std::ios_base::binary | std::ios_base::in);
      CerealPortableBinaryInputArchive iar(fin);
      test_in(iar);
      iar.close();
    }
    world.gop.barrier();

    if (is_write_proc) {
      cout << endl << "testing JSON Cereal archive" << endl;
      std::ofstream fout(f, std::ios_base::out | std::ios_base::trunc);
      CerealJSONOutputArchive oar(fout);
      test_out(oar);
      oar.close();
    }
    world.gop.barrier();

    if (is_read_proc) {
      std::ifstream fin(f, std::ios_base::in);
      CerealJSONInputArchive iar(fin);
      test_in(iar);
      iar.close();
    }
    world.gop.barrier();

    if (is_write_proc) {
      cout << endl << "testing XML Cereal archive" << endl;
      std::ofstream fout(f, std::ios_base::out | std::ios_base::trunc);
      CerealXMLOutputArchive oar(fout);
      test_out(oar);
      oar.close();
    }
    world.gop.barrier();

    if (is_read_proc) {
      std::ifstream fin(f, std::ios_base::in);
      CerealXMLInputArchive iar(fin);
      test_in(iar);
      iar.close();
    }
    world.gop.barrier();
  }
#endif // MADNESS_HAS_CEREAL

  madness::finalize();
  return 0;
}
