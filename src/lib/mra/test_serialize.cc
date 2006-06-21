#include<vector>
using namespace std;

#include <iostream>
using std::cout;
using std::endl;

#include <cmath>

#include <misc/communicator.h>
using madness::Communicator;

#include <mra/twoscale.h>
#include <mra/legendre.h>

/// \file serialize/test.cc
/// \brief Tests serialization by some of the archives

#define ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE

#include <serialize/textfsar.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;

#include <serialize/binfsar.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;

#include <serialize/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

#include <mra/mra.h>
using madness::Function;
//using madness::Function::save_local;
using madness::FunctionDefaults;
using madness::FunctionFactory;
using madness::FunctionOctTree;
using madness::FunctionNode;

#include <octtree/octtree.h>
using madness::OctTree;
using madness::OctTreeT;

#include <misc/misc.h>
using madness::redirectio;

// A is a class that provides a symmetric serialize method
class A {
public:
    float a;
    A(float a = 0.0) : a(a) {};

    template <class Archive>
    inline void serialize(const Archive& ar) {
        ar & a;
    }
};

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

// D is a class that provides a serialize method
class D {
private:
    int *a;

public:
    D(int value) {
        a = new int [value];
    };

    int get() {
        return *a;
    };

    void set(int value) {
        *a = value;
    };

    ~D() {
        delete [] a;
    };

//  template <class Archive>
//  inline void serialize(const Archive& ar) {ar & a;}
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
        } else cout << value << endl;
    };
};

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
                int value;
                bool flag;
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
using madness::Tensor;

template <class OutputArchive>
void test_out(const OutputArchive& oar) {
    const int n = 3;
    A a, an[n];
    B b, bn[n];
    C c, cn[n];
    int i, in[n];
    double *p = new double[n];
    A *q = new A[n];
    vector<int> v(3);
    pair<int,double> pp(33,99.0);
    map<short,double_complex> m;
    const char* teststr = "hello \n dude !";
    string str(teststr);
    linked_list list(0);
    Tensor<int> ti(n,n);
    Tensor<double> td(n,n);
    Tensor<double_complex> tc(n,n);
    double pi = atan(1.0)*4.0;
    double e = exp(1.0);

    // Initialize data
    a.a = b.b = c.c = i = 1;
    for (int k=0; k<n; k++) {
        p[k] = q[k].a = an[k].a = v[k] = cn[k].c = in[k] = k;
        bn[k].b = k&1;
        m[k] = double_complex(k,k);
        list.append(k+1);
    }
    ti.fillindex();
    td.fillindex();
    tc.fillindex();

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
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " double *p wrapped" << std::endl);
    oar << wrap(p,n);
    oar & wrap(p,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A *q wrapped" << std::endl);
    oar << wrap(q,n);
    oar & wrap(q,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<int>" << std::endl);
    oar << v;
    oar & v;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " pair<int,double>" << std::endl);
    oar << pp;
    oar & pp;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " map<short,double<complex>>" << std::endl);
    oar << m;
    oar & m;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " string" << std::endl);
    oar << str;
    oar & str;

    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " Tensors" << std::endl);
    oar << ti << td << tc;
    oar & ti & td & tc;

    oar & 1.0 & i & a & b & c & in & an & bn & cn & wrap(p,n) & wrap(q,n) & pp & m & str;
}

template <class InputArchive>
void test_in(const InputArchive& iar) {
    const int n = 3;
    A a, an[n];
    B b, bn[n];
    C c, cn[n];
    int i, in[n];
    double *p = new double[n];
    A *q = new A[n];
    vector<int> v(3);
    pair<int,double> pp(33,99.0);
    map<short,double_complex> m;
    const char* teststr = "hello \n dude !";
    string str(teststr);
    linked_list list;
    Tensor<int> ti;
    Tensor<double> td;
    Tensor<double_complex> tc;
    double pi, e;

    // Destroy in-core data
    a.a = b.b = c.c = i = 0;
    for (int k=0; k<n; k++) {
        p[k] = q[k].a = an[k].a = v[k] = cn[k].c = in[k] = -1;
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
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " double *p wrapped" << std::endl);
    iar & wrap(p,n);
    iar >> wrap(p,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A *q wrapped" << std::endl);
    iar & wrap(q,n);
    iar >> wrap(q,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<int>" << std::endl);
    iar & v;
    iar >> v;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " pair<int,double>" << std::endl);
    iar & pp;
    iar >> pp;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " map<short,double<complex>>" << std::endl);
    iar & m;
    iar >> m;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " string" << std::endl);
    iar & str;
    iar >> str;

    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " Tensors" << std::endl);
    iar >> ti >> td >> tc;
    iar & ti & td & tc;

    iar & 1.0 & i & a & b & c & in & an & bn & cn & wrap(p,n) & wrap(q,n) & pp & m & str;
    // Test data
    bool status = true;

#define TEST(cond) status &= cond;  \
                   if (!(cond)) std::cout << #cond << " failed" << std::endl
    TEST(a.a == 1);
    TEST(b.b == 1);
    TEST(c.c == 1);
    TEST(i == 1);
    for (int k=0; k<n; k++) {
        TEST(an[k].a == k);
        TEST(bn[k].b == (k&1));
        TEST(cn[k].c == k);
        TEST(in[k] == k);
        TEST(p[k] == k);
        TEST(q[k].a == k);
        TEST(v[k] == k);
        TEST(m[k] == double_complex(k,k));
    }
    TEST(pp.first==33 && pp.second==99.0);
    TEST(str == string(teststr));
    //  TEST((ti-Tensor<int>(n,n).fillindex).normf() == 0);
    //TEST((td-Tensor<double>(n,n).fillindex).normf() == 0);
    //TEST((tc-Tensor<double_complex>(n,n).fillindex).normf() == 0);


#undef TEST

    if (status)
        std::cout << "Serialization appears to work." << std::endl;
    else
        std::cout << "Sorry, back to the drawing board.\n";
}

const double PI = 3.1415926535897932384;

double fred(double x, double y, double z) {
    double fac = std::pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*std::exp(-65.0*(x*x+y*y+z*z));
}

/*
template <class Archive> 
struct ArchiveStoreImpl<Archive, D> {
//template <class OutputArchive, class Function>
void Function_Save(const OutputArchive& oar, Function& ofunc) {
  MAD_ARCHIVE_DEBUG(std::cout << std::endl << " Function.k " << std::endl);
  oar << ofunc.coeff(tree);
  oar & ofunc.coeff(tree);
}
 
template <class InputArchive, class Function>
void Function_Load(const InputArchive& iar, Function& ifunc) {
  ifunc.k = 0;
  MAD_ARCHIVE_DEBUG(std::cout << std::endl << " Function.k " << std::endl);
  iar & ifunc.data[0];
  iar >> ifunc.data[0];
  cout << " test2 in Function_load " << ifunc.k << endl;
}
*/

/*
template <class Archive> 
struct ArchiveStoreImpl<Archive, D> {
  static inline void store(const Archive& ar, const D& d) {
     ar & d.get();
  };
};
 
template <class Archive>
struct ArchiveLoadImpl<Archive, D> {
  static inline void load(const Archive& ar, D& d) {
     int value;
     ar & value;
     d.set(value);
  };
};
*/

/*
template <class Archive, class OctTreeT>
//struct ArchiveLoadImpl<Archive,OctTreeT> {
  void _load(const Archive& ar, OctTreeT *tree) {
//    typedef Tensor<T> TensorT;
    bool flag;
    ar & flag;
    if(flag) {
      ar & flag;
      if(flag) {
        Tensor t;
        ar & t;
        set_coeff(tree, t);
      }
      FORIJK(_load(ar, tree->child(i, j, k)););
    }
  }
//}
 
template <class Archive, class OctTreeT>
//struct ArchiveStoreImpl<Archive, OctTreeT> {
  void _save(Archive& ar, const OctTreeT *tree) {
    ar & isactive(tree);
    if(isactive(tree)) {
      const Tensor *t = coeff(tree);
      ar & (*t != 0);
      if(*t) ar & *t;
    }
    FORIJK(ar & (tree->child(i,j,k)!=0););
  }
//}
*/

int main(int argc, char* argv[]) {
  madness::archive::archive_initialize_type_names();
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(A);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(B);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(C);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(linked_list);
   ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(madness::archive::pair_int_double);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(madness::archive::pair_short_complex_double);
  ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(madness::archive::map_short_complex_double);
  
  {
    MPI::Init(argc, argv);
    Communicator comm;
    redirectio(comm);
    comm.print();
    load_coeffs(comm);
    load_quadrature(comm);
    FunctionDefaults::tree = new FunctionOctTree(OctTree<FunctionNode>::create_default(comm,2));

    Function<double> ftest = FunctionFactory<double>(fred).k(3).refine(1).compress(1).initial_level(2).thresh(1e-5);
    Function<double> ftest2 = FunctionFactory<double>();
    const char* f = "tserialize.dat";
    ftest.compress();
    ftest.truncate();
    //ftest.save_local(oar);
    ftest.save(f, comm);
    //ftest.load_local(iar);
    ftest2.load(f, comm);
    cout << " class subtraction test " << (ftest - ftest2).norm2sq() << endl;

    comm.close();
    MPI::Finalize();
    }

    return 0;
}
