/*!
  \file examples/dirac-hatom.cc
  \brief solve the hydrogen atom using the 4-component Dirac equation

  example usage:  dirac-hatom --thresh=1.e-5 --k=7 --dirac="charge=100"


*/
#include <iostream>
#include<madchem.h>

using namespace madness;


enum Uplo {upper, lower};
static bool debug=false;
static double alpha1=constants::fine_structure_constant;
static double clight=1.0/alpha1;
static double shift=0.0;
static double epsilon=1.e-12;
static bool use_ble=false;
static bool use_bsp=false;
static double lo=1.e-7;
static double nucrad=1e-8;
static double nucexpt=3.0/nucrad; // actually the squareroot of the Gaussian exponent

real_function_3d mask;

class DiracParameters : public QCCalculationParametersBase {
public:
    /// ctor reading out the input file
    DiracParameters() {
        initialize<double>("shift",0.0);
	initialize<double>("nucrad",1e-8);
        initialize<double>("charge",10.0);
        initialize<double>("lo",1.e-7,"lo for the operator apply");
        initialize<bool>("use_ble",false);
        initialize<bool>("use_bsp",false);
        initialize<long>("nstates",1);
        initialize<long>("ansatz",0);
    }

    DiracParameters(World& world, const commandlineparser& parser) : DiracParameters() {
        // read input file
        read_input_and_commandline_options(world,parser,"dirac");
    }

    double shift() const {return get<double>("shift");}
    double nucrad() const {return get<double>("nucrad");}
    double charge() const {return get<double>("charge");}
    double lo() const {return get<double>("lo");}
    bool use_ble() const {return get<bool>("use_ble");}
    bool use_bsp() const {return get<bool>("use_bsp");}
    long nstates() const {return get<long>("nstates");}
    long ansatz() const {return get<long>("ansatz");}


};

double compute_gamma(const double nuclear_charge) {
    return sqrt(1-nuclear_charge*nuclear_charge*alpha1*alpha1);
}

struct stepfunction : public FunctionFunctorInterface<double,3> {
    int axis=-1;
    stepfunction(const int axis) : axis(axis) {
        MADNESS_CHECK(axis>=0 && axis<3);
    }
    double operator()(const coord_3d& r) const final { return r[axis]/(r.normf()+epsilon); }
    Level special_level() const final {return 20;};
    std::vector<Vector<double, 3UL>> special_points() const override {
        coord_3d o={0.0,0.0,0.0};
        return {o};
    }
};

struct ncf_cusp {
    double a=1.3;
    double Z=-1;
    ncf_cusp(double a, double Z) :a(a),Z(Z) {}
    double operator()(const double& r) const {
        if (a<0.0) return 1.0;
        return 1.0+(1.0/(a-1.0))*exp(-a*Z*r);
    }
};

struct Sigma_ncf_cusp {
    double a=1.3;
    double Z=-1;
    Sigma_ncf_cusp(double a, double Z) :a(a),Z(Z) {}
    double operator()(const double& r) const {
        if (a<0.0) return 0.0;
        const double denominator=1.0 + 1.0 / (a - 1.0) * exp(-a * Z * r);
        const double numerator=-a*Z/(a-1.0) * exp(-a * Z * r);
        return numerator/denominator;
    }
};

struct ncf_singularity {
    double gamma=0.0;
    ncf_singularity(double gamma) : gamma(gamma) {}
    double operator()(const double& r) const {
        if (gamma<0.0) return 1.0;
        return 1.0/(std::pow(r,1.0-gamma)+epsilon);
    }
};

struct ncf : public FunctionFunctorInterface<double,3> {
    ncf_singularity ncfs;
    ncf_cusp ncfc;
    long power=1;
    ncf(double gamma, double a, double Z) : ncfs(gamma), ncfc(a,Z) {}
    double operator()(const coord_3d& r) const override {
        double rr=r.normf();
        return std::pow(ncfs(rr) * ncfc(rr),power);
    }

    Level special_level() const final {return 20;};
    std::vector<Vector<double, 3UL>> special_points() const final {
        coord_3d o={0.0,0.0,0.0};
        return {o};
    }
};

double generalized_laguerre(const double alpha, const long n, const double r) {
    if (n<0) return 0.0;
    else if (n==0) return 1.0;
    else if (n==1) return (1.0+alpha-r);
    else if (n==2) return 0.5*r*r - (alpha + 2.0)*r + 0.5*(alpha+1)*(alpha+2);
    else
    MADNESS_EXCEPTION("generalized Laguerre polynomial implemented only up to order n=2",1);
}


/// returns the complex value of a given spherical harmonic
struct SphericalHarmonics{
    long l, m;
    bool zero=false;
    SphericalHarmonics(const long l, const long m) : l(l), m(m) {
        if (l<0) this->l=-l-1;
        if (abs(this->m)>this->l) zero=true;
    }

    double_complex operator()(const coord_3d& xyz) const {
        if (zero) return {0.0,0.0};
        // const double r=xyz.normf();
        // const double r2=r*r;
        // const double x=xyz[0];
        // const double y=xyz[1];
        // const double z=xyz[2];
        const double_complex i=double_complex(0.0,1.0);
        auto stepx=stepfunction(0);
        auto stepy=stepfunction(1);
        auto stepz=stepfunction(2);
        double_complex step_x_m_iy=stepx(xyz) - i* stepy(xyz);
        double_complex step_x_p_iy=stepx(xyz) + i* stepy(xyz);
        if ((l==0) and (m== 0)) return 0.5*sqrt(1.0/constants::pi);

//        if ((l==1) and (m==-1)) return 0.5*sqrt(1.5/constants::pi) * (x - i*y)/r;
//        if ((l==1) and (m== 0)) return 0.5*sqrt(3.0/constants::pi) * z/r;
//        if ((l==1) and (m== 1)) return -0.5*sqrt(1.5/constants::pi) * (x + i*y)/r;
        if ((l==1) and (m==-1)) return 0.5*sqrt(1.5/constants::pi) * (stepx(xyz) - i *stepy(xyz));
        if ((l==1) and (m== 0)) return 0.5*sqrt(3.0/constants::pi) * stepz(xyz);
        if ((l==1) and (m== 1)) return -0.5*sqrt(1.5/constants::pi) * (stepx(xyz) + i *stepy(xyz));

        if ((l==2) and (m==-2)) return  0.25*sqrt(7.5/constants::pi) * step_x_m_iy * step_x_m_iy;
        if ((l==2) and (m==-1)) return  0.5 *sqrt(7.5/constants::pi) * step_x_m_iy * stepz(xyz);
        if ((l==2) and (m== 0)) return  0.25*sqrt(5.0/constants::pi) * (3.0*stepz(xyz)*stepz(xyz)- 1.0);
        if ((l==2) and (m== 1)) return -0.5 *sqrt(7.5/constants::pi) * step_x_p_iy * stepz(xyz);
        if ((l==2) and (m== 2)) return  0.25*sqrt(7.5/constants::pi) * step_x_p_iy * step_x_p_iy;

        if ((l==3) and (m==-3)) return  0.125*sqrt(35/constants::pi) * step_x_m_iy * step_x_m_iy * step_x_m_iy;
        if ((l==3) and (m==-2)) return  0.25*sqrt(105/constants::pi) * step_x_m_iy * step_x_m_iy * stepz(xyz);
        if ((l==3) and (m==-1)) return  0.125*sqrt(21/constants::pi) * step_x_m_iy * (5*stepz(xyz)*stepz(xyz) - 1.0);
        if ((l==3) and (m== 0)) return  0.25*sqrt(7.0/constants::pi) * (5*stepz(xyz)*stepz(xyz)*stepz(xyz) - 3.0*stepz(xyz));
        if ((l==3) and (m== 1)) return -0.125*sqrt(21/constants::pi) * step_x_p_iy * (5*stepz(xyz)*stepz(xyz) - 1.0);
        if ((l==3) and (m== 2)) return  0.25*sqrt(105/constants::pi) * step_x_p_iy * step_x_p_iy * stepz(xyz);
        if ((l==3) and (m== 3)) return -0.125*sqrt(35/constants::pi) * step_x_p_iy * step_x_p_iy * step_x_p_iy;
        MADNESS_EXCEPTION("out of range in SphericalHarmonics",1);

        return double_complex(0.0,0.0);
    }

};

struct sgl_guess {
    long n,l,m;
    double Z;
    sgl_guess(const long n, const long l, const long m, const double Z) : n(n), l(l), m(m), Z(Z) {}

    double_complex operator()(const coord_3d& coord) const {
        double_complex Y=SphericalHarmonics(l,m)(coord);
        double r=coord.normf();
        double rho=2.0*Z*r/n;
        double R= generalized_laguerre(2.*l+1.,n-l-1,rho);
        double e=exp(-0.5*rho);
        return e*R*std::pow(rho,double(l))*Y;
    }

    double energy() const {
        return 0.5*Z*Z/(n*n);
    }

    complex_function_3d get_wf(World& world) const {
        std::vector<coord_3d> special_points(1,coord_3d({0.0,0.0,0.0}));
        complex_function_3d wf = complex_factory_3d(world).functor(*this).special_points(special_points);
        double norm=wf.norm2();
        wf.scale(1.0/norm);
        return wf;
    }



};


/// following Dyall: p 104f
struct Xi {
    long k, component, a;
    double m, l;
    Xi(const long k, const double m, const double component, const double j, const long l)
      : k(k),  component(component), a(0), m(m), l(l) {
        MADNESS_CHECK(component==0 or component==1);
        a=compute_a(j,l);
    }

    static long compute_a(const double j, const long l) {
        if (j>l) return 1;
        if (j<l) return -1;
        throw;
    }

    double_complex operator()(const coord_3d& xyz) const {
        double_complex result;
        double denominator2=2.0*l+1;
        if (component==0) {
            double numerator=a*sqrt((l+0.5+a*m)/denominator2);
            auto Y=SphericalHarmonics(k,lround(m-0.5));
            result=numerator*Y(xyz);
        } else if (component==1) {
            double numerator=sqrt((l+0.5-a*m)/denominator2);
            auto Y=SphericalHarmonics(k,lround(m+0.5));
            result=numerator*Y(xyz);
        }
        return result;
    }
};

struct Omega {
    long k, component;
    double m;
    Omega(const long k, const double m, const double component) : k(k), component(component), m(m) {
        MADNESS_CHECK(component==0 or component==1);
    }
    double_complex operator()(const coord_3d& xyz) const {
        double sgnk= (k>0) ? 1.0 : -1.0;
        double_complex result;
        if (component==0) {
            double nn=sqrt((k+0.5-m)/(2*k+1));
            auto Y=SphericalHarmonics(k,lround(m-0.5));
            result=nn*Y(xyz);
        } else if (component==1) {
            double nn=sqrt((k+0.5+m)/(2*k+1));
            auto Y=SphericalHarmonics(k,lround(m+0.5));
            result=-sgnk*nn*Y(xyz);
        }
        return result;
    }
};


double compute_electronic_energy(const double energy) {
    double c=1.0/alpha1;
    return energy-c*c-shift;
}


double mask1d(double x) {
    /* Iterated first beta function to switch smoothly
           from 0->1 in [0,1].  n iterations produce 2*n-1
           zero derivatives at the end points. Order of polyn
           is 3^n.

           Currently use one iteration so that first deriv.
           is zero at interior boundary and is exactly representable
           by low order multiwavelet without refinement */

    x = (x * x * (3. - 2. * x));
    return x;
}

double mask3d(const Vector<double,3>& ruser) {
    coord_3d rsim;
    user_to_sim(ruser, rsim);
    double x = rsim[0], y = rsim[1], z = rsim[2];
    double lo = 0.0625, hi = 1.0 - lo, result = 1.0;
    double rlo = 1.0 / lo;

    if (x < lo)
        result *= mask1d(x * rlo);
    else if (x > hi)
        result *= mask1d((1.0 - x) * rlo);
    if (y < lo)
        result *= mask1d(y * rlo);
    else if (y > hi)
        result *= mask1d((1.0 - y) * rlo);
    if (z < lo)
        result *= mask1d(z * rlo);
    else if (z > hi)
        result *= mask1d((1.0 - z) * rlo);

    return result;
}

template<typename T, std::size_t NDIM>
class MyDerivativeOperator : public SCFOperatorBase<T,NDIM> {
    typedef Function<T,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef Tensor<T> tensorT;

public:

    MyDerivativeOperator(World& world, const int axis1) : world(world), axis(axis1) {}
    void set_ble1() {ble=true;};
    void set_bspline1() {bsp=true;};

    std::string info() const {return "D"+std::to_string(axis);}

    functionT operator()(const functionT& ket) const {
        vecfuncT vket(1,ket);
        return this->operator()(vket)[0];
    }

    vecfuncT operator()(const vecfuncT& vket) const {
        auto gradop = free_space_derivative<T,NDIM>(world, axis);
        if (ble) gradop.set_ble1();
        else if (bsp) gradop.set_bspline1();
        vecfuncT dvket=apply(world, gradop, vket, false);
        world.gop.fence();
        return dvket;
    }

    T operator()(const functionT& bra, const functionT& ket) const {
        vecfuncT vbra(1,bra), vket(1,ket);
        Tensor<T> tmat=this->operator()(vbra,vket); 
        return tmat(0l,0l);
    }

    tensorT operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        vecfuncT dvket=this->operator()(vket);
        return matrix_inner(world,vbra,dvket, bra_equiv_ket);
    }

private:
    World& world;
    int axis;
    bool ble=false;
    bool bsp=false;
};


/// defines a 4-spinor
class Spinor {
public:
    vector_complex_function_3d components;
    World& world() const {return components.front().world();}
    Spinor() {
        components.resize(4);
    }

    Spinor(World& world) {
        components=zero_functions_compressed<double_complex,3>(world,4);
    }

    Spinor(const vector_complex_function_3d& components) : components(components){}

    Spinor& operator+=(const Spinor& other) {
        components+=other.components;
        return *this;
    }

    Spinor& operator-=(const Spinor& other) {
        components-=other.components;
        return *this;
    }

    Spinor operator+(const Spinor& other) const {
        Spinor result;
        result.components=copy(world(),components); // deep copy
        result.components+=other.components;
        return result;
    }

    Spinor operator-(const Spinor& other) const {
        Spinor result;
        result.components=copy(world(),components); // deep copy
        result.components-=other.components;
        return result;
    }

    Spinor& truncate() {
        double thresh = components[0].thresh();
        components[0].truncate(thresh,false);
        components[1].truncate(thresh,false);
        components[2].truncate(thresh/clight,false);
	components[3].truncate(thresh/clight);
        return *this;
    }

    void print_norms(std::string text) {
        auto norms = norm2s(world(), components);
        auto realnorms = norm2s(world(), real(components));
        auto imagnorms = norm2s(world(), imag(components));
	size_t sizes[4];
	for (size_t i=0; i<4; i++) sizes[i] = components[i].size();
        print(text,norms);
        print("  -- real,imag",realnorms,imagnorms);
	print("  -- sizes", sizes);
	size_t k = components[0].k();
	for (auto& s : sizes) s /= k*k*k;
	print("  -- nodes", sizes);
	
//        print("  -- moments  ",x1,y1,z1,x2,y2,z2);
//        plot(text);
    }

    void plot(const std::string filename) const {
        plot_plane(world(), real(components), "Re_"+filename);
        plot_plane(world(), imag(components), "Im_"+filename);
    }

    friend double_complex inner(const Spinor& bra, const Spinor& ket) {
        return inner(bra.components,ket.components);
    }
};

template<typename T>
Spinor operator*(const T fac, const Spinor& arg) {
    return Spinor(fac*arg.components);
}

template<typename T>
Spinor operator*(const Spinor& arg, const T fac) {
    return Spinor(fac*arg.components);
}


Spinor copy(const Spinor& other) {
    return Spinor(copy(other.world(),other.components));
}

std::vector<Spinor> copy(const std::vector<Spinor>& other) {
    std::vector<Spinor> result;
    for (auto& oo : other) result.push_back(copy(oo.world(),oo.components));
    return result;
}


template<typename T>
std::vector<Spinor> operator*(const std::vector<Spinor>& arg, const T fac) {
    std::vector<Spinor> result=copy(arg);
    for (auto& r : result) r=r*fac;
    return result;
}
double_complex inner(const std::vector<Spinor>& bra, const std::vector<Spinor> ket) {
    double_complex result=0.0;
    for (size_t i=0; i<bra.size(); ++i) result+=inner(bra[i],ket[i]);
    return result;
}

Tensor<double_complex> matrix_inner(const std::vector<Spinor>& bra, const std::vector<Spinor> ket) {
    Tensor<double_complex> result(bra.size(),ket.size());
    for (size_t i=0; i<bra.size(); ++i) {
        for (size_t j=0; j<ket.size(); ++j) {
            result(i,j)=inner(bra[i],ket[j]);
        }
    }
    return result;
}

std::vector<Spinor> operator-=(std::vector<Spinor>& left, const std::vector<Spinor>& right) {
    for (size_t i=0; i<right.size(); ++i) left[i]-=right[i];
    return left;
}

std::vector<Spinor> operator+=(std::vector<Spinor>& left, const std::vector<Spinor>& right) {
    for (size_t i=0; i<right.size(); ++i) left[i]+=right[i];
    return left;
}

std::vector<Spinor> operator-(const std::vector<Spinor>& left, const std::vector<Spinor>& right) {
    std::vector<Spinor> result=copy(left);
    result-=right;
    return result;
}


std::vector<Spinor> truncate(std::vector<Spinor> arg) {
    for (auto& a : arg) a.truncate();
    return arg;
}

struct LProjector {
    long lmax=3;
    World& world;
    std::map<std::pair<long,long>,complex_function_3d> Ylm;
    LProjector(World& world) : world(world) {
        for (long l=0; l<lmax; ++l) {
            for (long m=-l; m<l+1; ++m) {
                SphericalHarmonics Y(l,m);
                auto lm=std::make_pair(l,m);
                auto functor=[&Y](const coord_3d& r){return Y(r)*100.0*exp(-10.0*r.normf());};
                Ylm[lm]=complex_factory_3d(world).functor(functor);
            }
        }
    }

    void analyze(const Spinor& f, const std::string text="") const {
        madness::print(text);
        for (int c=0; c<4; ++c) {
            auto bla=f.components[c];
            LProjector lproj(world);
            lproj(bla,"component "+std::to_string(c));
        }
    }

    void operator()(const complex_function_3d& f, const std::string text="") const {
        madness::print("  "+text);
        for (long l=0; l<lmax; ++l) {
            for (long m = -l; m < l+1; ++m) {
                auto lm=std::make_pair(l,m);
                const complex_function_3d& Y=Ylm.find(lm)->second;
                double_complex ovlp=inner(Y,f);
                if (std::abs(ovlp)>1.e-7) madness::print("< lm | f> ",l,m,ovlp);
            }
        }
    }
};



// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
struct spinorallocator {
    World& world;
    const int n;

    /// @param[in]	world	the world
    /// @param[in]	nn		the number of functions in a given vector
    spinorallocator(World& world, const int n) : world(world), n(n) {
    }

    /// allocate a vector of n empty functions
   std::vector<Spinor> operator()() {
        std::vector<Spinor> r(n);
        for (auto & rr : r) rr=Spinor(world);
        return r;
    }
};

/// class defining an operator in matrix form, fixed to size (4,4)
class MatrixOperator {
public:
    typedef std::vector<std::pair<double_complex,std::shared_ptr<SCFOperatorBase<double_complex,3>>>> opT;

    explicit MatrixOperator(const int i=4, const int j=4) {
        elements.resize(i);
        for (auto& e : elements) e.resize(j);
    }

    std::size_t nrow() const {return elements.size();}
    std::size_t ncol() const {
        if (elements.size()) return elements[0].size();
        return 0;
    }

    virtual Spinor operator()(const Spinor& arg) const {
        World& world=arg.components[0].world();
        double_complex norm1=inner(arg,arg);
        Spinor result;
        for (auto& c : result.components) c=complex_factory_3d(world).compressed();
        for (size_t i=0; i<ncol(); ++i) {
            for (size_t j=0; j<nrow(); ++j) {
                const opT& ops=elements[i][j];
                for (const auto& op : ops) {
                    auto fac=op.first;
                    result.components[i]+=fac * (*op.second)(arg.components[j]);
                }
            }
        }
        double_complex norm2=inner(arg,arg);
        result.truncate();

        if (std::abs(norm1-norm2)/std::abs(norm1)>1.e-10) throw;

        return result;
    }

    virtual std::vector<Spinor> operator()(const std::vector<Spinor>& arg) const {
        std::vector<Spinor> result;
        for (auto& a : arg) result.push_back((*this)(a));
        return result;
    }


    /// add a submatrix to this

    /// @param[in]  istart   row where to add the submatrix
    /// @param[in]  jstart   column where to add the submatrix
    void add_submatrix(int istart, int jstart, const MatrixOperator& submatrix) {
        if (istart+submatrix.ncol()>ncol()) throw std::runtime_error("submatrix too large: too many columns");
        if (jstart+submatrix.nrow()>nrow()) throw std::runtime_error("submatrix too large: too many rows");
        for (size_t i=istart; i<istart+submatrix.ncol(); ++i) {
            for (size_t j=jstart; j<jstart+submatrix.nrow(); ++j) {
                for (auto& elem : submatrix.elements[i-istart][j-jstart]) elements[i][j].push_back(elem);
            }
        }
    }

    MatrixOperator& operator+=(const MatrixOperator& other) {
        for (size_t i=0; i<ncol(); ++i) {
            for (size_t j=0; j<nrow(); ++j) {
                for (auto& elem : other.elements[i][j]) elements[i][j].push_back(elem);
            }
        }
        return *this;
    }

    MatrixOperator operator+(const MatrixOperator& other) const {
        MatrixOperator result;
        for (size_t i=0; i<ncol(); ++i) {
            for (size_t j=0; j<nrow(); ++j) {
                for (auto& elem : this->elements[i][j]) result.elements[i][j].push_back(elem);
                for (auto& elem : other.elements[i][j]) result.elements[i][j].push_back(elem);
            }
        }
        return result;
    }

    void add_operator(const int i, const int j, const double_complex& fac, const std::shared_ptr<SCFOperatorBase<double_complex,3>> op) {
        elements[i][j].push_back(std::make_pair(fac,op));
    }

    void print(std::string name="") const {
        madness::print(name);
        for (size_t i=0; i<ncol(); i++) {
	  for (size_t j=0; j<nrow(); j++) {
                const opT& ops=elements[i][j];
                for (auto op : ops) {
                    auto fac=op.first;
                    std::cout << " " << fac << " " << (op.second)->info();
                }
                std::cout << " ||| ";
            }
            std::cout << std::endl;
        }
    }

    /// matrix containing prefactor and operator
    std::vector<std::vector<opT>> elements;
};

class Metric : public MatrixOperator{
public:
    double c0=1.0, c1=1.0, c2=1.0, c3=1.0;
    virtual Spinor operator()(const Spinor& arg) const {
        auto result=copy(arg);
        result.components[0].scale(c0);
        result.components[1].scale(c1);
        result.components[2].scale(c2);
        result.components[3].scale(c3);
        return result;
    }

    void print() const {
        madness::print("metric ",c0,c1,c2,c3);
    }
};


void show_norms(const Spinor& bra, const Spinor& ket, const std::string& name) {
    Metric m;
    auto norms = inner(ket.world(), m(bra).components, ket.components);
    madness::print("norms of ",name,": ",norms);
}

MatrixOperator make_Hdiag(World& world, const LocalPotentialOperator<double_complex,3>& V1) {
    MatrixOperator Hv;
    Hv.add_operator(0,0, 1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(V1));
    Hv.add_operator(1,1, 1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(V1));
    Hv.add_operator(2,2, 1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(V1));
    Hv.add_operator(3,3, 1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(V1));
    return Hv;
}

// Potential due to gaussian with exponent
double potn(double r) {
  double nr = r*nucexpt;
  if (nr > 1e-8) return std::erf(nr)/r;
  else return 2.0*nucexpt/std::sqrt(constants::pi); // to avoid divide by zero at origin
}

/// this is the nuclear potential on the diagonal
MatrixOperator make_Hv(World& world, const double nuclear_charge) {
    complex_function_3d V=complex_factory_3d(world)
      .functor([&nuclear_charge](const coord_3d& r){return double_complex(-nuclear_charge*potn(r.normf()));});
          //.functor([&nuclear_charge](const coord_3d& r){return double_complex(-nuclear_charge/(r.normf()+1.e-13));});
    auto V1=LocalPotentialOperator<double_complex,3>(world,"V",V);
    return make_Hdiag(world,V1);
}



/// returns a (2,2) matrix
MatrixOperator make_sp(World& world) {
    MatrixOperator sp(2,2);
    const double_complex ii=double_complex(0.0,1.0);
    const double alpha=constants::fine_structure_constant;
    const double c =  1.0/alpha;
    auto Dz=MyDerivativeOperator<double_complex,3>(world,2);
    auto Dx=MyDerivativeOperator<double_complex,3>(world,0);
    auto Dy=MyDerivativeOperator<double_complex,3>(world,1);
    if (use_ble) {
        MADNESS_CHECK(!use_bsp);
        Dz.set_ble1();
        Dx.set_ble1();
        Dy.set_ble1();
    }
    if (use_bsp) {
        MADNESS_CHECK(!use_ble);
        Dz.set_bspline1();
        Dx.set_bspline1();
        Dy.set_bspline1();
    }

    sp.add_operator(0,0,-c*ii,std::make_shared<MyDerivativeOperator<double_complex,3>>(Dz));
    sp.add_operator(0,1,-c*ii,std::make_shared<MyDerivativeOperator<double_complex,3>>(Dx));
    sp.add_operator(0,1,-c,std::make_shared<MyDerivativeOperator<double_complex,3>>(Dy));

    sp.add_operator(1,0,-c*ii,std::make_shared<MyDerivativeOperator<double_complex,3>>(Dx));
    sp.add_operator(1,0,c,std::make_shared<MyDerivativeOperator<double_complex,3>>(Dy));
    sp.add_operator(1,1,c*ii,std::make_shared<MyDerivativeOperator<double_complex,3>>(Dz));
    return sp;
}


MatrixOperator make_alpha_p(World& world) {
    MatrixOperator Hd;
    MatrixOperator sp = make_sp(world);
    Hd.add_submatrix(0, 2, sp);
    Hd.add_submatrix(2, 0, sp);
    return Hd;
}

/// this is c sigma p + beta m c^2

/// @param[in]  ll  scalar on the LL (1,1) block
/// @param[in]  ss  scalar on the SS (2,2) block
MatrixOperator make_Hd(World& world, const std::pair<double_complex,std::string>& ll,
                       const std::pair<double_complex,std::string>& ss) {
    MatrixOperator Hd=make_alpha_p(world);

    complex_function_3d V=complex_factory_3d(world).functor([](const coord_3d& r){return double_complex(1.0,0.0);});
    auto V1=LocalPotentialOperator<double_complex,3>(world,ll.second,V);
    auto V2=LocalPotentialOperator<double_complex,3>(world,ss.second,V);

    if (ll.first!=0.0) {
        Hd.add_operator(0,0,ll.first,std::make_shared<LocalPotentialOperator<double_complex,3>>(V1));
        Hd.add_operator(1,1,ll.first,std::make_shared<LocalPotentialOperator<double_complex,3>>(V1));
    }
    if (ss.first!=0.0) {
        Hd.add_operator(2, 2, ss.first, std::make_shared<LocalPotentialOperator<double_complex, 3>>(V1));
        Hd.add_operator(3, 3, ss.first, std::make_shared<LocalPotentialOperator<double_complex, 3>>(V1));
    }
    return Hd;
}

template<typename ansatzT>
std::vector<Spinor> schrodinger2dirac(const std::vector<complex_function_3d> wf, const ansatzT& ansatz, const double nuclear_charge) {
    World& world=wf.front().world();
    std::vector<Spinor> sgl;
    for (auto& w : wf ){
        Spinor tmp(world);
        tmp.components[0]=w;
        sgl.push_back(tmp);
    }
    MatrixOperator sp=make_alpha_p(world);

    std::vector<Spinor> result1;
    for (auto& s : sgl)  {
        Spinor tmp=(sp(s)+s).truncate();
        tmp.components[2].scale(0.5*alpha1*alpha1);
        tmp.components[3].scale(0.5*alpha1*alpha1);
        result1.push_back(tmp);
    }
    MatrixOperator Rinv=ansatz.Rinv(world);
    std::vector<Spinor> result=Rinv(result1);

    for (size_t i=0; i<result.size(); ++i) result[i].print_norms("dirac"+std::to_string(i));
    return result;
}

Tensor<double_complex> get_fock_transformation(World& world, const std::vector<Spinor>& spinors,
                                               const Tensor<double_complex>& overlap, const Tensor<double_complex>& fock) {

    const double thresh_degenerate=1.e-6;
    Tensor<double_complex> U;
    Tensor<double> evals;
    sygvp(world, fock, overlap, 1, U, evals);

    Tensor<double> occ(spinors.size());
    occ=1.0;
    Localizer::undo_reordering(U, occ, evals);
    Localizer::undo_degenerate_rotations(U, evals, thresh_degenerate);

    world.gop.broadcast(U.ptr(), U.size(), 0);
    world.gop.broadcast(evals.ptr(), evals.size(), 0);
    return U;
}

/// Transforms a vector of functions according to new[i] = sum[j] old[j]*c[j,i]

/// Uses sparsity in the transformation matrix --- set small elements to
/// zero to take advantage of this.
std::vector<Spinor>
transform(World& world, const std::vector<Spinor>& v,
          const Tensor<double_complex>& c) {

    PROFILE_BLOCK(Vtransformsp);
    int n = v.size();  // n is the old dimension
    int m = c.dim(1);  // m is the new dimension
    MADNESS_ASSERT(n==c.dim(0));

    std::vector<Spinor> vc(m);
    for (size_t i=0; i<vc.size(); ++i) vc[i]=Spinor(world);

    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            if (c(j,i) != double_complex(0.0)) gaxpy(world,double_complex(1.0),vc[i].components,
                                                     c(j,i),v[j].components,false);
        }
    }

    world.gop.fence(); // must always fence here to ensure gaxpy's are complete
    return vc;
}

std::vector<Spinor> orthonormalize_fock(const std::vector<Spinor>& arg,
                                        const std::vector<Spinor>& bra,
                                        Tensor<double_complex>& fock) {
    World& world=arg.front().world();
    Tensor<double_complex> ovlp(arg.size(),arg.size());
    for (size_t i=0; i<arg.size(); ++i) {
        for (size_t j=i; j<arg.size(); ++j) {
            ovlp(i,j)=inner(bra[i],arg[j]);
            ovlp(j,i)=std::conj(ovlp(i,j));
        }
    }

    Tensor<double_complex> U= get_fock_transformation(world,arg,ovlp,fock);
    fock=inner(conj(U),inner(fock,U),0,0);
//    print("fock in orthonormalize_fock");
//    print(fock);

    std::vector<Spinor> result=transform(world,arg,U);
    return result;

}

struct AnsatzBase {
public:
    [[nodiscard]] virtual std::string filename() const {return "ansatz"+this->name(); }
    [[nodiscard]] virtual std::string name() const =0;

    AnsatzBase(const double Z, const double a) : nuclear_charge(Z), a(a) {}
    int iansatz=0;
    double nuclear_charge=0.0;
    double a=-1.3;

    virtual void normalize(Spinor& bra, Spinor& ket) const {
        Metric m;
        double_complex norm2=inner(bra,m(ket));
        double norm=sqrt(real(norm2));
        scale(ket.world(),ket.components,1.0/norm);
        if (&bra!=&ket) scale(bra.world(),bra.components,1.0/norm);
    }

    virtual void normalize(Spinor& ket) const {
        auto bra=make_bra(ket);
        normalize(bra,ket);
    }

    virtual void normalize(std::vector<Spinor>& ket) const {
        auto bra=make_vbra(ket);
        normalize(bra,ket);
    }

    virtual void normalize(std::vector<Spinor>& bra, std::vector<Spinor>& ket) const {
        for (size_t i=0; i<ket.size(); ++i) normalize(bra[i],ket[i]);
    }


    virtual Spinor make_guess(World& world) const = 0;
    virtual MatrixOperator make_Hd(World& world) const = 0;
    virtual MatrixOperator R(World& world) const {
        MADNESS_EXCEPTION("no R implemented in this ansatz",1);
    }
    virtual MatrixOperator Rinv(World& world) const {
        MADNESS_EXCEPTION("no Rinv implemented in this ansatz",1);
    }

    virtual std::vector<Spinor> make_vbra(const std::vector<Spinor>& ket) const {
        std::vector<Spinor> result;
        for (const auto& k : ket) result.push_back(this->make_bra(k));
        return truncate(result);
    }
    [[nodiscard]] virtual Spinor make_bra(const Spinor& ket) const = 0;

    [[nodiscard]] virtual double mu(const double energy) const {
        return sqrt(-energy*energy*alpha1*alpha1 + 1.0/(alpha1*alpha1));
    }
    [[nodiscard]] double get_cusp_a() const {return a;}
};

struct Ansatz0 : public AnsatzBase {
public:

    Ansatz0(const double nuclear_charge) : AnsatzBase(nuclear_charge,1.0) {
        this->a=-1.0;
        iansatz=0;
    };

    [[nodiscard]] std::string name() const {
        return "0";
    }
    Spinor make_guess(World& world) const {
        Spinor result;
        const double_complex ii(0.0,1.0);
        const double_complex one(1.0,0.0);
        const double n=1;
        const double Z=double(nuclear_charge);
        const double alpha=constants::fine_structure_constant;
        const double gamma=compute_gamma(nuclear_charge);
//        print("gamma-1",gamma-1.0);
        const double C=nuclear_charge/n;
        // m=0.5;
        result.components[0]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return std::pow(r.normf(),gamma-1.0)*one*(1+gamma)*exp(-C*r.normf());});
        result.components[1]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return std::pow(r.normf(),gamma-1.0)*0.0*one;});
        result.components[2]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return std::pow(r.normf(),gamma-1.0)*ii*Z*alpha*r[2]/r.normf()*exp(-C*r.normf());});
        result.components[3]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return std::pow(r.normf(),gamma-1.0)*ii*Z*alpha*(r[0] + ii*r[1])/r.normf()*exp(-C*r.normf());});
        // m=-0.5;
//        print("make_guess with m=-0.5");
//        result.components[0]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return std::pow(r.normf(),gamma-1.0)*0.0*one;});
//        result.components[1]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return std::pow(r.normf(),gamma-1.0)*one*(1+gamma)*exp(-C*r.normf());});
//        result.components[2]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return std::pow(r.normf(),gamma-1.0)*ii*Z*alpha*(r[0] - ii*r[1])/r.normf()*exp(-C*r.normf());});
//        result.components[3]=complex_factory_3d(world).functor([&Z,&gamma,&alpha,&C,&ii,&one](const coord_3d& r){return -std::pow(r.normf(),gamma-1.0)*ii*Z*alpha*r[2]/r.normf()*exp(-C*r.normf());});
        return result;
    }

    MatrixOperator make_Hv(World& world) const {
        return ::make_Hv(world,nuclear_charge);
    }

    Spinor make_bra(const Spinor& ket) const {
        return Spinor(copy(ket.world(),ket.components));
    }
    MatrixOperator make_Hd(World& world) const {
        double c2=1.0/(alpha1*alpha1);
        return ::make_Hd(world,{c2,"mc2"},{-c2,"-mc2"});
    }

//    double mu(const double energy) const {
//        return sqrt(-energy*energy*alpha1*alpha1 + 1.0/(alpha1*alpha1));
//    }
    double energy() const {
        return compute_gamma(nuclear_charge)/(alpha1*alpha1);
    }

    MatrixOperator R(World& world) const {
        complex_function_3d one1=complex_factory_3d(world).functor([](const coord_3d& r) {return double_complex(1.0,0.0);});
        auto one = LocalPotentialOperator<double_complex, 3>(world, "1" , one1);
        return make_Hdiag(world,one);
    }
    MatrixOperator Rinv(World& world) const {
        complex_function_3d one1=complex_factory_3d(world).functor([](const coord_3d& r) {return double_complex(1.0,0.0);});
        auto one = LocalPotentialOperator<double_complex, 3>(world, "1" , one1);
        return make_Hdiag(world,one);
    }

};

MatrixOperator moments(World& world, int axis, int order) {
    MatrixOperator result;
    //int a=axis;
    //int o=order;
    complex_function_3d m=complex_factory_3d(world)
            .functor([&axis, &order](const coord_3d& r){return double_complex(std::pow(r[axis],double(order)),0.0);});
    auto mm=LocalPotentialOperator<double_complex,3>(world,"moment",m);
    result.add_operator(0,0,1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(mm));
    result.add_operator(1,1,1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(mm));
    result.add_operator(2,2,1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(mm));
    result.add_operator(3,3,1.0,std::make_shared<LocalPotentialOperator<double_complex,3>>(mm));
    return result;
}

struct ExactSpinor : public FunctionFunctorInterface<double_complex,3> {
    long n, k, l;
    mutable int component=0;
    double E=0.0, C=0.0, gamma=0.0, Z, j, m;
    bool regularized=true;
    double cusp_a=-1.0;
    bool compute_F=false;
    ExactSpinor(const int n, const char lc, const double j, const int Z, const double m=0.0)
            : ExactSpinor(n, l_char_to_int(lc),j,Z,m) { }
    ExactSpinor(const int n, const int l, const double j, const int Z,const double m=0.0)
        : n(n), l(l), Z(Z), j(j), m(m) {
        if (std::abs(j-(l+0.5))<1.e-10) k=lround(-j-0.5);       // j = l+1/2
        else k=lround(j+0.5);

        if (m==0.0) this->m=j;
        gamma=sqrt(k*k - Z*Z*alpha1*alpha1);
        E= compute_energy();
        C=compute_C();
    }

    std::string l_to_string(const long l) const {
        if (l==0) return "S";
        if (l==1) return "P";
        if (l==2) return "D";
        return "failed";
    }

    std::string filename() const {
        return "es_"+std::to_string(n)+l_to_string(l)+std::to_string(j)+"_m"+std::to_string(m);
        return "es_n"+std::to_string(n)+"_k"+std::to_string(k)+"_j"+std::to_string(j)+"_m"+std::to_string(m);
    }

    void set_ansatz(const AnsatzBase& ansatz) {
        compute_F = (ansatz.iansatz == 3);
        cusp_a=ansatz.get_cusp_a();
        regularized= !(ansatz.iansatz == 0);
    }

    static int l_char_to_int(const char lc) {
        int ll=0;
        if (lc=='S') ll=0;
        else if (lc=='P') ll=1;
        else if (lc=='D') ll=2;
        else {
            MADNESS_EXCEPTION("confused L quantum in ExactSpinor",1);
        }
        return ll;
    }

    Level special_level() const final {return 20;};
    std::vector<Vector<double, 3UL>> special_points() const final {
        coord_3d o={0.0,0.0,0.0};
        return {o};
    }

    double compute_energy() const {
        double c=1.0/alpha1;
        MADNESS_CHECK(gamma!=0.0);
        double E=c * c * 1.0/sqrt(1.0 + Z*Z/(c*c)*std::pow(n-std::abs(k)+gamma,-2.0));
        return E;
    }
    double compute_C() const {
        MADNESS_CHECK(E!=0.0);
        const double c=1.0/alpha1;
        return sqrt(c*c - E*E/(c*c));
    }
    // the energy dependent term for the spinor
    double compute_en() const {
        const double c=1.0/alpha1;
        return (gamma * c*c - k*get_energy()) / (c*C);
    }

    double get_energy() const {
        return E;
    }

    void print() const {
        madness::print("exact solution for n=",n,", k=",k, "j=",j, "m=",m);
        char lc='S';
        if (l==1) lc='P';
        if (l==2) lc='D';
        MADNESS_CHECK(l<3);
        madness::print("term symbol",n,lc,j);
        madness::print("energy =    ",E, E-1.0/(alpha1*alpha1));
        madness::print("compute_F   ",compute_F);
        madness::print("regularized ",regularized);
        madness::print("cusp_a      ",cusp_a);
        madness::print("C        ",C);
    }

    double_complex operator()(const coord_3d& c) const override {
        if (compute_F) return Fvalue(c);
        else return psivalue(c);
    }

    double_complex Fvalue(const coord_3d& coord) const {
        double r=coord.normf();
        double rho=2*C*r;
        double gamma1= compute_gamma(Z);
        ncf_cusp cusp(cusp_a, Z);
        double radial=exp(-rho*0.5)/cusp(r);
        if (k*k>1) radial*=std::pow(rho,gamma-gamma1);

        double_complex i={0.0,1.0};

        const long absk=std::abs(k);
        const double c=1.0/alpha1;
        const double Lnk1= generalized_laguerre(2*gamma+1.0,n-absk-1,rho);
        const double Lnk = generalized_laguerre(2*gamma-1.0,n-absk  ,rho);
        const double En=compute_en();

        const double G=Z/c * (gamma1 + gamma + 1 -k)*rho*Lnk1 + (gamma1+1)* (gamma - gamma1 -k +1)*En *Lnk;
        const double_complex F=i*(gamma1+1)*(gamma1+gamma - 1 - k) * rho * Lnk1 + i*Z/c *(gamma1- gamma + 1 + k) *En *Lnk;

        if (component==0) {
            return radial * G * Omega(k,m,0)(coord);
        } else if (component==1) {
            return radial * G * Omega(k,m,1)(coord);
        } else if (component==2) {
            return radial * F * Omega(-k,m,0)(coord);
        } else if (component==3) {
            return radial * F * Omega(-k,m,1)(coord);
        }
        MADNESS_EXCEPTION("confused component in ExactSpinor::Fvalue",1);
        return 0.0;
    }

    double_complex psivalue(const coord_3d& coord) const {
        double r=coord.normf();
        double rho=2*C*r;
        double radial=1.0;
        ncf_cusp ncf(cusp_a,Z);
        radial*=std::pow(2*C,gamma)*exp(-0.5*rho)/ncf(r);
        // three cases for R^{-1} * r^{gamma_k-1}:
        // 1. regularization with R=r^{gamma1-1}, |k|==1: factor: 1
        // 2. regularization with R=r^{gamma1-1}, |k|>1: factor: r^{gamma_k-gamma1}
        // 3. no regularization, factor r^{gamma_k-1};
        double gamma1= compute_gamma(Z);
        long absk= std::abs(k);
        if (regularized) {
            if (absk>1) {
                radial *= std::pow(r,gamma-gamma1);          // exponent is positive
            }
        } else {
            ncf_singularity ncf_s(gamma);
            radial *= ncf_s(r);
        }


        const double Lnk1= generalized_laguerre(2*gamma+1.0,n-absk-1,rho);
        const double Lnk = generalized_laguerre(2*gamma-1.0,n-absk  ,rho);
        const double En=compute_en();
        const double c=1.0/alpha1;

        double_complex i={0.0,1.0};
        double g=radial * (Z/c*rho* Lnk1 + (gamma - k)*En * Lnk);
        double f=radial * ((gamma-k)*rho*Lnk1 + Z/c*En * Lnk);
//        double f=Z*alpha1*radial;
        //double sgnk= (k>0) ? 1.0 : -1.0;

//        return angular(coord,g,f);
//        if (component==0) {
//            return g * Xi(k,m,0,j,l)(coord);
//        } else if (component==1) {
//            return g * Xi(k,m,1,j,l)(coord);
//        } else if (component==2) {
//            return i * f * Xi(-k,m,0,j,l)(coord);
//        } else if (component==3) {
//            return i * f * Xi(-k,m,1,j,l)(coord);
//        }

        if (component==0) {
            return g * Omega(k,m,0)(coord);
        } else if (component==1) {
            return g * Omega(k,m,1)(coord);
        } else if (component==2) {
            return i * f * Omega(-k,m,0)(coord);
        } else if (component==3) {
            return i * f * Omega(-k,m,1)(coord);
        }

        MADNESS_EXCEPTION("confused component in ExactSpinor",1);
        return {0.0,0.0};
    }

    double_complex angular(const coord_3d& c, const double g, const double f) const {

        double_complex i={0.0,1.0};
        double_complex prefac =std::pow(i,l)*std::pow(-1.0,m+0.5);
        if (component==0) {   // j = l+1/2 : k=-j-0.5 == j=1/2 ; k=-1
            double nn = (l==lround(j-0.5)) ? -sqrt((j+m)/(2.0*j)) : sqrt((j-m+1)/(2*j+2));
            return prefac * g * nn *SphericalHarmonics(l,lround(m-0.5))(c);
//            return g/r * sqrt(double_complex((k + 0.5 - m)/(2.0*k + 1))) *SphericalHarmonics(k,lround(m-0.5))(c);
        } else if (component==1) {
            double nn = (l==lround(j-0.5)) ? sqrt((j-m)/(2.0*j)) : sqrt((j+m+1)/(2*j+2));
            return prefac * g * nn *SphericalHarmonics(l,lround(m+0.5))(c);
//            return -g/r * sgnk* sqrt(double_complex((k + 0.5 + m)/(2.0*k + 1))) *SphericalHarmonics(k,lround(m+0.5))(c);
        } else if (component==2) {
            double nn = (l==lround(j-0.5)) ? sqrt((j-m+1)/(2.0*j+2.0)) : sqrt((j+m)/(2.0*j));
            long ll=  (l==lround(j-0.5))  ? l+1 : l-1;
            return prefac * i*f * nn *SphericalHarmonics(ll,lround(m-0.5))(c);
//            return i*f/r * sqrt(double_complex((-k + 0.5 - m)/(-2.0*k + 1))) *SphericalHarmonics(-k,lround(m-0.5))(c);
        } else if (component==3) {
            double nn = (l==lround(j-0.5)) ? sqrt((j+m+1)/(2.0*j+2.0)) : -sqrt((j-m)/(2.0*j));
            long ll=  (l==lround(j-0.5))  ? l+1 : l-1;
            return prefac * i*f * nn *SphericalHarmonics(ll,lround(m+0.5))(c);
//            return -i*f/r * sgnk* sqrt(double_complex((-k + 0.5 - m)/(-2.0*k + 1))) *SphericalHarmonics(-k,lround(m+0.5))(c);
        }
        MADNESS_EXCEPTION("confused component in ExactSpinor::angular",1);
        return {0.0,0.0};
    }

    Spinor get_spinor(World& world) const {
        Spinor spinor;
        component=0;
        spinor.components[0]=complex_factory_3d(world).functor(*this);
        component=1;
        spinor.components[1]=complex_factory_3d(world).functor(*this);
        component=2;
        spinor.components[2]=complex_factory_3d(world).functor(*this);
        component=3;
        spinor.components[3]=complex_factory_3d(world).functor(*this);
        return spinor;
    }


};


template<typename AnsatzT>
void orthonormalize(std::vector<Spinor>& arg, const AnsatzT ansatz) {
    World& world=arg.front().world();

    double maxq;
    do {
        auto bra=ansatz.make_vbra(arg);
        ansatz.normalize(bra,arg);
        Tensor<double_complex> S=matrix_inner(bra,arg);
        Tensor<double_complex> Q = NemoBase::Q2(S);
        maxq=0.0;
        for (int i=0; i<Q.dim(0); ++i)
            for (int j=0; j<i; ++j)
                maxq = std::max(maxq,std::abs(Q(i,j)));
        arg = transform(world, arg, Q);
        truncate(arg);
    } while (maxq>0.01);

    auto bra=ansatz.make_vbra(arg);
}



template<typename ansatzT>
Spinor apply_bsh(ansatzT& ansatz, const MatrixOperator& Hd, const MatrixOperator& Hv, const MatrixOperator& metric,
                 const Spinor& spinor, const double energy) {
    World& world=spinor.world();
    timer t(world);

    const double mu=ansatz.mu(energy);
    double fac=  alpha1*alpha1;
    print("energy, mu, bsh_prefac in bsh: ",energy, mu,fac);

    auto bra=ansatz.make_bra(spinor);
    if (debug) show_norms(bra,spinor,"spinor before vpsi");
    auto vpsi=-2.0*Hv(spinor);

    vpsi.components[2] *= clight; // Scale small components so that integral operator maintains precision
    vpsi.components[3] *= clight;

    vpsi.truncate();
    if (debug) show_norms(bra,vpsi,"norms of vpsi");
    if (debug) show_norms(ansatz.make_bra(vpsi),vpsi,"<vpsi | vpsi>");
    t.tag("Vpsi");

    auto g=BSHOperator<3>(world,mu,lo,FunctionDefaults<3>::get_thresh());
    auto gvpsi1=apply(world,g,vpsi.components);
    gvpsi1 = truncate(gvpsi1); /// Truncate before undoing scaling
    
    vpsi.components[2] *= 1.0/clight; // Undo scaling
    vpsi.components[3] *= 1.0/clight;
    gvpsi1[2] *= 1.0/clight;
    gvpsi1[3] *= 1.0/clight;
    t.tag("GVpsi");

    auto gvpsi=Spinor(gvpsi1);
    
    if (debug) show_norms(ansatz.make_bra(gvpsi),gvpsi,"|| gvpsi ||");

    auto result1=Hd(gvpsi);
    if (debug) show_norms(ansatz.make_bra(result1),result1,"|| Hd(gvpsi) ||");

    auto result2=energy*metric(gvpsi);
    if (debug) show_norms(ansatz.make_bra(result2),result2,"|| energy*N*(gvpsi) ||");
    Spinor result=0.5*fac*(result1 + result2);

    if (debug) show_norms(ansatz.make_bra(result),result,"<result | result>");
    t.tag("HdGVpsi");

    return result.truncate();
}

void apply_mask(std::vector<Spinor>& v) {
  for (auto& s : v) {
    for (auto& f : s.components) {
      f = f*mask;
    }
  }
}

template<typename AnsatzT>
std::vector<Spinor> iterate(const std::vector<Spinor>& input, const std::vector<double> energy, const AnsatzT& ansatz, const int maxiter) {

    print_header2("starting iterations");
    World& world=input.front().world();
    spinorallocator alloc(world,input.size());
    XNonlinearSolver<std::vector<Spinor>,double_complex,spinorallocator> solver(alloc);
    solver.set_maxsub(5); // was 3
    solver.do_print=true;

    auto Hv=ansatz.make_Hv(world);
    auto Hd=ansatz.make_Hd(world);
    auto H=Hd+Hv;
    auto metric=  Metric();
    std::vector<Spinor> current=copy(input);
    orthonormalize(current,ansatz);
    for (int iter=0; iter<maxiter; ++iter) {
        if (iter<2) solver.clear_subspace();    // start KAIN after 1st iteration since initial guess might be not very smooth
        double wall0=wall_time();

        print_header3("Iteration "+std::to_string(iter));
        orthonormalize(current,ansatz);
        std::vector<Spinor> newpsi;
        for (auto& n : current) n.print_norms("norm of current spinor");
        for (size_t i=0; i<current.size(); ++i) newpsi.push_back(apply_bsh(ansatz,Hd,Hv,metric,current[i],energy[i]));
        for (auto& n : newpsi) n.print_norms("norm of updated spinor");
        auto residual=truncate(current-newpsi);
        double res=0.0;
        for (const auto& r : residual) res+=norm2(world,r.components);
        newpsi=truncate(solver.update(current,residual,1.e-4,3));
	apply_mask(newpsi);
        orthonormalize(newpsi,ansatz);
	auto bra=ansatz.make_vbra(newpsi);
        Tensor<double_complex> fock;
	{
	  auto Hpsi = H(newpsi);
	  auto Hvpsi = Hv(newpsi);
	  fock =matrix_inner(bra,Hpsi);
	
	  Vector<double,3> lo = {0.0, 0.0, 0.0};
	  Vector<double,3> hi = {0.0, 0.0, 0.5*FunctionDefaults<3>::get_cell()(2,1)};
	  Hpsi[0].components.insert(Hpsi[0].components.end(), Hvpsi[0].components.begin(), Hvpsi[0].components.end());
	  Hpsi[0].components.insert(Hpsi[0].components.end(), newpsi[0].components.begin(), newpsi[0].components.end());
	  plot_line("Hpsi.dat", 10001, lo, hi, Hpsi[0].components);
	}
	
        fock+=conj_transpose(fock);
        fock*=0.5;
        newpsi=truncate(orthonormalize_fock(newpsi,bra,fock));
        bra=ansatz.make_vbra(newpsi);
        std::vector<double> energy_differences;
        for (size_t i=0; i<current.size(); ++i) {
//            newpsi[i].plot("psi"+std::to_string(i)+"_iteration"+std::to_string(iter)+"_ansatz"+ansatz.filename());
            double en=real(inner(bra[i],H(newpsi[i])));
//            show_norms(bra,H(newpsi),"energy contributions");
            double el_energy=compute_electronic_energy(en);
            double exact_el_energy=compute_electronic_energy(energy[i]);
            double diff=(el_energy-exact_el_energy);
            energy_differences.push_back(diff);
            printf("energy, el. energy, exact el. energy, difference %12.8f %12.8f %12.8f %4.1e\n", en, el_energy,exact_el_energy,diff);
        }
        current=newpsi;
        double wall1=wall_time();
        printf("elapsed time in iteration %2d: %6.2f with error %4.1e \n",iter,wall1-wall0,res );
//        printf("elapsed time in iteration %2d: %6.2f with energy/diff %12.8f %.2e \n",iter,wall1-wall0,compute_electronic_energy(en),
//               compute_electronic_energy(en) - compute_electronic_energy(energy));
    }
    return current;
}

template<typename ansatzT>
void run(World& world, ansatzT ansatz, const int nuclear_charge, const commandlineparser& parser, const int nstates) {

  //double thresh=FunctionDefaults<3>::get_thresh();
  //long tmode=FunctionDefaults<3>::get_truncate_mode();
//    Nemo nemo(world,parser);
//    if (world.rank()==0) nemo.get_param().print("dft","end");

    mask = real_factory_3d(world).f(&mask3d);
    mask.truncate();
    mask.reconstruct();

    std::vector<Spinor> guesses;
    std::vector<double> energies;

//    guesses.push_back(guess);
//    energies.push_back(ansatz.energy());
    ExactSpinor psi1s_half=ExactSpinor(1,'S',0.5,nuclear_charge,0.5);
    ExactSpinor psi1s_mhalf=ExactSpinor(1,'S',0.5,nuclear_charge,-0.5);
    ExactSpinor psi2s_half=ExactSpinor(2,'S',0.5,nuclear_charge,0.5);
    ExactSpinor psi2s_mhalf=ExactSpinor(2,'S',0.5,nuclear_charge,-0.5);
    ExactSpinor psi2p1_half  =ExactSpinor(2,'P',0.5,nuclear_charge, 0.5);
    ExactSpinor psi2p1_mhalf =ExactSpinor(2,'P',0.5,nuclear_charge,-0.5);
    ExactSpinor psi2p2_thalf =ExactSpinor(2,'P',1.5,nuclear_charge, 1.5);
    ExactSpinor psi2p2_half  =ExactSpinor(2,'P',1.5,nuclear_charge, 0.5);
    ExactSpinor psi2p2_mhalf =ExactSpinor(2,'P',1.5,nuclear_charge,-0.5);
    ExactSpinor psi2p2_mthalf=ExactSpinor(2,'P',1.5,nuclear_charge,-1.5);
//    std::vector<ExactSpinor> states ={psi1s,psi2p};
    std::vector<ExactSpinor> states ={psi1s_half,psi1s_mhalf,       // 1S 1/2
                                      psi2s_half,psi2s_mhalf,       // 2S 1/2
                                      psi2p1_half,psi2p1_mhalf,     // 2P 1/2
                                      psi2p2_mhalf,psi2p2_thalf,psi2p2_mhalf,psi2p2_mthalf}; // 2P 3/2


    sgl_guess sgl_1s=sgl_guess(1,0,0,nuclear_charge);
    sgl_guess sgl_2s=sgl_guess(2,0,0,nuclear_charge);
    sgl_guess sgl_2p0=sgl_guess(2,1,0,nuclear_charge);
    sgl_guess sgl_2p1=sgl_guess(2,1,1,nuclear_charge);
    sgl_guess sgl_2pm1=sgl_guess(2,1,-1,nuclear_charge);
    sgl_guess sgl_3s=sgl_guess(3,0,0,nuclear_charge);
    std::vector<sgl_guess> sgl_states={sgl_1s,sgl_2s,sgl_2p0,sgl_2p1,sgl_2pm1,sgl_3s};
    std::vector<Spinor> guess;
    const bool sglguess=true;
    if (sglguess) {
        print("\nUsing Schroedinger guess\n");

        std::vector<complex_function_3d> wf;
        for (int i=0; i<nstates; ++i) wf.push_back(sgl_states[i].get_wf(world));
        guess= schrodinger2dirac(wf,ansatz,nuclear_charge);
    } else {
        print("\nUsing exact spinor guess\n");
        for (auto& state: states) {
            state.set_ansatz(ansatz);
            guess.push_back(state.get_spinor(world));
            state.print();
            guess.back().print_norms("guess");

        }
    }
    apply_mask(guess);
    orthonormalize(guess,ansatz);
    auto bra=ansatz.make_vbra(guess);
    Tensor<double_complex> S=matrix_inner(bra,guess);
    print("initial overlap after  orthonormalization");
    print(S);


//    for (ExactSpinor& state : states) {
    for (int i=0; i<nstates; ++i) {
        ExactSpinor& state=states[i];
        state.compute_F=true;
        state.cusp_a=ansatz.get_cusp_a();
//        Spinor guess=state.get_spinor(world);
        Spinor spinorguess=guess[i];
        guesses.push_back(spinorguess);
        energies.push_back(state.get_energy());
        ansatz.normalize(guesses.back());
    }


    //const double alpha=constants::fine_structure_constant;
    //const double c=1.0/alpha;
    //const double gamma= compute_gamma(nuclear_charge);
    //double electronic_energy=gamma*c*c - c*c;
    //double energy=ansatz.energy();

    auto Hv=ansatz.make_Hv(world);
    auto Hd=ansatz.make_Hd(world);
    auto H=Hv+Hd;
    print_header2("Hamilton matrices");
    Hv.print("Hamiltonian Hv");
    Hd.print("Hamiltonian Hd");
    H.print("Hamiltonian Hd+Hv");

    auto result=iterate(guesses,energies,ansatz,10);

}

template<typename ansatzT>
void eigenvector_test(World& world, const ansatzT ansatz, ExactSpinor es) {
    print("=============================================================");
    print("Ansatz", ansatz.name());
    LProjector lproj(world);

    es.set_ansatz(ansatz);
    es.print();
    auto exactF = es.get_spinor(world);
    ansatz.normalize(exactF);
    exactF.print_norms("exactf normalized");
    exactF.plot("exact"+es.filename()+ansatz.filename());
    lproj.analyze(exactF,"ExactSpinor");

    auto exactF1 = ansatz.make_guess(world);
    ansatz.normalize(exactF1);
    lproj.analyze(exactF1,"ansatz.make_guess");
    exactF1.print_norms("make_guess normalized");
    auto diff1=exactF-exactF1;
    diff1.print_norms("exactf-make_guess normalized");

    Spinor spinor = copy(exactF);

    auto Hv = ansatz.make_Hv(world);
    auto Hd = ansatz.make_Hd(world);
    auto H = Hv + Hd;
    H.print("H");
    Hd.print("Hd");
    Hv.print("Hv");



    Spinor bra = ansatz.make_bra(spinor);
    ansatz.normalize(bra, spinor);
    auto norms = norm2s(world, spinor.components);

    print("");
    auto Hdspinor = Hd(spinor);
    auto Hvspinor = Hv(spinor);
    auto Hspinor = H(spinor);
    auto hnorms = norm2s(world, Hspinor.components);
    auto energy_norms=norms;
    for (auto& c : energy_norms) c*=es.get_energy();
    print("E * component norms", energy_norms);
    print("H(spinor) component norms", hnorms);
    auto en = inner(bra, Hspinor);

    auto diff = Hspinor - en * spinor;
    spinor.print_norms("spinor");
    Hspinor.print_norms("Hspinor");
    Hdspinor.print_norms("Hdspinor");
    Hvspinor.print_norms("Hvspinor");
    diff.print_norms("diff_Hspinor_en_spinor");
    double c=1.0/alpha1;
    print("energy", en, real(en - c * c), "difference", real(en - c * c) - (es.get_energy() - c * c));print("");
    print("");
}

int main(int argc, char* argv[]) {
    World& world=initialize(argc,argv);
    if (world.rank()==0) {
        print("\n");
        print_centered("Dirac hydrogen atom");
    }
    startup(world,argc,argv,true);
    if (world.rank()==0) print(madness::info::print_revision_information());


    commandlineparser parser(argc,argv);
//    parser.set_keyval("dft","'k=8'");
    if (world.rank()==0) {
        print("\ncommand line parameters");
        parser.print_map();
    }

    FunctionDefaults<3>::set_cubic_cell(-20,20);
    FunctionDefaults<3>::set_k(12);
    FunctionDefaults<3>::set_thresh(1.e-10);
    int tmode=1;
    if (parser.key_exists("k")) FunctionDefaults<3>::set_k(atoi(parser.value("k").c_str()));
    if (parser.key_exists("thresh")) FunctionDefaults<3>::set_thresh(atof(parser.value("thresh").c_str()));
    if (parser.key_exists("L")) FunctionDefaults<3>::set_cubic_cell(-atof(parser.value("L").c_str()),atof(parser.value("L").c_str()));
    if (parser.key_exists("truncate_mode")) tmode=atoi(parser.value("truncate_mode").c_str());
    FunctionDefaults<3>::set_truncate_mode(tmode);

    DiracParameters parameters(world,parser);
    parameters.print("dirac");
    use_ble=parameters.use_ble();
    use_bsp=parameters.use_bsp();
    nucrad=parameters.nucrad();
    nucexpt=3.0/nucrad;
    print("square root of nuclear exponent", nucexpt);
    lo=parameters.lo();

    print("\nCalculation parameters");
    print("thresh      ",FunctionDefaults<3>::get_thresh());
    print("k           ",FunctionDefaults<3>::get_k());
    print("trunc mode  ",FunctionDefaults<3>::get_truncate_mode());
    print("cell        ",FunctionDefaults<3>::get_cell_width());


    const double alpha=constants::fine_structure_constant;
    const double c=1.0/alpha;
    const double gamma= compute_gamma(parameters.charge());
    print("speed of light",c);
    print("fine structure constant",alpha);
    //const int k=1;
    print("gamma",gamma);
    double energy_exact=gamma*c*c - c*c;
    print("1s energy for Z=",parameters.charge(),": ",energy_exact);

    ExactSpinor es1s(1,'S',0.5,parameters.charge());
    ExactSpinor es2s(2,'S',0.5,parameters.charge());
    ExactSpinor es2p1(2,'P',0.5,parameters.charge());
    ExactSpinor es2p2(2,'P',1.5,parameters.charge());

    print("exact energies: 1s, 2s, 2p1/2, 2p3/2",es1s.get_energy()-c*c,
          es2s.get_energy()-c*c,
          es2p1.get_energy()-c*c,
          es2p2.get_energy()-c*c);




//    eigenvector_test(world,Ansatz1(nuclear_charge,nemo_factor),ExactSpinor(1,'S',0.5,nuclear_charge, 0.5));
//    eigenvector_test(world,Ansatz1(nuclear_charge,nemo_factor),ExactSpinor(2,'S',0.5,nuclear_charge, 0.5));
//    eigenvector_test(world,Ansatz1(nuclear_charge,nemo_factor),ExactSpinor(2,'P',0.5,nuclear_charge, 0.5));
//    eigenvector_test(world,Ansatz1(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge, 1.5));
//    eigenvector_test(world,Ansatz0(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge, 0.5));
//    eigenvector_test(world,Ansatz0(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,-0.5));
//    eigenvector_test(world,Ansatz0(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,-1.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,nemo_factor),ExactSpinor(2,'S',0.5,nuclear_charge,-0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge, 1.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge, 0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,-0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,-1.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,nemo_factor),ExactSpinor(2,'P',0.5,nuclear_charge, 0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,nemo_factor),ExactSpinor(2,'P',0.5,nuclear_charge, 0.5));
//    eigenvector_test(world,Ansatz1(nuclear_charge,1),ExactSpinor(1,'S',0.5,nuclear_charge));
//    eigenvector_test(world,Ansatz2(nuclear_charge,1),ExactSpinor(1,'S',0.5,nuclear_charge));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,nemo_factor),ExactSpinor(1,'S',0.5,nuclear_charge));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,1.3),ExactSpinor(3,'D',2.5,nuclear_charge));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,nemo_factor),ExactSpinor(1,'S',0.5,nuclear_charge,0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,nemo_factor),ExactSpinor(1,'S',0.5,nuclear_charge,-0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,1.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,-0.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,nemo_factor),ExactSpinor(2,'P',1.5,nuclear_charge,-1.5));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,-1.2),ExactSpinor(2,'P',1.5,nuclear_charge));
//    eigenvector_test(world,Ansatz3(nuclear_charge,1,-1.2),ExactSpinor(1,'S',0.5,nuclear_charge));


    try {
        run(world,Ansatz0(parameters.charge()),parameters.charge(),parser,parameters.nstates());
    } catch (...) {
        std::cout << "caught an error " << std::endl;
    }
    finalize();
    return 0;


}
