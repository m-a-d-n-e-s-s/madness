//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness.h>
#include <chem/nemo.h>
#include <chem/projector.h>
#include <chem/SCFOperators.h>

using namespace madness;

namespace madness {


class GaussianGuess : public FunctionFunctorInterface<double,3> {
public:
    GaussianGuess(const Atom& atom, const double e, const int i, const int j,
            const int k) : x(atom.x), y(atom.y), z(atom.z),
            exponent(e), i(i), j(j), k(k) {
    }

    double x,y,z;       ///< location of the guess
    double exponent;    ///< exponent of the gues
    int i,j,k;          ///< cartesian exponents

    double operator()(const coord_3d& xyz) const {
        double xx=x-xyz[0];
        double yy=y-xyz[1];
        double zz=z-xyz[2];
        const double e=exponent*(xx*xx + yy*yy + zz*zz);
        return pow(xx,i)*pow(yy,j)*pow(zz,k)*exp(-e);
    }
};


/// make a set of Cartesian Gaussian functions located on atom

/// @param[in]  atom    where the Gaussian guess is located
/// @param[in]  l       the l-quantum number
/// @param[in]  exponent    exponent of the Gaussian function
std::vector<GaussianGuess> make_guess(const Atom& atom, int l, const double exponent) {
    std::vector<GaussianGuess> gg;
    const int maxlp1=l+1;

    // loop over all cartesian components of l
    for (int i=0; i<1000; ++i) {
        std::vector<int> ijk(3);
        ijk[0]=i%maxlp1;
        ijk[1]=(i/maxlp1)%maxlp1;
        ijk[2]=(i/maxlp1/maxlp1)%maxlp1;
        int current_l=ijk[0]+ijk[1]+ijk[2];
        if (current_l==l) {
            gg.push_back(GaussianGuess(atom,exponent,ijk[0],ijk[1],ijk[2]));
        }
        if (ijk[0]+ijk[1]+ijk[2]==3*l) break;
    }
    return gg;
}








class PNO {
public:
    typedef std::shared_ptr<operatorT> poperatorT;

    mutable double sss,ttt;
    void START_TIMER(World& world) const {
        world.gop.fence(); ttt=wall_time(); sss=cpu_time();
    }

    void END_TIMER(World& world, const char* msg) const {
        ttt=wall_time()-ttt; sss=cpu_time()-sss;
        if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
    }

    /// POD for PNO keywords
    struct Parameters {
        /// max l quantum number in the guess
        int lmax;

        /// type of partial wave guess function, e.g. 3s2p1d for a specific atom
        std::map<std::string,std::vector<int> > pw_guess;

        /// number of frozen orbitals; note the difference to the "pair" keyword where you
        /// request a specific orbital. Here you freeze lowest orbitals, i.e. if you find
        ///  freeze 1
        /// in the input file the 0th orbital is kept frozen, and orbital 1 is the first to
        /// be correlated.
        int freeze;

        int i, j;

        int maxiter;    ///< max number of iterations

        /// ctor reading out the input file
        Parameters(const std::string& input) : lmax(8),
                freeze(0), i(-1), j(-1), maxiter(20) {

            read_parameters(input);

        }

        /// read pno parameters from the input file
        void read_parameters(const std::string& input) {
            // get the parameters from the input file
            std::ifstream f(input.c_str());
            position_stream(f, "pno");
            std::string s;
            std::string str_pw_guess,symbol;

            while (f >> s) {
                if (s == "end") break;
                else if (s == "maxiter") f >> maxiter;
                else if (s == "freeze") f >> freeze;
                else if (s == "pair") f >> i >> j;
                else if (s == "pw_guess") {
                    if (f >> symbol >> str_pw_guess) {
                        pw_guess[symbol]=convert_pw_guess(str_pw_guess);
                    } else {
                        print("line",symbol,str_pw_guess);
                        MADNESS_EXCEPTION("failed to read partial wave guess",1);
                    }
                }
                else continue;
            }
        }

        /// convert "11s8p2d" into a vector of ints
        std::vector<int> convert_pw_guess(const std::string& input) const {
            std::vector<int> pw_guess(lmax);

            int ones=0;
            for (std::size_t i=0; i<input.size(); ++i) {
                if ((input[i]>='0') and (input[i]<='9')) {
                    if (ones==0) ones=atoi(&input[i]);  // converts all digits
                } else {
                    const char lquant=input[i];
                    if (lquant=='s') {pw_guess[0]=ones;}
                    if (lquant=='p') {pw_guess[1]=ones;}
                    if (lquant=='d') {pw_guess[2]=ones;}
                    if (lquant=='f') {pw_guess[3]=ones;}
                    if (lquant=='g') {pw_guess[4]=ones;}
                    if (lquant=='h') {pw_guess[5]=ones;}
                    if (lquant=='i') {pw_guess[6]=ones;}
                    if (lquant=='k') {pw_guess[7]=ones;}
                    ones=0;
                }
            }
            std::cout << "final partial wave guess" << std::endl;
            for (std::size_t i=0; i<pw_guess.size(); ++i) {
                if (pw_guess[i]>0) print("l ",i,pw_guess[i]);
            }
            return pw_guess;
        }


    };

    PNO(World& world, const Nemo& nemo, const std::string input) : sss(), ttt(),
        world(world), param(input), nemo(nemo),
        J(world,&nemo),
        K(world,&nemo,0),
        T(world),
        V(world,nemo.nuclear_correlation),
        F(world,nemo.get_calc().get(),nemo.nuclear_correlation),
        Q(world,nemo.get_calc()->amo) {

        poisson = std::shared_ptr<real_convolution_3d>(
              CoulombOperatorPtr(world, nemo.get_calc()->param.lo,
                        nemo.get_calc()->param.econv));
        bsh = std::shared_ptr<real_convolution_3d>(
                BSHOperatorPtr3D(world, 1.e-8, nemo.get_calc()->param.lo,
                        nemo.get_calc()->param.econv));
    }

    void solve() const {
        const vecfuncT& amo=nemo.get_calc()->amo;

        // plot mos
        for (std::size_t i=0; i<amo.size(); ++i) {
            std::string name="amo"+stringify(i);
            plot_plane(world,amo[i],name);
        }

        double energy=0.0;
        if ((param.i >= 0) and (param.j>=0)) {
            energy+=solve_pair(param.i,param.j);
        } else {
            for (std::size_t i=param.freeze; i<amo.size(); ++i) {
                for (std::size_t j=i; j<amo.size(); ++j) {
                    energy+=solve_pair(i,j);
                }
            }
        }

    }

    double solve_pair(const int i, const int j) const {
        const vecfuncT& amo=nemo.get_calc()->amo;
        print("solving pair ",i,j);

        vecfuncT virtuals=guess_virtuals();

        orthonormalize_cholesky(virtuals);
        vecfuncT virtuals_bar=copy(world,virtuals);

        // compute the energies of the occupied orbitals i and j
        const double f_ii=F(amo[i],amo[i]);
        const double f_jj=F(amo[j],amo[j]);
        print("diagonal fock matrix elements ",f_ii,f_jj);
        Tensor<double> amplitudes;

        for (int iter=0; iter<param.maxiter; iter++) {

            // recompute intermediates
            vecfuncT V_baraj_i=compute_V_aj_i(amo[i],amo[j],virtuals_bar);
            vecfuncT V_ai_j=compute_V_aj_i(amo[j],amo[i],virtuals);
            vecfuncT V_barai_j=compute_V_aj_i(amo[j],amo[i],virtuals_bar);
            vecfuncT V_aj_i=compute_V_aj_i(amo[i],amo[j],virtuals);
            Tensor<double> fmat=F(virtuals,virtuals);
            Tensor<double> fmat_bar=F(virtuals_bar,virtuals_bar);

            // compute energy
            Tensor<double> VBV=compute_hylleraas(i,j,V_aj_i,V_baraj_i,virtuals,virtuals_bar,
                    fmat,fmat_bar,f_ii+f_jj,amplitudes);
            double energy=VBV.sum();

            // sort the amplitudes and intermediates with increasing energy
            // and recompute the Fock matrices if necessary
            if (sort_virtuals(VBV,amplitudes,virtuals,virtuals_bar,
                    V_ai_j,V_baraj_i,V_barai_j,V_aj_i)) {
                fmat=F(virtuals,virtuals);
                fmat_bar=F(virtuals_bar,virtuals_bar);
            }


            if (world.rank() == 0) printf("in iteration %2d at time %6.1f: %12.8f\n",
                    iter, wall_time(),energy);
            for (std::size_t a=0; a<virtuals.size(); ++a) {
                std::string name="virtual"+stringify(a)+"_iteration"+stringify(iter);
                plot_plane(world,virtuals[a],name);
            }

            // will change virtuals and invalidate V_aj_i
            vecfuncT newvirtuals=update_virtuals(i,j,V_baraj_i,V_barai_j,virtuals,virtuals_bar,
                    f_ii+f_jj,fmat_bar,amplitudes);
            virtuals_bar=update_virtuals(i,j,V_ai_j,V_aj_i,virtuals_bar,virtuals,
                    f_ii+f_jj,fmat,amplitudes);
            virtuals=newvirtuals;

        }

        // compute the fock matrix of the virtuals
        vecfuncT V_baraj_i=compute_V_aj_i(amo[i],amo[j],virtuals_bar);
        vecfuncT V_aj_i=compute_V_aj_i(amo[i],amo[j],virtuals);
        Tensor<double> fmat=F(virtuals,virtuals);
        Tensor<double> fmat_bar=F(virtuals_bar,virtuals_bar);
        Tensor<double> VBV=compute_hylleraas(i,j,V_aj_i,V_baraj_i,virtuals,virtuals_bar,
                fmat,fmat_bar,f_ii+f_jj,amplitudes);
        double energy=VBV.sum();
        return energy;
    }

    /// update the virtual functions

    /// @param[in]  amplitudes  the amplitudes t_a
    vecfuncT update_virtuals(const int i, const int j,
            const vecfuncT& V_baraj_i, const vecfuncT& V_barai_j,
            const vecfuncT& virtuals, const vecfuncT& virtuals_bar,
            const double e0, const Tensor<double>& fmat,
            const Tensor<double> amplitudes) const {

        const int nvir=virtuals.size();

        START_TIMER(world);
        // compute the inhomogeneous term
        vecfuncT V_baraj_i_scaled=copy(world,V_baraj_i);
        vecfuncT V_barai_j_scaled=copy(world,V_barai_j);

        for (std::size_t a=0; a<V_baraj_i.size(); ++a) {
            V_baraj_i_scaled[a].scale(1./amplitudes(a));
            V_barai_j_scaled[a].scale(0.5/amplitudes(a));
        }
//        if (i!=j) {
            V_baraj_i_scaled=sub(world,V_baraj_i_scaled,V_barai_j_scaled);
//        }

        // compute the coupling term \sum_b |b> ...
        Tensor<double> transf(nvir,nvir);
        for (int a=0; a<fmat.dim(0); ++a) {
            for (int b=0; b<fmat.dim(1); ++b) {
                transf(a,b)=-fmat(a,b)*amplitudes(a)/amplitudes(b);
            }
            transf(a,a)=0.0;
        }

        // combine coupling and inhomogeneous term
        vecfuncT btilde=transform(world,virtuals,transf);
        btilde=sub(world,btilde,V_baraj_i_scaled);

        // compute (J-K+V) | a >
        vecfuncT Ja=J(virtuals);
        vecfuncT Ka=K(virtuals);
        vecfuncT Va=V(virtuals);
        vecfuncT JKVa=add(world,sub(world,Ja,Ka),Va);
        truncate(world,JKVa);
        END_TIMER(world,"prep update");
        vecfuncT Vpsi=sub(world,btilde,JKVa);
//        if (i!=j) {
            vecfuncT coupling=compute_coupling_terms(virtuals,virtuals_bar,
                    amplitudes);
            Vpsi=add(world,Vpsi,coupling);
//        }

        START_TIMER(world);
        scale(world,Vpsi,2.0);
        vecfuncT GVpsi;

        Tensor<double> evals(fmat.dim(0));
        for (long a=0; a<evals.size(); ++a) evals(a)=e0-fmat(a,a);
        std::vector<poperatorT> bsh3= nemo.get_calc()->make_bsh_operators(world,evals);
        GVpsi=apply(world,bsh3,Vpsi);

//        if (i!=j) {
            Tensor<double> S_bbara=matrix_inner(world,virtuals,virtuals_bar);
            for (int a=0; a<S_bbara.dim(1); ++a) {
                for (int b=0; b<S_bbara.dim(0); ++b) {
                    S_bbara(b,a)=S_bbara(b,a)*amplitudes(b)/amplitudes(a);
                }
            }
            vecfuncT tmp=transform(world,virtuals_bar,S_bbara);
            scale(world,tmp,0.5);
            GVpsi=add(world,GVpsi,tmp);
//        }
        END_TIMER(world,"BSH");


        START_TIMER(world);
        // post-processing: project out occ space, orthogonalize
        GVpsi=Q(GVpsi);
        orthonormalize_cholesky(GVpsi);
        check_orthonormality(GVpsi);
        END_TIMER(world,"post processing");

        return GVpsi;
    }

    vecfuncT compute_coupling_terms(const vecfuncT& virtuals,
            const vecfuncT& virtuals_bar, const Tensor<double>& amplitudes) const {
        START_TIMER(world);

        Tensor<double> S_bbara=matrix_inner(world,virtuals,virtuals_bar);
        Tensor<double> F_bbara=F(virtuals,virtuals_bar);

        for (int b=0; b<S_bbara.dim(0); ++b) {
            for (int a=0; a<S_bbara.dim(1); ++a) {
                S_bbara(b,a)=S_bbara(b,a)*amplitudes(b)/amplitudes(a);
                F_bbara(b,a)=F_bbara(b,a)*amplitudes(b)/amplitudes(a);
            }
        }

        vecfuncT Fb=transform(world,virtuals_bar,F_bbara);

        Tensor<double> F_barabarb=F(virtuals_bar,virtuals_bar);
        std::vector<double> F_aa(virtuals_bar.size());
        for (std::size_t i=0; i<F_aa.size(); ++i) F_aa[i]=F_barabarb(i,i);

        // terms \sum_b (J-K+V) |\bar b> S_{\bar b a}
        vecfuncT Jb=J(virtuals_bar);
        vecfuncT Kb=K(virtuals_bar);
        vecfuncT Vb=V(virtuals_bar);
        vecfuncT JKVa=add(world,sub(world,Jb,Kb),Vb);
        vecfuncT barc=transform(world,JKVa,S_bbara);

        // term F_{\bar a \bar a} \sum_b |\bar b>S_{\bar b a}
        vecfuncT Faab=transform(world,virtuals_bar,S_bbara);
        scale(world,Faab,F_aa);

        vecfuncT result=sub(world,add(world,barc,Fb),Faab);
        truncate(world,result);
        scale(world,result,0.5);
        END_TIMER(world,"compute coupling");

        return result;
    }

    /// guess a set up virtual orbitals -- currently from the minimal AO basis
    vecfuncT guess_virtuals() const {
        vecfuncT virtuals;
        for (int i=0; i<nemo.molecule().natom(); ++i) {

            // get the partial wave basis for the atom
            const Atom& atom=nemo.molecule().get_atom(i);
            const std::string symbol=atomic_number_to_symbol(atom.atomic_number);
            std::vector<int> pw_guess;
            if (param.pw_guess.find(symbol)!=param.pw_guess.end()) {
                pw_guess=param.pw_guess.find(symbol)->second;
            } else {
                if (world.rank()==0) print("symbol",symbol);
                MADNESS_EXCEPTION("could not find basis for atom",1);
            }

            print("working on atom",symbol);
            print("pw_guess",pw_guess);
            for (int l=0; l<param.lmax; ++l) {
                for (int i=0; i<pw_guess[l]; ++i) {
                    double e=double(pw_guess[l])/(double(i+1.0));
                    append(virtuals,guess_virtual_gaussian_shell(atom,l,e));
                    print("l,e",l,e);
                }
            }
        }
        print("number of guess virtuals: ",virtuals.size());
        virtuals=Q(virtuals);
        // plot guess
        for (std::size_t i=0; i<virtuals.size(); ++i) {
            std::string name="virtuals"+stringify(i);
            plot_plane(world,virtuals[i],name);
        }

        return virtuals;
    }

    /// return a shell of l-quantum l and exponent e, including all cartesian components
    vecfuncT guess_virtual_gaussian_shell(const Atom& atom, const int l,
            const double e) const {
        vecfuncT virtuals;
        std::vector<GaussianGuess> gg=make_guess(atom,l,e);
        for (std::size_t m=0; m<gg.size(); ++m) {
            virtuals.push_back(real_factory_3d(world).functor(gg[m]).truncate_on_project());
        }
        normalize(world,virtuals);
        return virtuals;
    }

    void append(vecfuncT& v, const vecfuncT& rhs) const {
        for (std::size_t i=0; i<rhs.size(); ++i) v.push_back(rhs[i]);
    }

    /// compute the function V_{\bar a j} | i> = - \int \dr' (-J12+K12+1/|r-r'|) a(r') j(r')

    /// the terms are expanded as follows:
    /// (-J1 +K1) | i(1) >  < a(2) | j(2) >
    ///  +  | i(1) > < a(2) | -J(2) + K(2) | j(2) >
    ///  +  i(1) * \int \dr2 1/|r12| a(2) j(2)
    /// the first line is zero due to orthogonality, the second line drops out
    /// due to the orthogonality projector.
    vecfuncT compute_V_aj_i(const real_function_3d& phi_i,
            const real_function_3d& phi_j, const vecfuncT& virtuals) const {

        const vecfuncT aj=mul(world,phi_j,virtuals);    // multiply \bar a j
        vecfuncT gaj=apply(world,*poisson,aj);        // \int \dr2 aj(2)/r12
        vecfuncT Vaj_i=mul(world,phi_i,gaj);
        vecfuncT Vaj_i1=Q(Vaj_i);
        truncate(world,Vaj_i1);
        return Vaj_i1;
    }


    /// compute the Hylleraas functional and the corresponding energy

    /// the virtuals are sorted according to their energy, so that some
    /// intermediates need to be resorted as well, in particular the V_aj_i
    /// intermediate and the Fock matrix. Its contents remain unchanged.
    // @param[in]  amo occupied orbitals
    /// @param[in,out]  virtuals    the optimal virtual orbitals; sorted upon exit
    /// @param[in]  e0  the zeroth-order energy (e_i + e_j)
    /// @param[out] t the optimized amplitudes
    /// @return the energy contributions of each virtual
    Tensor<double> compute_hylleraas(const int i, const int j,
            const vecfuncT& V_aj_i, const vecfuncT& V_baraj_i,
            const vecfuncT& virtuals, const vecfuncT& virtuals_bar,
            const Tensor<double>& f_aa, const Tensor<double>& f_barabara,
            const double e0, Tensor<double>& t) const {

        const std::size_t nvir=virtuals.size();
        Tensor<double> B(nvir,nvir),V;

//        if (i==j) {
//
//            V=inner(world,virtuals,V_ai_j);
//            for (std::size_t a=0; a<nvir; ++a) B(a,a)=f_aa(a,a) *2.0 - e0;
//
//
//        } else {
            const Tensor<double> f_abarb=F(virtuals,virtuals_bar);  // F_{a \bar b}
            const Tensor<double> f_barab=F(virtuals_bar,virtuals);  // F_{\bar a b}
            const Tensor<double> S_abarb=matrix_inner(world,virtuals,virtuals_bar);
            const Tensor<double> S_barab=matrix_inner(world,virtuals_bar,virtuals);
            Tensor<double> V_abara=inner(world,virtuals,V_baraj_i);
            Tensor<double> V_baraa=inner(world,virtuals_bar,V_aj_i);
            V=2.0*V_abara-V_baraa;

            for (std::size_t a=0; a<nvir; ++a) {
                B(a,a)=2.0*(f_aa(a,a)+f_barabara(a,a) - e0);
                for (std::size_t b=0; b<nvir; ++b) {
                    B(a,b)-=(f_abarb(a,b)*S_barab(a,b)+f_barab(a,b)*S_abarb(a,b)
                            -e0*S_abarb(a,b)*S_barab(a,b));
                }
            }
//        }
        double fac=2.0;
        if (i==j) fac-=1.0;
        B.scale(fac);
        V.scale(fac);
        Tensor<double> Binv=inverse(B);

        t=-inner(Binv,V);
        V.emul(t);
        return V;
    }

    struct sort_helper {
        sort_helper(double e, double a,
                const real_function_3d& v1, const real_function_3d& v2,
                const real_function_3d& v3, const real_function_3d& v4,
                const real_function_3d& v5, const real_function_3d& v6)
            : energy(e), amplitude(a), v1(v1), v2(v2), v3(v3), v4(v4),
              v5(v5), v6(v6) {}
        double energy;
        double amplitude;
        real_function_3d v1, v2, v3, v4, v5, v6;
    };

    static bool comp(const sort_helper& rhs, const sort_helper& lhs) {
        return rhs.energy<lhs.energy;
    }

    /// sort the virtuals according to their pair energies

    /// @return if the virtuals have been resorted
    bool sort_virtuals(Tensor<double>& VT, Tensor<double> amplitudes,
            vecfuncT& v1, vecfuncT& v2, vecfuncT& v3, vecfuncT& v4,
            vecfuncT& v5, vecfuncT& v6) const {
        std::vector<sort_helper> pairs;
        for (long i=0; i<VT.size(); ++i) {
            pairs.push_back(sort_helper(VT(i),amplitudes(i),v1[i],v2[i],v3[i],
                    v4[i],v5[i],v6[i]));
        }
        if (std::is_sorted(pairs.begin(),pairs.end(),PNO::comp)) {
            return false;
        } else {
            std::sort(pairs.begin(),pairs.end(),PNO::comp);
            for (long i=0; i<VT.size(); ++i) {
                VT(i)=pairs[i].energy;
                amplitudes(i)=pairs[i].amplitude;
                v1[i]=pairs[i].v1;
                v2[i]=pairs[i].v2;
                v3[i]=pairs[i].v3;
                v4[i]=pairs[i].v4;
                v5[i]=pairs[i].v5;
                v6[i]=pairs[i].v6;
            }
        }
        return true;
    }

    void check_orthonormality(const vecfuncT& v) const {
        Tensor<double> ovlp=matrix_inner(world,v,v);
        for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)-=1.0;
        double error=ovlp.normf()/ovlp.size();
        if (error>1.e-14) print("orthonormality error: ",error);
    }

    void orthonormalize_cholesky(vecfuncT& v) const {
        Tensor<double> ovlp=matrix_inner(world,v,v);
        cholesky(ovlp); // destroys ovlp
        Tensor<double> L=transpose(ovlp);
        Tensor<double> Linv=inverse(L);
        Tensor<double> U=transpose(Linv);
        v=transform(world,v,U);
        truncate(world,v);
    }

private:

    World& world;
    Parameters param;   ///< calculation parameters
    Nemo nemo;
    Coulomb J;
    Exchange K;
    Kinetic<double,3> T;
    Nuclear V;
    Fock F;
    QProjector<double,3> Q;
    std::shared_ptr<real_convolution_3d> poisson;
    std::shared_ptr<real_convolution_3d> bsh;

};
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
    startup(world,argc,argv);
    std::cout.precision(6);

    const std::string input="input";
    std::shared_ptr<SCF> calc(new SCF(world,input.c_str()));
    Nemo nemo(world,calc);
    const double energy=nemo.value();
    if (world.rank()==0) print("nemo energy: ",energy);
    if (world.rank() == 0) printf(" at time %.1f\n", wall_time());
    
    const vecfuncT nemos=nemo.get_calc()->amo;
    const vecfuncT R2nemos=mul(world,nemo.nuclear_correlation->square(),nemos);

    Fock F(world,nemo.get_calc().get(),nemo.nuclear_correlation);
    Tensor<double> fmat=F(nemos,nemos);
    print("Fock matrix");
    print(fmat);

    PNO pno(world,nemo,input);
    pno.solve();

    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}

