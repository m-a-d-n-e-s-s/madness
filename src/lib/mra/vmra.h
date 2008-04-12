#ifndef MADNESS_VMRA_H
#define MADNESS_VMRA_H

/// \file vmra.h

/// \brief Defines operations on vectors of Functions

#include <mra/mra.h>
#include <cstdio>

namespace madness {

    /// Compress a vector of functions
    template <typename T, int NDIM> 
    void compress(World& world,
                  const std::vector< Function<T,NDIM> >& v, 
                  bool fence=true) {
        bool must_fence = false;
        for (unsigned int i=0; i<v.size(); i++) {
            if (!v[i].is_compressed()) {
                v[i].compress(false);
                must_fence = true;
            }
        }
        
        if (fence && must_fence) world.gop.fence();
    }
    
    
    /// Reconstruct a vector of functions
    template <typename T, int NDIM> 
    void reconstruct(World& world,
                     const std::vector< Function<T,NDIM> >& v, 
                     bool fence=true) {
        bool must_fence = false;
        for (unsigned int i=0; i<v.size(); i++) {
            if (v[i].is_compressed()) {
                v[i].reconstruct(false);
                must_fence = true;
            }
        }
        
        if (fence && must_fence) world.gop.fence();
    }
    

    /// Generates non-standard form of a vector of functions
    template <typename T, int NDIM> 
    void nonstandard(World& world,
                     std::vector< Function<T,NDIM> >& v, 
                     bool fence=true) {
        reconstruct(world, v);
        for (unsigned int i=0; i<v.size(); i++) {
            v[i].nonstandard(false);
        }
        if (fence) world.gop.fence();
    }
    

    /// Generates standard form of a vector of functions
    template <typename T, int NDIM> 
    void standard(World& world,
                     std::vector< Function<T,NDIM> >& v, 
                     bool fence=true) {
        for (unsigned int i=0; i<v.size(); i++) {
            v[i].standard(false);
        }
        if (fence) world.gop.fence();
    }
    

    /// Truncates a vector of functions
    template <typename T, int NDIM> 
    void truncate(World& world,
                  std::vector< Function<T,NDIM> >& v, 
                  double tol=0.0,
                  bool fence=true) {

        compress(world, v);

        for (unsigned int i=0; i<v.size(); i++) {
            v[i].truncate(tol, false);
        }
        
        if (fence) world.gop.fence();
    }

    
    /// Differentiates a vector of functions
    template <typename T, int NDIM> 
    std::vector< Function<T,NDIM> > diff(World& world,
                                         const std::vector< Function<T,NDIM> >& v, 
                                         int axis,
                                         bool fence=true) {
        reconstruct(world, v);
        
        std::vector< Function<T,NDIM> > df(v.size());
        for (unsigned int i=0; i<v.size(); i++) {
            df[i] = diff(v[i],axis,false);
        }
        if (fence) world.gop.fence();
        return df;
    }
    
    /// Generates a vector of zero functions
    template <typename T, int NDIM>
    std::vector< Function<T,NDIM> >
    zero_functions(World& world, int n) {
        std::vector< Function<T,NDIM> > r(n);
        for (int i=0; i<n; i++) 
            r[i] = Function<T,NDIM>(FunctionFactory<T,NDIM>(world));

        return r;
    }


    /// Transforms a vector of functions according to new[i] = sum[j] old[j]*c[j,i]
    template <typename T, typename R, int NDIM> 
    std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > 
    transform(World& world,
              const std::vector< Function<T,NDIM> >& v, 
              const Tensor<R>& c,
              bool fence=true) {

        typedef TENSOR_RESULT_TYPE(T,R) resultT;
        int n = v.size();  // n is the old dimension
        int m = c.dim[1];  // m is the new dimension
        MADNESS_ASSERT(n==c.dim[0]);
        std::vector< Function<resultT,NDIM> > vc = zero_functions<resultT,NDIM>(world, m);
        
        compress(world, v);
        compress(world, vc); 
        
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                vc[i].gaxpy(1.0,v[j],c(j,i),false);
            }
        }
        if (fence) world.gop.fence();
        return vc;
    }
    
    
    /// Scales inplace a vector of functions by distinct values
    template <typename T, int NDIM> 
    void scale(World& world,
               std::vector< Function<T,NDIM> >& v, 
               const std::vector<double>& factors,
               bool fence=true) {
        for (unsigned int i=0; i<v.size(); i++) v[i].scale(factors[i],false);
        if (fence) world.gop.fence();
    }
    

    /// Computes the 2-norm of a vector of functions
    template <typename T, int NDIM> 
    std::vector<double> norm2(World& world, 
                              const std::vector< Function<T,NDIM> >& v) {
        std::vector<double> norms(v.size());
        for (unsigned int i=0; i<v.size(); i++) norms[i] = v[i].norm2sq_local();
        world.gop.sum(&norms[0], norms.size());
        for (unsigned int i=0; i<v.size(); i++) norms[i] = sqrt(norms[i]);
        return norms;
    }
    

    /// Computes the matrix inner product of two function vectors - q(i,j) = inner(f[i],g[j])
    template <typename T, typename R, int NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > matrix_inner(World& world,
                                                   const std::vector< Function<T,NDIM> >& f, 
                                                   const std::vector< Function<R,NDIM> >& g,
                                                   bool sym=false) {
        long n=f.size(), m=g.size();
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n,m);
        if (sym) MADNESS_ASSERT(n==m);
        
        compress(world, f);
        compress(world, g);

        for (long i=0; i<n; i++) {
            long jtop = m;
            if (sym) jtop = i+1;
            for (long j=0; j<jtop; j++) {
                r(i,j) = f[i].inner_local(g[j]);
                if (sym) r(j,i) = r(i,j);
            }
        }
        
        world.gop.sum(r.ptr(),n*m);
        return r;
    }
    
    /// Computes the inner product of two function vectors - q(i) = inner(f[i],g[i])
    template <typename T, typename R, int NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > inner(World& world,
                                            const std::vector< Function<T,NDIM> >& f, 
                                            const std::vector< Function<R,NDIM> >& g) {
        long n=f.size(), m=g.size();
        MADNESS_ASSERT(n==m);
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n);
        
        compress(world, f);
        compress(world, g);
        
        for (long i=0; i<n; i++) {
            r(i) = f[i].inner_local(g[i]);
        }
        
        world.gop.sum(r.ptr(),n);
        return r;
    }
    

    /// Multiplies a function against a vector of functions --- q[i] = a * v[i]
    template <typename T, typename R, int NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul(World& world,
        const Function<T,NDIM>& a, 
        const std::vector< Function<R,NDIM> >& v, 
        bool fence=true) 
    {
        reconstruct(world, v);
        a.reconstruct();
        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > av(v.size());
        for (unsigned int i=0; i<v.size(); i++) {
            av[i] = mul(a, v[i], false);
        }
        if (fence) world.gop.fence();
        return av;
    }

    /// Multiplies a function against a vector of functions using sparsity of v[i] --- q[i] = a * v[i]
    template <typename T, typename R, int NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul_sparse(World& world,
               const Function<T,NDIM>& a, 
               const std::vector< Function<R,NDIM> >& v, 
               double tol,
               bool fence=true) 
    {
        reconstruct(world, v);
        a.reconstruct();
        a.norm_tree();
        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > av(v.size());
        for (unsigned int i=0; i<v.size(); i++) {
            av[i] = mul_sparse<T,R,NDIM>(v[i], a, tol, false);
        }
        if (fence) world.gop.fence();
        return av;
    }

    /// Sets the threshold in a vector of functions
    template <typename T, int NDIM>
    void set_thresh(World& world, std::vector< Function<T,NDIM> >& v, double thresh, bool fence=true) {
        for (unsigned int j=0; j<v.size(); j++) {
            v[j].set_thresh(thresh,false);
        }
        if (fence) world.gop.fence();
    }

    /// Computes the square of a vector of functions --- q[i] = v[i]**2
    template <typename T, int NDIM>
    std::vector< Function<T,NDIM> >
    square(World& world,
        const std::vector< Function<T,NDIM> >& v, 
        bool fence=true) 
    {
        reconstruct(world, v);
        std::vector< Function<T,NDIM> > vsq(v.size());
        for (unsigned int i=0; i<v.size(); i++) {
            vsq[i] = square(v[i], false);
        }
        if (fence) world.gop.fence();
        return vsq;
    }

    /// Returns new vector of functions --- q[i] = a[i] + b[i]
    template <typename T, typename R, int NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    add(World& world,
        const std::vector< Function<T,NDIM> >& a, 
        const std::vector< Function<R,NDIM> >& b, 
        bool fence=true) 
    {
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(a.size());
        for (unsigned int i=0; i<a.size(); i++) {
            r[i] = add(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }


    /// Returns new vector of functions --- q[i] = a[i] - b[i]
    template <typename T, typename R, int NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    sub(World& world,
        const std::vector< Function<T,NDIM> >& a, 
        const std::vector< Function<R,NDIM> >& b, 
        bool fence=true) 
    {
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(a.size());
        for (unsigned int i=0; i<a.size(); i++) {
            r[i] = sub(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }

    
    /// Generalized A*X+Y for vectors of functions ---- q[i] = alpha*a[i] + beta*b[i]
    template <typename T, typename Q, typename R, int NDIM>
    void gaxpy(World& world,
               Q alpha,
               std::vector< Function<T,NDIM> >& a, 
               Q beta,
               const std::vector< Function<R,NDIM> >& b, 
               bool fence=true) 
    {
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        for (unsigned int i=0; i<a.size(); i++) {
            a[i].gaxpy(alpha, b[i], beta, false);
        }
        if (fence) world.gop.fence();
    }


    /// Applies a vector of operators to a vector of functions --- q[i] = apply(op[i],f[i])
    template <typename opT, typename R, int NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> >
    apply(World& world,
          const std::vector< SharedPtr<opT> >& op, 
          const std::vector< Function<R,NDIM> > f) {

        MADNESS_ASSERT(f.size()==op.size());

        std::vector< Function<R,NDIM> >& ncf = *const_cast< std::vector< Function<R,NDIM> >* >(&f);

        reconstruct(world, f);
        nonstandard(world, ncf);

        std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> > result(f.size());
        for (unsigned int i=0; i<f.size(); i++) {
            result[i] = apply_only(*op[i], f[i], false);
        }

        world.gop.fence();

        standard(world, ncf, false);  // restores promise of logical constness
	reconstruct(world, result);

        return result;
    }

    
    /// Applies an operator to a vector of functions --- q[i] = apply(op,f[i])
    template <typename opT, typename R, int NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> >
    apply(World& world,
          opT& op, 
          const std::vector< Function<R,NDIM> > f) {

        std::vector< Function<R,NDIM> >& ncf = *const_cast< std::vector< Function<R,NDIM> >* >(&f);

        reconstruct(world, f);
        nonstandard(world, ncf);

        std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> > result(f.size());
        for (unsigned int i=0; i<f.size(); i++) {
            result[i] = apply_only(op, f[i], false);
        }

        world.gop.fence();

        standard(world, ncf, false);  // restores promise of logical constness
	reconstruct(world, result);

        return result;
    }

    /// Normalizes a vector of functions --- v[i] = v[i].scale(1.0/v[i].norm2())
    template <typename T, int NDIM>
    void normalize(World& world, vector< Function<T,NDIM> >& v, bool fence=true) {
        vector<double> nn = norm2(world, v);
        for (unsigned int i=0; i<v.size(); i++) v[i].scale(1.0/nn[i],false);
        if (fence) world.gop.fence();
    }

}
#endif
