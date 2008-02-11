#ifndef MADNESS_VMRA_H
#define MADNESS_VMRA_H

namespace madness {

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
    
    
    template <typename T, int NDIM> 
    void scale(World& world,
               std::vector< Function<T,NDIM> >& v, 
               const std::vector<double>& factors,
               bool fence=true) {
        for (unsigned int i=0; i<v.size(); i++) v[i].scale(factors[i],false);
        if (fence) world.gop.fence();
    }
    
    template <typename T, int NDIM> 
    std::vector<double> norm2(World& world, 
                              const std::vector< Function<T,NDIM> >& v) {
        std::vector<double> norms(v.size());
        for (unsigned int i=0; i<v.size(); i++) norms[i] = v[i].norm2sq_local();
        world.gop.sum(&norms[0], norms.size());
        for (unsigned int i=0; i<v.size(); i++) norms[i] = sqrt(norms[i]);
        return norms;
    }
    
    template <typename T, typename R, int NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > inner(World& world,
                                            const std::vector< Function<T,NDIM> >& f, 
                                            const std::vector< Function<R,NDIM> >& g) {
        long n=f.size(), m=g.size();
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n,m);
        
        compress(world, g);
        if ((void *) &f != (void *) &g) compress(world, g);
        
        for (long i=0; i<n; i++) {
            for (long j=0; j<m; j++) {
                r(i,j) = f[i].inner_local(g[j]);
            }
        }
        
        world.gop.sum(r.ptr(),n*m);
        return r;
    }
    
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
}
#endif
