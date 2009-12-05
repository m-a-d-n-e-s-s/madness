#ifndef MADNESS_MRA_SIMPLECACHE_H__INCLUDED
#define MADNESS_MRA_SIMPLECACHE_H__INCLUDED

#include <mra/mra.h>

namespace madness {
    /// Simplified interface around hash_map to cache stuff for 1D 

    /// This is a write once cache --- subsequent writes of elements
    /// have no effect (so that pointers/references to cached data
    /// cannot be invalidated)
    template <typename Q, int NDIM>
    class SimpleCache {
    private:
        typedef ConcurrentHashMap< Key<NDIM>, Q, KeyHash<NDIM> > mapT;
        typedef std::pair<Key<NDIM>, Q> pairT;
        mapT cache;

    public:
        SimpleCache() : cache() {};

        SimpleCache(const SimpleCache& c) : cache(c.cache) {};

        SimpleCache& operator=(const SimpleCache& c) {
            if (this != &c) {
                cache.clear();
                cache = c.cache;
            }
            return *this;
        }

        /// If key is present return pointer to cached value, otherwise return NULL
        inline const Q* getptr(const Key<NDIM>& key) const {
            typename mapT::const_iterator test = cache.find(key);
            if (test == cache.end()) return 0;
            return &(test->second);
        }


        /// If key=(n,l) is present return pointer to cached value, otherwise return NULL

        /// This for the convenience (backward compatibility) of 1D routines
        inline const Q* getptr(Level n, Translation l) const {
            Key<NDIM> key(n,Vector<Translation,NDIM>(l));
            return getptr(key);
        }


        /// If key=(n,l) is present return pointer to cached value, otherwise return NULL

        /// This for the convenience (backward compatibility) of 1D routines
        inline const Q* getptr(Level n, const Key<NDIM>& disp) const {
            Key<NDIM> key(n,disp.translation());
            return getptr(key);
        }


        /// Set value associated with key ... gives ownership of a new copy to the container
        inline void set(const Key<NDIM>& key, const Q& val) {
            cache.insert(pairT(key,val));
        }

        inline void set(Level n, Translation l, const Q& val) {
            Key<NDIM> key(n,Vector<Translation,NDIM>(l));
            set(key, val);
        }

        inline void set(Level n, const Key<NDIM>& disp, const Q& val) {
            Key<NDIM> key(n,disp.translation());
            set(key, val);
        }
    };
}
#endif // MADNESS_MRA_SIMPLECACHE_H__INCLUDED
