#ifndef MAD_KEY_H
#define MAD_KEY_H

/// \file key.h
/// \brief Multidimension Key for MRA tree and associated iterators 

namespace madness {

    template <int NDIM> class KeyChildIterator;

    /// Key is the index for a node of the 2^NDIM-tree

    /// See KeyChildIterator for facile generation of children,
    /// and foreach_child(parent,op) for facile applicaiton of operators
    /// to child keys.
    template <int NDIM>
    class Key {
        friend class KeyChildIterator<NDIM>;
    private:
        Level n;
        Array<Translation,NDIM> l;
        hashT hashval;

        // Recomputes hashval
        void rehash() {
            hashval = madness::hash(n,madness::hash(l));
        };
    public:
        Key() {};

        Key(Level n, const Array<Translation,NDIM> l) 
            : n(n), l(l) {
            rehash();
        };

        bool operator==(const Key& other) const {
            if (hashval != other.hashval) return false;
            if (n != other.n) return false;
            return l == other.l;
        };

        bool operator<(const Key& other) const {
            if (hashval < other.hashval) return true;
            else if (hashval > other.hashval) return false;
            else if (n < other.n) return true;
            else if (n > other.n) return false;
            else return (l < other.l);
        };

        inline hashT hash() const {
            return hashval;
        };

        template <typename Archive>
        inline void serialize(Archive& ar) {
            ar & archive::wrap((unsigned char*) this, sizeof(this));
        }

        Level level() const {
            return n;
        };

        const Array<Translation,NDIM>& translation() const {
            return l;
        };
    };

    template <int NDIM>
    std::ostream& operator<<(std::ostream& s, const Key<NDIM>& key) {
        s << "(" << key.level() << "," << key.translation() << ")";
        return s;
    };

    /// Iterates in lexical order thru all children of a key

    /// Example usage:
    /// \code
    ///    for (KeyChildIterator<NDIM> it(key); it; ++it) print(it.key());
    /// \endcode
    template <int NDIM>
    class KeyChildIterator {
        Key<NDIM> parent;
        Key<NDIM> child;
        Array<Translation,NDIM> p;
        bool finished;

    public:
        KeyChildIterator(const Key<NDIM>& parent) 
            : parent(parent)
            , child(parent.n+1,parent.l*2)
            , p(0) 
            , finished(false)
        {};

        /// Pre-increment of an iterator (i.e., ++it)
        KeyChildIterator& operator++() {
            if (finished) return *this;
            int i;
            for (i=0; i<NDIM; i++) {
                if (p[i] == 0) {
                    p[i]++;
                    child.l[i]++;
                    for (int j=0; j<i; j++) {
                        p[j]--;
                        child.l[j]--;
                    }
                    break;
                }
            }
            finished = (i == NDIM);
            return *this;
        };

        /// True if iterator is not at end
        operator bool() const {
            return !finished;
        };

        /// Returns the key of the child
        inline const Key<NDIM>& key() const {
            return child;
        };
    };


    /// Applies op(key) to each child key of parent
    template <int NDIM, typename opT>
    inline void foreach_child(const Key<NDIM>& parent, opT& op) {
        for (KeyChildIterator<NDIM> it(parent); it; ++it) op(it.key());
    };

}

#endif

