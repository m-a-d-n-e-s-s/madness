#ifndef MAD_KEY_H
#define MAD_KEY_H

/// \file key.h
/// \brief Multidimension Key for MRA tree and associated iterators 

namespace madness {

    typedef unsigned long Translation;
    typedef long Level;

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

	Key(int k) : n(0) {
	    for (int i = 0; i < NDIM; i++)
		l[i] = 0;
            rehash();
	    
	};
        bool operator==(const Key& other) const {
            if (hashval != other.hashval) return false;
            if (n != other.n) return false;
            return l == other.l;
        };

/*
        bool operator<(const Key& other) const {
            if (hashval < other.hashval) return true;
            else if (hashval > other.hashval) return false;
            else if (n < other.n) return true;
            else if (n > other.n) return false;
            else return (l < other.l);
        };
*/
	// Sorry, but I need this really really ugly "less than" operator. -hqi

	bool operator<(const Key& other) const {
	    if (*this == other) return false; // I am not less than self
	    int ans;
	    if (this->n == other.n) {
		ans = ordering(*this, other);
	    }
	    else if (this->n > other.n) {
		Level dn = this->n - other.n;
		Key newthis = this->parent(dn);
		if (newthis == other) {
		    ans = 1;
		}
		else {
		    ans = ordering(newthis, other);
		}
	    }
	    else {
		Level dn = other.n - this->n;
		Key newother = other.parent(dn);
		if (newother == *this) {
		    ans = 0;
		}
		else {
		    ans = ordering(*this, newother);
		    if (ans < 0)
			ans = 1;
		    else
			ans = 0;
		}
	    }
	    if (ans > 0)
		return true;
	    else
		return false;
	};

	int ordering(const Key& k1, const Key& k2) const {
            bool egalite = true;
            Array<int,NDIM> dl;
//          print("ordering: comparing", k1, "and", k2);
            for (unsigned int i = 0; i < NDIM; i++) {
                dl[i] = k1.l[i] - k2.l[i];
                if (k1.l[i]/2 != k2.l[i]/2) {
                    egalite = false;
//                  cout << "ordering: k1 and k2 do not have same parent" << endl;
                }
            }
            if (!egalite) {
                return (ordering(k1.parent(), k2.parent()));
            }
            else {
                for (unsigned int i = 0; i < NDIM; i++) {
                    if (dl[i] > 0) {
//                      cout << "ordering: dl[" << i << "] > 0; returning -1" << endl;
                        return -1;
                    }
                    else if (dl[i] < 0) {
//                      cout << "ordering: dl[" << i << "] < 0; returning 1" << endl;
                        return 1;
                    }
                }
//              cout << "ordering: no dL greater or less than zero; returning 0" << endl;
                return 0;
            }
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

        /// Returns the key of the parent

        /// Default is the immediate parent (generation=1).  To get
        /// the grandparent use generation=2, and similarly for
        /// great-grandparents.
        ///
        /// !! If there is no such parent it quietly returns the
        /// closest match (which may be self if this is the top of the
        /// tree).
        Key parent(int generation=1) const {
            Array<Translation,NDIM> pl;
            if (generation > n) generation = n;
            for (int i=0; i<NDIM; i++) pl[i] = l[i]>>generation;
            return Key(n-generation,pl);
        };

	bool is_child_of(const Key key) const {
	    if (*this == key) {
		return true; // I am child of myself
	    }
	    else if (this->n <= key.n) {
            	return false; // I can't be child of something lower on the tree
            }
            else {
            	Level dn = this->n - key.n;
            	Key mama = this->parent(dn);
            	return (mama == key);
            }
	};

	bool is_parent_of(const Key key) const {
	    return (key.is_child_of(*this));
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

