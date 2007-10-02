/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
  $Id$
*/

  
#ifndef MAD_KEY_H
#define MAD_KEY_H

/// \file key.h
/// \brief Multidimension Key for MRA tree and associated iterators 

#include <mra/power.h>

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
        Vector<Translation,NDIM> l;
        hashT hashval;

        // Recomputes hashval
        void rehash() {
            hashval = madness::hash(n,madness::hash(l));
        };

	// Helper function for operator <
	int encode(int dig) const {
	    int retval = 0;
	    for (int j=0; j < NDIM; j++) {
		// retval += ((l[j]/2^{n-1-dig}) mod 2) * 2^j
		retval += ((l[j]>>(n-1-dig))%2)<<j;
	    }
	    return retval;
	};

	// Helper function for (Level, Translation) constructor
	Vector<Translation,NDIM> decode(Level level, Translation k) const {
	    Vector<Translation,NDIM> L(0);
	    int twotoD = power<NDIM>();
	    int powr=1, divisor=2;
	    for (Level i = 0; i < level; i++) {
		Translation r = k%twotoD;
		for (int j=0; j<NDIM; j++) {
		    L[NDIM-j-1]+=(r%divisor)*powr;
		    r/=divisor;
		}
		k/=twotoD;
		powr*=2;
	    }
	    return L;
	};
    public:
        /// Default constructor makes an \em uninitialized key
        Key() {};

        /// Constructor with given n, l
        Key(Level n, const Vector<Translation,NDIM> l) 
            : n(n), l(l) {
            rehash();
        };

        /// Constructor with given n and l=0
	Key(int n) : n(n), l(0) {
            rehash();
	};

        /// Constructor from lexical index in depth first order
	Key(Level n, Translation p) : n(n) {
	    l = decode(n,p);
	    rehash();
	};

        /// Equality test
        bool operator==(const Key& other) const {
            if (hashval != other.hashval) return false;
            if (n != other.n) return false;
            bool result = l == other.l;
	    if (result && hashval!=other.hashval) {
		print("!!  keys same but hash is different",hashval,other.hashval,*this,other);
                MADNESS_EXCEPTION("Tell HQI not RJ3!",0);
            }
            return result;
        };

	/// Comparison based upon depth first lexical order
	bool operator<(const Key& other) const {
	    if (*this == other) return false; // I am not less than self
	    Level nmin;
	    bool retval = false;

	    if (this->n > other.n) {
		nmin = other.n;
		retval = true;
	    }
	    else {
		nmin = this->n;
	    }
		
	    for (Level i = 0; i < nmin; i++) {
		int tthis = this->encode(i), tother = other.encode(i);
		if (tthis != tother) {
		    return (tthis < tother);
		}   
	    }
	    return retval;
	}

        inline hashT hash() const {
            return hashval;
        };

        template <typename Archive>
        inline void serialize(Archive& ar) {
            ar & archive::wrap((unsigned char*) this, sizeof(*this));
        }

        Level level() const {
            return n;
        };

        const Vector<Translation,NDIM>& translation() const {
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
            Vector<Translation,NDIM> pl;
            if (generation > n) generation = n;
            for (int i=0; i<NDIM; i++) pl[i] = l[i]>>generation;
            return Key(n-generation,pl);
        };

	bool is_child_of(const Key key) const {
            if (this->n < key.n) {
            	return false; // I can't be child of something lower on the tree
            }
	    else if (this->n == key.n) {
                return (*this == key); // I am child of myself
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
        Vector<Translation,NDIM> p;
        bool finished;

    public:
	KeyChildIterator()
	    : p(0)
	    , finished(true)
	{};

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
            child.rehash();
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


    /// Applies member function of obj to each child key of parent
    template <int NDIM, typename objT>
    inline void foreach_child(const Key<NDIM>& parent, objT* obj, void (objT::*memfun)(const Key<NDIM>&)) {
        for (KeyChildIterator<NDIM> it(parent); it; ++it) (obj->*memfun)(it.key());
    };

}

#endif

