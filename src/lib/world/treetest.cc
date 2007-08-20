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

  
  $Id: world.cc 289 2007-08-13 15:51:09Z HartmanBaker $
*/

  
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>

using namespace madness;
using namespace std;


struct Key {
    typedef unsigned long ulong;
    ulong n, i, j, k;
    hashT hashval;

    Key() {};  // Empty default constructor for speed - but is therefore junk

    Key(ulong n, ulong i, ulong j, ulong k)
        : n(n), i(i), j(j), k(k), hashval(madness::hash(&this->n,4,0)) {};

    hashT hash() const {
        return hashval;
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & n & i & j & k & hashval;
    }

    bool operator==(const Key& b) const {
        // It's extremely probable that different keys will have a different hash
        return hashval==b.hashval && n==b.n && i==b.i && j==b.j && k==b.k;
    }

    Key parent() const {
	return Key(n-1, i>>1, j>>1, k>>1);
    }
};

ostream& operator<<(ostream& s, const Key& key) {
    s << "Key(" << key.n << "," << key.i << "," << key.j << "," << key.k << ")";
    return s;
}

class KeyChildIterator {
    typedef unsigned long ulong;
    Key parent;
    Key child;
    ulong i,j,k;
public:
    KeyChildIterator(const Key& parent) 
	: parent(parent) 
	, child(parent.n+1,2*parent.i,2*parent.j,2*parent.k)
	, i(0)
	, j(0)
	, k(0) 
	{};

    KeyChildIterator& operator++() {
	if (k == 2) return *this;
	i++;
	if (i == 2) {
	    i = 0;
	    j++;
	    if (j == 2) {
		j = 0;
		k++;
	    }
	}
	// k==2 indicates end
	child = Key(parent.n+1, 2*parent.i+i, 2*parent.j+j, 2*parent.k+k);
	return *this;
    };

    const Key& key() {
	return child;
    };

    operator bool() const {
	return k!=2;
    };
};

struct Node;
ostream& operator<<(ostream& s, const Node& node);

struct Node {
    typedef WorldContainer<Key,Node> dcT;
    double value;
    bool isleaf;
    int nrecvd;

    Node() : value(0.0), isleaf(true), nrecvd(0) {};
    Node(double value) : value(value), isleaf(true), nrecvd(0) {};
    Node(const Node& node) : value(node.value), isleaf(node.isleaf) {};

    Void random_insert(const dcT& constd, const Key& key, double valin) {
        dcT& d = const_cast<dcT&>(constd);
        value = valin;
        isleaf = true;
        if (value > 0.25 && d.size()<9) {
            isleaf = false;
            World& world = d.world();
            double ran = 0.5; //world.mpi.drand();
	    for (KeyChildIterator it(key); it; ++it) {
		d.task(it.key(),&Node::random_insert, d, it.key(), value*ran);
	    }		
        }
        return None;
    };

    template <class Archive>
    void serialize(Archive& ar) {
        ar & value & isleaf;
    }

    bool is_leaf() const {return isleaf;};

    double get() const {return value;};
  
    Void set(double v) {value = v; return None;};

    void zero_nrecvd() {nrecvd = 0;};

    // Void recursive_print(dcT& d) const {  // d should be const
    Void recursive_print(const dcT& constd, const Key& key) const {
        dcT& d = const_cast<dcT&>(constd); // !!!!!!!!!!!!!!
	cout << d.world().rank() << ": RP: " << key << " " << *this;
	if (! is_leaf()) {
	    for (KeyChildIterator it(key); it; ++it) {
		d.send(it.key(),&Node::recursive_print, d, it.key());
	    }	
	}	
	return None;
    };

    Void receive_sum(const dcT& constd, const Key& key, double gift) {
        dcT& d = const_cast<dcT&>(constd); // !!!!!!!!!!!!!!
	value += gift;
	nrecvd++;
	if (nrecvd == 8 && key.n!=0) {
	    Key parent = key.parent();
	    d.send(parent, &Node::receive_sum, d, parent, get());
	}
	return None;
    };

    double do_sum(vector< Future<double> > v) {
	for (int i=0; i<8; i++) 
	    value += v[i].get();
	return value;
    };

    Future<double> do_sum_spawn(const dcT& constd, const Key& key) {
        dcT& d = const_cast<dcT&>(constd); // !!!!!!!!!!!!!!
	if (is_leaf()) {
	    return Future<double>(value);
	}
	else {
	    vector< Future<double> > v = future_vector_factory<double>(8);
	    int i=0;
	    for (KeyChildIterator it(key); it; ++it, ++i) {
		const Key& child = it.key();
		v[i] = d.send(child, &Node::do_sum_spawn, d, child);
	    }
	    return d.task(key, &Node::do_sum, v);
	}
    };
};

ostream& operator<<(ostream& s, const Node& node) {
  s << "Node(" << node.get() << "," << node.is_leaf() << ")" << endl;
  return s;
}

void doit(World& world) {
    ProcessID me = world.rank();
    WorldContainer<Key,Node> d(world);
    Key root(0,0,0,0);
    
    // First build an oct-tree with random depth
    world.mpi.srand();
    print("first ran#",world.mpi.drand());
    world.gop.fence();
    if (me == 0) d.send(root,&Node::random_insert,d,root,1.0);
    world.gop.fence();
    print("FINISHED INSERTING");

    // Now using AM walk the tree ... first in "non-process-centric" manner
    if (me == 0) d.send(root, &Node::recursive_print, d, root);
    world.gop.fence();
    print("FINISHED PRINTING");

//     // Now add up the tree
//     for (WorldContainer<Key,Node>::iterator it=d.begin(); it!=d.end(); ++it) {
// 	it->second.zero_nrecvd();
//     }
//     world.gop.fence();
//     print("FINISHED ZEROING COUNTERS");

//     for (WorldContainer<Key,Node>::iterator it=d.begin(); it!=d.end(); ++it) {
// 	const Key& key = it->first;
// 	Node& node = it->second;
// 	if (node.is_leaf()) {
// 	    Key parent = key.parent();
// 	    d.send(parent, &Node::receive_sum, d, parent, node.get());
// 	}
//     }
//     world.gop.fence();
//     print("FINISHED SUMMING");

    if (me == 0) d.send(root, &Node::recursive_print, d, root);
    world.gop.fence();

    if (me == 0) {
	print("AND THE FINAL ANSWER IS", d.send(root, &Node::do_sum_spawn, d, root).get());
    }
    world.gop.fence();
}


int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    redirectio(world);

    try {
	doit(world);
    } catch (MPI::Exception e) {
        error("caught an MPI exception");
    } catch (madness::MadnessException e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    print("done with final fence");
    MPI::Finalize();


    return 0;
}
