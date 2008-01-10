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

/// \file loadbal.h
/// \brief Declares and partially implements MyPmap, LoadBalImpl and associated load balancing classes.

  
#ifndef LOADBAL_H
#define LOADBAL_H

#include <tensor/random.h>

namespace madness {

    typedef int Cost;
    typedef double CompCost;

    /// Finds exponent k such that d^k <= me < d^{k+1}
    inline int nearest_power(int me, int d) {
        int k = 0;
        while (me != 0) {
            if (me%d == 0) {
                k++;
                me/=d;
            } else {
                break;
            }
        }
        return k;
    };

    template <int D> class LBNode;
    template <int D> struct TreeCoords;
    template <int D> class MyPmap;
    template <int D> class LBTree;
    class NodeData;
    template <int D> class LoadBalImpl;

    /// Convenient typedef shortcuts

    /// Makes it easier to handle these unwieldy templated types
    template <int D>
    struct DClass {
        typedef Key<D> KeyD;
        typedef const Key<D> KeyDConst;
        typedef TreeCoords<D> TreeCoords;
        typedef LBNode<D> NodeD;
        typedef const LBNode<D> NodeDConst;
        typedef MyPmap<D> MyPmap;
        typedef LBTree<D> treeT;
      typedef std::vector< std::vector< madness::TreeCoords<D> > > vvTreeCoords;
    };

    template <int D>
    class PartitionInfo {
    public:
	Cost maxcost;
	Cost cost_left;
	Cost skel_cost;
	unsigned int partition_number;
	unsigned int step_num;
	double facter;
	PartitionInfo(double f=1.1) :
	    maxcost(0)
	    , cost_left(0)
	    , skel_cost(0)
	    , partition_number(0)
	    , step_num(0)
	    , facter(f) { };

	void reset(unsigned int p=1) {
	  maxcost = 0;
	  cost_left = skel_cost;
	  partition_number = p;
	  step_num++;
	}

	template <typename Archive>
        void serialize(const Archive& ar) {
            ar & maxcost & cost_left & skel_cost & partition_number & step_num & facter;
        }
    };

    template <int D>
    std::ostream& operator<<(std::ostream& s, const PartitionInfo<D>& pi) {
        s << "maxcost = " << pi.maxcost << ", cost_left = " << pi.cost_left << 
	  ", skel_cost = " << pi.skel_cost << ", partition_number = " << 
	  pi.partition_number << ", step_num = " << pi.step_num << 
	  ", facter = " << pi.facter;
        return s;
    };



    /// Diagnostic data contained in fascimile tree
    /// Diagnostic data, including the cost of the node and the subtree headed by that node,
    /// along with a bool flag used during depth-first partitioning
    class NodeData {
        friend std::ostream& operator<<(std::ostream& s, const NodeData& nd);
    public:
        int cost;
        int subcost;
        bool is_taken;
        NodeData(int c = 1, int s = 1, bool i = false) : cost(c), subcost(s), is_taken(i) {};
        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & cost & subcost & is_taken;
        }
        void print() {
            cout << "cost = " << cost << ", subcost = " << subcost << ", is_taken = " << is_taken << endl;
        }
	template <typename functionT>
	void set_data(functionT function) {
	  functionT tmp = (*function);
	  this->cost = (int) tmp;
	}

	template <typename functionT, typename arg1T>
	void set_data(functionT function, const arg1T& arg1) {
	  functionT tmp = (*function)(arg1);
	  this->cost = (int) tmp;
	}
    };


    inline std::ostream& operator<<(std::ostream& s, const NodeData& nd) {
        s << "cost " << nd.cost << ", subcost " << nd.subcost << ", is_taken " << nd.is_taken;
        return s;
    };


    /// The node that is used in the fascimile copy of the tree to be load balanced

    /// The node used in the tree that is operated upon and load balanced in LoadBalImpl.
    template <int D>
    class LBNode {
    private:
        NodeData data;
        std::vector<bool> c; /// Existence of each child individually

        void all_children(bool status=false) {
            c.clear();
            c.assign(dim, status);
        };

    public:
	mutable KeyChildIterator<D> rpit;
        static int dim; /// Number of children in standard tree (e.g. 2^D)
	int nrecvd;

        LBNode() {
	    rpit = KeyChildIterator<D>();
            data = NodeData();
            all_children();
	    nrecvd = 0;
        };

        LBNode(const NodeData& d, bool children=false, int n=0) : data(d), nrecvd(n) {
	    rpit = KeyChildIterator<D>();
            all_children(children);
        };

	LBNode(const LBNode& node) : data(node.data), c(node.c), rpit(node.rpit), nrecvd(node.nrecvd) { };


	/// Determines whether node has any children at all
        bool has_children() const {
            for (int i = 0; i < dim; i++)
                if (c[i]) return true;
            return false;
        };

        bool has_child(unsigned int i) const {
            return c[i];
        };

        bool has_child(int i) const {
            return c[i];
        };

	int get_num_children() const {
	    int nkids = 0;
	    for (int i=0; i < dim; i++) {
		if (has_child(i)) nkids++;
	    }
	    return nkids;
	}

        void set_child(int i, bool setto = true) {
            c[i] = setto;
        };

	void set_all_children(bool setto = true) {
	  all_children(setto);
	}

        void set_data(const NodeData& d) {
            data = d;
        };

	template <typename functionT>
	void set_cost(functionT function) {
	  data.set_cost<functionT>(function);
	}

	template <typename functionT, typename arg1T>
	void set_cost(functionT function, const arg1T& arg1) {
	  data.set_cost<functionT, arg1T>(function, arg1);
	}

        NodeData get_data() const {
            return data;
        };

        vector<bool> get_c() const {
            return c;
        };

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & data & c & rpit;
        }
    };


    template <int D>
    std::ostream& operator<<(std::ostream& s, const LBNode<D>& node) {
        s << "data = " << node.get_data() << ", c = " << node.get_c();
	if (node.rpit) {
	  s  << ", key_iterator = " << node.rpit.key();
	} else {
	  s << ", key_iterator = <EMPTY>";
	}
        return s;
    };

    template <int D>
    std::ostream& operator<<(std::ostream& s, typename DClass<D>::NodeDConst& node) {
        s << "data = " << node.get_data() << ", c = " << node.get_c();
	if (node.rpit) {
	  s << ", key_iterator = " << node.rpit.key();
	} else {
	  s << ", key_iterator = <EMPTY>";
	}
        return s;
    };


    template <int D>
	int LBNode<D>::dim = power<D>();



    /// Key + owner, struct used to determine mapping of tree nodes
    template <int D>
    class TreeCoords {
    public:
        Key<D> key;
        ProcessID owner;

        TreeCoords(const Key<D>& k, ProcessID o) : key(Key<D>(k)), owner(o) {};
        TreeCoords(const TreeCoords& t) : key(Key<D>(t.key)), owner(t.owner) {};
        TreeCoords() : key(Key<D>()), owner(-1) {};
        void print() const {
            madness::print(key, "   owner =", owner);
        };

        bool operator< (const TreeCoords t) const {
            return (this->key < t.key);
        };

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & key & owner;
        }
    };

    template <int D>
    inline std::ostream& operator<<(std::ostream& s, const TreeCoords<D>& tc) {
        s << tc.key << "   owner = " << tc.owner;
        return s;
    };



    template<int D>
    class ProcMapImpl {
	public:
/*
#ifdef WORLDDC_USES_GNU_HASH_MAP
        template <typename T>
        struct PMLocalHash {
            std::size_t operator()(const T& t) const {
                return hash(t);
            };
        };
        typedef HASH_MAP_NAMESPACE::hash_map< typename DClass<D>::KeyDConst,ProcessID,PMLocalHash<typename DClass<D>::KeyDConst > > Mapinfo;
#else
*/
        typedef std::map<typename DClass<D>::KeyDConst,ProcessID> Mapinfo;
/*
#endif
*/
	typedef typename Mapinfo::iterator iterator;
	typedef const iterator iterator_const;
	typedef std::pair< typename DClass<D>::KeyDConst, ProcessID > pairT;

	ProcMapImpl() {};
	ProcMapImpl(std::vector< TreeCoords<D> > v) {
	    int vlen = v.size();
	    for (int i = 0; i < vlen; i++) {
	    	themap.insert(std::make_pair(v[i].key, v[i].owner));
	    }
	};

	ProcMapImpl(const TreeCoords<D>& t) {
	    themap.insert(std::make_pair(t.key, t.owner));
	};
	void insert(const TreeCoords<D>& t) {
	    themap.insert(std::make_pair(t.key, t.owner));
	};
	void erase(const TreeCoords<D>& t) {
	    themap.erase(t.key);
	};

        ProcessID find_owner(const Key<D>& key) const {
	    typename std::map<typename DClass<D>::KeyDConst,ProcessID>::const_iterator it = themap.find(key);
	    if (it != themap.end()) {
		return it->second;
	    } else if (key.level() == 0) {
		madness::print("find_owner: owner of ", key, "not found but returning 0");
		return 0;
	    } else {
		return this->find_owner(key.parent());
	    }
	};

	void print() {
	    for (iterator it = themap.begin(); it != themap.end(); ++it) {
		madness::print(it->first, "   ", it->second);
	    }
	}

	private:
	Mapinfo themap;
    };


    /// Procmap implemented using Tree of TreeCoords

    template <int D>
    class MyPmap : public WorldDCPmapInterface< Key<D> > {
    private:
        Tensor<ProcessID> simple_key_map; // map of keys at level n
        bool simplemap;
	const int nproc;
        const ProcessID me;
	const int n;
	SharedPtr< ProcMapImpl<D> > tree_map;
        typedef Key<D> KeyD;

	/// private method that builds the Tree underlying the procmap
	void build_tree_map(std::vector< TreeCoords<D> > v) {
	    tree_map = SharedPtr< ProcMapImpl<D> > (new ProcMapImpl<D>(v));
	};

	ProcessID simple_hash(const KeyD& key) const {
            KeyD parent = (key.level() > n) ? key.parent(key.level()-n) : key;
            return simple_key_map((const long *) &(parent.translation()[0]));
	};

    public:
        MyPmap() : simplemap(false) {};

        static bool costmapcmp(const std::pair<KeyD,double>& a, const std::pair<KeyD,double>& b) {
            return a.second > b.second;
        };

	MyPmap(World& world) 
            : simplemap(true)
            , nproc(world.nproc())
            , me(world.rank())
            , n((std::log(world.size())/std::log(2.0)+3)/D + 2)  // 16*nproc = 2^(nD)
        {
            // We set n to have about 16 tasks per processor and we try to
            // give each process a mix of large, medium, and small
            // tasks.  Currently estimate cost as inversely
            // proportional to distance from center but we could
            // enable the user to provide a function.

            //if (world.rank() == 0) madness::print("DIM",D,"N IN MAP IS",n);

            std::vector<long> vdim(D);
            for (int i=0; i<D; i++) vdim[i] = 1L<<n;
            simple_key_map = Tensor<ProcessID>(vdim);
            
            std::list< std::pair<KeyD,double> > costmap;
            Vector<Translation,D> l;
            long cent = (1L<<n) / 2;
            for (TensorIterator<ProcessID> iter=simple_key_map.unary_iterator(0,false,false); iter._p0; ++iter) { 
                double dist = 0.01;
                for (int i=0; i<D; i++) {
                    l[i] = iter.ind[i];
                    dist += (l[i] - cent)*(l[i] - cent);
                }
                double cost = 1.0/dist; // actually dist squared
                cost *= (1.0 + 0.001*RandomNumber<double>()); // To shuffle (nearly) equal values
                costmap.push_back(std::pair<KeyD,double>(KeyD(n,l),cost));
            }
            costmap.sort(costmapcmp);
//             if (world.rank() == 0) {
//                 for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
//                     madness::print("costmap", it->first, it->second);
//                 }
//             }
            ProcessID p = 0;
            for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
                const long *l = (const long *) &(it->first.translation()[0]);
                simple_key_map(l)  = p;
                p++;
                if (p == world.size()) p = 0;
            }
//             if (world.rank() == 0) {
//                 madness::print("SIMPLE MAP", D,"\n", simple_key_map);
//             }
        };

        MyPmap(World& world, vector<TreeCoords<D> > v) : simplemap(false), nproc(world.nproc()), me(world.rank()), n(0) {
            build_tree_map(v);
        };

        MyPmap(const MyPmap<D>& other) : simplemap(other.staticmap), nproc(other.nproc), me(other.me), n(other.n), tree_map(other.tree_map) {};

        MyPmap<D>& operator=(const MyPmap<D>& other) {
            if (this != &other) {
                simple_key_map = other.simple_key_map; // shallow copy
                simplemap = other.simplemap;
                nproc = other.nproc;
		me = other.me;
		n = other.n;
                tree_map = other.tree_map;
            }
            return *this;
        };

        void print() const {
	    if (!simplemap) {
		tree_map->print();
	    } else {
		madness::print("MyPmap: simple map with n =", n);
	    }
        };

	/// Find the owner of a given key
        ProcessID owner(const KeyD& key) const {
            if (simplemap)
                return simple_hash(key);
            else {
		return tree_map->find_owner(key);
            }
        };
    };


    /// The container in which the fascimile tree with its keys mapping to LBNodes is stored
    template <int D>
    class LBTree : public WorldObject< LBTree<D> > {
    public:
	typedef WorldObject<LBTree<D> > woT;
        typedef WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD> dcT;
	typedef typename dcT::iterator iterator;

	World& world;

	static typename DClass<D>::KeyDConst root;
	static typename DClass<D>::vvTreeCoords list_of_list;
	static std::vector<Cost> cost_list;

	PartitionInfo<D> partition_info;
	std::vector<typename DClass<D>::TreeCoords> temp_list;
	
    private:
        dcT impl;

    public:
        LBTree(World& world, const SharedPtr< WorldDCPmapInterface<typename DClass<D>::KeyD> >& pmap) : woT(world)
	    , world(world)
	    , impl(world,pmap) {
	    impl.process_pending();
	    this->process_pending();
        };
	/// Initialize the LBTree by converting a FunctionImpl to a LBTree
        template <typename T>
        inline void init_tree(const SharedPtr< FunctionImpl<T,D> >& f) {
            for (typename FunctionImpl<T,D>::dcT::iterator it = f->coeffs.begin(); it != f->coeffs.end(); ++it) {
            	// convert Node to LBNode
            	NodeData nd;
		typename DClass<D>::KeyD key = it->first;
            	if (!(it->second.has_children())) {
                	typename DClass<D>::NodeD lbnode(nd,false);
                	// insert into impl
                	impl.insert(key, lbnode);
            	} else {
                	typename DClass<D>::NodeD lbnode(nd,true);
                	// insert into impl
                	impl.insert(key, lbnode);
                }
            }
        }


        // Methods:

        template <typename T>
        inline void add_tree(const SharedPtr< FunctionImpl<T,D> >& f) {
            for (typename FunctionImpl<T,D>::dcT::iterator it = f->coeffs.begin(); it != f->coeffs.end(); ++it) {
            	// convert Node to LBNode
            	NodeData nd;
		typename DClass<D>::KeyD key = it->first;
		typename DClass<D>::treeT::iterator tree_it = impl.find(key);
		if (tree_it != impl.end()) {
		  typename DClass<D>::NodeD lbnode = tree_it->second;
		  if (it->second.has_children()) {
		    lbnode.set_all_children(true);
		  }
		  NodeData nd=lbnode.get_data();
		  nd.cost++;
		  nd.subcost++;
		  lbnode.set_data(nd);
		  impl.insert(key, lbnode);
		} else {
		    if (!(it->second.has_children())) {
		      typename DClass<D>::NodeD lbnode(nd,false);
		      // insert into impl
		      impl.insert(key, lbnode);
		    } else {
		      typename DClass<D>::NodeD lbnode(nd,true);
		      // insert into impl
		      impl.insert(key, lbnode);
		    }
		}
            }
        }

	ProcessID owner(typename DClass<D>::KeyDConst& key) {
	    return impl.owner(key);
	}

        void print(typename DClass<D>::KeyDConst& key) {
            typename DClass<D>::treeT::iterator it = impl.find(key);
            if (it == impl.end()) return;
            for (Level i = 0; i < key.level(); i++) cout << "  ";
            madness::print(key, it->second);
            for (KeyChildIterator<D> kit(key); kit; ++kit) {
                print(kit.key());
            }
        };

	void find_partitions(PartitionInfo<D>& pi);
	bool verify_partition(std::vector<TreeCoords<D> >& part_list);
	Void launch_make_partition(PartitionInfo<D> pi, bool first_time);
	Void meld_all(bool first_time);

        Cost fix_cost();


	void init_fix_cost();
	void fix_cost_spawn();
	Void fix_cost_sum(typename DClass<D>::KeyDConst& key, Cost c);


        void rollup();

	void reset(bool taken);

        void meld(typename DClass<D>::treeT::iterator it);


        Void make_partition(typename DClass<D>::KeyDConst& key, Cost partition_size,
			    Cost used_up, PartitionInfo<D> pi, bool downward);
	Void totally_reset(PartitionInfo<D> pi);
	Void add_to_partition(typename DClass<D>::TreeCoords p);

	typename DClass<D>::KeyD first_child(typename DClass<D>::KeyDConst& key, const typename DClass<D>::NodeD& node);
	typename DClass<D>::KeyD next_sibling(typename DClass<D>::KeyDConst& key);

	bool reset_partition(Cost& partition_size, Cost& used_up, PartitionInfo<D>& pi);


        MyPmap<D>& get_mypmap() {
            return *static_cast< MyPmap<D>* >(impl.get_pmap().get());
        };
       
        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & root & list_of_list & cost_list & impl;
        }


    };

    template <int D>
	typename DClass<D>::KeyDConst LBTree<D>::root(0);

    template <int D>
      typename DClass<D>::vvTreeCoords LBTree<D>::list_of_list;

    template <int D>
      typename std::vector<Cost> LBTree<D>::cost_list;

    /// Implementation of load balancing

    /// Implements the load balancing algorithm upon the tree underlying a function.
    template <int D>
    class LoadBalImpl {
    private:
	typedef MyPmap<D> Pmap;
	int k;
	double comm_bandw;
	double comm_latency;
	double flop_time;
        SharedPtr<typename DClass<D>::treeT> skeltree;
	World& world;

	template<typename T>
        void construct_skel(SharedPtr<FunctionImpl<T,D> > f) {
            skeltree = SharedPtr<typename DClass<D>::treeT>(new typename DClass<D>::treeT(f->world,
                       f->coeffs.get_pmap()));
//            madness::print("about to initialize tree");
	    skeltree->template init_tree<T>(f);
//            madness::print("just initialized tree");
	    pi.skel_cost = skeltree->fix_cost();
	    //	    madness::print("construct_skel: pi.skel_cost =", pi.skel_cost); 
	    pi.cost_left = pi.skel_cost;
        };

    public:
	PartitionInfo<D> pi;
        //Constructors
	template <typename T>
	LoadBalImpl(Function<T,D> f, double a=1e-8, double b=1e-5, double c=5e-10, double facter=1.1) : k(f.get_impl()->k) 
	    , comm_bandw(a)
	    , comm_latency(b)
	    , flop_time(c)
	    , world(f.get_impl()->world)
	    , pi(PartitionInfo<D>(facter)) {
            construct_skel(f.get_impl());
	    pi.partition_number = f.get_impl()->world.mpi.nproc()-1;
        };

        ~LoadBalImpl() {};

        //Methods

	/// Returns a shared pointer to a new process map, which can then be used to redistribute the function
        SharedPtr< WorldDCPmapInterface< Key<D> > > load_balance() {
            return SharedPtr< WorldDCPmapInterface< Key<D> > >(new MyPmap<D>(world, find_best_partition()));
        };

        std::vector<typename DClass<D>::TreeCoords> find_best_partition();

	CompCost compute_comp_cost(Cost c, int n);

	template <typename T>
	void add_tree(Function<T,D> f) {
	  skeltree->add_tree(f.get_impl());
	    pi.skel_cost = skeltree->fix_cost();
	    pi.cost_left = pi.skel_cost;

	}

    };


    Cost compute_partition_size(Cost cost, unsigned int parts);



}

#endif
