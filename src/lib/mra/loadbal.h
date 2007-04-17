#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#ifndef LOADBAL_H
#define LOADBAL_H

#include <world/world.h>
#include <mra/key.h>
#include <mra/mra.h>
//using namespace madness;
using namespace std;

namespace madness {

typedef int Cost;
typedef double CompCost;

inline int nearest_power(int me, int d) {
    int k = 0;
    while (me != 0) {
	if (me%d == 0) {
	    k++;
	    me/=d;
	}
	else {
	    break;
	}
    }
    return k;
};



template <typename Data, int D> class LBNode;
template <int D> struct TreeCoords;
template <int D> struct Tree;
template <int D> class MyProcmap;
template <int D, typename Pmap> class LBTree;
class NodeData;

template <int D>
struct DClass {
    typedef Key<D> KeyD;
    typedef const Key<D> KeyDConst;
    typedef TreeCoords<D> TreeCoords;
    typedef Tree<D> Tree;
    typedef LBNode<NodeData,D> NodeD;
    typedef const LBNode<NodeData,D> NodeDConst;
    typedef MyProcmap<D> MyProcMap;
//    typedef WorldContainer< KeyD,NodeD,MyProcMap > treeT;
    typedef LBTree<D,MyProcMap> treeT;
};

template <typename T, int D, typename Pmap>
void migrate(SharedPtr<FunctionImpl<T,D,Pmap> > tfrom, SharedPtr<FunctionImpl<T,D,Pmap> > tto);

template <typename T, int D, typename Pmap>
void migrate_data(SharedPtr<FunctionImpl<T,D,Pmap> > tfrom, SharedPtr<FunctionImpl<T,D,Pmap> > tto, 
	typename DClass<D>::KeyD key);

template <typename Data, int D>
class LBNode {
private:
    Data data;
    std::vector<bool> c;

    void allchildren(bool status=false) {
        c.clear();
        c.assign(dim, status);
    };

public:
    static int dim;

    LBNode() {
	data = Data();
	allchildren();
    };

    LBNode(Data d, bool children=false) : data(d) {
	allchildren(children);
    };

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

    void set_child(int i, bool setto = true) {
	c[i] = setto;
    };
    
    void set_data(Data d) {
	data = d;
    };

    Data get_data() const {
	return data;
    };

    vector<bool> get_c() const {
	return c;
    };

    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & data & c;
    }
};


template <typename Data, int D>
std::ostream& operator<<(std::ostream& s, const LBNode<Data, D>& node) {
    s << "data = " << node.get_data() << ", c = " << node.get_c();
    return s;
};

template <int D>
std::ostream& operator<<(std::ostream& s, typename DClass<D>::NodeDConst& node) {
    s << "data = " << node.get_data() << ", c = " << node.get_c();
    return s;
};


template <typename Data, int D>
int LBNode<Data,D>::dim = power<D>();


class NodeData {
    friend std::ostream& operator<<(std::ostream& s, const NodeData& nd);
public:
    int cost;
    int subcost;
    bool istaken;
    NodeData(int c = 1, int s = 1, bool i = false) : cost(c), subcost(s), istaken(i) {};
    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & cost & subcost & istaken;
    };
    void print() {
	cout << "cost = " << cost << ", subcost = " << subcost << ", istaken = " << istaken << endl;
    };
};


inline std::ostream& operator<<(std::ostream& s, const NodeData& nd) {
    s << "cost " << nd.cost << ", subcost " << nd.subcost << ", istaken " << nd.istaken;
    return s;
};



template <int D>
struct TreeCoords {
    Key<D> key;
    ProcessID owner;

    TreeCoords(const Key<D> k, ProcessID o) : key(Key<D>(k)), owner(o) {};
    TreeCoords(const TreeCoords& t) : key(Key<D>(t.key)), owner(t.owner) {};
    TreeCoords() : key(Key<D>()), owner(-1) {};
    void print() const {
	madness::print(key, "   owner =", owner);
    };

    bool operator< (const TreeCoords t) const {
	return (this->key < t.key);
    };
};



template <int D>
struct Tree {
    TreeCoords<D> data;
    vector<SharedPtr<Tree> > children;
    Tree* parent;

    Tree() {};
    Tree(TreeCoords<D> d) : data(d), parent(0) {};
    Tree(TreeCoords<D> d, Tree* p) : data(d), parent(p) {};

    Tree(const Tree<D>& tree) : data(tree.data), parent(0) {};
    Tree(const Tree<D>& tree, Tree* p) : data(tree.data), parent(p) {};

    Tree<D>& operator=(const Tree<D>& other) {
	if (this != &other) {
	    this->data = other.data;
	    this->parent = other.parent;
	    this->children = other.children;
	}
        return *this;
    };

    void insertChild(TreeCoords<D> d) {
	Tree* c = new Tree(d, this);
	children.insert(children.begin(),SharedPtr<Tree<D> > (c));
    };

    void insertChild(const Tree<D>& tree) {
	Tree* c = new Tree(tree, this);
	children.insert(children.begin(),SharedPtr<Tree<D> > (c));
    };

    void print() {
	data.print();
	int csize = children.size();
	for (int j = 0; j < csize; j++) {
	    children[j]->print();
	}
    };

    bool isForeparentOf(Key<D> key) const {
	return (this->data.key.is_parent_of(key));
    };

    void findOwner(const Key<D> key, ProcessID *ow) const {
	if (this->isForeparentOf(key)) {
	    *ow = this->data.owner;
	    int csize = children.size();
	    for (int j = 0; j < csize; j++) {
		children[j]->findOwner(key, ow);
	    }
	}
    };

    bool fill(TreeCoords<D> node) {
	bool success = false;
	if (this->isForeparentOf(node.key)) {
	    int csize = children.size();
	    for (int i = 0; i < csize; i++) {
		if (children[i]->isForeparentOf(node.key)) {
		    success = children[i]->fill(node);
		}
	    }
	    if (!success) {
		this->insertChild(node);
		success = true;
	    }
	}
	return success;
    }
};



template <int D>
class MyProcmap {
private:
    int whichmap;
    const ProcessID owner;
    Tree<D> treeMap;
    typedef Key<D> KeyD;


    ProcessID getOwner(const KeyD& key) const {
	ProcessID owner;
	treeMap.findOwner(key, &owner);
	return owner;
    };


    void buildTreeMap(vector<TreeCoords<D> > v) {
	sort(v.begin(), v.end());
	int vlen = v.size();

	if (vlen == 0) throw "empty map!!!";

	treeMap = Tree<D>(v[vlen-1]);
	for (int j = vlen-2; j >= 0; j--) {
	    treeMap.fill(v[j]);
	}
    };
	

public:
    MyProcmap() : whichmap(0), owner(0) {};
    MyProcmap(World& world) : whichmap(1), owner(0) {
	int NP = world.nproc();
	const int level = nearest_power(NP, D);
	int twotoD = power<D>();
	int NPin = (int) pow((double)twotoD,level);
	vector<TreeCoords<D> > v;
	
	for (Translation i=0; i < (Translation)NPin; i++) {
	    KeyD key(level,i);
	    if ((i%twotoD) == 0) {
		key = key.parent(nearest_power(i, twotoD));
	    }
	    v.push_back(TreeCoords<D>(key,i));
	    buildTreeMap(v);
	}
    }; 
    MyProcmap(World& world, ProcessID owner) : whichmap(0), owner(owner) {};

    MyProcmap(World& world, vector<TreeCoords<D> > v) : whichmap(1), owner(1) {
	buildTreeMap(v);
	treeMap.print();
    };

    MyProcmap(const MyProcmap<D>& other) : whichmap(other.whichmap), owner(other.owner), treeMap(other.treeMap) {};

    MyProcmap<D>& operator=(const MyProcmap<D>& other) {
	if (this != &other) {
	    whichmap = other.whichmap;
	    owner = other.owner;
	    treeMap = other.treeMap;
	}
        return *this;
    };

    ProcessID operator()(const KeyD& key) const {
	if (whichmap == 0)
	    return owner;
	else
	    return getOwner(key);
    };

};

template <int D, typename Pmap=MyProcmap<D> > 
class LBTree : public WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD,Pmap> {
    // No new variables necessary
    public:
	typedef WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD, Pmap> dcT;
	LBTree() {};
	LBTree(World& world, const Pmap& pmap) : dcT(world,pmap) {
	};
	template <typename T>
	inline void init_tree(SharedPtr<FunctionImpl<T,D,Pmap> > f, typename DClass<D>::KeyDConst key) {
madness::print("beginning of init_tree");
	    // find Node associated with key
	    typename FunctionImpl<T,D,Pmap>::iterator it = f->find(key);
	    if (it == f->end()) return;
	    // convert Node to LBNode
	    if (!(it->second.has_children())) {
		typename DClass<D>::NodeD lbnode(false);
	        // insert into "this"
		this->insert(key, lbnode);
	    }
	    else {
		typename DClass<D>::NodeD lbnode(true);
	        // insert into "this"
		this->insert(key, lbnode);
		// then, call for each child
		for (KeyChildIterator<D> kit(key); kit; ++kit) {
		    this->init_tree<T>(f, kit.key());
		}
	    }
madness::print("end of init_tree");
	};

	// Methods:
	Cost fixCost(typename DClass<D>::KeyDConst& key);
	Cost depthFirstPartition(typename DClass<D>::KeyDConst& key,
        	vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces,
        	Cost totalcost = 0, Cost *maxcost = 0);
	void rollup(typename DClass<D>::KeyDConst& key);
	void meld(typename DClass<D>::KeyDConst& key);
	Cost makePartition(typename DClass<D>::KeyDConst& key, 
		vector<typename DClass<D>::KeyD>* klist, Cost partitionSize, 
		bool lastPartition, Cost usedUp, bool *atleaf);
	void removeCost(typename DClass<D>::KeyDConst& key, Cost c);
	Cost computeCost(typename DClass<D>::KeyDConst& key);

	// inherited methods
	typename WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD, Pmap>::iterator 
		end() {
	    return WorldContainer<typename DClass<D>::KeyD, typename DClass<D>::NodeD, Pmap>::end();
	};
	typename WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD, Pmap>::iterator 
		find(typename DClass<D>::KeyDConst& key) {
	    return WorldContainer<typename DClass<D>::KeyD, typename DClass<D>::NodeD, Pmap>::find(key);
	};
	Pmap get_procmap() {
	    return WorldContainer<typename DClass<D>::KeyD, typename DClass<D>::NodeD, Pmap>::get_procmap();
	};
};

template <typename T, int D, typename Pmap=MyProcmap<D> >
class LoadBalImpl {
    private:
	Function<T,D,Pmap> f;
//	typename DClass<D>::treeT skeltree;
	SharedPtr<typename DClass<D>::treeT> skeltree;

	void construct_skel(SharedPtr<FunctionImpl<T,D,Pmap> > f) {
	    skeltree = SharedPtr<typename DClass<D>::treeT>(new typename DClass<D>::treeT(f->world,
		f->get_procmap()));
	    typename DClass<D>::KeyD root(0);
	    skeltree->template init_tree<T>(f,root);
	};

    public:
	//Constructors
	LoadBalImpl() {};
	LoadBalImpl(Function<T,D,Pmap> f) : f(f) {
	    madness::print("LoadBalImpl (Function) constructor: f.impl", &f.impl);
	    construct_skel(f.impl);
	};
	~LoadBalImpl() {};

	//Methods
	inline void loadBalance() {
	    partition(findBestPartition());
	};

	vector<typename DClass<D>::TreeCoords> findBestPartition();

	void partition(vector<typename DClass<D>::TreeCoords> v) {
	    // implement partition: copy to new FunctionImpl and replace within f
	    Pmap pmap(f.impl->world, v);
	    SharedPtr<FunctionImpl<T,D,Pmap> > newimpl(new FunctionImpl<T,D,Pmap>(*(f.impl.get()),pmap));
	    madness::migrate<T,D,Pmap>(f.impl, newimpl);
	    f.impl = newimpl;
	};

	World& world() {
	    return (f.impl->world);
	};
};

CompCost computeCompCost(Cost c, int n);

Cost computePartitionSize(Cost cost, unsigned int parts);

}

#endif
