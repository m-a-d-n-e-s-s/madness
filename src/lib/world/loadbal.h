#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#ifndef LOADBAL_H
#define LOADBAL_H

#include <world/world.h>
#include <mra/key.h>
//using namespace madness;
using namespace std;

namespace madness {

typedef int Cost;
typedef double CompCost;

template <typename Data, unsigned int D> class LBNode;
template <unsigned int D> class Key;
template <unsigned int D> struct TreeCoords;
template <unsigned int D> struct Tree;
template <unsigned int D> class MyProcmap;
class NodeData;

template <unsigned int D>
struct DClass {
    typedef Key<D> KeyD;
    typedef const Key<D> KeyDConst;
    typedef TreeCoords<D> TreeCoords;
    typedef Tree<D> Tree;
    typedef LBNode<NodeData,D> NodeD;
    typedef const LBNode<NodeData,D> NodeDConst;
    typedef MyProcmap<D> MyProcMap;
    typedef WorldContainer< KeyD,NodeD,MyProcMap > treeT;
};


template <typename Data, unsigned int D>
class LBNode {
private:
    Data data;
    std::vector<bool> c;

    void nochildren() {
        c.clear();
        c.assign(dim, false);
    };

public:
    static unsigned int dim;

    LBNode() {
	data = Data();
	nochildren();
    }
    LBNode(Data d) : data(d) {
	nochildren();
    };

    unsigned int compute_index(vector<unsigned int> v) {
	unsigned int vlen = v.size();
	unsigned int index = 0, twoD = 1;
	for (unsigned int i = vlen-1; i >= 0; i--) {
	    index+= (twoD*v[i]);
	    twoD*=2;
	}
	return index;
    };
	
    unsigned int compute_index(vector<int> v) {
	unsigned int vlen = v.size();
	unsigned int index = 0, twoD = 1;
	for (unsigned int i = vlen-1; i >= 0; i--) {
	    index+= (twoD*v[i]);
	    twoD*=2;
	}
	return index;
    };
	
    bool has_children() const {
	for (unsigned int i = 0; i < dim; i++)
	    if (c[i]) return true;
	return false;
    };

    bool has_child(unsigned int i) const {
	return c[i];
    };

    bool has_child(int i) const {
	return c[i];
    };

    bool has_child(vector<int> vi) const {
	return c[compute_index(vi)];
    };

    void set_child(int i, bool setto = true) {
	c[i] = setto;
    };
    
    void set_child(vector<int> vi, bool setto = true) {
	c[compute_index(vi)] = setto;
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


template <typename Data, unsigned int D>
std::ostream& operator<<(std::ostream& s, const LBNode<Data, D>& node) {
    s << "data = " << node.get_data() << ", c = " << node.get_c();
    return s;
};

template <unsigned int D>
std::ostream& operator<<(std::ostream& s, typename DClass<D>::NodeDConst& node) {
    s << "data = " << node.get_data() << ", c = " << node.get_c();
    return s;
};


template <typename Data, unsigned int D>
unsigned int LBNode<Data,D>::dim = (unsigned int) pow(2.0, (int) D);


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



/*
template <unsigned int D>
class Key {
public:
    unsigned int n;
    mutable Array<unsigned int, D> L;
    hashT hashval; // Saving this here is a big optimization

    // Recomputes hashval
    void rehash() {
        hashval = madness::hash(n,madness::hash(L));
    };

    Key() {};

    Key(const Key& k) : n(k.n), L(k.L) {
       hashval = k.hashval;
    };

    Key(int n, int i, int j) : n(n) {
	L[0] = i; L[1] = j;
	rehash();
    };

    Key(int n, vector<int> v) : n(n), L(Array<unsigned int, D>(v)) {
	rehash();
    };

    Key(unsigned int n, vector<unsigned int> v) : n(n), L(Array<unsigned int, D>(v)) {
	rehash();
    };

    Key(unsigned int n, Array<unsigned int, D> L) : n(n), L(L) {
	rehash();
    };

    hashT hash() const {
	return hashval;
    };


    Key myChild(int k) const {
	Array<unsigned int,D> LL;
	for (unsigned int i = 0; i < D; i++) {
	    LL[i] = 2*L[i] + k%2;
	    k/=2;
	}
	return Key(n+1, LL);
    };

    Key myParent(int k=1) const {
	if (k == 0) return Key(*this);
	Array<unsigned int,D> LL;
	unsigned int twotok = (unsigned int) pow(2.0, k);
	for (unsigned int i = 0; i < D; i++) {
	    LL[i] = (L[i]/twotok);
	}
	return Key(n-k, LL);
    };


    bool isChildOf(const Key& key) const {
	if (*this == key) {
	    return true; // I am child of myself
	}
	else if (this->n <= key.n) {
	    return false; // I can't be child of something lower on the tree
	}
	else {
	    unsigned int dn = this->n - key.n;
	    Key mama = this->myParent(dn);
	    return (mama == key);
	}
    };

    bool isParentOf(const Key& key) const {
	return key.isChildOf(*this);
    };
	    

    int ordering(const Key& k1, const Key& k2) const {
	bool egalite = true;
	Array<int,D> dL;
//	cout << "ordering: comparing ";
//	k1.print();
//	cout << " and ";
//	k2.print();
//	cout << endl;
	for (unsigned int i = 0; i < D; i++) {
	    dL[i] = k1.L[i] - k2.L[i];
	    if (k1.L[i]/2 != k2.L[i]/2) {
		egalite = false;
//		cout << "ordering: k1 and k2 do not have same parent" << endl;
	    }
	}
	if (!egalite) {
	    return (ordering(k1.myParent(), k2.myParent()));
	}
	else {
	    for (unsigned int i = 0; i < D; i++) {
		if (dL[i] > 0) {
//		    cout << "ordering: dL[" << i << "] > 0; returning -1" << endl;
		    return -1;
		}
		else if (dL[i] < 0) {
//		    cout << "ordering: dL[" << i << "] < 0; returning 1" << endl;
		    return 1;
		}
	    }
//	    cout << "ordering: no dL greater or less than zero; returning 0" << endl;
	    return 0;
	}
    };

    bool operator==(const Key& a) const {
	if (n != a.n) return false;
	for (unsigned int i = 0; i < D; i++) {
	    if (L[i] != a.L[i]) return false;
	}
	return true;
    };

    bool operator<(const Key& a) const {
	int ans;

	if (n == a.n) {
	    ans = ordering(*this, a);
	}
	else if (n > a.n) {
	    unsigned int dn = n - a.n;
	    Key newthis = this->myParent(dn);
	    if (newthis == a) {
		ans = 1;
	    }
	    else {
	    	ans = ordering(newthis, a);
	    }
	}
	else {
	    unsigned int dn = a.n - n;
	    Key newa = a.myParent(dn);
	    if (newa == *this) {
		ans = 0;
	    }
	    else {
	    	ans = ordering(*this, newa);
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


    void print() const {
	std::cout << "n = " << n << ", (";
	for (unsigned int i = 0; i < D-1; i++) {
	    std::cout << L[i] << ", ";
	}
	std::cout << L[D-1] << ")";
	std::cout << "  hash=" << hashval;
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & n & L & hashval;
    }
};
*/


template <unsigned int D>
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


template <unsigned int D>
struct Tree {
    TreeCoords<D> data;
    vector<Tree*> children;
    Tree* parent;

    Tree() {};
    Tree(TreeCoords<D> d) : data(d), parent(0) {};
    Tree(TreeCoords<D> d, Tree* p) : data(d), parent(p) {};

    void insertChild(TreeCoords<D> d) {
	Tree* c = new Tree(d, this);
	children.insert(children.begin(),c);
    };

    void print() {
	data.print();
	int csize = children.size();
	for (int j = 0; j < csize; j++) {
	    children[j]->print();
	}
    };

    bool isForeparentOf(Key<D> key) const {
//	return (this->data.key.isParentOf(key));
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

    bool fill(vector<TreeCoords<D> > v, int j) {
	bool success = false;
	if (this->isForeparentOf(v[j].key)) {
	    int csize = children.size();
	    for (int i = 0; i < csize; i++) {
		if (children[i]->isForeparentOf(v[j].key)) {
		    success = children[i]->fill(v, j);
		}
	    }
	    if (!success) {
		this->insertChild(v[j]);
		success = true;
	    }
	}
	return success;
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



template <unsigned int D>
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
//	    treeMap.fill(v, j);
	    treeMap.fill(v[j]);
	}
    };
	

public:
    MyProcmap() : whichmap(0), owner(0) {};
    MyProcmap(ProcessID owner) : whichmap(0), owner(owner) {};

    MyProcmap(vector<TreeCoords<D> > v) : whichmap(1), owner(1) {
	buildTreeMap(v);
	treeMap.print();
    };

    ProcessID operator()(const KeyD& key) const {
	if (whichmap == 0)
	    return owner;
	else
	    return getOwner(key);
    };

};


template <unsigned int D>
void build_tree(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key);

template <unsigned int D>
void print_tree(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key);

template <unsigned int D>
Cost computeCost(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key);

template <unsigned int D>
void meld(typename DClass<D>::treeT& tree, typename DClass<D>::KeyDConst& key);

template <unsigned int D>
void rollup(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key);

template <unsigned int D>
Cost fixCost(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key); 

Cost computePartitionSize(Cost cost, unsigned int parts);

template <unsigned int D>
Cost makePartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, 
	vector<typename DClass<D>::KeyD >* klist, Cost partitionSize, bool lastPartition, 
	Cost usedUp, bool *atleaf);

template <unsigned int D>
Cost depthFirstPartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, 
	vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces, Cost totalcost, Cost *maxcost); 

template <unsigned int D>
void removeCost(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, Cost c); 

CompCost computeCompCost(Cost c, int n);

template <unsigned int D>
void findBestPartition(typename DClass<D>::treeT tree, typename DClass<D>::KeyD key, 
	vector<typename DClass<D>::TreeCoords>* klist, unsigned int npieces);

template <unsigned int D>
void migrate_data(typename DClass<D>::treeT tfrom, typename DClass<D>::treeT tto, 
	typename DClass<D>::KeyD key);

template <unsigned int D>
void migrate(typename DClass<D>::treeT tfrom, typename DClass<D>::treeT tto); 

}

#endif
