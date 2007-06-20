#ifndef LOADBAL_H
#define LOADBAL_H

namespace madness {

    typedef int Cost;
    typedef double CompCost;

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

    template <typename Data, int D> class LBNode;
    template <int D> struct TreeCoords;
    template <int D> struct Tree;
    template <int D> class MyPmap;
    template <int D> class LBTree;
    class NodeData;

    template <int D>
    struct DClass {
        typedef Key<D> KeyD;
        typedef const Key<D> KeyDConst;
        typedef TreeCoords<D> TreeCoords;
        typedef Tree<D> Tree;
        typedef LBNode<NodeData,D> NodeD;
        typedef const LBNode<NodeData,D> NodeDConst;
        typedef MyPmap<D> MyPmap;
        typedef LBTree<D> treeT;
    };

//     template <typename T, int D>
//     void migrate(SharedPtr<FunctionImpl<T,D> > tfrom, SharedPtr<FunctionImpl<T,D> > tto);

//     template <typename T, int D>
//     void migrate_data(SharedPtr<FunctionImpl<T,D> > tfrom, SharedPtr<FunctionImpl<T,D> > tto,
//                       typename DClass<D>::KeyD key);

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

        void print() const {
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
//madness::print("findOwner: at node", this->data.key);
            if (this->isForeparentOf(key)) {
//madness::print("findOwner: node", this->data.key, "is foreparent of", key, "so owner =", this->data.owner);
                *ow = this->data.owner;
                if (this->data.key.level() < key.level()) {
                    int csize = children.size();
                    for (int j = 0; j < csize; j++) {
//madness::print("findOwner: recursively call on ", this->children[j]->data.key);
                        children[j]->findOwner(key, ow);
                    }
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
    class MyPmap : public WorldDCPmapInterface< Key<D> > {
    private:
        bool staticmap;
        const ProcessID staticmap_owner;
        Tree<D>* treeMap;
        typedef Key<D> KeyD;

        void buildTreeMap(vector<TreeCoords<D> > v) {
            sort(v.begin(), v.end());
            int vlen = v.size();

            if (vlen == 0) throw "empty map!!!";

            treeMap = new Tree<D>(v[vlen-1]);
            for (int j = vlen-2; j >= 0; j--) {
                treeMap->fill(v[j]);
            }
        };


    public:
        MyPmap() : staticmap(false), staticmap_owner(0) {};

        MyPmap(World& world) : staticmap(false), staticmap_owner(0) {
            int NP = world.nproc();
            int twotoD = power<D>();
            const int level = nearest_power(NP, twotoD);
            int NPin = (int) pow((double)twotoD,level);
            vector<TreeCoords<D> > v;

            for (Translation i=0; i < (Translation)NPin; i++) {
                KeyD key(level,i);
                if ((i%twotoD) == 0) {
                    key = key.parent(nearest_power(NPin-i, twotoD));
                }
                v.push_back(TreeCoords<D>(key,i));
            }
            buildTreeMap(v);
            madness::print("MyPmap constructor");
            treeMap->print();
        };

        MyPmap(World& world, ProcessID owner) : staticmap(true), staticmap_owner(owner) {};

        MyPmap(World& world, vector<TreeCoords<D> > v) : staticmap(false), staticmap_owner(0) {
            buildTreeMap(v);
            madness::print("");
            treeMap->print();
        };

        MyPmap(const MyPmap<D>& other) : staticmap(other.staticmap), staticmap_owner(other.staticmap_owner), treeMap(other.treeMap) {};

        MyPmap<D>& operator=(const MyPmap<D>& other) {
            if (this != &other) {
                staticmap = other.staticmap;
                owner = other.owner;
                treeMap = other.treeMap;
            }
            return *this;
        };

        void print() const {
            treeMap->print();
        };

        ProcessID owner(const KeyD& key) const {
            if (staticmap)
                return staticmap_owner;
            else {
                ProcessID owner;
                treeMap->findOwner(key, &owner);
                return owner;
            }
        };
    };

    template <int D>
    class LBTree : public WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD> {
        // No new variables necessary
    public:
        typedef WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD> dcT;
        LBTree() {};
        LBTree(World& world, const SharedPtr< WorldDCPmapInterface<typename DClass<D>::KeyD> >& pmap) : dcT(world,pmap) {
            madness::print("LBTree(world, pmap) constructor");
            const MyPmap<D>* ppp = &(this->get_mypmap());
	    ppp->print();
            madness::print("LBTree(world, pmap) constructor (goodbye)");
        };
        template <typename T>
        inline void init_tree(SharedPtr< FunctionImpl<T,D> > f, typename DClass<D>::KeyDConst key) {
            // find Node associated with key
            typename FunctionImpl<T,D>::dcT::iterator it = f->coeffs.find(key);
            if (it == f->coeffs.end()) return;
            // convert Node to LBNode
            NodeData nd;
            if (!(it->second.has_children())) {
                typename DClass<D>::NodeD lbnode(nd,false);
                // insert into "this"
                this->insert(key, lbnode);
            } else {
                typename DClass<D>::NodeD lbnode(nd,true);
                // insert into "this"
                this->insert(key, lbnode);
                // then, call for each child
                for (KeyChildIterator<D> kit(key); kit; ++kit) {
                    this->init_tree<T>(f, kit.key());
                }
            }
        };

        // Methods:
        void print(typename DClass<D>::KeyDConst& key) {
            typename DClass<D>::treeT::iterator it = this->find(key);
            if (it == this->end()) return;
            for (Level i = 0; i < key.level(); i++) cout << "  ";
            madness::print(key, it->second);
            for (KeyChildIterator<D> kit(key); kit; ++kit) {
                print(kit.key());
            }
        };

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
        typename WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD>::iterator 
        end() {
            return WorldContainer<typename DClass<D>::KeyD, typename DClass<D>::NodeD>::end();
        };

        typename WorldContainer<typename DClass<D>::KeyD,typename DClass<D>::NodeD>::iterator
        find(typename DClass<D>::KeyDConst& key) {
            return WorldContainer<typename DClass<D>::KeyD, typename DClass<D>::NodeD>::find(key);
        };

//         const SharedPtr<WorldDCPmapInterface< typename DClass<D>::KeyD >& get_pmap() {
//             return WorldContainer<typename DClass<D>::KeyD, typename DClass<D>::NodeD>::get_pmap();
//         };

        MyPmap<D>& get_mypmap() {
            return *static_cast< MyPmap<D>* >(this->get_pmap().get());
        };

    };

    template <typename T, int D>
    class LoadBalImpl {
    private:
	typedef MyPmap<D> Pmap;
        Function<T,D> f;
        SharedPtr<typename DClass<D>::treeT> skeltree;

        void construct_skel(SharedPtr<FunctionImpl<T,D> > f) {
            skeltree = SharedPtr<typename DClass<D>::treeT>(new typename DClass<D>::treeT(f->world,
                       f->coeffs.get_pmap()));
            typename DClass<D>::KeyD root(0);
            madness::print("about to initialize tree");
            if (f->world.mpi.rank() == 0) {
                skeltree->template init_tree<T>(f,root);
            }
            madness::print("just initialized tree");
        };

    public:
        //Constructors
        LoadBalImpl() {};

        LoadBalImpl(Function<T,D> f) : f(f) {
            madness::print("LoadBalImpl (Function) constructor: f.impl", &f.get_impl());
            construct_skel(f.get_impl());
        };

        ~LoadBalImpl() {};

        //Methods

        SharedPtr< WorldDCPmapInterface< Key<D> > > loadBalance() {
            return SharedPtr< WorldDCPmapInterface< Key<D> > >(new MyPmap<D>(f.get_impl()->world, findBestPartition()));
        };

        vector<typename DClass<D>::TreeCoords> findBestPartition();
    };

    CompCost computeCompCost(Cost c, int n);

    Cost computePartitionSize(Cost cost, unsigned int parts);

}

#endif
