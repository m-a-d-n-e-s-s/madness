// #include <madness/world/MADworld.h>
// #include <madness/misc/misc.h>
// #include <madness/tensor/tensor.h>
// #include <madness/misc/ran.h>
// #include <madness/mra/key.h>
// #include <madness/mra/funcimpl.h>


This is not presently in use but is left here since it is actually useful.  It provides a process map that can use a cost function for partitioning subtrees.  Originally written by Rebecca Hartman-Baker.

namespace madness {
    /// Procmap implemented using Tree of TreeCoords

    template <int D>
    class MyPmap : public WorldDCPmapInterface< Key<D> > {
    private:
        unsigned int map_type; // 0 = simple map, 1 = gaussian distributed map, 2 = treecoords list
        const int nproc;
        const int n;
        std::shared_ptr< ProcMapImpl<D> > tree_map; // for map_type 2
        Tensor<ProcessID> simple_key_map; // map of keys at level n, for map_type 1
        typedef Key<D> KeyD;

        /// private method that builds the Tree underlying the procmap
        void build_tree_map(std::vector< TreeCoords<D> > v) {
            tree_map = std::shared_ptr< ProcMapImpl<D> > (new ProcMapImpl<D>(v));
        }

        ProcessID simple_hash(const KeyD& key) const {
	    if (key.level() == 0) return 0;
            KeyD parent = (key.level() > n) ? key.parent(key.level()-n) : key;
            return (parent.hash()%nproc);
        }
        ProcessID not_so_simple_hash(const KeyD& key) const {
            KeyD parent = (key.level() > n) ? key.parent(key.level()-n) : key;
            return simple_key_map((const long *) &(parent.translation()[0]));
        }

        void prepare_not_so_simple_map(World& world) {
	    std::vector<long> vdim(D);
            for (int i=0; i<D; ++i) vdim[i] = 1L<<n;
            simple_key_map = Tensor<ProcessID>(vdim);

            std::list< std::pair<KeyD,double> > costmap;
            Vector<Translation,D> l;
            long cent = (1L<<n) / 2;
            for (TensorIterator<ProcessID> iter=simple_key_map.unary_iterator(0,false,false); iter._p0; ++iter) {
                double dist = 0.01;
                for (int i=0; i<D; ++i) {
                    l[i] = iter.ind[i];
                    dist += (l[i] - cent)*(l[i] - cent);
                }
                double cost = 1.0/dist; // actually dist squared
                cost *= (1.0 + 0.001*RandomValue<double>()); // To shuffle (nearly) equal values
                costmap.push_back(std::pair<KeyD,double>(KeyD(n,l),cost));
            }
            costmap.sort(costmapcmp);
//            if (world.rank() == 0) {
//                for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
//                    madness::print("costmap", it->first, it->second);
//                }
//            }
            ProcessID p = 0;
            for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
                const long *l = (const long *) &(it->first.translation()[0]);
                simple_key_map(l)  = p;
                ++p;
                if (p == world.size()) p = 0;
            }
//            if (world.rank() == 0) {
//                madness::print("SIMPLE MAP", D,"\n", simple_key_map);
//            }
        }

    public:
        MyPmap() : map_type(2) {};

        static bool costmapcmp(const std::pair<KeyD,double>& a, const std::pair<KeyD,double>& b) {
            return a.second > b.second;
        }
        MyPmap(World& world)
                : map_type(1)
                , nproc(world.nproc())
                , n(int((std::log((double)world.size())/std::log(2.0)+3)/D) + 2) { // 16*nproc = 2^(nD)
            // We set n to have about 16 tasks per processor and we try to
            // give each process a mix of large, medium, and small
            // tasks.  Currently estimate cost as inversely
            // proportional to distance from center but we could
            // enable the user to provide a function.

            //if (world.rank() == 0) madness::print("DIM",D,"N IN MAP IS",n);
            prepare_not_so_simple_map(world);
        }

        MyPmap(World& world, unsigned int map_type, int n=100)
                : map_type(map_type)
                , nproc(world.nproc())
                , n(n) {
            if (map_type==1) {
                n =int((std::log((double)world.size())/std::log(2.0)+3)/D) + 2; // 16*nproc = 2^(nD)
                prepare_not_so_simple_map(world);
            }
        }

        MyPmap(World& world, std::vector<TreeCoords<D> > v) : map_type(2), nproc(world.nproc()), n(0) {
            build_tree_map(v);
        }

        MyPmap(const MyPmap<D>& other) : map_type(other.map_type), nproc(other.nproc), n(other.n), tree_map(other.tree_map) {};

        MyPmap<D>& operator=(const MyPmap<D>& other) {
            if (this != &other) {
                map_type = other.map_type;
		simple_key_map = other.simple_key_map; // shallow copy
                nproc = other.nproc;
                n = other.n;
                tree_map = other.tree_map;
            }
            return *this;
        }

        void print() const {
            if (map_type == 2) {
                tree_map->print();
            } else if (map_type == 1) {
	        madness::print("MyPmap: gaussian distributed map with n =", n);
	    } else {
                madness::print("MyPmap: simple map with n =", n);
            }
        }

        /// Find the owner of a given key
        ProcessID owner(const KeyD& key) const {
	    if (map_type == 0) {
                return simple_hash(key);
	    } else if (map_type == 1) {
	        return not_so_simple_hash(key);
            } else {
                return tree_map->find_owner(key);
            }
        }
    };
}
