#ifndef WORLD_LOADBAL_DEUX
#define WORLD_LOADBAL_DEUX

#include <mra/mra.h>

namespace madness {

    template <int NDIM>
    class LBNodeDeux {
        typedef Key<NDIM> keyT;
        typedef LBNodeDeux<NDIM> nodeT;
        typedef WorldContainer<keyT,nodeT> treeT;
        volatile double child_cost[1<<NDIM];
        volatile double my_cost;
        volatile double total_cost;
        volatile bool gotkids;
        volatile int nsummed;
        
        /// Computes index of child key in this node using last bit of translations
        int index(const keyT& key) {
            int ind = 0;
            for (int d=0; d<NDIM; d++) ind += ((key.translation()[d])&0x1) << d;
            return ind;
        }

    public:
        LBNodeDeux() 
            : my_cost(0.0), total_cost(0.0), gotkids(false), nsummed(0) {
            for (int i=0; i<(1<<NDIM); i++) child_cost[i] = 0.0;
        }

        bool has_children() const {return gotkids;}

        double get_total_cost() const {return total_cost;}

        /// Accumulates cost into this node
        Void add(double cost, bool got_kids) {
            total_cost = (my_cost += cost);
            gotkids = gotkids || got_kids;
        }

        /// Accumulates cost up the tree from children
        Void sum(treeT* tree, const keyT& child, double value) {
            child_cost[index(child)] = value;
            nsummed++;
            if (nsummed == (1<<NDIM)) {
                for (int i=0; i<(1<<NDIM); i++) total_cost += child_cost[i];
                if (child.level() > 1) {
                    keyT key = child.parent();
                    keyT parent = key.parent();
                    tree->send(parent, &nodeT::sum, tree, key, double(total_cost));
                }
            }
        }


        /// Logically deletes this node by setting cost to -1

        /// Cannot actually erase this node from the container since the send() handler
        /// is holding an accessor to it.
        Void deleter(treeT* tree, const keyT& key) {
            total_cost = my_cost = -1.0;
            if (has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT child = kit.key();
                    tree->send(child, &nodeT::deleter, tree, child);
                }
            }
        }
        

        /// Descends tree deleting all except internal nodes and sub-tree parents
        Void partition(treeT* tree, const keyT& key, double avg) {
            if (has_children()) {
                // Sort children in descending cost order
                keyT keys[1<<NDIM];
                double vals[1<<NDIM];
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT child = kit.key();
                    int ind = index(child);
                    keys[ind] = child;
                    vals[ind] = child_cost[ind];
                }
                for (int i=0; i<(1<<NDIM); i++) {
                    for (int j=i+1; j<(1<<NDIM); j++) {
                        if (vals[i] < vals[j]) {
                            std::swap(vals[i],vals[j]);
                            std::swap(keys[i],keys[j]);
                        }
                    }
                }
                
                // Split off subtrees in decreasing cost order
                for (int i=0; i<(1<<NDIM); i++) {
                    if (total_cost <= avg) {
                        tree->send(keys[i], &nodeT::deleter, tree, keys[i]);
                    }
                    else {
                        total_cost -= vals[i];
                        tree->send(keys[i], &nodeT::partition, tree, keys[i], avg);
                    }
                }
            }
        }

        /// Printing for the curious
        Void print(treeT* tree, const keyT& key) const {
            for(int i=0; i<key.level(); i++) std::cout << "  ";
            madness::print(key, my_cost, total_cost);
            if (gotkids) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    tree->send(kit.key(), &nodeT::print, tree, kit.key());
                }
            }
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & archive::wrap_opaque(this,sizeof(nodeT));
        }
    };
        

    template <int NDIM>
    class LoadBalanceDeux {
        typedef Key<NDIM> keyT;
        typedef LBNodeDeux<NDIM> nodeT;
        typedef WorldContainer<keyT,nodeT> treeT;
        typedef typename treeT::iterator iteratorT;
        World& world;
        treeT tree;
        

        template <typename T, typename costT>
        struct add_op {
            LoadBalanceDeux* lb;
            const costT& costfn;
            add_op(LoadBalanceDeux* lb, const costT& costfn) : lb(lb), costfn(costfn) {}
            void operator()(const keyT& key, const FunctionNode<T,NDIM>& node) const {
                lb->tree.send(key, &nodeT::add, costfn(key,node), node.has_children());
            }
        };

        /// Sums costs up the tree returning to everyone the total cost
        double sum() {
            world.gop.fence();
            for (iteratorT it=tree.begin(); it!=tree.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (!node.has_children() && key.level() > 0) {
                    tree.send(key.parent(), &nodeT::sum, &tree, key, node.get_total_cost());
                }
            }
            world.gop.fence();
            double total;
            keyT key0(0);
            if (world.rank() == tree.owner(key0)) {
                total = tree.find(key0).get()->second.get_total_cost();
            }
            world.gop.broadcast(total, tree.owner(key0));
            world.gop.fence();

            return total;
        }

        void doprint() {
            keyT key0(0);
            if (world.rank() == tree.owner(key0)) {
                tree.send(key0, &nodeT::print, &tree, key0);
            }
            world.gop.fence();
        }


    public:
        LoadBalanceDeux(World& world) 
            : world(world)
            , tree(world, FunctionDefaults<NDIM>::get_pmap())
        {};

        /// Accumulates cost from a function
        template <typename T, typename costT>
        void add(const Function<T,NDIM>& f, const costT& costfn, bool fence=true) {
            const_cast<Function<T,NDIM>&>(f).unaryop_node(add_op<T,costT>(this,costfn), fence);
        }

        /// Performs the partitioning of the tree
        void partition() {
            double avg = sum()/world.size();
            avg = avg/4.0;
            print("The average cost is", avg);
            
            keyT key0(0);
            if (world.rank() == tree.owner(key0)) {
                tree.send(key0, &nodeT::partition, &tree, key0, avg*1.1);
            }
            world.gop.fence();

            // Delete dead nodes
            for (iteratorT it=tree.begin(); it!=tree.end();) {
                const keyT key = it->first;
                const nodeT& node = it->second;
                ++it;
                if (node.get_total_cost() < 0) tree.erase(key);
            }
            world.gop.fence();
            doprint();
        }
        
        
    };
}


#endif

