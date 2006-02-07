#ifndef OCTTREE_H
#define OCTTREE_H

/// \file octtree.h
/// \brief Defines OctTreeLayout and OctTree

#include <iostream>
#include <algorithm>
#include <cmath>
#include <list>

#include <mad_types.h>
#include <misc/print.h>
#include <misc/communicator.h>
#include <misc/shared_ptr.h>
#include <serialize/archive.h>


#define FORIJK(expr) \
    do { \
       for (int i=0; i<2; i++) \
          for (int j=0; j<2; j++) \
   	     for (int k=0; k<2; k++) {expr} \
    } while(0)

#define FOREACH_CHILD(OctTreeT, t, expr) \
    do { \
      for (int i=0; i<2; i++) \
        for (int j=0; j<2; j++) \
	  for (int k=0; k<2; k++) { \
            OctTreeT *child = ((t)->child(i,j,k)); \
            if (child) {expr} \
          } \
    } while(0)

#define FOREACH_LOCAL_CHILD(OctTreeT, t, expr) \
    do { \
      for (int i=0; i<2; i++) \
        for (int j=0; j<2; j++) \
	  for (int k=0; k<2; k++) { \
            OctTreeT *child = ((t)->child(i,j,k)); \
            if (child && child->islocal()) {expr} \
          } \
    } while(0)

#define FOREACH_REMOTE_CHILD(OctTreeT, t, expr) \
    do { \
      for (int i=0; i<2; i++) \
        for (int j=0; j<2; j++) \
	  for (int k=0; k<2; k++) { \
            OctTreeT *child = ((t)->child(i,j,k)); \
            if (child && child->isremote()) {expr} \
          } \
    } while(0)

namespace madness {
    
    // declare it up here so that it will be recognized by OctTree:

    template <class T> class OctTreeFlat;

    /// Raise integer two to a positive power 
    inline Translation two_to_power(long n) {return Translation(1)<<n;}


    /// NULL operation to help generalize & combine tree recursions
    template <typename T> static inline void no_op(const T *t) {};

    /// A distributed OctTree that contains connectivity and data
    
    /// A local node will have \c _remote=false.  If any of its parent or
    /// child pointers are null, it means those nodes do not exist anywhere.
    ///
    /// A node with a null parent pointer is the parent of the local
    /// subtree.  If such a node is local, it is the parent of the
    /// entire tree.
    ///
    /// A remote node will have \c _remote=true and \c _remote_proc will
    /// contain the rank of the remote process in the default
    /// communicator.   Any non-null child or parent pointers will point to
    /// local nodes.  Null pointers only indicate that the information for
    /// that node is stored remotely ... i.e., only local connection info
    /// is stored and the remote node may, or may not, have children.
    template <class T>
    class OctTree {
    private:
        typedef OctTree<T> OctTreeT;
        Translation _x;             ///< x translation (0,...,2**n-1)
        Translation _y;             ///< y translation (0,...,2**n-1)
        Translation _z;             ///< z translation (0,...,2**n-1)
        Level _n;                     ///< Level in tree
        
        bool _remote;               ///< True if this node is remote
        ProcessID _rank;     ///< if (_remote) The rank of the remote process

        OctTreeT* _p;          ///< Points to parent node, or null if no parent
        shared_ptr<OctTreeT> _c[2][2][2]; ///< Null or Points to child nodes, or null if not there

        // (why not shared ptr ?)
        const Communicator* _comm; ///< Info on parallel processes

//        T _data;			///< The payload stored by value

	Cost _cost;			///< Cost associated with node
	Cost _localSubtreeCost;		///< Cost associated with local parts of node's subtree
	bool _visited;			///< Whether a node has been visited and assigned to a partition

        /// Bit reversal of integer of given length
        static unsigned long bit_reverse(unsigned long i, int nbits) {
            unsigned long r = 0;
            while (nbits--) {
                r = (r<<1) | (i&1);
                i = i>>1;
            }
            return r;
        };
        
        /// Map level + translation to process in default mapping
        static ProcessID default_map_to_rank(const Communicator& comm, Level n, 
                                             Translation x, Translation y, Translation z) {
            long npx, npy, npz;
            comm.mesh(npx, npy, npz);
            long twon = two_to_power(n);
            long px = bit_reverse(x,n)&(std::min(twon,npx)-1);
            long py = bit_reverse(y,n)&(std::min(twon,npy)-1);
            long pz = bit_reverse(z,n)&(std::min(twon,npz)-1);
            return comm.rank_from_coords(px,py,pz);
        };
        
        /// Recursively add missing remote children
        void add_remote_children(Level nmax) {
            if (n() == nmax) return;
            
            FORIJK(OctTreeT *c = child(i,j,k);
                   if (c) 
                   c->add_remote_children(nmax);
                   else if (islocal()) {
                       Level in = n()+1;
                       Translation ix = ((x())<<1) + i;
                       Translation iy = ((y())<<1) + j;
                       Translation iz = ((z())<<1) + k;
                       ProcessID rank = default_map_to_rank(*_comm,in,ix,iy,iz);
                       insert_remote_child(i,j,k,rank);
                   });
        };
        
    public:
	// can't figure out another way to get to the FunctionNode without allowing copy
        T _data;			///< The payload stored by value
        /// Default constructor makes empty node with n=-1 to indicate invalid.
        OctTree() :
            _x(0), _y(0), _z(0), _n(-1),
            _remote(false),  _rank(-1), 
            _p(0), _comm(0)
        {};

        /// Constructor makes node with most info provided
        OctTree(Level n, Translation x, Translation y, Translation z,
                      bool remote, OctTreeT* parent, ProcessID remote_proc,
                      const Communicator* comm) :
            _x(x), _y(y), _z(z), _n(n),
            _remote(remote),  _rank(remote_proc), 
	    _p(parent), _comm(comm), _cost(1), _localSubtreeCost(1),
	    _visited(false)
            {
                FORIJK(_c[i][j][k] = 0;);
            };

        /// Constructor makes node with all info provided
        OctTree(Level n, Translation x, Translation y, Translation z,
                      bool remote, OctTreeT* parent, ProcessID remote_proc,
                      const Communicator* comm, Cost cost, 
		      Cost localSubtreeCost, bool visited) :
            _x(x), _y(y), _z(z), _n(n),
            _remote(remote),  _rank(remote_proc), 
	    _p(parent), _comm(comm), _cost(cost), 
	    _localSubtreeCost(localSubtreeCost), _visited(visited)
            {
                FORIJK(_c[i][j][k] = 0;);
            };

	/// Constructor makes node from flattened data type, with parent
	OctTree(OctTreeFlat<T> t, OctTreeT* parent)
	{
	    _x = t.x(); _y = t.y(); _z = t.z(); _n = t.n();
	    _remote = t.remote(); _rank = t.rank();
	    _cost = t.cost(); _localSubtreeCost = t.localSubtreeCost();
	    _data = t.data(); _visited = t.visited();
	    FORIJK(_c[i][j][k] = 0;);
	    _p = parent;
	};

	/// Constructor makes node from flattened data type, with parent (except data)
	OctTree(OctTreeFlat<T> t, OctTreeT* parent, bool makeData)
	{
	    _x = t.x(); _y = t.y(); _z = t.z(); _n = t.n();
	    _remote = t.remote(); _rank = t.rank();
	    _cost = t.cost(); _localSubtreeCost = t.localSubtreeCost();
	    _visited = t.visited();
	    if (makeData)
		_data = T(t.data());
	    FORIJK(_c[i][j][k] = 0;);
	    _p = parent;
	};

	/// Constructor makes node from flattened data type  
	OctTree(OctTreeFlat<T> t)
	{
	    _x = t.x(); _y = t.y(); _z = t.z(); _n = t.n();
	    _remote = t.remote(); _rank = t.rank();
	    _cost = t.cost(); _localSubtreeCost = t.localSubtreeCost();
	    _data = t.data(); _visited = t.visited();
	    FORIJK(_c[i][j][k] = 0;);
	    _p = 0;
	};

	/// Constructor makes node from flattened data type (except data)
	OctTree(OctTreeFlat<T> *t, bool makeData)
	{
	    _x = t->x(); _y = t->y(); _z = t->z(); _n = t->n();
	    _remote = t->remote(); _rank = t->rank();
	    _cost = t->cost(); _localSubtreeCost = t->localSubtreeCost();
	    _visited = t->visited();
	    if (makeData)
		_data = T(t->data());
std::cout << "_data" << std::endl;
std::cout << "   sample from it: " << _data[1](1,1,1) << std::endl;;
	    FORIJK(_c[i][j][k] = 0;);
	    _p = 0;
	};

        ~OctTree() {
        }


        /// Returns a reference to the data
        inline T& data() {return _data;};

        /// Returns a const reference to the data
        inline const T& data() const {return _data;};

	/// Set data
	inline void setData(T data) {_data = data;};

	inline void setData(T *data)
	{
	    _data = *data;
	}

        /// Get cost of node
        inline Cost getCost() {return _cost;};

        /// Set cost of node
        inline void setCost(Cost cost) {_cost = cost;};

        /// Get cost for local subtree
        inline Cost getLocalSubtreeCost() {return _localSubtreeCost;};

        /// Level of node in tree (0 is highest)
        inline long n() const {return _n;};
        
	/// Set x index of node (0,...,2**n-1)
	inline void setx(Translation x) {_x = x;};

        /// x index of node (0,...,2**n-1)
        inline Translation x() const {return _x;};
        
	/// Set y index of node (0,...,2**n-1)
	inline void sety(Translation y) {_y = y;};

        /// y index of node (0,...,2**n-1)
        inline Translation y() const {return _y;};
        
	/// Set z index of node (0,...,2**n-1)
	inline void setz(Translation z) {_z = z;};

        /// z index of node (0,...,2**n-1)
        inline Translation z() const {return _z;};
        
	/// Set remote (true if node is remote, false otherwise)
	inline void setRemote(bool remote) {_remote = remote;};

        /// returns true if node is local, false otherwise
        inline bool islocal() const {return !_remote;};
        
        /// returns true if node is remote, false otherwise
        inline bool isremote() const {return _remote;};

	/// returns true if parent to node is remote, false otherwise
	inline bool isLocalRoot() {return parent().isremote();};

	/// returns true if node has children
	inline bool isParent() {return (!(child(0,0,0) == 0));}
        
        /// returns true if node is the parent of the local subtree
        inline bool islocalsubtreeparent() const {return (_p==0);};
        
	/// returns true if node has been visited by partitioner
	inline bool isVisited() {return _visited;};

        /// returns the communicator pointer
        inline const Communicator* comm() const {return _comm;};
        
	/// Set ID of owning remote process (-1 if not remote)
	inline void setRank(ProcessID rank) {_rank = rank;};

        /// if (isremote()) returns id of owning remote process (-1 if not remote)
        inline ProcessID rank() const {return _rank;};
        
        /// returns pointer to parent (null if absent or if both this node & parent are remote)
        inline OctTreeT* parent() const {return _p;};
        
        /// returns pointer to child (null if absent or if both this node & child are remote)
        inline OctTreeT* child(int x, int y, int z) const {
            return (OctTreeT*) _c[x][y][z];
        };

	/// insert child from flattened data structure, returning pointer to child
	OctTreeT* insert_flat_child(OctTreeFlat<T> t)
	{
	    Translation x, y, z, px, py, pz;
	    x = t.x();
	    y = t.y();
	    z = t.z();
	    px = this->x();
	    py = this->y();
	    pz = this->z();

	    x -= (2*px);
	    y -= (2*py);
	    z -= (2*pz);

	    _c[x][y][z] = shared_ptr<OctTreeT>(new OctTreeT(t, this));
	    return _c[x][y][z];
	};
        
	/// insert child (except data) from flattened data structure, returning pointer to child
	OctTreeT* insert_flat_child(OctTreeFlat<T> t, bool makeData)
	{
	    Translation x, y, z, px, py, pz;
	    x = t.x();
	    y = t.y();
	    z = t.z();
	    px = this->x();
	    py = this->y();
	    pz = this->z();

	    x -= (2*px);
	    y -= (2*py);
	    z -= (2*pz);

	    _c[x][y][z] = shared_ptr<OctTreeT>(new OctTreeT(t, this, makeData));
	    return _c[x][y][z];
	};
        
        /// insert local child (x, y, z in {0,1}) returning pointer to child
        OctTreeT* insert_local_child(int x, int y, int z) {
            _c[x][y][z] = shared_ptr<OctTreeT>(new OctTreeT(_n + 1,
                                                       2*this->_x + x,
                                                       2*this->_y + y,
                                                       2*this->_z + z,
                                                       false,
                                                       this,
                                                       -1,
                                                       _comm));
            return _c[x][y][z];
        };
        
        /// insert remote child (x, y, z in {0,1}) returning pointer to child
        OctTreeT* insert_remote_child(int x, int y, int z, ProcessID remote_proc) {
            _c[x][y][z] = shared_ptr<OctTreeT>(new OctTreeT(_n + 1,
                                                       2*this->_x + x,
                                                       2*this->_y + y,
                                                       2*this->_z + z,
                                                       true,
                                                       this,
                                                       remote_proc,
                                                       _comm));
            return _c[x][y][z];
        };

        /// Starting from the current node, walk the tree to find the requested node.
        
        /// If the requested node exists and is local, return a pointer to
        /// it.  Otherwise, return a pointer to the deepest (largest n)
        /// existing node of which the requested node would be a
        /// descendent (i.e., its parent, or grand parent, etc.).  If it
        /// is not known in the local sub-tree, return 0.  The returned
        /// node might be remote.
        OctTreeT* find(Level n, Translation x, Translation y, Translation z) const {
            if (n >= _n) {
                /// Determine if the desired box is a descendent of this node
                Level nn = n - _n;
                Translation xx = (x>>nn), yy = (y>>nn), zz = (z>>nn);
                if (xx==_x && yy==_y && zz==_z) {
                    if (nn == 0) return (OctTreeT *) this;	// Found it!
                    nn--;
                    xx = (x>>nn)&1; 
                    yy = (y>>nn)&1;
                    zz = (z>>nn)&1;
                    OctTreeT* c = child(xx,yy,zz);
                    if (c) 
                        return c->find(n,x,y,z); // Look below
                    else {
                        return (OctTreeT* ) this; // Does not exist, return closest parent
                    }
                }
            }
            if (_p) return _p->find(n,x,y,z); // Pass up the tree
            return 0;			// Not in the (sub)tree we are connected to
        };

	/// Set _visited to true for this node and all its subnodes, and _rank to p
	void setVisited(ProcessID p)
	{
	    _visited = true;
	    this->setRank(p);

	    FOREACH_LOCAL_CHILD(OctTreeT, this,
		child->setVisited(p);
	    );
	}

	/// Set _visited to true for this node and all its subnodes
	void setVisited()
	{
	    _visited = true;

	    FOREACH_LOCAL_CHILD(OctTreeT, this,
		child->setVisited();
	    );
	}

	/// Depth-first traversal of tree (prints out diagnostic info)
	void depthFirstTraverse()
	{
	    std::cout << "layer " << n() << ", (x,y,z) = " << x() << ", "
			<< y() << ", " << z() << std::endl;
	    std::cout << "      hasChildren? " << this->isParent()
		 << "    isRemote? " << this->isremote()
		 << "    rank? " << this->rank() << std::endl;
	    FOREACH_LOCAL_CHILD(OctTreeT, this,
		child->depthFirstTraverse();
		);
	}

	/// Tally number of nodes of tree (probably different than cost)
	int tallyNodes()
	{
	    int nnodes = 1;
	    if (isParent())
		FOREACH_CHILD(OctTreeT, this,
		    if (child->isremote())
		    {
			nnodes++;
		    }
		    else
		    {
			nnodes += child->tallyNodes();
		    }
		);
	    return nnodes;
	}
        
	/// Compute cost of subtree
	int computeCost()
	{
	    int total = 0;
	    if (isParent())
	    {
		FOREACH_LOCAL_CHILD(OctTreeT, this,
		    total += child->computeCost();
		);
	    }
	    total += getCost();
	    _localSubtreeCost = total;
	    return total;
	}

	/// Compute cost of subtree and set all nodes' visited flag
	int computeCost(bool visited)
	{
	    int total = 0;
	    if (isParent())
	    {
		FOREACH_LOCAL_CHILD(OctTreeT, this,
		    total += child->computeCost();
		);
	    }
	    total += getCost();
	    _localSubtreeCost = total;
	    _visited = visited;
	    return total;
	}

        /// Apply a user-provided function (or functor) to all nodes
        template <class Op1, class Op2>
        void recur(Op1 op1, Op2 op2) {
            op1(this);
            FOREACH_CHILD(OctTreeT, this, child->recur(op1,op2););
            op2(this);
        }
        
        /// Apply the operation at each node starting at the top of each sub-tree
        template <class Operation>
        void recur_down(Operation op) {
            recur<Operation, void (*) (const OctTreeT *)>(op, no_op<OctTreeT>);
        }
        
        /// Apply the operation at each node starting at the bottom of each sub-tree
        template <class Operation>
        void recur_up(Operation op) {
            recur<void (*) (const OctTreeT *), Operation>(no_op<OctTreeT>,op);
        }

        /// Example of how to walk down the tree
        
        /// Can be made more efficient by bundling the communication from
        /// parent to a node (this can be done in the function class by using
        /// ghost data in the parent)
        void walk_down() const {
            int buf;
            if (parent() && parent()->isremote()) {
                print_coords();
                madness::print(" ",_comm->rank()," receiving from",parent()->rank());
                _comm->Recv(buf, parent()->rank(), 1);
            }
            
            // do work here if not remote
            
            FORIJK(const OctTreeT* c = child(i,j,k);
                   if (c) {
                       if (c->isremote()) {
                           c->print_coords();
                           madness::print(" ",_comm->rank()," sending to",c->rank());
                           _comm->Send(buf, c->rank(), 1);
                       }
                       else {
                           c->walk_down();
                       }
                   });
        };
        
        
        /// Example of how to walk up the tree
        void walk_up() const {
            int buf;
            FORIJK(const OctTreeT* c = child(i,j,k);
                   if (c) {
                       c->walk_up();
                       if (c->isremote()) {
                           c->print_coords();
                           madness::print(" ",_comm->rank()," receiving from",c->rank());
                           _comm->Recv(buf, c->rank(), 1);
                       }
                   });
            
            // do work here if not remote
            
            if (parent() && parent()->isremote()) {
                print_coords();
                madness::print(" ",_comm->rank()," sending to",parent()->rank());
                _comm->Send(buf, parent()->rank(), 1);
            }
        };
        
        /// Print the coords to std::cout (mostly for debugging)
        void print_coords() const {
            std::cout << n() << "(" << x() << ", " << y() << ", " << z() << ")";
            if (isremote()) {
                if (islocalsubtreeparent()) 
                    std::cout << " <--- " << rank();
                else
                    std::cout << " ---> " << rank();
            }
        };
        
        void print() const {
            std::cout << _comm->rank() << ": ";
            for (int i=0; i<n(); i++) std::cout << "  ";
            print_coords(); 
            std::cout << std::endl;
            FORIJK(if (_c[i][j][k]) _c[i][j][k]->print(););
        };
        

        /// Create a default or initial parallel decomposition for an OctTree
        
        /// This is currently uniform partitioning down to level nmax but
        /// needs improving.  This call involves no communication.  The
        /// utility of bit reversed process ids in each dimension should be
        /// clear from the below picture of a tree with the binary value of
        /// the process number.  Note that at each level if the binary string
        /// is reversed, the processes are numbered sequentially going across
        /// a row.  Also note that we can mask off the highest bits (or the
        /// lowest bits of the bit reversed string) to identify parents.
        ///
        /// \verbatim
        ///                                        0
        ///                                        0
        ///                    0                                       1
        ///                    0                                       1
        ///          0                   2                   1                    3
        ///         00                  10                  01                   11
        ///
        ///      0        4         2         6         1         5         3         7
        ///     000      100       010       110       001       101       011       111
        ///
        ///   0    8    4   12    2   10    6   14    1    9    5   13    3   11    7   15
        /// 0000 1000 0100 1100 0010 1010 0110 1110 0001 1001 0101 1101 0011 1011 0111 1111
        /// \endverbatim
        static OctTreeT* create_default(const Communicator& comm, Level nmax) {
            OctTreeT* p = 0;
            
            long npx, npy, npz;         // dimensions of 3D process mesh
            comm.mesh(npx, npy, npz);
            
            // Ensure all processes are included by using a fine enough level
            Translation nproc_max = std::max(std::max(npx, npy),npz);
            while (two_to_power(nmax) < nproc_max) nmax++;
            
            // Create the sub-tree of nodes that I own.  Must also add a
            // remote parent if necessary.
            for (Level n=0; n<=nmax; n++) {
                Translation twon=two_to_power(n);
                for (Translation ix=0; ix<twon; ix++) {
                    for (Translation iy=0; iy<twon; iy++) {
                        for (Translation iz=0; iz<twon; iz++) {
                            ProcessID rank = OctTreeT::default_map_to_rank(comm,n,ix,iy,iz);
                            if (rank == comm.rank()) {
                                if (p==0 && rank!=0) { 
                                    // Add remote parent before adding any children
                                    ProcessID prank = OctTreeT::default_map_to_rank(comm, n-1,
                                                                                         ix>>1,iy>>1,iz>>1);
                                    p = new OctTreeT(n-1,ix>>1,iy>>1,iz>>1,true,0,prank,&comm);
                                }
                                
                                if (p) {
                                    // Add child to my sub-tree
                                    OctTreeT* c = p->find(n, ix, iy, iz);
                                    if (!c) {
                                        madness::print("looking for",n,ix,iy,iz);
                                        p->print();
                                        throw "find failed while adding?";
                                    }
                                    c->insert_local_child(ix&1, iy&1, iz&1);
                                }
                                else {
                                    // This node is (0,0,0,0) and I am process 0
                                    if (rank != 0) throw "confusion adding parent?";
                                    p = new OctTreeT(n,ix,iy,iz,false,0,-1,&comm);
                                }
                            }
                        }
                    }
                }
            }
            
            // Missing children above the lowest level correspond
            // to connections to remote nodes which we now make.
            p->add_remote_children(nmax);
            
            return p;
        };
        
    };

    template <class T>
    class OctTreeFlat {
    private:
	Translation _x;		///< x translation (0,...,2**n-1)
	Translation _y;		///< y translation (0,...,2**n-1)
	Translation _z;		///< z translation (0,...,2**n-1)
	Level _n; 		///< Level in tree

	// NOTE: All bools changed to ints because bool is not C type

	int _remote; 		///< True if this node is remote

	Cost _cost;		///< Cost associated with node
	Cost _localSubtreeCost;	///< Cost associated with local parts of node's subtree
	int _visited;		///< Whether a node has been visited and assigned to a partition
	int _hasChildren;	///< New variables, used in reconstruction
	int _hasParent;

	ProcessID _rank;	///< if (_remote) The rank of remote process
//	T _data;		///< The payload stored by value

	typedef OctTreeFlat<T> OctTreeFlatT;
	typedef typename std::list<OctTreeFlatT,std::allocator<OctTreeFlatT> >::iterator LI;

    public:
	T _data;		///< The payload stored by value
	/// Default constructor
	OctTreeFlat():
	    _x(0), _y(0), _z(0), _n(-1), _remote(0), _rank(-1),
	    _cost(1), _localSubtreeCost(1), _visited(0)
	{
//	    _x = 0; _y = 0; _z = 0; _n = -1; _remote = 0; _rank = -1;
//	    _cost = 1; _localSubtreeCost = 1; _visited = 0;
	};

	/// Constructor with some info
	OctTreeFlat(Level n, Translation x, Translation y, Translation z,
		int remote, ProcessID remote_proc):
	    _x(x), _y(y), _z(z), _n(n), _remote(remote_proc), _rank(-1),
	    _cost(1), _localSubtreeCost(1), _visited(0)
	{
//	    _x = x; _y = y; _z = z; _n = n; _remote = remote_proc; _rank = -1;
//	    _cost = 1; _localSubtreeCost = 1; _visited = 0;
	};
	
	/// Constructor made with OctTree (the crucial constructor)
	OctTreeFlat(OctTree<T> t)
	{
	    _x = t.x(); _y = t.y(); _z = t.z(); _n = t.n();
	    _remote = (int) t.isremote(); _rank = t.rank(); _cost = t.getCost();
	    _localSubtreeCost = t.getLocalSubtreeCost();
	    _visited = (int) t.isVisited();
	    if (t.isParent() && t.islocal())
	    {
		_hasChildren = 1;
	    }
	    else
	    {
		_hasChildren = 0;
	    }
	    if (t.parent() == 0)
	    {
		_hasParent = 0;
	    }
	    else
	    {
		_hasParent = 1;
	    }
	    _rank = t.rank();
	    _data = t.data();
	};
	
	/// Constructor made with OctTree (the crucial constructor)
	OctTreeFlat(OctTree<T> *t, bool makeData)
	{
std::cout << "making OctTreeFlat in constructor" << std::endl;
	    _x = t->x(); _y = t->y(); _z = t->z(); _n = t->n();
	    _remote = (int) t->isremote(); _rank = t->rank(); _cost = t->getCost();
	    _localSubtreeCost = t->getLocalSubtreeCost();
	    _visited = (int) t->isVisited();
//	    _x = t.x(); _y = t.y(); _z = t.z(); _n = t.n();
//	    _remote = (int) t.isremote(); _rank = t.rank(); _cost = t.getCost();
//	    _localSubtreeCost = t.getLocalSubtreeCost();
//	    _visited = (int) t.isVisited();
//	    if (t.isParent() && t.islocal())
	    if (t->isParent() && t->islocal())
	    {
		_hasChildren = 1;
	    }
	    else
	    {
		_hasChildren = 0;
	    }
//	    if (t.parent() == 0)
	    if (t->parent() == 0)
	    {
		_hasParent = 0;
	    }
	    else
	    {
		_hasParent = 1;
	    }
//	    _rank = t.rank();
	    _rank = t->rank();
	    if (makeData)
	    	_data = T(t->data());
////	    	_data = t.data();
//std::cout << "_data" << std::endl;
//std::cout << "   sample from it: " << std::endl;
	};

	~OctTreeFlat() {};

	template <class Archive>
        inline void serialize(const Archive& ar) {
//	    ar & _x & _y & _z & _n & _remote & _rank & _cost & _localSubtreeCost & _visited & _hasChildren & _hasParent & _rank;} 
	    ar & _x & _y & _z & _n & _remote & _rank & _cost & _localSubtreeCost & _visited & _hasChildren & _hasParent & _rank & _data;} 

	inline Translation x() {return _x;};
	inline Translation y() {return _y;};
	inline Translation z() {return _z;};
	inline Level n() {return _n;};
	inline int remote() {return _remote;};
	inline ProcessID rank() {return _rank;};
	inline T data() {return _data;};
	inline Cost cost() {return _cost;};
	inline Cost localSubtreeCost() {return _localSubtreeCost;};
	inline int visited() {return _visited;};
	inline bool hasChildren() {return _hasChildren;};
	inline bool hasParent() {return _hasParent;};

	inline void setData(T *data)
	{
	    &_data = data;
	}
    };
}

#endif
