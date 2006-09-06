#ifndef SENDRECV_H
#define SENDRECV_H

/// \file sendrecv.h
/// \brief Defines send/recv for OctTree

#include <iostream>
#include <algorithm>
#include <cmath>
#include <list>

//#include <mad_types.h>
//#include <misc/print.h>
//#include <misc/communicator.h>
//#include <misc/shared_ptr.h>
//#include <misc/misc.h>
#include <serialize/mpiar.h>
#include <octtree/octtree.h>
using namespace madness;

#define EXIT_CODE -1
#define COMPUTE_LOCAL_COST 1
#define GLOBAL_COST 2
#define GLOB_PART_SUB_SIGNAL 3

namespace madness {

    template <class T>
    struct less {
        bool operator() (const T p1, const T p2) {
            if (!p1)
                return true;
            if (!p2)
                return false;
            return (!(*p1 < *p2));
        }
    };

    /****************************************************************************************/
    /* class RootList: shorthand, minimal way of communicating that an OctTree<T> needs to  */
    /* be moved, and from current_owner to future_owner.  All public for ease of use.	*/
    /****************************************************************************************/

    class RootList {
    public:
        Translation x, y, z;
        Level n;
        ProcessID current_owner, future_owner;
        RootList(): x(0), y(0), z(0), n(-1), current_owner(0), future_owner(0) {};
        RootList(Translation x0, Translation y0, Translation z0, Level n0, ProcessID current,
                 ProcessID future): x(x0), y(y0), z(z0), n(n0), current_owner(current),
        future_owner(future) {};
        template <class T>
        RootList(OctTree<T> *t, ProcessID current, ProcessID future): x(t->x()), y(t->y()),
        z(t->z()), n(t->n()), current_owner(current), future_owner(future) {};

        template <class Archive>
        inline void serialize(const Archive& ar) {
            ar & x & y & z & n & current_owner & future_owner;
        }

        bool equals(RootList r) {
            if ((n == r.n) && (x == r.x) && (y == r.y) && (z == r.z))
                return true;
            else
                return false;
        }

	template <class T>
        bool equals(OctTree<T>* r) {
//std::cout << "bool equals: at beginning" << std::endl;
            if ((n == r->n()) && (x == r->x()) && (y == r->y()) && (z == r->z()))
                return true;
            else
                return false;
        }

        bool isDescendant(RootList r) {
            if (r.n >= n)
                return false;

            Level dn = (Level) pow(2, n-r.n);
            if ((x/dn == r.x) && (y/dn == r.y) && (z/dn == r.z))
                return true;
            else
                return false;
        }


        friend bool operator < (const RootList& t1, const RootList& t2) {
            std::cout << "operator <: beginning:" << std::endl;
            std::cout << "operator <: t1: n = " << t1.n << " (" << t1.x << "," << t1.y << "," << t1.z << ")"<< std::endl;
            std::cout << "operator <: t2: n = " << t1.n << " (" << t2.x << "," << t2.y << "," << t2.z << ")"<< std::endl;
            if (t1.n > t2.n) {
                std::cout << "operator <: t1.n > t2.n" << std::endl;
                return true;
            } else if (t1.n < t2.n) {
                std::cout << "operator <: t1.n < t2.n" << std::endl;
                return false;
            } else {
                std::cout << "operator <: computing n1 and n2" << std::endl;
                Translation s1, s2, n1 = (Translation) pow(2,t1.n-1), n2 = (Translation) pow(2,t2.n-1);
                std::cout << "operator <: n1 = " << n1 << ", n2 = " << n2 << std::endl;
                if (n1 == 0) n1 = 1;
                if (n2 == 0) n2 = 1;
                s1 = (t1.x/n1)*4 + (t1.y/n1)*2 + t1.z/n1;
                s2 = (t2.x/n2)*4 + (t2.y/n2)*2 + t2.z/n2;
                std::cout << "operator <: computed s1 and s2: " << s1 << " " << s2 << std::endl;
                if (s1 < s2)
                    return true;
                else if (s1 > s2)
                    return false;
                else {
                    if (t1.x < t2.x)
                        return true;
                    else if (t1.x > t2.x)
                        return false;
                    else {
                        if (t1.y < t2.y)
                            return true;
                        else if (t1.y > t2.y)
                            return false;
                        else {
                            if (t1.z < t2.z)
                                return true;
                            else
                                return false;
                        }
                    }
                }
            }
        }
    };


    template <class T>
    void parallelPartition(OctTree<T> *root, std::vector<RootList> *globalList,
                           std::vector<OctTree<T>* > *treeList); 

    template <class T>
    void partitionManager(OctTree<T> *root, std::vector<RootList> *globalList,
                          std::vector<OctTree<T>* > *treeList); 

    template <class T>
    std::vector<RootList> globalPartitionSubtree(OctTree<T> *root, Cost *spaceleft, int p); 

    template <class T>
    void partitionWorker(std::vector<OctTree<T>* > *treeList); 

    template <class T>
    void partitionWorkerSub(OctTree<T> *root, std::vector<RootList> *localList, Cost *spaceleft,
                            std::vector<Cost> *costList, std::vector<RootList> *remoteList, int p); 

    template <class T>
    Cost computeGlobalCost(RootList root, OctTree<T> *tree, std::vector<OctTree<T>* > const treeList); 

    template <class T>
    Cost globalCostManager(RootList root, std::vector<OctTree<T>* > const treeList);

    template <class T>
    Cost globalCostWorker(RootList root, std::vector<OctTree<T>* > const treeList);

    template <class T>
    std::vector<RootList> findRemoteChildren(OctTree<T> *t, int p); 

    template <class T>
    void exchangeTrees(std::vector<RootList> *globalList, std::vector<SharedPtr<OctTree<T> > > *treeList,
                       bool glue);

    template <class T>
    void sendSubtree(OctTree<T> *tree, ProcessID me, ProcessID dest); 

    template <class T>
    void sendSubtree(SharedPtr<OctTree<T> > tree, ProcessID me, ProcessID dest);

    template <class T>
    void recvSubtree(OctTree<T> *t, ProcessID me, ProcessID source);

    template <class T>
    void recvSubtree(SharedPtr<OctTree<T> > p, ProcessID me, ProcessID source);

    template <class T>
    void sendMsg(T msg, ProcessID me, ProcessID dest);

    template <class T>
    void recvMsg(T *msg, ProcessID me, ProcessID source);

    template <class T>
    void glueTrees(std::vector<SharedPtr<OctTree<T> > > *treeList);


    template <class T>
    SharedPtr<OctTree<char> > createGhostTree(std::vector<SharedPtr<OctTree<T> > > *treeList, 
	std::vector<RootList> *ownerList);
    
    std::vector<RootList> createMergeList(std::vector<RootList> globalList, 
	std::vector<RootList> ownerList);

    template <class T>
    void serialLoadBalance(std::vector<SharedPtr<OctTree<T> > > *treeList);


    template <class T>
    SharedPtr<OctTree<T> > findPtr(SharedPtr<OctTree<T> > tree, Level n, Translation x, 
	Translation y, Translation z);

    template <class T>
    class depthFirstOrder {
	public:
	    bool operator() (const SharedPtr<OctTree<T> > t1, const SharedPtr<OctTree<T> > t2)
	    {
		int ans;
		std::cout << "depthFirstOrder: t1: n = " << t1->n() << ", (" << t1->x() << ","
			<< t1->y() << "," << t1->z() << ")" << std::endl;
		std::cout << "depthFirstOrder: t2: n = " << t2->n() << ", (" << t2->x() << ","
			<< t2->y() << "," << t2->z() << ")" << std::endl;
	        if (t1->n() == t2->n())
		{	
		    bool retval;
		    ans = ordering(RootList(t1->x(), t1->y(), t1->z(), t1->n(), 0, 0), 
				RootList(t2->x(), t2->y(), t2->z(), t2->n(), 0, 0));
		    if (ans == 1)
		    {
		        std::cout << "t1 < t2 so return true" << std::endl;
			retval = true;
		    }
		    else if (ans == 0)
		    {
			std::cout << "t1 !< t2 so return false" << std::endl;
			retval = false;
		    }
		    else
		    {
			std::cout << "t1 == t2 so false" << std::endl;
			retval = false;
		    }
		    return retval;
		}

	        if (t1->n() > t2->n())
       		{
		    std::cout << "depthFirstOrder: n1 > n2" << std::endl;
	            Translation x1 = t1->x(), y1 = t1->y(), z1 = t1->z();
	            Level dn = t1->n() - t2->n();
		    std::cout << "depthFirstOrder: dn = " << dn << std::endl;
		    Translation twotodn = (Translation) pow(2.0, dn);
        	    x1 /= twotodn; y1 /= twotodn; z1 /= twotodn;
               	    RootList r1 = RootList(x1, y1, z1, t2->n(), 0, 0);
               	    RootList r2 = RootList(t2->x(), t2->y(), t2->z(), t2->n(), 0, 0);
               	    ans = ordering(r1, r2);
		    bool retval;
		    if (ans == 1)
		    {
		        std::cout << "r1 < r2 so return true" << std::endl;
			retval = true;
		    }
		    else if (ans == 0)
		    {
		        std::cout << "r1 !< r2 so return false" << std::endl;
			retval = false;
		    }
		    else
		    {
			std::cout << "r1 == r2 but n1 > n2 so compare t1 with t2's first child" << std::endl;
			retval = tiebreaking(t1, t2);
		    }
		    return retval;
        	}
        	else
        	{
		    std::cout << "depthFirstOrder: n1 < n2" << std::endl;
            	    Translation x2 = t2->x(), y2 = t2->y(), z2 = t2->z();
                    Level dn = t2->n() - t1->n();
		    std::cout << "depthFirstOrder: dn = " << dn << std::endl;
		    Translation twotodn = (Translation) pow(2.0, dn);
                    x2 /= twotodn; y2 /= twotodn; z2 /= twotodn;
                    RootList r1 = RootList(t1->x(), t1->y(), t1->z(), t1->n(), 0, 0);
                    RootList r2 = RootList(x2, y2, z2, t1->n(), 0, 0);
                    ans = ordering(r1, r2);
		    bool retval;
		    if (ans == 1)
		    {
			std::cout << "r1 < r2 so return true" << std::endl;
			retval = true;
		    }
		    else if (ans == 0)
		    {
			std::cout << "r1 !< r2 so return false" << std::endl;
			retval = false;
		    }
		    else
		    {
			std::cout << "r1 == r2 but n1 < n2 so compare t1's first child with t2" << std::endl;
			retval = !(tiebreaking(t2, t1));
		    }
		    return retval;
                }
	    }

	    int ordering(RootList r1, RootList r2)
	    {
		long dx, dy, dz;
		dx = r1.x - r2.x; dy = r1.y - r2.y; dz = r1.z - r2.z;
		std::cout << "ordering: r1 = (" << r1.x << "," << r1.y << "," << r1.z << "), n = "
			<< r1.n << std::endl;
		std::cout << "ordering: r2 = (" << r2.x << "," << r2.y << "," << r2.z << "), n = "
			<< r2.n << std::endl;
		std::cout << "ordering: dx = " << dx << ", dy = " << dy << ", dz = " << dz 
			<< ", n = " << r1.n << std::endl;
		if ((r1.x/2 == r2.x/2) && (r1.y/2 == r2.y/2) && (r1.z/2 == r2.z/2))
		{
		    if (dx == 0)
		    {
			if (dy == 0)
			{
			    std::cout << "ordering: about to return and dz = " << dz << std::endl;
			    if (dz < 0)
				return 1;
			    else if (dz > 0)
				return 0;
			    else return -1;
			}
			else
			{
			    std::cout << "ordering: about to return and dy = " << dy << std::endl;
			    if (dy < 0)
				return 1;
			    else
				return -1;
			}
		    }
		    else
		    {
			std::cout << "ordering: about to return and dx = " << dx << std::endl;
			if (dx < 0)
			    return 1;
			else
			    return 0;
		    }
		}
		else
		{
		    std::cout << "ordering: recursively calling self " << std::endl;
		    return (ordering(RootList(r1.x/2, r1.y/2, r1.z/2, r1.n-1, 0, 0),
				RootList(r2.x/2, r2.y/2, r2.z/2, r2.n-1, 0, 0)));
		}
	    }

	    bool tiebreaking(SharedPtr<OctTree<T> > t1, SharedPtr<OctTree<T> > t2)
	    {
		// for the sake of argument, assume t1->n() > t2->n()

		bool ans, haskids = t2->isParent();
		if (haskids)
		{
		    SharedPtr<OctTree<T> > c = SharedPtr<OctTree<T> >();
		    for (int i = 0; i < 2 ;i++)
		    {
		        for (int j = 0; j < 2; j++)
		        {
			    for (int k = 0; k < 2; k++)
			    {
			        c = t2->child(i,j,k);
			        if (c.get()) break;
			    }
			    if (c.get()) break;
		        }
		        if (c.get()) break;
		    }
		    int ord;
		    if (c->n() == t1->n())
		    {
			ord = ordering(RootList(t1->x(), t1->y(), t1->z(), t1->n(), 0, 0),
				RootList(c->x(), c->y(), c->z(), c->n(), 0, 0));
			if (ord == 1)
			    ans = true;
			else
			    ans = false;
		    }
		    else
		    {
			ans = tiebreaking(t1, c);
		    }
		}
		else
		{
		    // no children, so compare ordering of t1 and t2
		    
	            Translation x1 = t1->x(), y1 = t1->y(), z1 = t1->z();
	            Level dn = t1->n() - t2->n();
		    std::cout << "depthFirstOrder: dn = " << dn << std::endl;
		    Translation twotodn = (Translation) pow(2.0, dn);
        	    x1 /= twotodn; y1 /= twotodn; z1 /= twotodn;
               	    RootList r1 = RootList(x1, y1, z1, t2->n(), 0, 0);
               	    RootList r2 = RootList(t2->x(), t2->y(), t2->z(), t2->n(), 0, 0);
		    int ord = ordering(r1, r2);
		    if (ord == 1)
			ans = true;
		    else
			ans = false;
		}
		return ans;
	    }
    };
}

#endif
