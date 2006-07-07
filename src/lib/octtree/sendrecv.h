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
            //std::cout << "operator <: beginning" << std::endl;
            if (t1.n > t2.n) {
                //std::cout << "operator <: t1.n > t2.n" << std::endl;
                return true;
            } else if (t1.n < t2.n) {
                //std::cout << "operator <: t1.n < t2.n" << std::endl;
                return false;
            } else {
                //std::cout << "operator <: computing n1 and n2" << std::endl;
                Translation s1, s2, n1 = (Translation) pow(2,t1.n-1), n2 = (Translation) pow(2,t2.n-1);
                //std::cout << "operator <: n1 = " << n1 << ", n2 = " << n2 << std::endl;
                if (n1 == 0) n1 = 1;
                if (n2 == 0) n2 = 1;
                s1 = (t1.x/n1)*4 + (t1.y/n1)*2 + t1.z/n1;
                s2 = (t2.x/n2)*4 + (t2.y/n2)*2 + t2.z/n2;
                //std::cout << "operator <: computed s1 and s2" << std::endl;
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
}

#endif
