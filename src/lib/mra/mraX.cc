#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
using std::abs;

#include <algorithm>
using std::max;
using std::sort;

#include <complex>

#include "mraX.h"
#include "twoscale.h"
#include "legendre.h"
#include <octtree/sendrecv.cc>

/// \file mra/mra.cc
/// \brief Implements Function and friends

/// Eventually, these notes provide a basic user's
/// programmer's manual.

namespace madness {

    // Definition and initialization of FunctionDefaults static members
    // It cannot be an instance of FunctionFactory since we want to
    // set the defaults independent of the data type.  The defaults
    // are only used in the FunctionFactory constructor except
    // for cell[][] which is used in FunctionData
    SharedPtr<FunctionOctTree> FunctionDefaults::tree = SharedPtr<FunctionOctTree> ();
//    SharedPtr<FunctionOctTree> FunctionDefaults::tree = 0;
    int FunctionDefaults::k = 7;
    double FunctionDefaults::thresh = 1e-5;
    int FunctionDefaults::initial_level = 2;
    int FunctionDefaults::max_refine_level = 30;
    int FunctionDefaults::truncate_method = 0;
    bool FunctionDefaults::compress = true;
    bool FunctionDefaults::refine = true;
    bool FunctionDefaults::autorefine = false;
    bool FunctionDefaults::debug = false;
    double FunctionDefaults::cell[3][2] = {{0.0,1.0}, {0.0,1.0}, {0.0,1.0}};

    template <typename T>
    FunctionCommonData<T> FunctionData<T>::commondata[FunctionData::MAXK+1];

    template <typename T>
    bool FunctionData<T>::initialized = false;

    template <typename T>
    void FunctionCommonData<T>::_make_dc_periodic() {
        // See ABGV for details
        r0 = Tensor<double>(k,k);
        rp = Tensor<double>(k,k);
        rm = Tensor<double>(k,k);

        double iphase = 1.0;
        for (int i=0; i<k; i++) {
            double jphase = 1.0;
            for (int j=0; j<k; j++) {
                double gammaij = sqrt(double((2*i+1)*(2*j+1)));
                double Kij;
                if (((i-j)>0) && (((i-j)%2)==1))
                    Kij = 2.0;
                else
                    Kij = 0.0;

                r0(i,j) = 0.5*(1.0 - iphase*jphase - 2.0*Kij)*gammaij;
                rm(i,j) = 0.5*jphase*gammaij;
                rp(i,j) =-0.5*iphase*gammaij;
                jphase = -jphase;
            }
            iphase = -iphase;
        }

        // Make the rank-1 forms of rm and rp
        rm_left = Tensor<double>(k);
        rm_right = Tensor<double>(k);
        rp_left = Tensor<double>(k);
        rp_right = Tensor<double>(k);

        iphase = 1.0;
        for (int i=0; i<k; i++) {
            double gamma = sqrt(0.5*(2*i+1));
            rm_left(i)  = rp_right(i) = gamma;
            rm_right(i) = rp_left(i)  = gamma*iphase;
            iphase *= -1.0;
        }
        rp_left.scale(-1.0);

//         Tensor<double> rm_test = outer(rm_left,rm_right);
//         Tensor<double> rp_test = outer(rp_left,rp_right);
    }

    template <typename T>
    void FunctionCommonData<T>::_init_twoscale() {
        if (! two_scale_hg(k, &hg)) throw "failed to get twoscale coefficients";
        hgT = transpose(hg);
        hgsonly = copy(hg(Slice(0,k-1),_));
    }

    template <typename T>
    void FunctionCommonData<T>::_init_quadrature() {
        quad_x = Tensor<double>(npt);
        quad_w = Tensor<double>(npt);
        quad_phi = Tensor<double>(npt,k);
        quad_phiw = Tensor<double>(npt,k);

        gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
        for (int i=0; i<npt; i++) {
            double phi[200];
            legendre_scaling_functions(quad_x(i),k,phi);
            for (int j=0; j<k; j++) {
                quad_phi(i,j) = phi[j];
                quad_phiw(i,j) = quad_w(i)*phi[j];
            }
        }
        quad_phit = transpose(quad_phi);
    }

    template <typename T>
    void Function<T>::_init(const FunctionFactory<T>& factory) {
        bool compress = factory._compress;
        bool refine = factory._refine;
        bool empty = factory._empty;

        if (data->f || data->vf) {
//std::cout << "_init: data->f or data->vf " << data->f << " " << data->vf << std::endl;
	    int tlen = data->trees->nTrees();
//std::cout << "_init: tlen = " << tlen << std::endl;
	    for (int i = 0; i < tlen; i++)
	    {
//std::cout << "_init: about to _fine_scale_project" << std::endl;
            	_fine_scale_projection(tree(i), data->initial_level);
//std::cout << "_init: done with _fine_scale_project" << std::endl;
            	if (refine) _refine(tree(i));
//std::cout << "_init: done with _refine" << std::endl;
	    }
            if (compress) this->compress();
        } else if (empty) {             // Do not set any coefficients at all
        } else {                        // Set coefficients so have zero function
	    int tlen = data->trees->nTrees();
	    for (int i = 0; i < tlen; i++)
	    {
            	if (tree(i)->n() == 0) set_coeff(tree(i),TensorT(2*k,2*k,2*k));
	    }
        }
        data->compressed = compress; // Just to be sure we are consistent.

        //tree()->print();
    }

    template <typename T>
    void Function<T>::_fine_scale_projection(OctTreeT* tree, Level initial_level) {
        // Recur down oct_tree to initial_level.  Refine locally
        // if necessary.  Project when get to desired level
//std::cout << "_fine_scale_projection: at beginning of function " << tree << std::endl;
        if (tree->n() < initial_level) {
//std::cout << "_fine_scale_projection: n < initial_level" << std::endl;
            set_active(tree);
//std::cout << "_fine_scale_projection: just set tree active" << std::endl;
            FORIJK(OctTreeT* c = tree->child(i,j,k);
//std::cout << "_fine_scale_projection: for child(" << i << "," << j << "," << k << ")" << std::endl;
                   if (!c && islocal(tree)) c = tree->insert_local_child(i,j,k);
                   if (c) _fine_scale_projection(c,initial_level););
        } else if (tree->n()==initial_level) {
//std::cout << "_fine_scale_projection: n == initial_level" << std::endl;
            set_active(tree);
//std::cout << "_fine_scale_projection: just set tree active" << std::endl;
            if (islocal(tree)) _project(tree);
//std::cout << "_fine_scale_projection: just projected the tree" << std::endl;
        }
//std::cout << "_fine_scale_projection: at end of function" << std::endl;
    }


    /// Private:  Projects function in given box
    template <typename T>
    void Function<T>::_project(OctTreeT* tree) {
        // We are at level n in the of wavelet coefficients, which
        // corresponds to level n+1 in the scaling function
        // coefficients.  Evaluate all 8 blocks of the scaling
        // function coefficients which are stored together in the tree.
        int npt = data->cdata->npt;
        Tensor<T> r(2*k,2*k,2*k);
        Tensor<T> fval(npt,npt,npt);
        Level n = tree->n();
        Translation x=2*tree->x(), y=2*tree->y(), z=2*tree->z();
        const Slice* s = data->cdata->s;
        for (int ix=0; ix<2; ix++) {
            for (int iy=0; iy<2; iy++) {
                for (int iz=0; iz<2; iz++) {
                    if (data->vf) {
                        _vfcube(n+1, x+ix, y+iy, z+iz, data->vf, fval);
                    } else if (data->f) {
                        _fcube(n+1, x+ix, y+iy, z+iz, data->f, fval);
                    }
                    // Can optimize away the tensor construction in transform3d
                    r(s[ix],s[iy],s[iz]) = transform3d(fval,data->cdata->quad_phiw);
                }
            }
        }
        set_coeff(tree,r.scale(pow(0.125,0.5*n)));
    }

    template <typename T>
    void Function<T>::_compress(OctTreeT* tree) {
//std::cout << "_compress: at beginning of function" << std::endl;
        FOREACH_CHILD(OctTreeT, tree,
                      if (isactive(child)) _compress(child););

        if (isremote(tree)) {
//std::cout << "_compress: tree n = " << tree->n() << ", (" << tree->x() << "," << tree->y() << 
//	"," << tree->z() << "), rank " << tree->rank() << " is remote" << std::endl;
            long k3 = k*k*k;
            if (tree->islocalsubtreeparent()) {
                // Send data from active children to remote parent
//std::cout << "_compress: this tree is localsubtreeparent" << std::endl;
                std::vector<Slice>& s0 = data->cdata->s0;
                TensorT& buf = data->cdata->work1;
                for (int i=0; i<2; i++) {
		    for (int j=0; j<2; j++) {
			for (int k=0; k<2; k++) {
			    OctTreeT *child = tree->child(i,j,k);
			    if (child) {
//                FOREACH_CHILD(OctTreeT, tree,
                              if (isactive(child)) {
//std::cout << "_compress: about to make c, coeffs of child" << std::endl;
                              const TensorT* c = coeff(child);
                                  if (!c) throw "compress: yo! show me the data(2)!";
//std::cout << "_compress: just made c, coeffs of child" << std::endl;
                                  buf(s0) = (*c)(s0);
//std::cout << "_compress: just did the buf(s0) thingy" << std::endl;
              //                                  print("compress: sending",buf.normf(),"to",tree->rank());
//std::cout << "_compress: want to send data on n = " << child->n() << ", (" << child->x() << "," <<
//	child->y() << "," << child->z() << ") to processor " << tree->rank() << std::endl;
//                                  comm()->Send(buf.ptr(), k3, tree->rank(), 2);
				MPI::COMM_WORLD.Send(buf.ptr(), k3, MPI_DOUBLE, tree->rank(), 2);
//std::cout << "_compress: sent data on n = " << child->n() << ", (" << child->x() << "," <<
//	child->y() << "," << child->z() << ") to processor " << tree->rank() << std::endl;
              //                                  print("child norm before",c->normf());
                                  (*c)(s0) = T(0.0);
              //                                  print("child norm after",c->normf());
                              }
			    }
			}
		    }
		}
		//            );
            } else {
                // Receive data from a remote child and store as
                // ghost ... parent will eventually delete it.
//std::cout << "_compress: before setting coeff TensorT* t business" << std::endl;
//                TensorT* t = set_coeff(tree,TensorT(k,k,k));
		TensorT* t = new TensorT();
//std::cout << "_compress: created TensorT* t" << std::endl;
                t = set_coeff(tree,TensorT(k,k,k));
//std::cout << "_compress: want to receive data on n = " << tree->n() << ", (" << tree->x() << "," <<
//	tree->y() << "," << tree->z() << ") from processor " << tree->rank() << std::endl;
//                comm()->Recv(t->ptr(), k3, tree->rank(), 2);
		MPI::COMM_WORLD.Recv(t->ptr(), k3, MPI_DOUBLE, tree->rank(), 2);
//                print("compress: received",t->normf(),"from",tree->rank());
            }
        } else {
//std::cout << "_compress: tree n = " << tree->n() << ", (" << tree->x() << "," <<
//	tree->y() << "," << tree->z() << ") is local" << std::endl;
            // Node is local.  May or may not have data.
            Slice *s = data->cdata->s;
            TensorT* t = coeff(tree);
//std::cout << "_compress: tensor t exists? " << t << std::endl;
            if (!t) t = set_coeff(tree,TensorT(2*k,2*k,2*k));
            FOREACH_CHILD(OctTreeT, tree,
                          if (isactive(child)) {
//std::cout << "_compress: child(" << i << "," << j << "," << k << ") is active" << std::endl;
                          TensorT* c = coeff(child);
//std::cout << "_compress: made c " << c << std::endl;
                              if (!c) throw "compress: yo! show me the data!";
//std::cout << "_compress: about to do crazy slicing of t += c thingy" << std::endl;
                              (*t)(s[i],s[j],s[k]) += (*c)(s[0],s[0],s[0]);
//std::cout << "_compress: did this crazy slicing of t += c thingy" << std::endl;
                              if (isremote(child)) {
//std::cout << "_compress: child(" << i << "," << j << "," << k << ") is remote" << std::endl;
                                  unset_coeff(child);
                              } else {
//std::cout << "_compress: child(" << i << "," << j << "," << k << ") is local" << std::endl;
                                  (*c)(s[0],s[0],s[0]) = T(0.0);
//std::cout << "_compress: did weird (*c) thingy" << std::endl;
                              }
                          }
			  else
			  {
//std::cout << "_compress: child(" << i << "," << j << "," << k << ") is not active" << std::endl;
			  }
                         );
            //set_coeff(tree,filter(*t));
//std::cout << "_compress: about to filter_inplace" << std::endl;
            filter_inplace(*t);
//std::cout << "_compress: just did filter_inplace" << std::endl;
        }
//std::cout << "_compress: at end of function" << std::endl;
    }


    template <typename T>
    void Function<T>::_reconstruct(OctTreeT* tree) {
        std::vector<Slice>& s0 = data->cdata->s0;
        Slice* s = data->cdata->s;
        TensorT& buf = data->cdata->work1;
        long k3 = k*k*k;

        if (isremote(tree)) {
            if (tree->islocalsubtreeparent()) {
                FOREACH_CHILD(OctTreeT, tree,
                              if (isactive(child)) {
                              //comm()->Recv(buf.ptr(),k3,tree->rank(),3);
		MPI::COMM_WORLD.Recv(buf.ptr(), k3, MPI_DOUBLE, tree->rank(), 3);
              //                                  print("reconstruct: received",buf.normf(),"from",tree->rank());
                                  (*coeff(child))(s0) = buf;
                              }
                             );
            } else {
                throw "_reconstruct: should not be here?";
            }
        } else {
            int nchild = 0;
            TensorT* t = coeff(tree);
            unfilter_inplace(*t);
            FOREACH_CHILD(OctTreeT, tree,
                          if (isactive(child)) {
                          nchild++;
                          if (islocal(child)) {
                                  (*coeff(child))(s0) = (*t)(s[i],s[j],s[k]);
                              } else {
                                  buf(s0) = (*t)(s[i],s[j],s[k]);
              //                                  print("reconstruct: sending",buf.normf(),"to",child->rank());
//                                  comm()->Send(buf.ptr(),k3,child->rank(),3);
				MPI::COMM_WORLD.Send(buf.ptr(), k3, MPI_DOUBLE, child->rank(), 3);
                              }
                          }
                         );
            if (nchild) unset_coeff(tree); // NEED TO KEEP IF WANT REDUNDANT TREE
        }

        FOREACH_CHILD(OctTreeT, tree,
                      if (isactive(child) && islocal(child)) _reconstruct(child););
    }

    /// Private: Refine an existing block of coefficients
    template <typename T>
    void Function<T>::_refine(OctTreeT* tree) {
        if (tree->islocalsubtreeparent() && isremote(tree)) {
            // This remote node is the parent of the local subtree.
            // Loop thru the local children and receive refinement
            // info from above.  The messages will come in order.
            FOREACH_CHILD(OctTreeT, tree,
                          bool dorefine;
                          comm()->Recv(dorefine, tree->rank(), 1);
              //                          print("   refine: got message",dorefine,"for",child->n(),child->x(),child->y(),child->z());
                          if (dorefine) {
                          set_active(child);
                              set_active(tree);
                              _project(child);
                          }
                          _refine(child);
                         );
        } else {
            TensorT* t = coeff(tree);
            if (t) {
                // Local node with data
                if (isremote(tree)) throw "refine: remote node with data?";
                TensorT d = filter(*t);
                d(data->cdata->s0) = 0.0;
                bool dorefine = (d.normf() > truncate_tol(data->thresh,tree->n()));
                if (dorefine)
//                    print("refine:",tree->n(),tree->x(),tree->y(),tree->z(),d.normf(),"dorefine =",dorefine);
                    if (dorefine) unset_coeff(tree);
                // First, send messages to remote children in order to get them working ASAP
                FOREACH_REMOTE_CHILD(OctTreeT, tree,
                                     comm()->Send(dorefine, child->rank(), 1);
                                     set_active(child);
                     //                                     print("sent",dorefine,"message to",child->n(),child->x(),child->y(),child->z());
                                    );
                // Next, handle local stuff
                FORIJK(OctTreeT* child = tree->child(i,j,k);
                       // If child does not exist, OK to refine locally
                       if (!child && dorefine) child = tree->insert_local_child(i,j,k);
                       if ( child && islocal(child)) {
                           if (dorefine) {
                                   _project(child);
                                   set_active(child);
                               }
                               _refine(child);
                           }
                      );
            } else {
                // This node does not have data and is not a remote subtree parent.
                // It must be either an empty local node or a remote leaf.
                // Again, first send msgs to remote nodes, then treat local guys.
                FOREACH_REMOTE_CHILD(OctTreeT, tree,
                                     comm()->Send(false, child->rank(), 1);
                     //                                     print("    sent false to",child->n(),child->x(),child->y(),child->z());
                                    );
                FOREACH_LOCAL_CHILD(OctTreeT, tree, _refine(child););
            }
        }
    }

    template <typename T>
    void Function<T>::_fcube(long n, long lx, long ly, long lz,
                             T (*f)(double,double,double), Tensor<T>& fcube) {
        double h = 1.0/two_to_power(n);
        double xlo = lx*h, ylo = ly*h, zlo = lz*h;

        int npt = data->cdata->npt;
        const Tensor<double>& quad_x = data->cdata->quad_x;

        for (int i=0; i<npt; i++) {
            double x = xlo + h*quad_x(i);
            for (int j=0; j<npt; j++) {
                double y = ylo + h*quad_x(j);
                for (int k=0; k<npt; k++) {
                    double z = zlo + h*quad_x(k);
                    fcube(i,j,k) = f(x,y,z);
                }
            }
        }
    }

    template <typename T>
    void Function<T>::_vfcube(long n, long lx, long ly, long lz,
                              void (*vf)(long l, const double*, const double*, const double*,
                                         T* RESTRICT),
                              Tensor<T>& fcube) {
        double h = 1.0/two_to_power(n);
        double xlo = lx*h, ylo = ly*h, zlo = lz*h;
        int npt = data->cdata->npt;
        const Tensor<double>& quad_x = data->cdata->quad_x;
#ifdef IBMXLC
        if (npt > 35) throw("hard dimension failure");
#endif
        double* x = new double[npt*npt*npt];
        double* y = new double[npt*npt*npt];
        double* z = new double[npt*npt*npt];

        int count = 0;
        for (int i=0; i<npt; i++) {
            double xx = xlo + h*quad_x(i);
            for (int j=0; j<npt; j++) {
                double yy = ylo + h*quad_x(j);
                for (int k=0; k<npt; k++) {
                    double zz = zlo + h*quad_x(k);
                    x[count] = xx;
                    y[count] = yy;
                    z[count] = zz;
                    count++;
                }
            }
        }

        vf(npt*npt*npt, x, y, z, fcube.ptr());

        delete[] x;
        delete[] y;
        delete[] z;
    }

    void balanceFunctionOctTree(SharedPtr< FunctionOctTree> trees)
    {
//	std::cout << "balanceFunctionOctTree: at very very beginning" << std::endl;
	int als = 0;
        std::vector< SharedPtr< OctTree< FunctionNode> > > treeList;

//	bool debug = true;
	bool debug = false;

	if (trees.get())
	    als = trees->getNalloc();

	if (debug)
	{
	    std::cout << "balanceFunctionOctTree: at very beginning, num allocated = " << als << std::endl;
	}

	treeList = trees->treeList();
	if (debug)
	{
	    std::cout << "balanceFunctionOctTree: just cleared trees" << std::endl;
	}

        serialLoadBalance(&treeList);

	if (debug)
	{
	    std::cout << "balanceFunctionOctTree: back from serialLoadBalance" << std::endl;
	}

	setRemoteActive(&treeList);
	trees->setTreeList(treeList);
	if (trees.get())
	    als = trees->getNalloc();
	if (debug)
	{
	    std::cout << "balanceFunctionOctTree: over and out, num allocated = " << als << std::endl;
	}
    }


    void setRemoteActive(std::vector<SharedPtr<OctTree<FunctionNode > > > *treeList)
    {
	std::vector<ActiveRootList> activeList, remoteParentList;
	std::vector<RootList> remoteChildList, parentRecvList;
	std::vector<bool> a;
	int tlen = treeList->size();
	Communicator comm;

//	bool debug = true;
	bool debug = false;

	if (debug)
	{
	    std::cout << "setRemoteActive: at beginning" << std::endl;
	}

	for (int i = 0; i < tlen; i++)
	{
	    makeRemoteParentList(&remoteParentList, (*treeList)[i]);
		
	    if (debug)
	    {
		std::cout << "setRemoteActive, beginning of i loop, i = " << i << std::endl;
	    }
	    if ((*treeList)[i]->isremote())
	    {
		if (debug)
		{
		    std::cout << "setRemoteActive, (*treeList)[" << i << "] is remote" << std::endl;
		}

		OctTree<FunctionNode> *p = new OctTree<FunctionNode> ();
		p = (*treeList)[i].get();
		RootList prl = RootList(p, p->rank(), p->rank());
		parentRecvList.push_back(prl);
//		ProcessID parentOwner = (*treeList)[i]->rank();
		ProcessID parentOwner = p->rank();
		if (debug)
		{
		    std::cout << "setRemoteActive, (*treeList)[" << i << "] has rank " << parentOwner
			 << std::endl;
		}
		FOREACH_LOCAL_CHILD(OctTree<FunctionNode >, p,
		    if (debug)
		    {
			std::cout << "setRemoteActive: about to get a vector from child (" << i <<
				"," << j << "," << k << ")" << std::endl;
		    }
		    a = child->data().getActiveList();
		    if (debug)
		    {
			std::cout << "setRemoteActive: got a vector" << std::endl;
		    }
		    activeList.push_back(ActiveRootList(child, parentOwner, a));
		    if (debug)
		    {
			std::cout << "setRemoteActive: pushed back child(" << i << "," << j <<
				"," << k << ") onto activeList" << std::endl;
		    }
		);
	    }
	    if (debug)
	    {
		std::cout << "setRemoteActive: before findRemoteChildrenList" << std::endl;
	    }
	    findRemoteChildrenList(&remoteChildList, (*treeList)[i]);
	    if (debug)
	    {
		std::cout << "setRemoteActive: after findRemoteChildrenList" << std::endl;
	    }
	}

	if (debug)
	{
	    std::cout << "setRemoteActive: after for loop, activeList size = " << activeList.size()
		<< std::endl;
	}

	int rplen = remoteParentList.size();
	for (int i = 0; i < rplen; i++)
	{
	    activeList.push_back(remoteParentList[i]);
	}

	// sort activeList by send, and take everybody who is sending to one processor, put 
	// on one vector, and send it
	sortBySend(&activeList);

	if (debug)
	{
	    int alen = activeList.size();
	    std::cout << "setRemoteActive: after sortBySend, size of activeList = " << alen 
		<< std::endl;
	    for (int i = 0; i < alen; i++)
	    {
		std::cout << "    send " << activeList[i].r.current_owner << ", n = " << 
			activeList[i].r.n << ", (" << activeList[i].r.x << "," << 
			activeList[i].r.y << "," << activeList[i].r.z << ")" << std::endl;
	    }
	}

	ProcessID np = comm.nproc();
	if (!activeList.empty())
	{
	    for (ProcessID i = 0; i < np; i++)
	    {
	    	std::vector<ActiveRootList> tmpList;
	    	bool keepon = true;
	    	while (keepon)
	    	{
		    if (activeList[0].r.current_owner == i)
		    {
		    	tmpList.push_back(activeList[0]);
		    	if (debug)
		    	{
			std::cout << "setRemoteActive: about to erase top of activeList" << std::endl;
		    	}
		    	activeList.erase(activeList.begin());
		    	if (debug)
		    	{
			std::cout << "setRemoteActive: erased top of activeList" << std::endl;
		    	}
		    }
		    else
		    {
		        keepon = false;
		    }
		    if (activeList.empty())
		    {
		        keepon = false;
		    }
	    	}
	    	if (!tmpList.empty())
	    	{
		    if (debug)
		    {
		    	std::cout << "setRemoteActive: about to send tmpList" << std::endl;
		    }
		    archive::MPIOutputArchive arsend(comm, i);
		    arsend & tmpList;
	    	}
	    }
	}

	int prl = parentRecvList.size();
	for (int i = 0; i < prl; i++)
	{
	    remoteChildList.push_back(parentRecvList[i]);
	}
	if (debug)
	{
	    std::cout << "setRemoteActive: about to sortBySend remoteChildList" << std::endl;
	    std::cout << "setRemoteActive: size of remoteChildList = " << remoteChildList.size()
		<< std::endl;
	}
	sortBySend(&remoteChildList);
	int k = 0;
	if (debug)
	{
	    std::cout << "setRemoteActive: about to receive list of active nodes" << std::endl;
	}
	if (!remoteChildList.empty())
	{
	    if (debug)
	    {
		std::cout << "setRemoteActive: remoteChildList is not empty; it's of size " 
			<< remoteChildList.size() << std::endl;
	    }

	    for (ProcessID i = 0; i < np; i++)
	    {
		if (k >= remoteChildList.size())
		{
		    break;
		}
	    	if (remoteChildList[k].current_owner == i)
	    	{
		    std::vector<ActiveRootList> tmpList;
		    archive::MPIInputArchive arrecv(comm, i);
		    if (debug)
		    {
			std::cout << "setRemoteActive: k = " << k << 
				"; about to receive tmpList from proc " << i << std::endl;
		    }
		    arrecv & tmpList;
		    if (debug)
		    {
			std::cout << "setRemoteActive: received tmpList from proc " << i << std::endl;
		    }
		    int len = tmpList.size();
		    for (int j = 0; j < len; j++)
		    {
		    	activeList.push_back(tmpList[j]);
		    }
		    while (remoteChildList[k].current_owner == i)
		    {
		    	k++;
		    }
	    	}
	    }

	    // now, sort activeList into depth-first traversal order
	    sort(activeList.begin(), activeList.end());
	    if (debug)
	    {
	    	std::cout << "setRemoteActive: sorted list of active nodes" << std::endl;
	    }
	
	    // then, depth-first traverse my list of trees, looking for activeChildrenList[k]
	    // in the list of trees, and filling in the appropriate "a" vector
	    int aclen = activeList.size();
	    for (int j = 0; j < aclen; j++)
	    {
		if (debug)
		{
		    std::cout << "setRemoteActive: at beginning of j loop, j = " << j << std::endl;
		}
	    	for (int i = 0; i < tlen; i++)
	    	{
		    if (debug)
		    {
			std::cout << "setRemoteActive: at beginning of inner i loop, i = " << i 
				<< std::endl;
		    }
	    	    OctTree<FunctionNode> *t = new OctTree<FunctionNode>();
	    	    t = (*treeList)[i]->find(activeList[j].r.n, activeList[j].r.x,
			activeList[j].r.y, activeList[j].r.z);
		    if ((t) && activeList[j].r.equals(t))
		    {
			if (debug)
			{
			    std::cout << "setRemoteActive: found tree n = " << t->n() << ", (" <<
				t->x() << "," << t->y() << "," << t->z() << ")" << std::endl;
			}
		    	// fill in the active a vector
			t->setData(new FunctionNode(activeList[j].activeList));
			if (debug)
			{
			    std::vector<bool> abool = t->data().getActiveList();
			    std::cout << "setRemoteActive: a = ";
			    for (int q = 0; q < abool.size(); q++)
			    {
			    	std::cout << abool[q] << " ";
			    }
			    std::cout << std::endl;
			    std::cout << "setRemoteActive: v.size() = " << t->data().getTensorListSize()
				<< std::endl;
			    std::cout << "setRemoteActive: set t->data" << std::endl;
			}
		    	break;
		    }
	    	}
	    }	
	}
	else
	{
	    if (debug)
	    {
	    	std::cout << "setRemoteActive: since remoteChildList is empty, I'm all done" 
			<< std::endl;
	    }
	}
    }


    void findRemoteChildrenList(std::vector<RootList> *childList, SharedPtr<OctTree<FunctionNode > > tree) 
    {
//	bool debug = true;
	bool debug = false;

	FOREACH_CHILD(OctTree<FunctionNode>, tree,
	    if (child->islocal())
	    {
		if (debug)
		{
		    std::cout << "findRemoteChildrenList: child n = " << child->n() << 
			", ( " << child->x() << "," << child->y() << "," << child->z() << 
			") is local" << std::endl;
		}
		findRemoteChildrenList(childList, tree->childPtr(i,j,k));
	    }
	    else
	    {
		childList->push_back(RootList(child, child->rank(), child->rank()));
		if (debug)
		{
		    std::cout << "findRemoteChildrenList: added n = " << child->n() << 
			", ( " << child->x() << "," << child->y() << "," << child->z() << ")"
			<< std::endl;
		}
	    }
	);
    }

    void makeRemoteParentList(std::vector<ActiveRootList> *parentList, SharedPtr<OctTree<FunctionNode> > tree)
    {
	std::vector<RootList> childList, tmpList;
	std::vector<bool> a;

//	bool debug = true;
	bool debug = false;

	findRemoteChildrenList(&childList, tree);

	if (debug)
	{
	    std::cout << "makeRemoteParentList: just returned from findRemoteChildrenList with list len =" 
		<< childList.size() << std::endl;
	}
	int tlen = childList.size();
	if (tlen > 0)
	{
	    tmpList.push_back(childList[0]);
	    if (debug)
	    {
	    	std::cout << "makeRemoteParentList: pushed back child 0, current_owner = " <<
			childList[0].current_owner << std::endl;
	    }
	}
	for (int i = 1; i < tlen; i++)
	{
	    // if they are siblings and belong to the same processor, do nothing
	    if ((childList[i].n == childList[i-1].n) && (childList[i].x/2 == childList[i-1].x/2) && 
		(childList[i].y/2 == childList[i-1].y/2) && (childList[i].z/2 == childList[i-1].z/2) 
		&& (childList[i].current_owner == childList[i-1].current_owner))
	    {
		if (debug)
		{
		    std::cout << "makeRemoteParentList: child[" << i << "] is a duplicate" << std::endl;
		}
	    }
	    else
	    {
		tmpList.push_back(childList[i]);
		if (debug)
		{
		    std::cout << "makeRemoteParentList: pushed back child[" << i << "], current_owner = " 
			<< childList[i].current_owner << ", tmpList.size() = " << tmpList.size() 
			<< std::endl;
		}
	    }
	}
	if (debug)
	{
	    std::cout << "makeRemoteParentList: just finished pushing non-duplicates onto tmpList" << 
		std::endl;
	}
	tlen = tmpList.size();
	for (int i = 0; i < tlen; i++)
	{
	    OctTree<FunctionNode> *t = new OctTree<FunctionNode>();
	    t = tree->findDown(tmpList[i].n-1, tmpList[i].x/2, tmpList[i].y/2, tmpList[i].z/2);
	    if (t)
	    {
		a = t->data().getActiveList();
		parentList->push_back(ActiveRootList(t, tmpList[i].current_owner, a));
	    }
	    else
	    {
		std::cout << "makeRemoteParentList: error: couldn't find parent!" << std::endl;
	    }
	}
	if (debug)
	{
	    std::cout << "makeRemoteParentList: just about to return, with list longer by " <<
		parentList->size() << std::endl;
	}
    }

    // Explicit instantiations for double and complex<double>

    template class Function<double>;
    template class Function< std::complex<double> >;
    template class FunctionData<double>;
    template class FunctionData< std::complex<double> >;
    template class FunctionCommonData<double>;

    template void serialLoadBalance<double>(std::vector< SharedPtr<OctTree<double> > >*);
    template void serialLoadBalance<std::complex<double> >(std::vector< SharedPtr<OctTree<std::complex<double> > > >*);
}

