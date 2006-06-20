#include <iostream>
using std::cout;
using std::endl;

#include <complex>
#include <octtree/octtree.h>
#include <mra/mra.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tasks/sav.h>
#include <tasks/tasks.h>
#include <misc/print.h>


/// \file mra/taskstuff.cc
/// \brief Implements Function stuff that uses tasks


namespace madness {
	
#define FOREACH_ACTIVE_CHILD(OctTreeT, t, expr) \
    do { \
      for (int i=0; i<2; i++) \
        for (int j=0; j<2; j++) \
	  for (int k=0; k<2; k++) { \
            OctTreeT *child = ((t)->child(i,j,k)); \
            if (child && isactive(child)) {expr} \
          } \
    } while(0)
	
    /// A task tailored to the needs of compress2
	template <typename T>
	class TaskCompress : public TaskInterface {
	public:
		typedef Tensor<T> TensorT;
		typedef SAV<TensorT> ArgT;
		typedef OctTree<FunctionNode> OctTreeT;
		typedef Function<T> FunctionT;
	private:
	    FunctionT* that;
	    OctTreeT* tree;
	    ArgT outarg;
	    ArgT inarg[2][2][2];
	public:
		TaskCompress(FunctionT* that, OctTreeT* tree, ArgT in[2][2][2], ArgT& out) 
		: that(that)
		, tree(tree)
		, outarg(out) 
		{
			FORIJK(inarg[i][j][k] = in[i][j][k];);
		};
		
		bool probe() const {
			FOREACH_CHILD(OctTreeT, tree, 
						  if (that->isactive(child) && !inarg[i][j][k].probe()) 
						  	  return false;);
			return true;
		};
		
		void run() {that->_compress2op(tree,inarg,outarg);};
		
		virtual ~TaskCompress() {};
	};
	
	static inline int taghash(Level n, Translation x, Translation y, Translation z) {
		return 1;
		int n4  = (n&15)<<27;
		int x9 = (x&511)<<18;
		int y9 = (y&511)<<9;
		int z9 = (z&511);
		return n4|x9|y9|z9;
	}; 
		
	/// Called by consumer to make a variable that will be set by producer
	template <typename T>
	typename Function<T>::ArgT Function<T>::input_arg(const OctTreeT* consumer, 
													  const OctTreeT* producer) {
		int tag = taghash(producer->n(),producer->x(),producer->y(),producer->z());
		if (consumer->islocal() && producer->islocal()) 
			return ArgT();
		else if (consumer->islocal() && producer->isremote())
			return ArgT(producer->rank(), tag, true, this->data->cdata->vk);
		else if (consumer->isremote() && producer->islocal())
			return ArgT(consumer->rank(), tag, false, this->data->cdata->vk);
		else
			throw "input_arg: should not happen?";
	};


	/// Compress function (scaling function to wavelet)

    /// Communication streams up the tree.
    /// Returns self for chaining.
    template <typename T>
    Function<T>& Function<T>::compress2() {
        if (!data->compressed) {
            if (isactive(tree())) {
            	ArgT dummy;
            	_compress2(tree(),dummy);
            	globalq.wait();
            }
            data->compressed = true;
        }
        return *this;
    };
	    
	template <typename T>
    void Function<T>::_compress2(OctTreeT* tree, Function<T>::ArgT& parent) {
    	ArgT args[2][2][2];
        FOREACH_ACTIVE_CHILD(OctTreeT, tree,
                           	 args[i][j][k] = this->input_arg(tree,child);
                      		 if (child->islocal()) _compress2(child,args[i][j][k]););
 	
 	    if (tree->islocal()) {
 	    	print(comm()->rank(),"adding task",tree->n(),tree->x(),tree->y(),tree->z());
 	    	globalq.add(new TaskCompress<T>(this,tree,args,parent));
 	    }
    };
    
    template <typename T>
    void Function<T>::_compress2op(OctTreeT* tree, Function<T>::ArgT args[2][2][2], Function<T>::ArgT& parent) {
 	    print(comm()->rank(),"executing task",tree->n(),tree->x(),tree->y(),tree->z());
        Slice *s = data->cdata->s;
        if (!coeff(tree)) set_coeff(tree,TensorT(2*k,2*k,2*k));
        TensorT& t = *coeff(tree);
        FOREACH_ACTIVE_CHILD(OctTreeT, tree,t(s[i],s[j],s[k]) += args[i][j][k].get(););
        filter_inplace(t);
        parent.set(madness::copy(t(data->cdata->s0)));
        if (tree->n() > 0) t(data->cdata->s0)=0.0;
    };

    template class Function<double>;
    template class Function< std::complex<double> >;
    template class FunctionData<double>;
    template class FunctionData< std::complex<double> >;
    template class FunctionCommonData<double>;
    template class TaskCompress<double>;
    template class TaskCompress< std::complex<double> >;
    
}
