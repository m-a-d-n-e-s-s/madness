#ifndef BASETNESOR_H
#define BASETENSOR_H

/// \file basetensor.h
/// \brief Declares BaseTensor

#include <tensor/slice.h>
#include <tensor/tensor_macros.h>
    
namespace madness {
    /// The base class for tensors defines generic capabilities.
    
    /// The base class manages the size, dimension and
    /// stride information, and provides operations to manipulate
    /// them.
    ///
    /// It also provides methods for type-safe operation on tensors using
    /// just the base class pointers. This interface is primarily useful
    /// only to the interface to Python, since Python is largely neutral
    /// to (and ignorant of) the type.  These are still being
    /// re-implemented after the big clean up.
    ///
    /// Since the base tensor class is virtual, you cannot have an
    /// instance of it.  Thus, in addition to methods that return information
    /// or perform checks, there are two types of base tensor
    /// operations.  
    /// - Inplace operations change \c *this , and return \c void .
    /// - Operations that must return a new tensor return a pointer to a tensor
    ///   allocated with \c new on the heap.  The caller is responsible for
    ///   eventually freeing the memory using \c delete .
    class BaseTensor {
    private:
#ifdef TENSOR_INSTANCE_COUNT
        static long instance_count;	///< For debug, count total# instances 
                                /// ... initialized in basetensor.cc ... not thread safe!
#endif
        
    protected:
        void set_dims_and_size(long nd, const long d[]);
        
    public:
        long size;			///< Number of elements in the tensor
        long id; 			///< Id from TensorTypeData<T> in type_data.h
        long ndim;			///< Number of dimensions (-1=invalid; 0=scalar; >0=tensor)
        long dim[TENSOR_MAXDIM];	///< Size of each dimension
        long stride[TENSOR_MAXDIM];   ///< Increment between elements in each dimension
        
#ifdef TENSOR_INSTANCE_COUNT
        BaseTensor() {instance_count++;};
        virtual ~BaseTensor() {instance_count--;};
#else
        BaseTensor() {};
        virtual ~BaseTensor() {};
#endif
        virtual BaseTensor* new_shallow_copy_base() const = 0;
        virtual BaseTensor* new_deep_copy_base() const = 0;
        virtual BaseTensor* slice_base(const std::vector<Slice>& s) const = 0;
        
        /// Returns the count of all current instances of tensors & slice tensors of all types.
#ifdef TENSOR_INSTANCE_COUNT
        static inline int get_instance_count() {return instance_count;};
#else
        static inline int get_instance_count() {return 0;};
#endif
        
        bool conforms(const BaseTensor *t) const;
        inline bool iscontiguous() const {
            /// Returns true if the tensor refers to contiguous memory locations.
            long size = 1;
            for (long i=ndim-1; i>=0; i--) {
                if (stride[i] != size) return false;
                size *= dim[i];
            }
            return true;
        };
        void reshape_inplace_base(const std::vector<long>& d);
        BaseTensor* reshape_base(const std::vector<long>& d) const;
        void flat_inplace_base();
        BaseTensor* flat_base() const;
        void splitdim_inplace_base(long i, long dimi0, long dimi1);
        BaseTensor* splitdim_base(long i, long dimi0, long dimi1) const;
        void fusedim_inplace_base(long i);
        BaseTensor* fusedim_base(long i) const;
        void swapdim_inplace_base(long i, long j);
        BaseTensor* swapdim_base(long idim, long jdim) const;
        void cycledim_inplace_base(long shift, long start, long end);
        BaseTensor* cycledim_base(long shift, long start = 0, long end = -1) const;
        void mapdim_inplace_base(const std::vector<long>& map);
        BaseTensor* mapdim_base(const std::vector<long>& map) const;
    };
    
}

#endif
