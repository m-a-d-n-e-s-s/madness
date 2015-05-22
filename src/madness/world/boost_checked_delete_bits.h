#ifndef MADNESS_BOOST_CHECKED_DELETE_BITS_H
#define MADNESS_BOOST_CHECKED_DELETE_BITS_H
namespace madness {
    namespace detail {

        // These checked delete and deleters are copied from Boost.
        // They ensure that compilers issue warnings if T is an incomplete type.

        /// Checked pointer delete function

        /// This function ensures that the pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        /// \param p The pointer to be deleted.
        /// \todo This can be replaced by std::default_delete<T>?
        template <typename T>
        inline void checked_delete(T* p) {
            // intentionally complex - simplification causes regressions
            typedef char type_must_be_complete[sizeof(T) ? 1 : -1];
            (void)sizeof(type_must_be_complete);
            delete p;
        }

        /// Checked array pointer delete function

        /// This function ensures that the pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        /// \param a The array pointer to be deleted.
        /// \todo This can be replaced by std::default_delete<T[]>?
        template <typename T>
        inline void checked_array_delete(T* a) {
            typedef char type_must_be_complete[sizeof(T) ? 1 : -1];
            (void)sizeof(type_must_be_complete);
            delete[] a;
        }

        /// Function to free memory for a shared_ptr using free()

        /// Checks the pointer to make sure it is a complete type, you will get
        /// a compiler error if it is not.
        template <typename T>
        inline void checked_free(T* p) {
            typedef char type_must_be_complete[sizeof(T) ? 1 : -1];
            (void)sizeof(type_must_be_complete);
            free(p);
        }

        /// Use this function with shared_ptr to do nothing for the pointer cleanup
        template <typename T>
        inline void no_delete(T*) { }

        inline void no_delete(void*) { }

    } // namespace detail
} // namespace madness
#endif
