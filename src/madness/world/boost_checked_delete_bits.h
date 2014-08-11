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

        /// Checked pointer delete functor

        /// This functor is used to delete a pointer. It ensures that the
        /// pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        template <typename T>
        struct CheckedDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* p) const {
                checked_delete(p);
            }
        };

        /// Checked array pointer delete functor

        /// This functor is used to delete an array pointer. It ensures that the
        /// pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        template <typename T>
        struct CheckedArrayDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* a) const {
                checked_array_delete(a);
            }
        };

        /// Deleter to free memory for a shared_ptr using free()

        /// Checks the pointer to make sure it is a complete type, you will get
        /// a compiler error if it is not.
        template <typename T>
        struct CheckedFree {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* p) const {
                checked_free(p);
            }
        };

        /// Use this deleter with shared_ptr to do nothing for the pointer cleanup
        template <typename T>
        struct NoDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T*) const { }
        };

    } // namespace detail
} // namespace madness
#endif
