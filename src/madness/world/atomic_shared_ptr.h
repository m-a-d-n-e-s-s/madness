/**
 \file atomic_shared_ptr.h
 \brief Defines \c atomic_shared_ptr — a std::shared_ptr-compatible
        wrapper backed by atomic storage.
 \ingroup world
*/

#ifndef MADNESS_WORLD_ATOMIC_SHARED_PTR_H__INCLUDED
#define MADNESS_WORLD_ATOMIC_SHARED_PTR_H__INCLUDED

#include <atomic>
#include <memory>
#include <madness/madness_config.h>

namespace madness {

/// A \c std::shared_ptr-like handle backed by atomic storage.
///
/// Uses \c std::atomic<std::shared_ptr<T>> when available (C++20), otherwise
/// falls back to the pre-C++20 atomic free functions on \c std::shared_ptr.
/// The public API mirrors \c std::shared_ptr (copy construction/assignment,
/// \c operator->, \c operator bool, implicit conversion, \c reset,
/// \c use_count, equality comparison) plus a single extension, \c exchange(),
/// needed for the move-optimisation in \c Future::get()&&.
///
/// \note User-declared copy constructor and copy-assignment operator suppress
/// the implicit move constructor and move-assignment operator (and
/// \c std::atomic is itself non-movable in the C++20 backend), so move
/// operations silently fall back to copy.
///
/// \note All loads are \c memory_order_relaxed, relying on the caller to provide necessary
/// ordering guarantees. The exception is \c exchange(), which defaults to \c memory_order_seq_cst to match the pre-C++20
/// \c std::atomic_exchange behaviour it replaces.
template<typename T>
class atomic_shared_ptr {
#if defined(__cpp_lib_atomic_shared_ptr)
    std::atomic<std::shared_ptr<T>> ptr_;

    std::shared_ptr<T> do_load() const noexcept {
        return ptr_.load(std::memory_order_relaxed);
    }
    void do_store(std::shared_ptr<T> p) noexcept {
        ptr_.store(std::move(p), std::memory_order_relaxed);
    }
#else
    std::shared_ptr<T> ptr_;

    // The atomic_*_explicit free functions on std::shared_ptr are deprecated
    // in C++20 (superseded by std::atomic<shared_ptr<T>>). We reach this
    // branch only when the library lacks __cpp_lib_atomic_shared_ptr, which
    // can happen even on a C++20 compiler (old stdlib), so suppress the
    // warning specifically for C++20 builds where it would fire.
#if __cplusplus >= 202002L
    MADNESS_PRAGMA_CLANG(diagnostic push)
    MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
    MADNESS_PRAGMA_GCC(diagnostic push)
    MADNESS_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")
#endif
    // Payload ordering supplied by FutureImpl's Spinlock; relaxed for identity.
    std::shared_ptr<T> do_load() const noexcept {
        return std::atomic_load_explicit(&ptr_, std::memory_order_relaxed);
    }
    void do_store(std::shared_ptr<T> p) noexcept {
        std::atomic_store_explicit(&ptr_, std::move(p), std::memory_order_relaxed);
    }
#if __cplusplus >= 202002L
    MADNESS_PRAGMA_CLANG(diagnostic pop)
    MADNESS_PRAGMA_GCC(diagnostic pop)
#endif
#endif

public:
    constexpr atomic_shared_ptr() noexcept = default;

    atomic_shared_ptr(std::shared_ptr<T> p) noexcept : ptr_(std::move(p)) {}

    /// Copy-constructs by atomically loading the other pointer.
    atomic_shared_ptr(const atomic_shared_ptr& other) noexcept : ptr_(other.do_load()) {}

    /// Copy-assigns by atomically loading from \p other and storing.
    atomic_shared_ptr& operator=(const atomic_shared_ptr& other) noexcept {
        if (this != &other) do_store(other.do_load());
        return *this;
    }

    /// Assigns a new shared_ptr value atomically.
    atomic_shared_ptr& operator=(std::shared_ptr<T> p) noexcept {
        do_store(std::move(p));
        return *this;
    }

    /// Replaces the stored pointer, adopting \p p (may be null).
    void reset(T* p) { do_store(std::shared_ptr<T>(p)); }

    /// Replaces the stored pointer with an empty shared pointer.
    void reset() { do_store(std::shared_ptr<T>()); }

    explicit operator bool() const noexcept { return bool(do_load()); }

    /// Dereferences the stored pointer.
    /// Safe because \c do_load() returns a \c shared_ptr temporary whose
    /// lifetime extends to the end of the full expression containing the
    /// \c -> call, keeping the referent alive across the dereference.
    T* operator->() noexcept { return do_load().get(); }

    /// Dereferences the stored pointer, const overload.
    /// Safe because \c do_load() returns a \c shared_ptr temporary whose
    /// lifetime extends to the end of the full expression containing the
    /// \c -> call, keeping the referent alive across the dereference.
    const T* operator->() const noexcept { return do_load().get(); }

    /// Implicit conversion to \c shared_ptr, used for lifetime extension and
    /// APIs that require a \c shared_ptr (e.g. \c RemoteReference).
    operator std::shared_ptr<T>() const noexcept { return do_load(); }

    long use_count() const noexcept { return do_load().use_count(); }

    bool operator==(const atomic_shared_ptr& r) const noexcept { return do_load() == r.do_load(); }
    bool operator!=(const atomic_shared_ptr& r) const noexcept { return do_load() != r.do_load(); }

    /// Atomically replaces the stored pointer with \p p and returns the old
    /// value. This is the only operation that differs from \c std::shared_ptr.
    /// Defaults to \c memory_order_seq_cst to match the pre-C++20
    /// \c std::atomic_exchange behaviour this replaces.
    std::shared_ptr<T> exchange(std::shared_ptr<T> p = {},
                                std::memory_order order = std::memory_order_seq_cst) noexcept {
#if defined(__cpp_lib_atomic_shared_ptr)
        return ptr_.exchange(std::move(p), order);
#else
#if __cplusplus >= 202002L
        MADNESS_PRAGMA_CLANG(diagnostic push)
        MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
        MADNESS_PRAGMA_GCC(diagnostic push)
        MADNESS_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")
#endif
        return std::atomic_exchange_explicit(&ptr_, std::move(p), order);
#if __cplusplus >= 202002L
        MADNESS_PRAGMA_CLANG(diagnostic pop)
        MADNESS_PRAGMA_GCC(diagnostic pop)
#endif
#endif
    }
};

} // namespace madness

#endif // MADNESS_WORLD_ATOMIC_SHARED_PTR_H__INCLUDED
