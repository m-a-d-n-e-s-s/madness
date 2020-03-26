/*
 * meta.h
 *
 *  Created on: Apr 15, 2017
 *      Author: evaleev
 */

#ifndef SRC_MADNESS_WORLD_META_H_
#define SRC_MADNESS_WORLD_META_H_

#include <type_traits>

namespace madness {
namespace meta {  // to make it easier importing another MP library

///////////////////////////////////////////////////////////////////////////////

template <typename... P>
struct typelist {};

///////////////////////////////////////////////////////////////////////////////

template <typename... Ts>
struct last_type {};

template <typename T0>
struct last_type<T0> {
  typedef T0 type;
};

template <typename T0, typename T1, typename... Ts>
struct last_type<T0, T1, Ts...> : last_type<T1, Ts...> {};

///////////////////////////////////////////////////////////////////////////////

template <template <typename...> class MetaFn, typename CurrentTypelist,
          typename... RestOfTypes>
struct drop_last_arg_and_apply_impl;

template <template <typename...> class MetaFn, typename... UpToT, typename T,
          typename... Rest>
struct drop_last_arg_and_apply_impl<MetaFn, typelist<UpToT...>, T, Rest...> {
  using type =
      typename drop_last_arg_and_apply_impl<MetaFn, typelist<UpToT..., T>,
                                            Rest...>::type;
};

template <template <typename...> class MetaFn, typename... UpToLast,
          typename Last>
struct drop_last_arg_and_apply_impl<MetaFn, typelist<UpToLast...>, Last> {
  using type = MetaFn<UpToLast...>;
};

template <template <typename...> class MetaFn, typename... Args>
struct drop_last_arg_and_apply {
  using type =
      typename drop_last_arg_and_apply_impl<MetaFn, typelist<>, Args...>::type;
};

///////////////////////////////////////////////////////////////////////////////

template <template <typename...> class MetaFn, typename Callable,
          typename CurrentTypelist, typename... RestOfTypes>
struct drop_last_arg_and_apply_callable_impl;

template <template <typename...> class MetaFn, typename Callable,
          typename... UpToT, typename T, typename... Rest>
struct drop_last_arg_and_apply_callable_impl<MetaFn, Callable,
                                             typelist<UpToT...>, T, Rest...> {
  using type = typename drop_last_arg_and_apply_callable_impl<
      MetaFn, Callable, typelist<UpToT..., T>, Rest...>::type;
};

template <template <typename...> class MetaFn, typename Callable,
          typename... UpToLast, typename Last>
struct drop_last_arg_and_apply_callable_impl<MetaFn, Callable,
                                             typelist<UpToLast...>, Last> {
  using type = MetaFn<Callable(UpToLast...)>;
};

template <template <typename...> class MetaFn, typename Callable,
          typename... Args>
struct drop_last_arg_and_apply_callable {
  using type =
      typename drop_last_arg_and_apply_callable_impl<MetaFn, Callable,
                                                     typelist<>, Args...>::type;
};

///////////////////////////////////////////////////////////////////////////////

// nonesuch struct from Library Fundamentals V2, source from https://en.cppreference.com/w/cpp/experimental/nonesuch

struct nonesuch {
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

// is_detected family from Library Fundamentals V2, source from https://en.cppreference.com/w/cpp/experimental/is_detected

namespace detail {

template <class Default, class Enabler, template <class...> class Op,
          class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

template <class Default, template <class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type = Op<Args...>;
};

} // namespace detail

template <template<class...> class Op, class... Args>
using is_detected = typename detail::detector<nonesuch, void, Op, Args...>::value_t;

template <template<class...> class Op, class... Args>
using detected_t = typename detail::detector<nonesuch, void, Op, Args...>::type;

template <class Default, template<class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

template< template<class...> class Op, class... Args >
constexpr bool is_detected_v = is_detected<Op, Args...>::value;

template< class Default, template<class...> class Op, class... Args >
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

template <class Expected, template<class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

template <class Expected, template<class...> class Op, class... Args>
constexpr bool is_detected_exact_v = is_detected_exact<Expected, Op, Args...>::value;

template <class To, template<class...> class Op, class... Args>
using is_detected_convertible = std::is_convertible<detected_t<Op, Args...>, To>;

template <class To, template<class...> class Op, class... Args>
constexpr bool is_detected_convertible_v = is_detected_convertible<To, Op, Args...>::value;

}  // namespace meta
}  // namespace madness

#endif /* SRC_MADNESS_WORLD_META_H_ */
