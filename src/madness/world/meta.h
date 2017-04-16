/*
 * meta.h
 *
 *  Created on: Apr 15, 2017
 *      Author: evaleev
 */

#ifndef SRC_MADNESS_WORLD_META_H_
#define SRC_MADNESS_WORLD_META_H_

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

}  // namespace meta
}  // namespace madness

#endif /* SRC_MADNESS_WORLD_META_H_ */
