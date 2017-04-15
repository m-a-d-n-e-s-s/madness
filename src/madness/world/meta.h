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

template <template <typename...> class MetaFn, typename T, typename... P>
struct drop_last_param_and_apply_impl;

template <template <typename...> class MetaFn, typename... P1, typename T,
          typename... P2>
struct drop_last_param_and_apply_impl<MetaFn, typelist<P1...>, T, P2...> {
  using type =
      typename drop_last_param_and_apply_impl<MetaFn, typelist<P1..., T>,
                                              P2...>::type;
};

template <template <typename...> class MetaFn, typename... P1, typename T,
          typename L>
struct drop_last_param_and_apply_impl<MetaFn, typelist<P1...>, T, L> {
  using type = MetaFn<P1..., T>;
};

template <template <typename...> class MetaFn, typename... P>
struct drop_last_param_and_apply {
  using type =
      typename drop_last_param_and_apply_impl<MetaFn, typelist<>, P...>::type;
};

///////////////////////////////////////////////////////////////////////////////

template <template <typename...> class MetaFn, typename T, typename... P>
struct drop_last_param_and_apply_callable_impl;

template <template <typename...> class MetaFn, typename Callable,
          typename... P1, typename T, typename... P2>
struct drop_last_param_and_apply_callable_impl<MetaFn, Callable,
                                               typelist<P1...>, T, P2...> {
  using type = typename drop_last_param_and_apply_callable_impl<
      MetaFn, Callable, typelist<P1..., T>, P2...>::type;
};

template <template <typename...> class MetaFn, typename Callable,
          typename... P1, typename T, typename L>
struct drop_last_param_and_apply_callable_impl<MetaFn, Callable,
                                               typelist<P1...>, T, L> {
  using type = MetaFn<Callable(P1..., T)>;
};

template <template <typename...> class MetaFn, typename Callable, typename... P>
struct drop_last_param_and_apply_callable {
  using type =
      typename drop_last_param_and_apply_callable_impl<MetaFn, Callable,
                                                       typelist<>, P...>::type;
};

}  // namespace meta
}  // namespace madness

#endif /* SRC_MADNESS_WORLD_META_H_ */
