#ifndef LION_FOUNDATION_TYPE_TRAITS_H
#define LION_FOUNDATION_TYPE_TRAITS_H
#pragma once


#include <type_traits>
#include <vector>
#include <array>

#include "cppad/cppad.hpp"


//! This file defines additional type traits to be used for SFINAE
//! - is_vector<T>
//! - is_matrix<T>

//! Class is_vector: contains a bool value which is only true for 
//!    - std::vector
//!    - std::array
//! Is not true for vector of vectors
template<typename T>
struct is_vector
{
    static constexpr auto value = false;
};

template<typename T>
struct is_vector<std::vector<T>> 
{
    using type = std::vector<T>;
    static constexpr auto value = true;
};

template<typename T>
struct is_vector<std::vector<std::vector<T>>>
{
    static constexpr auto value = false;
};

template<typename T, std::size_t N>
struct is_vector<std::array<T,N>> 
{
    using type = std::array<T,N>;
    static constexpr auto value = true;
};

template<typename T, std::size_t N, std::size_t M>
struct is_vector<std::array<std::array<T,M>,N>>
{
    static constexpr auto value = false;
};

//! Class is_matrix: contains a bool value which is only true for 
//!    - std::vector<std::vector>
//!    - std::array<std::array>
template<typename T>
struct is_matrix
{
    static constexpr auto value = false;
};

template<typename T>
struct is_matrix<std::vector<std::vector<T>>> 
{
    using type = std::vector<std::vector<T>>;
    static constexpr auto value = true;
};

template<typename T, std::size_t N, std::size_t M>
struct is_matrix<std::array<std::array<T,M>,N>> 
{
    using type = std::array<std::array<T,M>,N>;
    static constexpr auto value = true;
};

template<typename T>
struct is_matrix<std::vector<std::vector<std::vector<T>>>>
{
    static constexpr auto value = false;
};

template<typename T, std::size_t N, std::size_t M, std::size_t L>
struct is_matrix<std::array<std::array<std::array<T,L>,M>,N>>
{
    static constexpr auto value = false;
};

template<typename U, typename V>
struct combine_types
{};

template<typename T>
struct combine_types<T,T>
{
    using type = T;
};

template<typename T>
struct combine_types<CppAD::AD<T>,T>
{
    using type = CppAD::AD<T>;
};

template<typename T>
struct combine_types<T,CppAD::AD<T>>
{
    using type = CppAD::AD<T>;
};

template<typename T, size_t N>
struct combine_types<T, std::array<T, N>>
{
    using type = std::array<T, N>;
};

template<typename T, typename U>
using combine_types_t = typename combine_types<T, U>::type;


/// Type-trait style function to eliminate copy constructors (via std::enable_if) from a class constructor template that
/// only has variadic arguments (i.e., "template<typename... Args> ClassName(Args&&...) { ... }").
//!
//! @tparam Class The constructed class.
//! @tparam FirstArg The first parameter of the class's ctor.
//! @tparam RestOfArgs The rest of parameters of the class's ctor.
//! @tparam std::enable_if_t<...> To declare that a ctor with this signature can never be a copy ctor.
//! @return True.
//!
template<class Class, typename FirstArg, typename... RestOfArgs,
    typename std::enable_if_t<sizeof...(RestOfArgs) != 0u || !std::is_same_v<Class, std::decay_t<FirstArg> > >* = nullptr>
constexpr auto is_not_copy_ctor() { return true; }

/// Type-trait style function to eliminate copy constructors (via std::enable_if) from a class constructor template that
/// only has variadic arguments (i.e., "template<typename... Args> Class(Args&&...) { ... }").
//!
//! @tparam Class The constructed class.
//! @tparam Arg The single parameter of the class's ctor.
//! @tparam std::enable_if_t<...> To declare that a ctor with this signature must be a copy ctor.
//! @return False.
//!
template<class Class, typename SingleArg,
    typename std::enable_if_t<std::is_same_v<Class, std::decay_t<SingleArg> > >* = nullptr>
constexpr auto is_not_copy_ctor() { return false; }

#ifdef _MSC_VER
    // In VS2017-2019 the "is_not_copy_ctor()" constexpr function does not
    // work, we'll use a macro until the compiler catches up!
#define WINDOWS_IS_NOT_COPY_CTOR(CLASSNAME)                                                                                                         \
    template<typename... Args,                                                                                                                      \
             typename std::enable_if_t<(sizeof...(Args) > 1u) ||                                                                                    \
                                       !std::is_same_v<std::decay_t<std::tuple_element_t<0u, std::tuple<Args..., void> > >, CLASSNAME> >* = nullptr>
#endif


/// Helper template to declare a type that may be a dummy,
/// empty type ("nothing") or an actual type, based on the
/// value of a constexpr boolean flag. This is the primary
/// template that gets used when the mentioned flag is false.
//!
//! @tparam ShouldBeType The constexpr boolean flag indicating
//!   if the declared type should be an actual type (if true)
//!   or an empty, dummy type (if false).
//! @tparam (empty) To be used by its specialization.
//!
template<bool ShouldBeType, typename = void>
struct type_or_nothing
{
  /// Struct representing an empty, dummy type.
  struct nothing
  {
    /// Class ctor (does nothing, admits anything).
    /**
     * @tparam Args Parameter pack that won't get used.
     */
    template<typename... Args>
    constexpr nothing(Args&&...) {}
  };

  /// Member type pointing to the empty, dummy type.
  using type = nothing;
};

/// Specialization of the above template that gets used
/// when the constexpr boolean flag is true.
//!
//! @tparam T The actual type we desire to declare.
//!
template<typename T>
struct type_or_nothing<true, T>
{
  /// Member type pointing to the actual type we desire
  /// to declare.
  using type = T;
};

template<bool ShouldBeType, typename T>
using type_or_nothing_t = typename type_or_nothing<ShouldBeType, T>::type;


//
// std::true_type if the tested type is an
// specialization of a template class.
//

template<typename T, template<typename...> class>
struct is_specialization : std::false_type {};

template<template<typename...> class T, typename... Ts>
struct is_specialization<T<Ts...>, T> : std::true_type {};

template<typename T, template<typename...> class C>
constexpr bool is_specialization_v = is_specialization<T, C>::value;


//
// std::true_type if the tested type is
// an std::array.
//

template<typename Arr, typename = void>
struct is_std_array : std::false_type {};

template<typename T, std::size_t N>
struct is_std_array<std::array<T, N> > : std::true_type {};

template<typename Arr>
constexpr bool is_std_array_v = is_std_array<Arr>::value;


//
// Resolves to the type that we get when applying the
// bracket operator (with an input index) to an object
// (e.g., "decltype(v[i])").
//

template<class C, typename IndexType = std::size_t>
using typeof_bracket_operator = std::decay_t<decltype(std::declval<C &>()[std::declval<IndexType>()])>;


//
// Constexpr version of the ternary operator.
//

template<bool condition, typename T, typename F>
constexpr auto&& ternary_constexpr(T &&val_if_true, F &&val_if_false)
{
    //
    // Returns "val_if_true" when "condition" is true,
    // and "val_if_false" otherwise.
    //

    if constexpr (condition) {
        return std::forward<T>(val_if_true);
    }
    else {
        return std::forward<F>(val_if_false);
    }
}

#endif
