#ifndef __TYPE_TRAITS_H__
#define __TYPE_TRAITS_H__

#include <type_traits>
#include <vector>
#include <array>
#include "lion/thirdparty/include/cppad/cppad.hpp"

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
    const static inline bool value = false;
};

template<typename T>
struct is_vector<std::vector<T>> 
{
    using type = std::vector<T>;
    const static inline bool value = true;
};

template<typename T>
struct is_vector<std::vector<std::vector<T>>>
{
    const static inline bool value = false;
};

template<typename T, std::size_t N>
struct is_vector<std::array<T,N>> 
{
    using type = std::array<T,N>;
    const static inline bool value = true;
};

template<typename T, std::size_t N, std::size_t M>
struct is_vector<std::array<std::array<T,M>,N>>
{
    const static inline bool value = false;
};

//! Class is_matrix: contains a bool value which is only true for 
//!    - std::vector<std::vector>
//!    - std::array<std::array>
template<typename T>
struct is_matrix
{
    const static inline bool value = false;
};

template<typename T>
struct is_matrix<std::vector<std::vector<T>>> 
{
    using type = std::vector<std::vector<T>>;
    const static inline bool value = true;
};

template<typename T, std::size_t N, std::size_t M>
struct is_matrix<std::array<std::array<T,M>,N>> 
{
    using type = std::array<std::array<T,M>,N>;
    const static inline bool value = true;
};

template<typename T>
struct is_matrix<std::vector<std::vector<std::vector<T>>>>
{
    const static inline bool value = false;
};

template<typename T, std::size_t N, std::size_t M, std::size_t L>
struct is_matrix<std::array<std::array<std::array<T,L>,M>,N>>
{
    const static inline bool value = false;
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
using combine_types_t = combine_types<T, U>::type;


/// Type-trait style function to eliminate copy constructors (via std::enable_if) from a class constructor template that
/// only has variadic arguments (i.e., "template<typename... Args> ClassName(Args&&...) { ... }").
/**
* @tparam Class The constructed class.
* @tparam FirstArg The first parameter of the class's ctor.
* @tparam RestOfArgs The rest of parameters of the class's ctor.
* @tparam std::enable_if_t<...> To declare that a ctor with this signature can never be a copy ctor.
* @return True.
*/
template<class Class, typename FirstArg, typename... RestOfArgs,
    typename std::enable_if_t<sizeof...(RestOfArgs) != 0u || !std::is_same_v<Class, std::decay_t<FirstArg> > >* = nullptr>
constexpr auto is_not_copy_ctor() { return true; }

/// Type-trait style function to eliminate copy constructors (via std::enable_if) from a class constructor template that
/// only has variadic arguments (i.e., "template<typename... Args> Class(Args&&...) { ... }").
/**
* @tparam Class The constructed class.
* @tparam Arg The single parameter of the class's ctor.
* @tparam std::enable_if_t<...> To declare that a ctor with this signature must be a copy ctor.
* @return False.
*/
template<class Class, typename SingleArg,
         typename std::enable_if_t<std::is_same_v<Class, std::decay_t<SingleArg> > >* = nullptr>
constexpr auto is_not_copy_ctor() { return false; }

#endif