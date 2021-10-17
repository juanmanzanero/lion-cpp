#ifndef __TYPE_TRAITS_H__
#define __TYPE_TRAITS_H__

#include <type_traits>
#include <vector>
#include <array>

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

#endif
