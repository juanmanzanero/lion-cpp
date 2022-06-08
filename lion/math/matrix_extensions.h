#ifndef __MATRIX_EXTENSIONS_H__
#define __MATRIX_EXTENSIONS_H__

#include <vector>
#include <assert.h>
#include <array>
#include "lion/foundation/type_traits.h"

//! Sum two vectors
//! @param[in] lhs: left hand side vector
//! @param[in] rhs: right hand side vector
template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type* = nullptr, typename is_vector<VEC2>::type* = nullptr>
constexpr VEC1 operator+(const VEC1& lhs, const VEC2& rhs );

//! Increment vector by vector
//! @param[in] lhs: left hand side vector
//! @param[in] rhs: right hand side vector
template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type* = nullptr, typename is_vector<VEC2>::type* = nullptr>
constexpr VEC1& operator+=(VEC1& lhs, const VEC2& rhs );

//! Subtract two vectors
//! @param[in] lhs: left hand side vector
//! @param[in] rhs: right hand side vector
template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type* = nullptr, typename is_vector<VEC2>::type* = nullptr>
constexpr VEC1 operator-(const VEC1& lhs, const VEC2& rhs );

//! Unitary minus operator for vector
//! @param[in] rhs: original vector
template<class T, size_t N>
constexpr std::array<T,N> operator-(const std::array<T,N>& rhs );

//! Sum two matrices of std::array
//! @param[in] lhs: left hand side vector
//! @param[in] rhs: right hand side vector
template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,M>,N> operator+(const std::array<std::array<T,M>,N>& lhs, 
                                        const std::array<std::array<T,M>,N>& rhs );

//! Subtract two matrices of std::array
//! @param[in] lhs: left hand side vector
//! @param[in] rhs: right hand side vector
template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,M>,N> operator-(const std::array<std::array<T,M>,N>& lhs, 
                                        const std::array<std::array<T,M>,N>& rhs );

//! Unitary minus of matrices of std::array
//! @param[in] rhs: right hand side vector
template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,M>,N> operator-(const std::array<std::array<T,M>,N>& rhs);

//! Matrix product of two std::array matrices
//! @param[in] lhs: left hand side vector
//! @param[in] rhs: right hand side vector
template<class T, size_t N, size_t L, size_t M>
constexpr std::array<std::array<T,M>,N> operator*(const std::array<std::array<T,L>,N>& lhs,
                                        const std::array<std::array<T,M>,L>& rhs );

// Implementation of max, and abs functions for vectors 
namespace ext 
{
    //! Get a vector with the component-wise maximum of a vector v and a constant k
    //! u_i = std::max(v_i, k)
    //! @param[in] v: the vector
    //! @param[in] k: the constant
    template<typename VEC, typename is_vector<VEC>::type* = nullptr>
    constexpr VEC max(const VEC& v, const double k);

    //! Get a vector with the component-wise maximum of two vectors
    //! v_i = std::max(lhs_i, rhs_i)
    //! @param[in] lhs: the first vector
    //! @param[in] rhs: the second vector
    template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type* = nullptr, typename is_vector<VEC2>::type* = nullptr>
    constexpr VEC1 max(const VEC1& lhs, const VEC2& rhs);

    //! Get the component-wise absolute value of a vector
    //! u_i = std::abs(v_i)
    //! @param[in] v: the vector
    template<typename VEC, typename is_vector<VEC>::type* = nullptr>
    constexpr VEC abs(const VEC& v);
}


//! Transpose a std::array matrix
//! @param[in] m: the matrix
template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,N>,M> transpose(const std::array<std::array<T,M>,N>& m);

//! Compute the matrix-vector product, based on std::arrays
//! @param[in] lhs: the matrix
//! @param[in] rhs: the vector
template<class T, size_t N, size_t M>
constexpr std::array<T,N> operator*(const std::array<std::array<T,M>,N>& lhs, const std::array<T,M>& rhs);

//! Compute the vector-matrix product, based on std::arrays
//! @param[in] lhs: the vector
//! @param[in] rhs: the matrix
template<class T, size_t N, size_t M>
constexpr std::array<T,M> operator*(const std::array<T,N>& lhs, const std::array<std::array<T,M>,N>& rhs);

//! Multiply a vector by a constant
//! @param[in] lhs: the vector
//! @param[in] rhs: the constant
template<typename VEC, typename T, typename is_vector<VEC>::type* = nullptr>
constexpr VEC operator*(const VEC& lhs, const T& rhs);

//! Multiply a vector by a constant
//! @param[in] lhs: the vector
//! @param[in] rhs: the constant
template<typename VEC, typename is_vector<VEC>::type* = nullptr>
constexpr VEC& operator*=(VEC& lhs, const double rhs);

//! Multiply a constant by a vector
//! @param[in] lhs: the constant
//! @param[in] rhs: the vector
template<typename VEC, typename T, typename is_vector<VEC>::type* = nullptr>
constexpr VEC operator*(const T& lhs, const VEC& rhs) { return rhs*lhs; }

//! Compute a component wise vector division
//! @param[in] lhs: dividend vector
//! @param[in] rhs: divisor vector
template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type* = nullptr, typename is_vector<VEC2>::type* = nullptr>
constexpr VEC1 operator/(const VEC1& lhs, const VEC2& rhs);

//! Solve a linear system of order N, with matrix A and vector b
//! @param[in] N: order of the system
//! @param[in] A: coefficients matrix
//! @param[in] b: independent terms vector
template<typename MAT, typename VEC, typename is_matrix<MAT>::type* = nullptr, typename is_vector<VEC>::type* = nullptr>
inline VEC linsolve(int N, MAT A, const VEC& b);

//! Print a std::array matrix
//! @param[in] m: the matrix
template <typename T,size_t N,size_t M>
std::ostream& operator<<(std::ostream &os, const std::array<std::array<T,M>,N>& m);

//! Print a std::array
//! @param[in] v: the array
template <typename T,size_t N>
std::ostream& operator<<(std::ostream &os, const std::array<T,N>& v);

//! Print a std::vector
//! @param[in] v: the vector
template<typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T>& v);

#include "matrix_extensions.hpp"
#endif
