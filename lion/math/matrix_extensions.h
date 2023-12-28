#ifndef __MATRIX_EXTENSIONS_H__
#define __MATRIX_EXTENSIONS_H__

#include <vector>
#include <assert.h>
#include <array>
#include <iostream>
#include <cmath>
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

//! Add a std::vector and a scalar
template<typename T>
constexpr std::vector<T> operator+(const std::vector<T>& lhs, const T& rhs);

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



template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type*, typename is_vector<VEC2>::type*>
constexpr VEC1 operator+(const VEC1& lhs, const VEC2& rhs )
{
    assert(lhs.size() == rhs.size());

    VEC1 result = lhs;

    for (size_t i = 0; i < lhs.size(); ++i)
        result[i] += rhs[i];

    return result;
}

template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type*, typename is_vector<VEC2>::type*>
constexpr VEC1& operator+=(VEC1& lhs, const VEC2& rhs )
{
    assert(lhs.size() == rhs.size());

    for (size_t i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];

    return lhs;
}



template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type*, typename is_vector<VEC2>::type*>
constexpr VEC1 operator-(const VEC1& lhs, const VEC2& rhs )
{
    assert(lhs.size() == rhs.size());

    VEC1 result = lhs;

    for (size_t i = 0; i < lhs.size(); ++i)
        result[i] -= rhs[i];

    return result;
}


template<class T, size_t N>
constexpr std::array<T,N> operator-(const std::array<T,N>& rhs )
{
    std::array<T,N> result = rhs;

    for (size_t i = 0; i < N; ++i)
        result[i] *= -1.0;

    return result;
}


template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,M>,N> operator+(const std::array<std::array<T,M>,N>& lhs, 
                                        const std::array<std::array<T,M>,N>& rhs )
{
    std::array<std::array<T,M>,N> result = {{{0.0}}};
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            result[i][j] = lhs[i][j]+rhs[i][j];

    return result;
}

template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,M>,N> operator-(const std::array<std::array<T,M>,N>& lhs, 
                                        const std::array<std::array<T,M>,N>& rhs )
{
    std::array<std::array<T,M>,N> result = {{{0.0}}};
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            result[i][j] = lhs[i][j]-rhs[i][j];

    return result;
}

template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,M>,N> operator-(const std::array<std::array<T,M>,N>& rhs)
{
    std::array<std::array<T,M>,N> result = {{{0.0}}};
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            result[i][j] = -rhs[i][j];

    return result;
}


template<class T, size_t N, size_t L, size_t M>
constexpr std::array<std::array<T,M>,N> operator*(const std::array<std::array<T,L>,N>& lhs,
                                        const std::array<std::array<T,M>,L>& rhs )
{
    std::array<std::array<T,M>,N> result = {{{0.0}}};

    for (size_t i = 0; i < N; ++i)
        for (size_t l = 0; l < L; ++l)
            for (size_t j = 0; j < M; ++j)
                result[i][j] += lhs[i][l]*rhs[l][j];

    return result;
}

template<typename VEC, typename is_vector<VEC>::type*>
constexpr VEC ext::max(const VEC& v, const double k)
{
    auto result = v;

    for (auto it = result.begin(); it != result.end(); ++it)
        *it = std::max(*it,k);

    return result;
}


template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type*, typename is_vector<VEC2>::type*>
constexpr VEC1 ext::max(const VEC1& lhs, const VEC2& rhs)
{
    assert(lhs.size() == rhs.size());

    auto result = lhs;

    for (size_t i = 0; i < lhs.size(); ++i)
        (rhs[i] > result[i] ? result[i] = rhs[i] : true);

    return result;
}

template<typename VEC, typename is_vector<VEC>::type*>
constexpr VEC ext::abs(const VEC& v)
{
    auto result = v;

    for (auto it = result.begin(); it != result.end(); ++it)
        *it = std::abs(*it);

    return result;
}


template<class T, size_t N, size_t M>
constexpr std::array<std::array<T,N>,M> transpose(const std::array<std::array<T,M>,N>& m)
{
    std::array<std::array<T,N>,M> result = {{{0.0}}};

    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            result[j][i] = m[i][j];

    return result;
}


template<class T, size_t N, size_t M>
constexpr std::array<T,N> operator*(const std::array<std::array<T,M>,N>& lhs, const std::array<T,M>& rhs)
{
    std::array<T,N> result = {{0.0}};

    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            result[i] += lhs[i][j]*rhs[j];

    return result;
}

template<class T, size_t N, size_t M>
constexpr std::array<T,M> operator*(const std::array<T,N>& lhs, const std::array<std::array<T,M>,N>& rhs)
{
    std::array<T,M> result = {{0.0}};

    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            result[j] += lhs[i]*rhs[i][j];

    return result;
}


template<typename VEC, typename T, typename is_vector<VEC>::type*>
constexpr VEC operator*(const VEC& lhs, const T& rhs)
{
    auto result = lhs;
    
    for (auto it = result.begin(); it != result.end(); ++it)
        *it *= rhs;

    return result;
}

template<typename VEC, typename is_vector<VEC>::type*>
constexpr VEC& operator*=(VEC& lhs, const double rhs)
{
    for (auto it = lhs.begin(); it != lhs.end(); ++it)
        *it *= rhs;

    return lhs;
}



template<typename VEC1, typename VEC2, typename is_vector<VEC1>::type*, typename is_vector<VEC2>::type*>
constexpr VEC1 operator/(const VEC1& lhs, const VEC2& rhs)
{
    assert( lhs.size() == rhs.size() );

    auto result = lhs;
    for (size_t i = 0; i < lhs.size(); ++i)
        result[i] *= (1.0/rhs[i]);

    return result;
}



template <typename T,size_t N,size_t M>
std::ostream& operator<<(std::ostream &os, const std::array<std::array<T,M>,N>& m)
{
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
            os << m[i][j] << " " ;
        os << '\n';
    }

    return os;
}


template<typename MAT, typename VEC, typename is_matrix<MAT>::type*, typename is_vector<VEC>::type*>
inline VEC linsolve(int N, MAT A, const VEC& b)
{
  int i;
  int j;
  int k;
  int p;
  double t;
//
//  Create X and set it to B.
//
  VEC x = b;

  for ( k = 0; k < N; k++ )
  {
//
//  Find the maximum element in column I.
//
    p = k;

    for ( i = k + 1; i < N; i++ )
    {
      if ( std::fabs ( A[p][k] ) < std::fabs ( A[i][k] ) )
      {
        p = i;
      }
    }

    if ( A[p][k] == 0.0 )
    {
      std::cerr << "\n";
      std::cerr << "R8RMAT_FS_NEW - Fatal error!\n";
      std::cerr << "  Zero pivot on step " << k << "\n";
      throw;
    }
//
//  Switch rows K and P.
//
    if ( k != p )
    {
      for ( j = 0; j < N; j++ )
      {
        t        = A[k][j];
        A[k][j] = A[p][j];
        A[p][j] = t;
      }
      t    = x[k];
      x[k] = x[p];
      x[p] = t;
    }
//
//  Scale the pivot row.
//
    t = A[k][k];
    A[k][k] = 1.0;
    for ( j = k + 1; j < N; j++ )
    {
      A[k][j] = A[k][j] / t;
    }
    x[k] = x[k] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = k + 1; i < N; i++ )
    {
      if ( A[i][k] != 0.0 )
      {
        t = - A[i][k];
        A[i][k] = 0.0;
        for ( j = k + 1; j < N; j++ )
        {
          A[i][j] = A[i][j] + t * A[k][j];
        }
        x[i] = x[i] + t * x[k];
      }
    }
  }
//
//  Back solve.
//
  for ( j = N - 1; 1 <= j; j-- )
  {
    for ( i = 0; i < j; i++ )
    {
      x[i] = x[i] - A[i][j] * x[j];
    }
  }

  return x;
}


template<typename T>
constexpr std::vector<T> operator+(const std::vector<T>& lhs, const T& rhs)
{
    auto result = lhs;
    std::for_each(result.begin(), result.end(), [&](auto& result_i) { result_i += rhs; });
    
    return result;
}


template <typename T,size_t N>
std::ostream& operator<<(std::ostream &os, const std::array<T,N>& v)
{
    for (size_t j = 0; j < N; ++j)
        os << v[j] << " " ;

    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T>& v)
{
    if (v.empty()) {
        os << "{}";
    }
    else {
        os << "{" ;
        for (auto it = v.cbegin(); it < std::prev(v.cend()); ++it) {
            os << *it << ", " ;
        }
        os << v.back() << "}" ;
    }

    return os;
}


#endif
