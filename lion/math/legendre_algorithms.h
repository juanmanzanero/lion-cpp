#ifndef __LEGENDRE_ALGORITHMS_H__
#define __LEGENDRE_ALGORITHMS_H__

#include "lion/foundation/types.h"
#include "lion/foundation/constants.h"

//! Compute the  Legendre Polynomial by the three point recursion
//!
//!          L (x) = 2k-1 xL     - k-1 L
//!           k      ----   k-1    ---  k-2
//!                    k            k
//!
//! @param[in] N: polynomial order
//! @param[in] x: evaluation point
template<typename T>
constexpr T legendre_polynomial(size_t N, const T& x );

//! Get the legendre polynomials up to order N evaluated at x
//! @param[in] N: final polynomial order
//! @param[in] x: evaluation point
template<typename T>
std::vector<T> legendre_polynomials(size_t N, const T& x );

//! Compute the Legendre Polynomial of degree N and its derivative
//! @param[in] N: polynomial order
//! @param[in] x: evaluation point
template<typename T>
constexpr std::pair<T,T> legendre_poly_and_derivative(size_t N, const T& x);

//! Compute the Gauss-Legendre quadrature nodes and weights
//! @param[in] N: polynomial order
std::pair<std::vector<scalar>,std::vector<scalar>> gauss_legendre_nodes_and_weights(size_t N);

//! Compute the j-th Lagrange polynomial with nodes {xj} evaluated at x
//! @param[in] x: evaluation point
//! @param[in] j: index of the Lagrange polynomial (starts at 0)
//! @param[in] N: polynomial order
//! @param[in] xj: interpolation nodes
scalar lagrange_polynomial(scalar x, size_t j, size_t N, const std::vector<scalar>& xj);

//! Compute the Gauss-Legendre-Lobatto quadrature nodes and weights
//! @param[in] N: polynomial order
std::pair<std::vector<scalar>,std::vector<scalar>> gauss_legendre_lobatto_nodes_and_weights(const size_t N);

//! Compute the function Q_N = L_{N+1} - L_{N-1} and
//! its derivative to find the Gauss-Lobatto points and
//! weights.
//! @param[in] N: polynomial order
//! @param[in] x: evaluation point
std::tuple<scalar,scalar,scalar> q_and_L_evaluation(const size_t N, const scalar x);

#include "legendre_algorithms.hpp"

#endif
