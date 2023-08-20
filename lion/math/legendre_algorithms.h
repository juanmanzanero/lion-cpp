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


//! Compute the barycentric weights of a given set of nodes
template<typename scalar_type = scalar>
inline std::vector<scalar_type> compute_barycentric_weights(const std::vector<scalar_type>& nodes)
{
    assert(nodes.size() > 0u);

    const auto polynomial_order = nodes.size() - 1;
    std::vector<scalar_type> weights(polynomial_order + 1, 1.0);

    for (size_t j = 1; j <= polynomial_order; ++j) {
        for (size_t k = 0; k < j; ++k) {
            weights[k] *= nodes[k] - nodes[j];
            weights[j] *= nodes[j] - nodes[k];
        }
    }

    std::for_each(weights.begin(), weights.end(), [](auto& w) { w = 1.0 / w; });

    return weights;
}

//! Compute the polynomial derivative matrix Dij = lj'(xi)
//! We store it rowmajor, it is easier to use later this way: Dij = D[(N+1)*i + j]
//! @param[in] nodes
template<typename scalar_type = scalar>
inline std::vector<scalar_type> compute_derivative_matrix(const std::vector<scalar_type>& nodes)
{
    assert(nodes.size() > 0u);

    const auto barycentric_weights = compute_barycentric_weights(nodes);
    const auto polynomial_order = nodes.size() - 1;
    std::vector<scalar_type> derivative_matrix((polynomial_order + 1) * (polynomial_order + 1), 0.0);

    for (size_t i = 0; i <= polynomial_order; ++i) {
        const auto diagonal_element = (polynomial_order + 2) * i;
        for (size_t j = 0; j <= polynomial_order; ++j) {
            if (j != i)
            {
                const auto element = (polynomial_order + 1) * i + j;

                derivative_matrix[element] = barycentric_weights[j] / (barycentric_weights[i] * (nodes[i] - nodes[j]));
                derivative_matrix[diagonal_element] -= derivative_matrix[element];
            }
        }
    }

    return derivative_matrix;

}

#include "legendre_algorithms.hpp"

#endif
