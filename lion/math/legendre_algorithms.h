#ifndef __LEGENDRE_ALGORITHMS_H__
#define __LEGENDRE_ALGORITHMS_H__

#include "lion/foundation/types.h"
#include "lion/foundation/constants.h"
#include "lion/foundation/lion_exception.h"


namespace lioncpp {
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




//
//     ----------------------------------------------------------------------
//     Compute the  Legendre Polynomial by the three point recursion
//
//           L (x) = 2k-1 xL     - k-1 L
//            k      ----   k-1    ---  k-2
//                     k            k
//
//     Compute the Legendre Polynomial of degree k and its derivative
//     ------------------------------------------------------------------------
//
template<typename T>
inline constexpr T legendre_polynomial(size_t N, const T& x )
{      
    if ( N == 0 ) return 1.0;
    else if ( N == 1 ) return x;
    else
    {
        T Lnm2 = 1.0;
        T Lnm1 = x;
        T Ln   = 0.0;

        for (size_t k = 2; k <= N; ++k )
        {
            Ln = ((2.0*k-1.0)*x*Lnm1 - (k-1.0)*Lnm2)/k;
            Lnm2 = Lnm1;
            Lnm1 = Ln;    
        } 

        return Ln;
    }
}

template<typename T>
inline std::vector<T> legendre_polynomials(size_t N, const T& x )
{      
    if ( N == 0 ) return {1.0};
    else if ( N == 1 ) return {1.0,x};
    else
    {
        std::vector<T> result(N+1);
        result[0] = 1.0;
        result[1] = x;

        for (size_t k = 2; k <= N; ++k )
        {
            result[k] = ((2.0*k-1.0)*x*result[k-1] - (k-1.0)*result[k-2])/k;
        } 

        return result;
    }
}

//     Compute the Legendre Polynomial of degree k and its derivative
template<typename T>
inline constexpr std::pair<T,T> legendre_poly_and_derivative(size_t N, const T& x)
{
    if ( N == 0 )  return {1.0, 0.0};
    else if ( N == 1 ) return {x,1.0};
    else
    {
        T Lnm2 = 1.0;   T dLnm2 = 0.0;
        T Lnm1 = x ;    T dLnm1 = 1.0;
        T Ln   = 0.0;   T dLn   = 0.0;
        for (size_t k = 2; k <= N; ++k)
        {
            Ln = ((2.0*k-1.0)*x*Lnm1 - (k-1.0)*Lnm2)/k;
            dLn = dLnm2 + (2.0*k-1.0)*Lnm1;

            Lnm2 = Lnm1;  dLnm2 = dLnm1;
            Lnm1 = Ln  ;  dLnm1 = dLn;
        }

        return {Ln,dLn};
    }
}


inline std::pair<std::vector<scalar>,std::vector<scalar>> gauss_legendre_nodes_and_weights(size_t N)
{
    const size_t n_newton_iterations = 10;
    const scalar tolerance_factor = 4.0;

    const scalar tolerance = tolerance_factor*eps;

    if ( N == 0 )
        return  {std::vector<scalar>{0.0}, std::vector<scalar>{2.0}} ;

    else if ( N == 1 )
        return {std::vector<scalar>{-std::sqrt(1.0/3.0),std::sqrt(1.0/3.0)}, std::vector<scalar>{1.0, 1.0}};

    else
    {
      
        const size_t n_div_2 = (N+1)/2;

        std::vector<scalar> nodes(N+1);
        std::vector<scalar> weights(N+1);

//      Iterate on half the interior nodes
//      ----------------------------------
        for (size_t j = 0; j < n_div_2; ++j)
        {
            scalar xj = -cos( (2.0*j+1.0)*pi/(2.0*N+2.0) );
            std::pair<scalar,scalar> Ln({0.0,0.0});

            for (size_t k = 0; k < n_newton_iterations; ++k)
            {
                Ln =  legendre_poly_and_derivative( N+1, xj);
                const scalar delta = -Ln.first/Ln.second;
                xj += delta;

                if ( std::abs(delta) <= tolerance*std::abs(xj) ) break;
            }
    
            nodes[j] = xj;
            weights[j] = 2.0/( (1.0-xj*xj)*Ln.second*Ln.second);
            nodes[N-j] = -nodes[j]; 
            weights[N-j] = weights[j];
        }
//
//      ---------------------------
//      Fill in middle if necessary
//      ---------------------------
//
        if ( (N % 2) == 0 )
        {
            std::pair<scalar,scalar> Ln = legendre_poly_and_derivative(N+1,0.0);
            nodes[N/2]   = 0.0;
            weights[N/2] = 2.0/(Ln.second*Ln.second);
        }

        return {nodes,weights};
    }
}


inline scalar lagrange_polynomial(scalar x, size_t j, size_t N, const std::vector<scalar>& xj)
{
    if ( j == 0 )
    {
        scalar p = (x - xj[1])/(xj[0] - xj[1]);
  
        for (size_t k = 2; k <= N; ++k)
            p *= (x - xj[k])/(xj[0] - xj[k]);
  
        return p;
    }
    else
    {
        scalar p = (x - xj[0])/(xj[j] - xj[0]);

        for ( size_t k = 1; k < j; ++k )
            p *= (x - xj[k])/(xj[j] - xj[k]);

        for ( size_t k = j+1; k <= N; ++k)
            p *= (x - xj[k])/(xj[j] - xj[k]);

        return p;
    }
}


inline std::pair<std::vector<scalar>,std::vector<scalar>> gauss_legendre_lobatto_nodes_and_weights(const size_t N)
{
    const size_t n_newton_iterations = 10;
    const scalar tolerance_factor = 4.0;
      
    const scalar tolerance = tolerance_factor*eps;

    if ( N == 0 )
        throw lion_exception("Order must be N>0 for GL points");

    else if ( N == 1 ) 
        return std::pair(std::vector<scalar>{-1.0,1.0}, std::vector<scalar>{1.0,1.0});

    else
    {
        std::vector<scalar> x(N+1), w(N+1);
    
        x[0] = -1.0;
        w[0] = 2.0/(N*(N+1.0));

        x[N] = 1.0;
        w[N] = w[0];

        const size_t n_div_2 = (N+1.0)/2.0;

        for (size_t j = 1; j < n_div_2; ++j)
        {
            scalar xj = -cos( (j+0.25)*pi/N - 3.0/(8.0*N*pi*(j+0.25)));

            for (size_t k = 0; k < n_newton_iterations; ++k)
            {
                auto q = q_and_L_evaluation(N,xj);
                const scalar delta = -std::get<0>(q)/std::get<1>(q);
                xj += delta;
    
                if ( std::abs(delta) <= tolerance*std::abs(xj) )
                    break;
            }

            auto l_n = legendre_polynomial(N, xj);

            x[j] = xj;
            w[j] = 2.0/(N*(N+1)*l_n*l_n);

            x[N-j] = -xj;
            w[N-j] = w[j];
        }
        // Fill in the middle if necessary
        if ( (N % 2) == 0 )
        {
            scalar l_n = legendre_polynomial(N, 0.0);
            x[N/2] = 0.0;
            w[N/2] = 2.0/(N*(N+1)*l_n*l_n);
        }

        return {x,w};
    }
}


inline std::tuple<scalar,scalar,scalar> q_and_L_evaluation(const size_t N, const scalar x)
{
    if ( N == 0 )
        throw lion_exception("q_and_L_evaluation cannot be called for N=0");

    else if ( N == 1 )
        throw lion_exception("q_and_L_evaluation cannot be called for N=1");

    else
    {
        scalar L_kM2 = 1.0;
        scalar L_prime_kM2 = 0.0;

        scalar L_kM1 = x;
        scalar L_prime_kM1 = 1.0;

        scalar L_k;
        scalar L_prime_k;

        for (size_t k = 2; k <= N; ++k)
        {
            L_k = ((2.0*k-1.0)*x*L_kM1 - (k-1.0)*L_kM2)/k;    
            L_prime_k = L_prime_kM2 + (2.0*k-1.0)*L_kM1;
            L_kM2 = L_kM1;
            L_kM1 = L_k;
            L_prime_kM2 = L_prime_kM1;
            L_prime_kM1 = L_prime_k;
        }

        const size_t k = N+1;
        L_k = ((2.0*k-1.0)*x*L_kM1 - (k-1.0)*L_kM2)/k;
        L_prime_k = L_prime_kM2 + (2.0*k-1.0)*L_kM1;

               /*   Q              Qprime             L_N */
        return { L_k-L_kM2, L_prime_k - L_prime_kM2, L_kM1 };
    }
}

}

#endif
