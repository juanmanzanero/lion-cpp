#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include "lion/foundation/types.h"
#include "legendre_algorithms.h"
#include "vector3d.hpp"

//!         A polynomial class
//!         ------------------
//!
//!  This class defines spectral like polynomials. The class can contain a series of _n_blocks 
//! polynomial blocks, where the polynomial is evaluated from with Legendre polynomials.
//! Each block expands from _ai[i] to _bi[i], and the polynomial expands from _a = _ai.front() 
//! to b = _bi.back().
//! 
//! Two modes are implemented to instantiate a polynomial from (x0,y0) pairs
//!     1) The "smooth" version: the (x0,y0) points are projected into a given order single polynomial
//!        --------------------  block. For this projection, the (x0,y0) pairs are divided in chunks
//!                              of 5 points and this sub-polynomials are projected onto the LGL points
//!                              to compute the final polynomial
//!
//!     2) The non-smooth version: the (x0,y0) are divided in polynomial blocks of N+1 (given) points: 
//!        ----------------------  P1 = (x0,x1,...,xN), P2 = (xN+1,...), ... which are 
//!                                used as sub-polynomials. When an evaluation at x is requested,
//!                                first the block that contains this point is found, and then
//!                                this polynomial is evaluated at x
//!                                 
template<class T>
class Polynomial
{
 public:
    using value_type = T;

    //! Default constructor (evaluates at 0)
    Polynomial();

    //! Constructor from pairs (x0,y0) and given order. If smooth == false, it is computed as a piecewise
    //! order N polynomial, taking order N chunks of (x0,y0). If smooth == true, the polynomial is 
    //! projected to an order N Legendre-Gauss-Lobatto basis. This projection uses an order 4 piecewise
    //! interpolation to avoid Runge's phenomenon
    //! @param[in] x0: x-coordinates
    //! @param[in] y0: y-coordinates
    //! @param[in] N: desired interpolation order
    //! @param[in] smooth: whether it is projected to an order N LGL basis or kept as order N piecewise
    Polynomial(const std::vector<scalar>& x0, 
               const std::vector<T>& y0, 
               const size_t N, 
               bool smooth = false
              );
    
    //! Constructor from polynomial orders, block bounds and legendre coefficients
    //! @param[in] N: polynomial orderorder
    //! @param[in] ai: left bounds for each blocks
    //! @param[in] bi: right bounds for each block
    //! @param[in] coeffs: Legendre coefficients of each block
    Polynomial(const std::vector<size_t>& N, const std::vector<scalar>& ai, const std::vector<scalar>& bi, 
               const std::vector<std::vector<T>>& coeffs);

    //! Constructor from length of the blocks and values at the LGL points for each block
    //! @param[in] x0: left bound of the polynomial
    //! @param[in] L: array that contains the length of the blocks
    //! @param[in] y0: values of the polynomials at the Legendre-Gauss-Lobatto points for each block
    Polynomial(const scalar x0, const std::vector<scalar>& L, const std::vector<std::vector<T>>& y0);

    //! Constructor from (x0,y0), providing different blocks: x0.size() polynomials are constructed as:
    //! p[i] = Polynomial(x0[i],y0[i],smooth)
    //! and are then merged using Polynomial(p)
    //! @param[in] x0: x0 coordinates for each block
    //! @param[in] y0: y0 coordinates for each block
    //! @param[in] continuous: if the blocks are expected to be continuous (it only enables the check)
    //! @param[in] smooth: whether each block is projected to an order N LGL basis or kept as 
    //!                    order N piecewise
    Polynomial(const std::vector<std::vector<scalar>>& x0, 
               const std::vector<std::vector<T>>& y0, 
               bool continuous = true, 
               bool smooth = false
              );

    //! Constructor from vector of polynomials: each polynomial is a block of the new poly
    //! @param[in] p: vector of polynomials
    Polynomial(const std::vector<Polynomial>& p);

    //! Non-safe evaluation of the polynomial at x 
    //! @param[in] x: evaluation point 
    template<typename U>
    typename combine_types<U,T>::type operator[](const U& x) const;

    //! Safe evaluation of the polynomial at x (exception if x is out of bounds)
    //! @param[in] x: evaluation point 
    template<typename U>
    typename combine_types<U,T>::type operator()(const U& x) const;

    //! Computes the derivative of p
    Polynomial derivative() const;

    //! Computes the integral of p such that int_p(a) = 0
    //! The integral is always a continuous polynomial
    Polynomial integral() const;

    //! Get the number of blocks
    constexpr const size_t& get_n_blocks() const { return _n_blocks; }

    //! Get the coefficients for a given block
    //! @param[in] n: index of the block
    const std::vector<T>& get_coeffs(size_t n = 0) const { return _coeffs[n]; }

    //! Get the polynomial left bound (a)
    constexpr const scalar& get_left_bound() const { return _a; }

    //! Get the polynomial right bound (b)
    constexpr const scalar& get_right_bound() const { return _b; }
    
 private:
    //! Compute the Legendre polynomial coefficients of the polynomial that contains the 
    //! (x0,y0) points
    //! @param[in] ix0: iterator to x0 points
    //! @param[in] iy0: iterator to y0 points
    //! @param[in] N: number of (x0,y0) points - 1 (order of the final polynomial)
    static std::vector<T> compute_coefficients(const typename std::vector<scalar>::const_iterator& ix0, 
                                               const typename std::vector<T>::const_iterator& iy0,
                                               const size_t N);

    //! Compute the Legendre coefficients from values at the LGL points (aka nodal to modal)
    //! @param[in] y0: polynomial values at the Legendre-Gauss-Lobatto points
    //! @param[in] N: polynomial order
    static std::vector<T> compute_coefficients(const std::vector<T>& y0, const size_t N);

    //! Compute the order N Legendre coefficients of a given polynomial 
    //! @param[in] p: the polynomial
    //! @param[in] N: the new polynomial order
    static std::vector<T> compute_coefficients(const Polynomial<T>& p, const size_t N);

    inline static size_t N_smooth = 4;  //! Polynomial order used to construct and intermediate 
                                        //! piecewise polynomial in the smooth polynomial version

    size_t _n_blocks;       //! The number of blocks
    std::vector<size_t> _n; //! The polynomial order of each block
    scalar _a;              //! The left bound
    scalar _b;              //! The right bound
    
    std::vector<scalar> _ai;    //! The left bound of each block
    std::vector<scalar> _bi;    //! The right bound of each block
    
    std::vector<std::vector<T>> _coeffs;  //! The Legendre polynomials coefficients for each block
};


//! An auxiliary class to handle a vector of polynomials. p_vector(x) = [p1(x),p2(x),p3(x),...]
template<class T>
class Polynomial_vector
{
 public: 
    //! Constructor from vector of polynomials 
    //! @param[in] p: vector of polynomials
    Polynomial_vector(const std::vector<Polynomial<T>>& p): _p(p) {};

    //! Constructor from number of polynomials: constructs n empty polynomials
    //! @param[in] n: number of polynomials
    Polynomial_vector(size_t n): _p(n, Polynomial<T>()) {};

    //! Return the evaluation of the polynomials at x as a vector
    //! @param[in] x: evaluation point
    template<typename U>
    std::vector<typename combine_types<U,T>::type> operator()(const U& x) const
    {
	using result_type = typename combine_types<U,T>::type;
        std::vector<result_type> result(_p.size());

        for (size_t i = 0; i < _p.size(); ++i)
            result[i] = _p[i](x);
    
        return result;
    }

    //! Return the polynomial with index i
    //! @param[in] i: index
    const Polynomial<T>& at(size_t i) const { return _p.at(i); }
          Polynomial<T>& at(size_t i)       { return _p.at(i); }


 private:
    std::vector<Polynomial<T>> _p;  //! The vector of polynomials
};

//! An auxiliary class to handle an array of polynomials. p_array(x) = [p1(x),p2(x),p3(x),...]
template<class T,size_t N>
class Polynomial_array
{
 public: 
    //! Default constructor
    Polynomial_array() = default;

    //! Constructor from number of polynomials: constructs n empty polynomials
    //! @param[in] n: number of polynomials
    Polynomial_array(const std::array<Polynomial<T>,N>& p): _p(p) {};

    //! Return the evaluation of the polynomials at x as a vector
    //! @param[in] x: evaluation point
    template<typename U>
    std::array<typename combine_types<U,T>::type,N> operator()(const U& x) const
    {
	using result_type = typename combine_types<U,T>::type;
        std::array<result_type,N> result;

        for (size_t i = 0; i < _p.size(); ++i)
            result[i] = _p[i](x);
    
        return result;
    }

    //! Return the polynomial with index i
    //! @param[in] i: index
    const Polynomial<T>& at(size_t i) const { return _p.at(i); }
          Polynomial<T>& at(size_t i)       { return _p.at(i); }

 private:
    std::array<Polynomial<T>,N> _p; //! The array of polynomials
};

using sPolynomial = Polynomial<scalar>;
using vPolynomial = Polynomial<sVector3d>;
using sPolynomial_vector = Polynomial_vector<scalar>;

#include "polynomial.hpp"

#endif
