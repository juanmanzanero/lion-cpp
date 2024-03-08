#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <numeric>

#include "lion/foundation/types.h"
#include "legendre_algorithms.h"
#include "vector3d.hpp"

namespace lioncpp {
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
template<typename T>
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

    //! Returns a block as polynomial
    Polynomial get_block(const size_t block_id) const;

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

    //! Check if the polynomial is zero
    constexpr bool is_zero() const;

    //! Check if the polynomial is zero, comparison functor provided by user
    template<typename UnaryOperation>
    constexpr bool is_zero(UnaryOperation unary_op) const;
    
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

    size_t _N_smooth = 4;  //! Polynomial order used to construct and intermediate 
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
template<typename T>
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
template<typename T,size_t N>
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




template<typename T>
inline Polynomial<T>::Polynomial()
: _n_blocks(0),
  _n(0),
  _a(0.0),
  _b(0.0),
  _ai(),
  _bi(),
  _coeffs()
{}


template<typename T>
inline Polynomial<T>::Polynomial(const std::vector<scalar>& x0, const std::vector<T>& y0, const size_t N, bool smooth)
: _n_blocks(smooth ? 1 : (x0.size()-1)/(N) + ( ((x0.size()-1 % (N))!=0) ? 1 : 0 )),
  _n(_n_blocks,N),
  _a(x0.front()),
  _b(x0.back()),
  _ai(_n_blocks),
  _bi(_n_blocks), 
  _coeffs(_n_blocks)
{
    if ( x0.size() != y0.size() )
        throw lion_exception("x0 and y0 must have the same length");

    if ( x0.size() == 0 )
    {
        *this = Polynomial();
        return;
    }

    if ( !smooth )
    {
        // Non-smooth version: construct piecewise order N polynomials

        auto ix0 = x0.cbegin();
        auto iy0 = y0.cbegin();
        
        // Fill all blocks except the last one
        for (int n = 0; n < (int)_n_blocks - 1; ++n)
        {
            _coeffs[n] = compute_coefficients(ix0,iy0,N);
    
            _ai[n] = *ix0;
            _bi[n] = *(ix0+N);
    
            ix0 += N;
            iy0 += N;
        }
    
        // Add the last block
        if ( x0.size() > N+1 )
        {
            // Get the remaining points
            ix0 = x0.cend() - N - 1 ;
            iy0 = y0.cend() - N - 1 ;
    
            _coeffs.back() = compute_coefficients(ix0,iy0,N);
    
            _ai.back() = *ix0;
            _bi.back() = x0.back();
        }
        else
        {
            // Put all points if they are less than the desired order
            ix0 = x0.cbegin();
            iy0 = y0.cbegin();

            _coeffs.back() = compute_coefficients(ix0,iy0,x0.size()-1);
            _ai.back() = x0.front();
            _bi.back() = x0.back();
            _n.back()  = x0.size()-1;
        }
    }
    else
    {
        _ai.front() = _a;
        _bi.front() = _b;

        // Smooth version: construct an auxiliar piecewise polynomial of order N_smooth
        // and interpolate it to LGL points afterwards
        Polynomial p(x0, y0, _N_smooth, false);

        _coeffs.front() = compute_coefficients(p, N);
    }
}


template<typename T>
inline Polynomial<T>::Polynomial(const std::vector<size_t>& N, const std::vector<scalar>& ai, const std::vector<scalar>& bi, 
    const std::vector<std::vector<T>>& coeffs)
: _n_blocks(coeffs.size()),
  _n(N),
  _a(ai.front()),
  _b(bi.back()),
  _ai(ai),
  _bi(bi),
  _coeffs(coeffs)
{
    if ( _coeffs.size() == 0 )
    {
        *this = Polynomial();
        return;
    }
}


template<typename T>
inline Polynomial<T>::Polynomial(const double x0, const std::vector<scalar>& L, const std::vector<std::vector<T>>& y0)
: _n_blocks(L.size()),
  _n(_n_blocks),
  _a(x0),
  _b(std::accumulate(L.begin(), L.end(), x0)),
  _ai(_n_blocks),
  _bi(_n_blocks),
  _coeffs(_n_blocks)
{
    // Check consistency
    if ( L.size() != y0.size() )
        throw lion_exception("L and y0 must have the same length");

    if ( L.size() == 0 )
    {
        *this = Polynomial();
        return;
    }

    // Set blocks bounds
    _ai.front() = _a;
    _bi.back()  = _b;
    for (size_t i = 0; i < _n_blocks-1; ++i)
    {
        _ai[i+1] = _ai[i] + L[i];
        _bi[i]   = _ai[i+1];
    }

    // Set blocks order and compute coefficients
    for (size_t i = 0; i < _n_blocks; ++i)
    {
        _n[i] = y0[i].size()-1;
        _coeffs[i] = compute_coefficients(y0[i],_n[i]);
    }
}


template<typename T>
inline Polynomial<T>::Polynomial(const std::vector<std::vector<scalar>>& x0, const std::vector<std::vector<T>>& y0, 
        bool continuous, bool smooth)
: _n_blocks(x0.size()),
  _n(_n_blocks),
  _a(x0.front().front()),
  _b(x0.back().back()),
  _ai(_n_blocks),
  _bi(_n_blocks),
  _coeffs(_n_blocks)
{
    // Check sizes consistency
    if ( x0.size() != y0.size() )
        throw lion_exception("x0 and y0 must have the same length");

    if ( x0.size() == 0 )
    {
        *this = Polynomial();
        return;
    }

    // Check x0 consistency
    for (auto it = x0.cbegin(); it != x0.cend()-1; ++it)
        if ( (*(it+1)).back() < ((*it).back() + eps) )
            throw lion_exception("blocks are not contiguous");

    // Check y0 consistency if continuous
    if ( continuous ) 
        for (auto it = y0.cbegin(); it != y0.cend()-1; ++it)
            if ( std::abs( (*it).back() - ((*(it+1)).front()) ) > eps )
                throw lion_exception("blocks are not continuous, and continuous = true");


    std::vector<Polynomial> p(x0.size());

    for (size_t i = 0; i < _n_blocks; ++i)
        p[i] = Polynomial(x0[i], y0[i], smooth);

    *this = Polynomial(p);
}


template<typename T>
inline Polynomial<T>::Polynomial(const std::vector<Polynomial>& p)
: _n_blocks(0)
{
    for (const Polynomial& pol : p)
    {
        _n_blocks += pol._n_blocks;

        for (size_t i = 0; i < pol._n_blocks; ++i)
        {
            _n.push_back(pol._n[i]);
            _ai.push_back(pol._ai[i]);
            _bi.push_back(pol._bi[i]);
            _coeffs.push_back(pol._coeffs[i]);
        }
    }
    _a = _ai.front();
    _b = _bi.back();
}


template<typename T>
inline Polynomial<T> Polynomial<T>::get_block(const size_t block_id) const
{
    if (block_id >= _n_blocks)
        throw lion_exception("[ERROR] Polynomial::get_block -> block_id is out of bounds");

    return { {_n[block_id]},{_ai[block_id]},{_bi[block_id]},{_coeffs[block_id]} };
}


template<class T>
template<typename U>
inline typename combine_types<U,T>::type Polynomial<T>::operator[](const U& x) const
{
    using return_type = typename combine_types<U,T>::type;

    if ( _n_blocks == 0 )
        return return_type{};

    size_t block = 0;

    while(x > _bi[block])
        ++block;

    const auto xi = 2.0*(x-_ai[block])/(_bi[block]-_ai[block]) - 1.0;

    const auto Lk = legendre_polynomials(_n[block],xi);

    auto result = return_type{};

    if constexpr (is_specialization_v<T, std::vector>) {
        if (!_coeffs.empty() && !_coeffs.front().empty()) {
            result.resize(_coeffs.front().front().size(), typename T::value_type{ 0 });
        }
    }


    for (size_t i = 0; i <= _n[block]; ++i)
        result += _coeffs[block][i]*Lk[i];

    return result;
}


template<class T>
template<typename U>
inline typename combine_types<U,T>::type Polynomial<T>::operator()(const U& x) const
{
    using return_type = typename combine_types<U,T>::type;

    if ( _n_blocks == 0 )
        return return_type{};

    // Check bounds
    if ( (x < _a-100*eps) || (x > _b+100.0*eps) )
    {
        std::ostringstream s_out;
        s_out << std::setprecision(17);
        s_out << "x is out of bounds" << std::endl;
        s_out << "  a: " << _a << std::endl;
        s_out << "  x: " <<  x << std::endl;
        s_out << "  b: " << _b ;
        throw lion_exception(s_out.str());
    }

    // Call operator[]
    return (*this)[x];
}


template<typename T>
inline Polynomial<T> Polynomial<T>::derivative() const
{
    std::vector<std::vector<T>> coeffs(_n_blocks);

    for (size_t n = 0; n < _n_blocks; ++n)
    {
        coeffs[n] = std::vector<T>(_n[n],T());
        for (size_t i = 0; i <= _n[n]; ++i)
            for (std::ptrdiff_t j = static_cast<std::ptrdiff_t>(i) - 1 ; j >= 0; j -= 2)
                coeffs[n][j] += _coeffs[n][i]*(2.0*(j+1.0)-1.0);
    
        for (size_t i = 0; i < _n[n]; ++i)
            coeffs[n][i] *= 2.0/(_bi[n]-_ai[n]);
    }

    std::vector<size_t> N_minus_one = _n;

    for (auto it = N_minus_one.begin(); it != N_minus_one.end(); ++it)
        --(*it);

    return Polynomial(N_minus_one,_ai,_bi,coeffs);
}


template<typename T>
inline Polynomial<T> Polynomial<T>::integral() const
{
    std::vector<std::vector<T>> coeffs(_n_blocks);

    T left_value = 0.0;

    for (size_t n = 0; n < _n_blocks; ++n)
    {
        coeffs[n] = std::vector<T>(_n[n]+2,0.0);
        // integral of L0 = L1
        coeffs[n][1] = _coeffs[n][0];
        for (size_t i = 1; i <= _n[n]; ++i)
        {
            coeffs[n][i+1] += _coeffs[n][i]/(2.0*i+1.0);
            coeffs[n][i-1] -= _coeffs[n][i]/(2.0*i+1.0);
        }

        for (size_t i = 0; i < _n[n]+2; ++i)
            coeffs[n][i] *= 0.5*(_bi[n]-_ai[n]);

        // Set the L0 coefficient to match left_value at p(-1)
        std::vector<scalar> Lk = legendre_polynomials(_n[n]+1,-1.0);

        T result = 0.0;

        for (size_t i = 1; i <= _n[n]+1; ++i)
            result += coeffs[n][i]*Lk[i];

        coeffs[n][0] = -result;

        // Update left_value to p(1.0) for the next block
        Lk = legendre_polynomials(_n[n]+1,1.0);

        left_value = 0.0;

        for (size_t i = 0; i <= _n[n]+1; ++i)
            left_value += coeffs[n][i]*Lk[i];
    }

    std::vector<size_t> N_plus_one = _n;

    for (auto it = N_plus_one.begin(); it != N_plus_one.end(); ++it)
        ++(*it);

    return Polynomial(N_plus_one,_ai,_bi,coeffs);
}


template<typename T>
constexpr bool Polynomial<T>::is_zero() const
{
    return is_zero([](const auto& coeff_i) { return std::abs(coeff_i) < 1.0e-12; });
}

template<typename T>
template<typename UnaryOperation>
constexpr bool Polynomial<T>::is_zero(UnaryOperation unary_op) const
{
    for (const auto& coeff : _coeffs)
    {
        if (!std::accumulate(coeff.cbegin(), coeff.cend(), true,
            [&](const bool& is_zero, const auto& coeff_i) { return is_zero && unary_op(coeff_i); }))
            return false;
    }

    return true;
}



template<typename T>
inline std::vector<T> Polynomial<T>::compute_coefficients(const typename std::vector<scalar>::const_iterator& ix0, 
        const typename std::vector<T>::const_iterator& iy0, const size_t N)
{
    std::vector<T> coeffs(N+1,T());


    // support T as std::vector<U>
    if constexpr (is_specialization_v<T, std::vector>) {
        std::for_each(coeffs.begin(), coeffs.end(), [&](auto& coeff) { coeff.resize(iy0->size(), typename T::value_type{ 0 }); });
    }

    // Get quadrature nodes and weights
    auto [xj, wj] = gauss_legendre_lobatto_nodes_and_weights(N);

    // Get points in local coordinates xi in [-1,1]
    std::vector<scalar> xi(N+1);
    
    for (size_t i = 0; i <= N; ++i)
        xi[i] = 2.0*(*(ix0+i)-(*ix0))/(*(ix0+N)-(*ix0)) - 1.0;

    // Get y evaluated at the interpolation points
    std::vector<T> yj(N+1,T());
    if constexpr (is_specialization_v<T, std::vector>) {
        std::for_each(yj.begin(), yj.end(), [&](auto& coeff) { coeff.resize(iy0->size(), typename T::value_type{ 0 }); });
    }

    for (size_t i = 0; i <= N; ++i)
        for (size_t j = 0; j <= N; ++j)
            yj[i] += *(iy0+j)*lagrange_polynomial(xj[i],j,N,xi);

    // Get the coefficients 
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j <= N; ++j)
            coeffs[i] += yj[j]*(wj[j]*legendre_polynomial(i,xj[j]));
 
        // |Ln|2 = 2/(2n+1) 
        coeffs[i] *= (i + 0.5);
    }

    // The last Ln norm is not exact and must be computed with the quadratures to get a good interpolant
    scalar Ln_norm = 0.0;
    for (size_t j = 0; j <= N; ++j)
    {
        const scalar Ln = legendre_polynomial(N,xj[j]);
        coeffs[N] += yj[j]*(wj[j]*Ln);
        Ln_norm += wj[j]*Ln*Ln;
    }

    coeffs[N] *= (1.0/Ln_norm);
 
    return coeffs;
}


template<typename T>
inline std::vector<T> Polynomial<T>::compute_coefficients(const std::vector<T>& y0, const size_t N)
{
    std::vector<T> coeffs(N+1,T());

    // Get quadrature nodes and weights
    auto [xj, wj] = gauss_legendre_lobatto_nodes_and_weights(N);

    // Get the coefficients 
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j <= N; ++j)
            coeffs[i] += y0[j]*(wj[j]*legendre_polynomial(i,xj[j]));
 
        // |Ln|2 = 2/(2n+1) 
        coeffs[i] *= (i + 0.5);
    }

    // The last Ln norm is not exact and must be computed with the quadratures to get a good interpolant
    scalar Ln_norm = 0.0;
    for (size_t j = 0; j <= N; ++j)
    {
        const scalar Ln = legendre_polynomial(N,xj[j]);
        coeffs[N] += y0[j]*(wj[j]*Ln);
        Ln_norm += wj[j]*Ln*Ln;
    }

    coeffs[N] *= (1.0/Ln_norm);
 
    return coeffs;
}


template<typename T>
inline std::vector<T> Polynomial<T>::compute_coefficients(const Polynomial& p, const size_t N)
{
    // Get quadrature nodes and weights
    auto xj = std::get<0>(gauss_legendre_lobatto_nodes_and_weights(N));

    const scalar& a = p.get_left_bound();
    const scalar& b = p.get_right_bound();

    std::vector<T> y0(N+1,T());

    for (size_t i = 0; i <= N; ++i)
        y0[i] = p(a + 0.5*(b-a)*(xj[i]+1.0));

    return compute_coefficients(y0, N);
}


}
#endif
