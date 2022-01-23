#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <numeric>

template<class T>
inline Polynomial<T>::Polynomial()
: _n_blocks(0),
  _n(0),
  _a(0.0),
  _b(0.0),
  _ai(),
  _bi(),
  _coeffs()
{}


template<class T>
inline Polynomial<T>::Polynomial(const std::vector<scalar>& x0, const std::vector<T>& y0, const size_t N, bool smooth)
: _n_blocks(smooth ? 1 : (x0.size()-1)/(N) + ( ((x0.size()-1 % (N))!=0) ? 1.0 : 0.0 )),
  _n(_n_blocks,N),
  _a(x0.front()),
  _b(x0.back()),
  _ai(_n_blocks),
  _bi(_n_blocks), 
  _coeffs(_n_blocks)
{
    if ( x0.size() != y0.size() )
        throw std::runtime_error("x0 and y0 must have the same length");

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
        Polynomial p(x0, y0, N_smooth, false);

        _coeffs.front() = compute_coefficients(p, N);
    }
}


template<class T>
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


template<class T>
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
        throw std::runtime_error("L and y0 must have the same length");

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


template<class T>
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
        throw std::runtime_error("x0 and y0 must have the same length");

    if ( x0.size() == 0 )
    {
        *this = Polynomial();
        return;
    }

    // Check x0 consistency
    for (auto it = x0.cbegin(); it != x0.cend()-1; ++it)
        if ( (*(it+1)).back() < ((*it).back() + eps) )
            throw std::runtime_error("blocks are not contiguous");

    // Check y0 consistency if continuous
    if ( continuous ) 
        for (auto it = y0.cbegin(); it != y0.cend()-1; ++it)
            if ( std::abs( (*it).back() - ((*(it+1)).front()) ) > eps )
                throw std::runtime_error("blocks are not continuous, and continuous = true");


    std::vector<Polynomial> p(x0.size());

    for (size_t i = 0; i < _n_blocks; ++i)
        p[i] = Polynomial(x0[i], y0[i], smooth);

    *this = Polynomial(p);
}


template<class T>
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


template<class T>
inline T Polynomial<T>::operator[](scalar x) const
{
    if ( _n_blocks == 0 )
        return T();

    size_t block = 0;

    while(x > _bi[block])
        ++block;

    scalar xi = 2.0*(x-_ai[block])/(_bi[block]-_ai[block]) - 1.0;

    std::vector<scalar> Lk = legendre_polynomials(_n[block],xi);

    T result = T();

    for (size_t i = 0; i <= _n[block]; ++i)
        result += _coeffs[block][i]*Lk[i];

    return result;
}


template<class T>
inline T Polynomial<T>::operator()(scalar x) const
{
    if ( _n_blocks == 0 )
        return T();

    // Check bounds
    if ( (x < _a) || (x > _b) )
    {
        std::ostringstream s_out;
        s_out << "x is out of bounds" << std::endl;
        s_out << "  a: " << _a << std::endl;
        s_out << "  x: " <<  x << std::endl;
        s_out << "  b: " << _b ;
        throw std::runtime_error(s_out.str());
    }

    // Call operator[]
    return operator[](x);
}


template<class T>
inline Polynomial<T> Polynomial<T>::derivative() const
{
    std::vector<std::vector<T>> coeffs(_n_blocks);

    for (size_t n = 0; n < _n_blocks; ++n)
    {
        coeffs[n] = std::vector<T>(_n[n],T());
        for (size_t i = 0; i <= _n[n]; ++i)
            for (int j = i-1 ; j >= 0; j -= 2)
                coeffs[n][j] += _coeffs[n][i]*(2.0*(j+1.0)-1.0);
    
        for (size_t i = 0; i < _n[n]; ++i)
            coeffs[n][i] *= 2.0/(_bi[n]-_ai[n]);
    }

    std::vector<size_t> N_minus_one = _n;

    for (auto it = N_minus_one.begin(); it != N_minus_one.end(); ++it)
        --(*it);

    return Polynomial(N_minus_one,_ai,_bi,coeffs);
}


template<class T>
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


template<class T>
inline std::vector<T> Polynomial<T>::compute_coefficients(const typename std::vector<scalar>::const_iterator& ix0, 
        const typename std::vector<T>::const_iterator& iy0, const size_t N)
{
    std::vector<T> coeffs(N+1,T());

    // Get quadrature nodes and weights
    auto [xj, wj] = gauss_legendre_lobatto_nodes_and_weights(N);

    // Get points in local coordinates xi in [-1,1]
    std::vector<scalar> xi(N+1);
    
    for (size_t i = 0; i <= N; ++i)
        xi[i] = 2.0*(*(ix0+i)-(*ix0))/(*(ix0+N)-(*ix0)) - 1.0;

    // Get y evaluated at the interpolation points
    std::vector<T> yj(N+1,T());

    for (size_t i = 0; i <= N; ++i)
        for (size_t j = 0; j <= N; ++j)
            yj[i] += *(iy0+j)*lagrange_polynomial(xj[i],j,N,xi);

    // Get the coefficients 
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j <= N; ++j)
            coeffs[i] += wj[j]*legendre_polynomial(i,xj[j])*yj[j];
 
        // |Ln|2 = 2/(2n+1) 
        coeffs[i] *= (i + 0.5);
    }

    // The last Ln norm is not exact and must be computed with the quadratures to get a good interpolant
    scalar Ln_norm = 0.0;
    for (size_t j = 0; j <= N; ++j)
    {
        const scalar Ln = legendre_polynomial(N,xj[j]);
        coeffs[N] += wj[j]*Ln*yj[j];
        Ln_norm += wj[j]*Ln*Ln;
    }

    coeffs[N] *= (1.0/Ln_norm);
 
    return coeffs;
}


template<class T>
inline std::vector<T> Polynomial<T>::compute_coefficients(const std::vector<T>& y0, const size_t N)
{
    std::vector<T> coeffs(N+1,T());

    // Get quadrature nodes and weights
    auto [xj, wj] = gauss_legendre_lobatto_nodes_and_weights(N);

    // Get the coefficients 
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j <= N; ++j)
            coeffs[i] += wj[j]*legendre_polynomial(i,xj[j])*y0[j];
 
        // |Ln|2 = 2/(2n+1) 
        coeffs[i] *= (i + 0.5);
    }

    // The last Ln norm is not exact and must be computed with the quadratures to get a good interpolant
    scalar Ln_norm = 0.0;
    for (size_t j = 0; j <= N; ++j)
    {
        const scalar Ln = legendre_polynomial(N,xj[j]);
        coeffs[N] += wj[j]*Ln*y0[j];
        Ln_norm += wj[j]*Ln*Ln;
    }

    coeffs[N] *= (1.0/Ln_norm);
 
    return coeffs;
}


template<class T>
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
#endif
