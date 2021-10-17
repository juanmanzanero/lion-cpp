#ifndef __POLYNOMIAL_PIECEWISE_HPP__
#define __POLYNOMIAL_PIECEWISE_HPP__

template<size_t Np>
Polynomial_piecewise<Np>::Polynomial_piecewise(const std::vector<double>& x, const std::vector<double>& y)
: _n_subpoly(std::max(1u,x.size() - 2*Np + 1)),
  _x(x),
  _y(_n_subpoly),
  _coeffs(x.size()-1)
{
    PRINTVARIABLE(JMT,_n_subpoly);
    PRINTVARIABLE(JMT, _x);
    // Check consistency
    assert(x.size() == y.size());

    // Fill y
    for (size_t i = 0; i < _n_subpoly; ++i )
        for (size_t n = 0; n < Npoints; ++n)
            _y[i][n] = y[i + n];

    PRINTVARIABLE(JMT, _y);

}

#endif
