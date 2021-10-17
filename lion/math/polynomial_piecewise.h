#ifndef __POLYNOMIAL_PIECEWISE_H__
#define __POLYNOMIAL_PIECEWISE_H__

#include <vector>
#include <array>
#include <iostream>
#include "polynomial.h"
#include "matrix_extensions.h"

template<size_t Np>
class Polynomial_piecewise
{
 public:
    constexpr static size_t Nhalf = Np;
    constexpr static size_t Norder = 2*Np-1;
    constexpr static size_t Npoints = 2*Np;

    //! Constructor from (x,y)
    Polynomial_piecewise(const std::vector<double>& x, const std::vector<double>& y);

            

    
 private:
    size_t _n_subpoly;                               //! number of sub-polynomials
    std::vector<double> _x;                          //! x-coordinates
    std::vector<std::array<double,Npoints>> _y;      //! y-coordinates of each piecewise sub-polynomial
    std::vector<std::array<double,Npoints>> _coeffs; //! Legendre coefficients of each sub-polynomial

};

#include "polynomial_piecewise.hpp"

#endif
