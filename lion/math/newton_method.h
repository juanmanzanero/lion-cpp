#ifndef __NEWTON_METHOD_H__
#define __NEWTON_METHOD_H__

#include "lion/math/matrix_extensions.h"

namespace lioncpp{

//!     An implementation of Newton's method to solve nonlinear systems
//!     ---------------------------------------------------------------
//!
//! @param F: the type of the functor
template<typename F, size_t N>
class Newton_method
{
 public:
    using Jacobian_type = std::array<std::array<timeseries,N>,N>;
    using Solution_type = std::array<timeseries,N>;

    //! An auxiliary struct to hold the solutions
    struct Newton_method_solution
    {
        Solution_type x;     //! Solution vector
        Solution_type err;   //! Error vector
        double max_err;                 //! Maximum error
        size_t n_iter;                  //! Number of iterations
        bool valid_solution;            //! Whether solution is acceptable or not
    };

    //! Solve
    //! @param[in] f: Functor of the form y=f(x)
    //! @param[in] x: initial guess
    //! @return solution struct
    static Newton_method_solution solve(F& f, const std::array<timeseries,N>& x0);

 private:
    static Solution_type iteration(F& f, const Solution_type& xi, const Solution_type& yi);

    inline static double _relaxation = 0.3;
    inline static size_t _n_max_iter = 1000;
    inline static double _err_tol = 1.0e-6;
    inline static double _eps_jac = 1.0e-4;
};


template<typename F, size_t N>
typename Newton_method<F,N>::Newton_method_solution Newton_method<F,N>::solve(F& f, const std::array<timeseries,N>& x0)
{
    auto x(x0);
    std::array<timeseries,N> err;
    timeseries max_err = 0.0;
    bool valid_solution = false;
    size_t iter;

    for (iter = 0; iter < _n_max_iter; ++iter)
    {
        err = f(x);
        auto abs_err = ext::abs(err);
        max_err = *std::max_element(abs_err.cbegin(), abs_err.cend());

        if ( max_err < _err_tol )
        {
            valid_solution = true; 
            break;
        }
         
        x = iteration(f, x, err);
    }

    return {x, err, max_err, iter, valid_solution};
}


template<typename F, size_t N>
std::array<timeseries,N> Newton_method<F,N>::iteration(F& f, const std::array<timeseries,N>& xi, const std::array<timeseries,N>& yi)
{
    // Compute the Jacobian
    Jacobian_type J;

    for (size_t icol = 0; icol < N; ++icol)
    {
        auto xi_eps = xi;
        xi_eps[icol] += _eps_jac;

        auto yi_peps = f(xi_eps);
        
        xi_eps[icol] -= 2.0*_eps_jac;
        auto yi_meps = f(xi_eps);

        auto Jcol = (yi_peps - yi_meps)*(0.5/_eps_jac);

        for (size_t irow = 0; irow < N; ++irow)
            J[irow][icol] = Jcol[irow];
    }

    // Perform the iteration
    return xi + _relaxation*linsolve(N,J,-yi);
}

}

#endif
