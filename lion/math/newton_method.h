#ifndef __NEWTON_METHOD_H__
#define __NEWTON_METHOD_H__

#include "lion/math/matrix_extensions.h"

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

#include "newton_method.hpp"

#endif
