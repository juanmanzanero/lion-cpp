#ifndef __SOLVE_NONLINEAR_SYSTEM_H__
#define __SOLVE_NONLINEAR_SYSTEM_H__

#include "lion/foundation/types.h"
#include <math.h>

struct Solve_nonlinear_system_options
{
    std::string hessian_approximation = "limited-memory";
    int print_level = 0;
    std::string mu_strategy = "monotone";
    std::string nlp_scaling_method = "none";    
    double constr_viol_tol = 1.0e-09;
    double acceptable_tol = 1.0e-6;
    double tol = 1.0;
    std::string derivative_test = "none";
    std::string jacobian_approximation = "finite-difference-values";
    bool throw_if_fail = true;
    double findiff_perturbation = 1.0e-6;
};

template<typename C>
class Solve_nonlinear_system
{
 public:
    struct Nonlinear_system_solution
    {
        bool solved;
        std::vector<scalar> x;
    };

    static Nonlinear_system_solution solve(const size_t n, const size_t nc, const std::vector<scalar>& x0, C& c, const std::vector<scalar>& x_lb, 
                                        const std::vector<scalar>& x_ub, const std::vector<scalar>& c_lb, 
                                        const std::vector<scalar>& c_ub, const Solve_nonlinear_system_options& = {});

 private:
    class Fitness_function
    {
     public:
        using argument_type = typename C::argument_type;   
        
        double operator()(const argument_type& x) { return 1.0; }
    };
};

#include "solve_nonlinear_system.hpp" 

#endif
