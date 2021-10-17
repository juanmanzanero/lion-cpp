#ifndef __OPTIMISE_H__
#define __OPTIMISE_H__

#include "lion/foundation/types.h"
#include <math.h>
#include <string>
#include <vector>

struct Optimise_options
{
    std::string hessian_approximation = "limited-memory";
    int print_level = 0;
    std::string mu_strategy = "monotone";
    std::string nlp_scaling_method = "none";    
    double constr_viol_tol = 1.0e-09;
    double acceptable_tol = 1.0e-6;
    double tol = 1.0e-8; 
    std::string derivative_test = "none";
    std::string jacobian_approximation = "finite-difference-values";
    double findiff_perturbation = 1.0e-6;
};

template<typename F, typename C>
class Optimise
{
 public:
    struct Optimisation_result
    {
        bool solved;
        std::vector<scalar> x;
        scalar f;
        scalar cons_err;
    };

    static Optimisation_result optimise(const size_t n, const size_t nc, const std::vector<scalar>& x0, F& f, C& c, const std::vector<scalar>& x_lb, 
                                        const std::vector<scalar>& x_ub, const std::vector<scalar>& c_lb, 
                                        const std::vector<scalar>& c_ub, const Optimise_options& = {});
};

#include "optimise.hpp" 

#endif
