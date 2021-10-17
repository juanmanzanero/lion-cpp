#ifndef __OPTIMISE_HPP__
#define __OPTIMISE_HPP__

#include "ipopt_handler.h"
#include "lion/thirdparty/include/coin-or/IpSolveStatistics.hpp"
#include "lion/math/matrix_extensions.h"
#include <math.h>
#include <iomanip>
#include "lion/foundation/utils.h"
#include "lion/thirdparty/include/logger.hpp"


template<typename F, typename C>
inline typename Optimise<F,C>::Optimisation_result Optimise<F,C>::optimise(const size_t n, const size_t nc,
    const std::vector<scalar>& x0, F& f, C& c, 
    const std::vector<scalar>& x_lb, const std::vector<scalar>& x_ub, const std::vector<scalar>& c_lb, 
    const std::vector<scalar>& c_ub, const Optimise_options& options)
{
    assert(x0.size() == n);
    assert(x_lb.size() == n);
    assert(x_ub.size() == n);
    assert(c_lb.size() == nc);
    assert(c_ub.size() == nc);

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    Ipopt::ApplicationReturnStatus status = app->Initialize();

    if ( status != Ipopt::Solve_Succeeded )
        throw std::runtime_error("Ipopt: error during initialization");

    // Configure options
    app->Options()->SetStringValue ("hessian_approximation", options.hessian_approximation);
    app->Options()->SetIntegerValue("print_level"          , options.print_level);
    app->Options()->SetStringValue ("mu_strategy"          , options.mu_strategy);
    app->Options()->SetStringValue ("nlp_scaling_method"   , options.nlp_scaling_method);
    app->Options()->SetNumericValue("constr_viol_tol"      , options.constr_viol_tol);
    app->Options()->SetNumericValue("acceptable_tol"       , options.acceptable_tol);
    app->Options()->SetNumericValue("tol"                  , options.tol);
    app->Options()->SetStringValue ("derivative_test"      , options.derivative_test);
    app->Options()->SetStringValue ("jacobian_approximation"   ,options.jacobian_approximation);
    app->Options()->SetNumericValue("findiff_perturbation"     ,options.findiff_perturbation);
    

    Ipopt::SmartPtr<Ipopt::TNLP> nlp = new Ipopt_handler<F,C>(n, nc, x0, f, c, x_lb, x_ub, c_lb, c_ub);
    status = app->OptimizeTNLP(nlp);
    
    Ipopt::Number dual_inf, constr_viol, complementarity, kkt_error, varbounds_viol;
    app->Statistics()->Infeasibilities(dual_inf, constr_viol, varbounds_viol, complementarity, kkt_error);

    bool succeed = (status == Ipopt::Solve_Succeeded) || (status == Ipopt::Solved_To_Acceptable_Level);

    if ( !succeed ) 
        throw std::runtime_error("Ipopt did not succeed");

    // Get solution
    const std::vector<scalar> x = static_cast<Ipopt_handler<F,C>*>(GetRawPtr(nlp))->x();

    // Check that all variables are within bounds
    for (size_t i = 0; i < n; ++i)
    {
        if ( x[i] < (x_lb[i] - 10.0*options.constr_viol_tol) ) {out(2) << "1" << std::endl; succeed = false; }
        if ( x[i] > (x_ub[i] + 10.0*options.constr_viol_tol) ) {out(2) << "2" << std::endl; succeed = false; }
    }

    // Check that constraints are satisfied
    typename C::argument_type x_c;
    if constexpr (std::is_same<typename C::argument_type, std::vector<double>>::value)
        x_c = std::vector<double>(n);    

    for (size_t i = 0; i < n; ++i)
        x_c[i] = x[i];

    auto constraints = c(x_c);

    for (size_t i = 0; i < nc; ++i)
    {
        if ( constraints[i] < (c_lb[i] - 10.0*options.constr_viol_tol) ) {out(2) << "3" << std::endl; succeed = false;}
        if ( constraints[i] > (c_ub[i] + 10.0*options.constr_viol_tol) ) {out(2) << "4" << std::endl; succeed = false;}
    }

    if ( !succeed )
    {
        out(2) << std::setprecision(16);
        for (size_t i = 0; i < n; ++i)
        {
            PRINTVARIABLE(JMT,x[i]);
        }
        for (size_t i = 0; i < nc; ++i)
        {
            PRINTVARIABLE(JMT,constraints[i]);
        }

        throw std::runtime_error("Even though Ipopt succeeded, the constraints are not satisfied to the desired tolerance");
    }


    return 
    {
        succeed, 
        x,
        app->Statistics()->FinalObjective(), 
        constr_viol 
    };
}

#endif
