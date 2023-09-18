#ifndef LIONCPP_MATH_IPOPT_OPTIMIZE_FINITE_DIFFERENCES_H
#define LIONCPP_MATH_IPOPT_OPTIMIZE_FINITE_DIFFERENCES_H

#include <math.h>
#include <string>
#include <vector>
#include <iomanip>

#include "lion/thirdparty/include/logger.hpp"
#include "lion/thirdparty/include/coin-or/IpSolveStatistics.hpp"

#include "lion/foundation/types.h"
#include "lion/foundation/utils.hpp"

#include "matrix_extensions.h"
#include "detail/ipopt_nlp_finite_differences.h"


namespace lioncpp
{

    struct Ipopt_optimize_NLP_finite_differences_options
    {
        int print_level = 0;
        std::string linear_solver = "mumps";
        std::string mu_strategy = "monotone";
        std::string nlp_scaling_method = "none";
        double constr_viol_tol = 1.0e-09;
        double acceptable_tol = 1.0e-6;
        double tol = 1.0e-8;
        std::string derivative_test = "none";
        std::string gradient_approximation = "exact";
        std::string jacobian_approximation = "finite-difference-values";
        std::string hessian_approximation = "limited-memory";
        double ipopt_findiff_perturbation = 1.0e-6;
        double exact_jac_findiff_perturbation = std::sqrt(2.0e-16);
        double exact_hess_findiff_perturbation = std::cbrt(2.0e-16);
    };


    template<typename Fitness_type, typename Constraints_type>
    class Ipopt_optimize_NLP_finite_differences
    {
    public:

        struct Result
        {
            bool solved;
            std::vector<scalar> x;
            scalar f;
            scalar cons_err;
        };

        static auto optimize(const size_t n, const size_t nc, const std::vector<scalar>& x0, Fitness_type& f, Constraints_type& c, const std::vector<scalar>& x_lb,
            const std::vector<scalar>& x_ub, const std::vector<scalar>& c_lb,
            const std::vector<scalar>& c_ub, const Ipopt_optimize_NLP_finite_differences_options & = {}) -> Result;
    };


    template<typename Fitness_type, typename Constraints_type>
    inline auto Ipopt_optimize_NLP_finite_differences<Fitness_type, Constraints_type>::optimize(const size_t n, const size_t nc,
        const std::vector<scalar>& x0, Fitness_type& f, Constraints_type& c,
        const std::vector<scalar>& x_lb, const std::vector<scalar>& x_ub, const std::vector<scalar>& c_lb,
        const std::vector<scalar>& c_ub, const Ipopt_optimize_NLP_finite_differences_options& options) -> Result
    {
        assert(x0.size() == n);
        assert(x_lb.size() == n);
        assert(x_ub.size() == n);
        assert(c_lb.size() == nc);
        assert(c_ub.size() == nc);

        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        Ipopt::ApplicationReturnStatus status = app->Initialize();

        if (status != Ipopt::Solve_Succeeded)
            throw lion_exception("Ipopt: error during initialization");

        // Configure options
        app->Options()->SetStringValue("linear_solver", options.linear_solver);
        app->Options()->SetStringValue("gradient_approximation", options.gradient_approximation);
        app->Options()->SetStringValue("jacobian_approximation", options.jacobian_approximation);
        app->Options()->SetNumericValue("findiff_perturbation", options.ipopt_findiff_perturbation);
        app->Options()->SetStringValue("hessian_approximation", options.hessian_approximation);
        app->Options()->SetIntegerValue("print_level", options.print_level);
        app->Options()->SetStringValue("mu_strategy", options.mu_strategy);
        app->Options()->SetStringValue("nlp_scaling_method", options.nlp_scaling_method);
        app->Options()->SetNumericValue("constr_viol_tol", options.constr_viol_tol);
        app->Options()->SetNumericValue("acceptable_tol", options.acceptable_tol);
        app->Options()->SetNumericValue("tol", options.tol);
        app->Options()->SetStringValue("derivative_test", options.derivative_test);


        Ipopt::SmartPtr<Ipopt::TNLP> nlp =
            new lioncpp::detail::Ipopt_NLP_finite_differences<Fitness_type, Constraints_type>(
                n, nc, x0, f, c, x_lb, x_ub, c_lb, c_ub,
                options.exact_jac_findiff_perturbation, options.exact_hess_findiff_perturbation);

        status = app->OptimizeTNLP(nlp);

        Ipopt::Number dual_inf, constr_viol, complementarity, kkt_error, varbounds_viol;
        app->Statistics()->Infeasibilities(dual_inf, constr_viol, varbounds_viol, complementarity, kkt_error);

        bool succeed = (status == Ipopt::Solve_Succeeded) || (status == Ipopt::Solved_To_Acceptable_Level);

        if (!succeed)
            throw lion_exception("Ipopt did not succeed");

        // Get solution
        const std::vector<scalar> x = static_cast<detail::Ipopt_NLP_finite_differences<Fitness_type, Constraints_type>*>(GetRawPtr(nlp))->x();

        // Check that all variables are within bounds
        for (size_t i = 0; i < n; ++i)
        {
            if (x[i] < (x_lb[i] - 10.0 * options.constr_viol_tol)) { out(2) << "1" << std::endl; succeed = false; }
            if (x[i] > (x_ub[i] + 10.0 * options.constr_viol_tol)) { out(2) << "2" << std::endl; succeed = false; }
        }

        // Check that constraints are satisfied
        typename Constraints_type::argument_type x_c;
        if constexpr (std::is_same<typename Constraints_type::argument_type, std::vector<double>>::value)
            x_c = std::vector<double>(n);

        for (size_t i = 0; i < n; ++i)
            x_c[i] = x[i];

        auto constraints = c(x_c);

        for (size_t i = 0; i < nc; ++i)
        {
            if (constraints[i] < (c_lb[i] - 10.0 * options.constr_viol_tol)) { out(2) << "3" << std::endl; succeed = false; }
            if (constraints[i] > (c_ub[i] + 10.0 * options.constr_viol_tol)) { out(2) << "4" << std::endl; succeed = false; }
        }

        if (!succeed)
        {
            out(2) << std::setprecision(16);
            for (size_t i = 0; i < n; ++i)
            {
                std::cout << "x[" << i << "]: " << x[i] << std::endl;
            }
            for (size_t i = 0; i < nc; ++i)
            {
                std::cout << "constraints[" << i << "]: " << constraints[i] << std::endl;
            }

            throw lion_exception("Even though Ipopt succeeded, the constraints are not satisfied to the desired tolerance");
        }


        return
        {
            succeed,
            x,
            app->Statistics()->FinalObjective(),
            constr_viol
        };
    }
}

#endif
