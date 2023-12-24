#ifndef LIONCPP_MATH_IPOPT_OPTIMIZE_FINITE_DIFFERENCES_H
#define LIONCPP_MATH_IPOPT_OPTIMIZE_FINITE_DIFFERENCES_H


#include <math.h>
#include <string>
#include <vector>
#include <iomanip>

#include "coin-or/IpSolveStatistics.hpp"
#include "lion/thirdparty/include/logger.hpp"

#include "lion/foundation/types.h"
#include "lion/foundation/utils.hpp"

#include "matrix_extensions.h"
#include "detail/ipopt_nlp_finite_differences.h"


namespace lioncpp {

struct Ipopt_optimize_NLP_finite_differences
{
    struct options
    {
        std::string linear_solver = "mumps";
        std::string nlp_scaling_method = "none";
        std::string gradient_approximation = "exact";         // if == "exact", it uses our finite differences; if == "finite-difference-values", it uses ipopt's ones
        std::string jacobian_approximation = "exact";         // if == "exact", it uses our finite differences; if == "finite-difference-values", it uses ipopt's ones
        std::string hessian_approximation = "limited-memory"; // if == "limited-memory", it uses ipopt's quasi-newton second derivatives; if == "exact", it uses our finite-difference centered hessian
        std::string mu_strategy = "monotone";

        bool strict_bounds = true;

        int max_iter = 3000;
        scalar tol = 1e-8;
        scalar constr_viol_tol = 1.e-09; // ipopt's default is 1e-4;
        scalar compl_inf_tol = 1e-4;
        scalar dual_inf_tol = 1.;

        int acceptable_iter = 100; // ipopt's default is 15;
        scalar acceptable_tol = 1e-20; // ipopt's default is 1e-6;
        scalar acceptable_constr_viol_tol = 1e-2;
        scalar acceptable_compl_inf_tol = 1e-2;
        scalar acceptable_dual_inf_tol = 1e10;
        scalar acceptable_obj_change_tol = 1e20;

        scalar ipopt_findiff_perturbation = 1.0e-6;
        scalar exact_jac_findiff_perturbation = std::sqrt(2.0e-16);
        scalar exact_hess_findiff_perturbation = std::cbrt(2.0e-16);
        bool exact_central_finite_differences = false; // if "false", takes forward finite-differences for the "exact" gradient and jacobian aproximations (the "exact" hessian approximation is always central)

        int print_level = 0;
        std::string derivative_test = "none"; // "none", "first-order" (checks 1st order), "second-order" (checks 1st and 2nd order). Requires "print_level = 5" to show the result of the test.
        bool throw_if_fail = true;
    };


    struct result
    {
        bool solved;
        std::vector<scalar> x;
        scalar f;
        scalar cons_err;
    };


    template<typename ArgumentType>
    struct no_fitness_T
    {
        using argument_type = ArgumentType;

        scalar operator()(const argument_type &) const
        {
            return 0.;
        }
    };

    using no_fitness = no_fitness_T<std::vector<scalar> >;


    template<typename ArgumentType>
    struct no_constraints_T
    {
        using argument_type = ArgumentType;

        const argument_type& operator()(const argument_type &) const
        {
            static const auto emptyarg = argument_type{};
            return emptyarg;
        }
    };

    using no_constraints = no_constraints_T<std::vector<scalar> >;


    template<typename ConstraintsType>
    static result solve_nonlinear_system(const size_t n, const size_t nc, const std::vector<scalar> &x0,
                                         ConstraintsType &&c,
                                         const std::vector<scalar> &x_lb, const std::vector<scalar> &x_ub,
                                         const std::vector<scalar> &c_lb, const std::vector<scalar> &c_ub,
                                         const options &opts)
    {
        //
        // Solves a nonlinear system of "nc" equations (a.k.a "constraints")
        // with "n" variables.
        //

        return optimize(n, nc, x0,
                        no_fitness_T<typename std::decay_t<ConstraintsType>::argument_type>{}, c,
                        x_lb, x_ub,
                        c_lb, c_ub,
                        opts);
    }


    template<typename FitnessType, typename ConstraintsType>
    static result optimize(const size_t n, const size_t nc, const std::vector<scalar>& x0,
                           FitnessType &&f, ConstraintsType &&c,
                           const std::vector<scalar>& x_lb, const std::vector<scalar>& x_ub,
                           const std::vector<scalar>& c_lb, const std::vector<scalar>& c_ub,
                           const options &opts)
    {
        //
        // Optimizes a nonlinear problem with "n" variables and "nc" constraints.
        //

        // validate inputs
        constexpr auto without_fitness = is_specialization_v<std::decay_t<FitnessType>, no_fitness_T>;
        constexpr auto without_constraints = is_specialization_v<std::decay_t<ConstraintsType>, no_constraints_T>;

        static_assert(!(without_fitness && without_constraints),
            "Ipopt_optimize_NLP_finite_differences::solve: will do nothing"
            " with no_fitness and no_constraints!");

        if constexpr (without_constraints) {
            if (nc != 0u || !c_lb.empty() || !c_ub.empty()) {
                throw std::runtime_error(
                    "Ipopt_optimize_NLP_finite_differences::solve:"
                    " inconsistent constraints.");
            }
        }

        assert(x0.size() == n);
        assert(x_lb.size() == n);
        assert(x_ub.size() == n);
        assert(c_lb.size() == nc);
        assert(c_ub.size() == nc);

        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        Ipopt::ApplicationReturnStatus status = app->Initialize();

        if (status != Ipopt::Solve_Succeeded) {
            throw lion_exception("Ipopt: error during initialization");
        }


        // configure ipopt's options
        app->Options()->SetStringValue("sb", "yes");
        app->Options()->SetStringValue("linear_solver", opts.linear_solver);
        if (opts.linear_solver == "ma97") {
            app->Options()->SetStringValue("ma97_order", "metis");
            app->Options()->SetStringValue("ma97_scaling", "mc77");
        }
        else if (opts.linear_solver == "ma27" || opts.linear_solver == "ma57") {
            app->Options()->SetStringValue("linear_system_scaling", "mc19");
        }

        app->Options()->SetStringValue("nlp_scaling_method", opts.nlp_scaling_method);

        app->Options()->SetStringValue("gradient_approximation", opts.gradient_approximation);
        app->Options()->SetStringValue("jacobian_approximation", opts.jacobian_approximation);
        app->Options()->SetStringValue("hessian_approximation", opts.hessian_approximation);

        app->Options()->SetStringValue("mu_strategy", opts.mu_strategy);
        if (opts.mu_strategy == "adaptive") {
            app->Options()->SetStringValue("adaptive_mu_globalization", "never-monotone-mode");
            app->Options()->SetStringValue("mu_oracle", "probing");
            app->Options()->SetStringValue("fixed_mu_oracle", "probing");
        }

        if (opts.strict_bounds) {
            app->Options()->SetNumericValue("bound_relax_factor", 0.);
            app->Options()->SetStringValue("honor_original_bounds", "yes");
        }

        app->Options()->SetIntegerValue("max_iter", opts.max_iter);
        app->Options()->SetNumericValue("tol", opts.tol);
        app->Options()->SetNumericValue("constr_viol_tol", opts.constr_viol_tol);
        app->Options()->SetNumericValue("compl_inf_tol", opts.compl_inf_tol);
        app->Options()->SetNumericValue("dual_inf_tol", opts.dual_inf_tol);

        app->Options()->SetIntegerValue("acceptable_iter", opts.acceptable_iter);
        app->Options()->SetNumericValue("acceptable_tol", opts.acceptable_tol);
        app->Options()->SetNumericValue("acceptable_constr_viol_tol", opts.acceptable_constr_viol_tol);
        app->Options()->SetNumericValue("acceptable_compl_inf_tol", opts.acceptable_compl_inf_tol);
        app->Options()->SetNumericValue("acceptable_dual_inf_tol", opts.acceptable_dual_inf_tol);
        app->Options()->SetNumericValue("acceptable_obj_change_tol", opts.acceptable_obj_change_tol);

        app->Options()->SetNumericValue("findiff_perturbation", opts.ipopt_findiff_perturbation);

        app->Options()->SetStringValue("derivative_test", opts.derivative_test);

        if (opts.derivative_test != "none" && opts.print_level < 5) {
            std::cout << "Ipopt_optimize_NLP_finite_differences::optimize:"
                         " WARNING, the selected print level ("
                      << opts.print_level
                      << ") will be overriden so that the derivative test's"
                         " results can be displayed..."
                      << std::endl;

            app->Options()->SetIntegerValue("print_level", 5);
        }
        else {
            app->Options()->SetIntegerValue("print_level", opts.print_level);
        }


        // create our TNLP
        Ipopt::SmartPtr<Ipopt::TNLP> nlp =
            new lioncpp::detail::Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>(
                n, nc, x0, std::forward<FitnessType>(f), std::forward<ConstraintsType>(c), x_lb, x_ub, c_lb, c_ub,
                opts.exact_jac_findiff_perturbation, opts.exact_hess_findiff_perturbation,
                opts.exact_central_finite_differences);

        status = app->OptimizeTNLP(nlp);

        Ipopt::Number dual_inf, constr_viol, complementarity, kkt_error, varbounds_viol;
        app->Statistics()->Infeasibilities(dual_inf, constr_viol, varbounds_viol, complementarity, kkt_error);

        bool succeed = (status == Ipopt::Solve_Succeeded) || (status == Ipopt::Solved_To_Acceptable_Level);

        if (!succeed && opts.throw_if_fail) {
            throw lion_exception("Ipopt did not succeed");
        }

        // Get solution
        const std::vector<scalar> x =
            static_cast<detail::Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>* >(GetRawPtr(nlp))->x();

        // Check that all variables are within bounds
        for (size_t i = 0; i < n; ++i)
        {
            if (x[i] < (x_lb[i] - 10.0 * opts.constr_viol_tol)) { out(2) << "1" << std::endl; succeed = false; }
            if (x[i] > (x_ub[i] + 10.0 * opts.constr_viol_tol)) { out(2) << "2" << std::endl; succeed = false; }
        }

        // Check that constraints are satisfied
        using constraints_argument_type = typename std::decay_t<ConstraintsType>::argument_type;

        constraints_argument_type x_c;
        if constexpr (is_specialization_v<constraints_argument_type, std::vector>) {
            x_c.resize(n);
        }

        for (size_t i = 0; i < n; ++i)
            x_c[i] = x[i];

        auto constraints = c(x_c);

        for (size_t i = 0; i < nc; ++i)
        {
            if (constraints[i] < (c_lb[i] - 10.0 * opts.constr_viol_tol)) { out(2) << "3" << std::endl; succeed = false; }
            if (constraints[i] > (c_ub[i] + 10.0 * opts.constr_viol_tol)) { out(2) << "4" << std::endl; succeed = false; }
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

            if (opts.throw_if_fail) {
                throw lion_exception("Even though Ipopt succeeded, the constraints are not satisfied to the desired tolerance");
            }
        }


        return
        {
            succeed,
            x,
            app->Statistics()->FinalObjective(),
            constr_viol
        };
    }
};

} // end namespace lioncpp

#endif
