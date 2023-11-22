#ifndef LION_MATH_IPOPT_OPTIMIZE_NLP_CPPAD_H
#define LION_MATH_IPOPT_OPTIMIZE_NLP_CPPAD_H
#pragma once

#include <string>

#include "coin-or/IpSmartPtr.hpp"
#include "coin-or/IpIpoptApplication.hpp"
#include "cppad/cppad.hpp"

#include "lion/foundation/types.h"
#include "lion/foundation/lion_exception.h"
#include "lion/math/matrix_extensions.h"
#include "lion/math/optimization_result.h"
#include "lion/math/ipopt_nlp_complete_output.h"
#include "lion/math/detail/ipopt_nlp_cppad.h"


//
// Defines functions to optimize, via Ipopt + CppAD, NLP's represented by
// either:
//
//  a) Function objects that must define a member operator
//    "void operator()(ADvector &fg, const ADvector &x)", in which "x"
//    represents the optimization variables and "fg" is the resulting
//    concatenation of "[scalar_cost; vector_of_constraints]".
//
//  b) Instances of class "CppAD::ADFun<Base>", which have been
//    taped adequately so that their "Forward(0, x)" call returns
//    the said "[scalar_cost; vector_of_constraints]".
//

namespace lioncpp {

namespace detail::ipopt_optimize_nlp_cppad {

//
// A type_trait to deduce the "ADvector" type that "ipopt_optimize_nlp_cppad"
// will use to solve the NLP represented by the "FG_eval" functor. If "FG_eval"
// declares a member "ADvector" typedef, then it will use that one. When the
// "FG_eval" is a "CppAD::ADFun<Base>", it will use "std::vector<CppAD::AD<Base> >".
// Else, it will use "std::vector<CppAD::AD<scalar> >".
//

template<class FG_eval, typename = void>
struct ADvector_type
{
    using type = std::vector<CppAD::AD<scalar> >;
};

template<class FG_eval>
struct ADvector_type<FG_eval, std::void_t<typename FG_eval::ADvector> >
{
    using type = typename FG_eval::ADvector;
};

template<typename Base>
struct ADvector_type<CppAD::ADFun<Base> >
{
    using type = std::vector<CppAD::AD<Base> >;
};

template<class FG_eval>
using ADvector_type_t = typename ADvector_type<FG_eval>::type;

} // end namespace detail::ipopt_optimize_nlp_cppad


// fwd. decl. of a helper function to set several options from a string into an Ipopt Application
void set_ipopt_app_options_from_string(Ipopt::SmartPtr<Ipopt::IpoptApplication> &app,
                                       const std::string &options,
                                       bool &retape,
                                       bool &sparse_forward,
                                       bool &sparse_reverse);


template<typename Dvector,
         class FG_eval,
         typename ADvector = detail::ipopt_optimize_nlp_cppad::ADvector_type_t<FG_eval> >
void ipopt_optimize_nlp_cppad(const std::string &options,
                              const Dvector &xi,
                              const Dvector &xl,
                              const Dvector &xu,
                              const Dvector &gl,
                              const Dvector &gu,
                              FG_eval &fg_eval,
                              Optimization_result<Dvector> &solution)
{
    //
    // Optimizes an NLP represented by functor "fg_eval".
    //

    bool ok = true;

    CPPAD_ASSERT_KNOWN(
        xi.size() == xl.size() && xi.size() == xu.size(),
        "ipopt::solve: size of xi, xl, and xu are not all equal."
    );

    CPPAD_ASSERT_KNOWN(
        gl.size() == gu.size(),
        "ipopt::solve: size of gl and gu are not equal."
    );

    size_t nx = xi.size();
    size_t ng = gl.size();

    // Create an IpoptApplication
    using Ipopt::IpoptApplication;
    Ipopt::SmartPtr<IpoptApplication> app = new IpoptApplication();

    bool retape = false, sparse_forward = false, sparse_reverse = false;
    set_ipopt_app_options_from_string(app, options, retape, sparse_forward, sparse_reverse);

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    ok &= status == Ipopt::Solve_Succeeded;
    if (!ok)
    {
        solution.status = Optimization_result<Dvector>::status_type::unknown;
        return;
    }

    // Create an interface from Ipopt to this specific problem.
    size_t nf = 1;
    Ipopt::SmartPtr<Ipopt::TNLP> cppad_nlp =
        new detail::Ipopt_NLP_CppAD<Dvector, ADvector, FG_eval>(
            nf,
            nx,
            ng,
            xi,
            xl,
            xu,
            gl,
            gu,
            fg_eval,
            retape,
            sparse_forward,
            sparse_reverse,
            solution);

    // Run the IpoptApplication
    app->OptimizeTNLP(cppad_nlp);

    solution.iter_count = app->Statistics()->IterationCount();

    return;
}


template<typename Dvector,
         class FG_eval,
         typename ADvector = detail::ipopt_optimize_nlp_cppad::ADvector_type_t<FG_eval> >
void ipopt_optimize_nlp_cppad(const std::string &options,
                              const Dvector &xi,
                              const Dvector &xl,
                              const Dvector &xu,
                              const Dvector &gl,
                              const Dvector &gu,
                              const Dvector &lambda_i,
                              const Dvector &zl_i,
                              const Dvector &zu_i,
                              FG_eval &fg_eval,
                              Optimization_result<Dvector> &solution)
{
    //
    // Optimizes an NLP represented by functor "fg_eval",
    // with warm-start.
    //

    bool ok = true;

    CPPAD_ASSERT_KNOWN(
        xi.size() == xl.size() && xi.size() == xu.size(),
        "ipopt::solve: size of xi, xl, and xu are not all equal."
    );

    CPPAD_ASSERT_KNOWN(
        gl.size() == gu.size(),
        "ipopt::solve: size of gl and gu are not equal."
    );

    size_t nx = xi.size();
    size_t ng = gl.size();

    // Create an IpoptApplication
    using Ipopt::IpoptApplication;
    Ipopt::SmartPtr<IpoptApplication> app = new IpoptApplication();

    // Set warm initialization
    app->Options()->SetStringValue("warm_start_init_point", "yes");
    app->Options()->SetNumericValue("warm_start_bound_frac", 1.0e-16);
    app->Options()->SetNumericValue("warm_start_mult_bound_push", 1.0e-16);
    app->Options()->SetNumericValue("warm_start_slack_bound_frac", 1.0e-16);
    app->Options()->SetNumericValue("warm_start_slack_bound_push", 1.0e-16);

    bool retape = false, sparse_forward = false, sparse_reverse = false;
    set_ipopt_app_options_from_string(app, options, retape, sparse_forward, sparse_reverse);

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    ok &= status == Ipopt::Solve_Succeeded;
    if (!ok)
    {
        solution.status = Optimization_result<Dvector>::status_type::unknown;
        return;
    }

    // Create an interface from Ipopt to this specific problem.
    // Note the assumption here that ADvector is same as cppd_ipopt::ADvector
    size_t nf = 1;
    Ipopt::SmartPtr<Ipopt::TNLP> cppad_nlp =
        new detail::Ipopt_NLP_CppAD<Dvector, ADvector, FG_eval>(
            nf,
            nx,
            ng,
            xi,
            xl,
            xu,
            gl,
            gu,
            lambda_i,
            zl_i,
            zu_i,
            fg_eval,
            retape,
            sparse_forward,
            sparse_reverse,
            solution
            );

    // Run the IpoptApplication
    app->OptimizeTNLP(cppad_nlp);

    solution.iter_count = app->Statistics()->IterationCount();

    return;
}


inline void set_ipopt_app_options_from_string(Ipopt::SmartPtr<Ipopt::IpoptApplication> &app,
                                              const std::string &options,
                                              bool &retape,
                                              bool &sparse_forward,
                                              bool &sparse_reverse)
{
    // process the options argument
    size_t begin_1, end_1, begin_2, end_2, begin_3, end_3;
    begin_1 = 0;
    retape = false;
    sparse_forward = false;
    sparse_reverse = false;

    while (begin_1 < options.size())
    {   
        // split this line into tokens
        while (options[begin_1] == ' ')
            begin_1++;
        end_1 = options.find_first_of(" \n", begin_1);
        begin_2 = end_1;
        while (options[begin_2] == ' ')
            begin_2++;
        end_2 = options.find_first_of(" \n", begin_2);
        begin_3 = end_2;
        while (options[begin_3] == ' ')
            begin_3++;
        end_3 = options.find_first_of(" \n", begin_3);

        // check for errors
        CPPAD_ASSERT_KNOWN(
            (end_1 != std::string::npos) &
            (end_2 != std::string::npos) &
            (end_3 != std::string::npos),
            "ipopt::solve: missing '\\n' at end of an option line"
        );
        CPPAD_ASSERT_KNOWN(
            (end_1 > begin_1) & (end_2 > begin_2),
            "ipopt::solve: an option line does not have two tokens"
        );

        // get first two tokens
        std::string tok_1 = options.substr(begin_1, end_1 - begin_1);
        std::string tok_2 = options.substr(begin_2, end_2 - begin_2);

        // get third token
        std::string tok_3;
        bool three_tok = false;
        three_tok |= tok_1 == "Sparse";
        three_tok |= tok_1 == "String";
        three_tok |= tok_1 == "Numeric";
        three_tok |= tok_1 == "Integer";
        if (three_tok)
        {
            CPPAD_ASSERT_KNOWN(
                (end_3 > begin_3),
                "ipopt::solve: a Sparse, String, Numeric, or Integer\n"
                "option line does not have three tokens."
            );
            tok_3 = options.substr(begin_3, end_3 - begin_3);
        }

        // switch on option type
        if (tok_1 == "Retape")
        {
            CPPAD_ASSERT_KNOWN(
                (tok_2 == "true") | (tok_2 == "false"),
                "ipopt::solve: Retape value is not true or false"
            );
            retape = (tok_2 == "true");
        }
        else if (tok_1 == "Sparse")
        {
            CPPAD_ASSERT_KNOWN(
                (tok_2 == "true") | (tok_2 == "false"),
                "ipopt::solve: Sparse value is not true or false"
            );
            CPPAD_ASSERT_KNOWN(
                (tok_3 == "forward") | (tok_3 == "reverse"),
                "ipopt::solve: Sparse direction is not forward or reverse"
            );
            if (tok_2 == "false")
            {
                sparse_forward = false;
                sparse_reverse = false;
            }
            else
            {
                sparse_forward = tok_3 == "forward";
                sparse_reverse = tok_3 == "reverse";
            }
        }
        else if (tok_1 == "String")
            app->Options()->SetStringValue(tok_2.c_str(), tok_3.c_str());
        else if (tok_1 == "Numeric")
        {
            Ipopt::Number value = std::atof(tok_3.c_str());
            app->Options()->SetNumericValue(tok_2.c_str(), value);
        }
        else if (tok_1 == "Integer")
        {
            Ipopt::Index value = std::atoi(tok_3.c_str());
            app->Options()->SetIntegerValue(tok_2.c_str(), value);
        }
        else
            CPPAD_ASSERT_KNOWN(
                false,
                "ipopt::solve: First token is not one of\n"
                "Retape, Sparse, String, Numeric, Integer"
            );

        begin_1 = end_3;
        while (options[begin_1] == ' ')
            begin_1++;
        if (options[begin_1] != '\n') CPPAD_ASSERT_KNOWN(
            false,
            "ipopt::solve: either more than three tokens "
            "or no '\\n' at end of a line"
        );
        begin_1++;
    }
}


template<typename Dvector,
         class FG_eval,
         typename ADvector = detail::ipopt_optimize_nlp_cppad::ADvector_type_t<FG_eval> >
NLP_complete_output ipopt_compute_complete_nlp_cppad(std::size_t nx,
                                                     std::size_t ng,
                                                     FG_eval &fg_eval,
                                                     const std::vector<scalar> &x,
                                                     const std::vector<scalar> &lambda)
{
    //
    // Completely characterizes the NLP represented by functor "fg_eval",
    // at point "x" with Lagrange multipliers "lambda"
    //

    bool ok = true;

    Optimization_result<Dvector> solution;

    bool retape = false, sparse_forward = true, sparse_reverse = false;

    // Create an interface from Ipopt to this specific problem.
    // Note the assumption here that ADvector is same as cppd_ipopt::ADvector
    size_t nf = 1;
    auto cppad_nlp = detail::Ipopt_NLP_CppAD<Dvector, ADvector, FG_eval>(
        nf,
        nx,
        ng,
        x,
        Dvector(nx, 0.0),
        Dvector(nx, 0.0),
        Dvector(ng, 0.0),
        Dvector(ng, 0.0),
        fg_eval,
        retape,
        sparse_forward,
        sparse_reverse,
        solution);

    NLP_complete_output output;

    Ipopt::TNLP::IndexStyleEnum index_style;

    // (1) Get NLP info (problem dimensions)
    cppad_nlp.get_nlp_info(output.number_of_variables, output.number_of_constraints, output.jacobian_nnz, output.hessian_nnz, index_style);

    // (2) Get Jacobian sparsity pattern
    output.col_jac.resize(output.jacobian_nnz);
    output.row_jac.resize(output.jacobian_nnz);

    cppad_nlp.eval_jac_g(output.number_of_variables, nullptr, false, output.number_of_constraints, output.jacobian_nnz, output.row_jac.data(), output.col_jac.data(), nullptr);

    // (3) Get Hessian sparsity pattern 
    output.col_hes.resize(output.hessian_nnz);
    output.row_hes.resize(output.hessian_nnz);

    cppad_nlp.eval_h(output.number_of_variables, nullptr, false, 1.0, output.number_of_constraints, nullptr, false, output.hessian_nnz, output.row_hes.data(), output.col_hes.data(), nullptr);

    // (4) Evaluate fitness function
    cppad_nlp.eval_f(output.number_of_variables, x.data(), true, output.fitness_function);

    // (5) Evaluate constraints
    output.constraints.resize(output.number_of_constraints);
    cppad_nlp.eval_g(output.number_of_variables, x.data(), false, output.number_of_constraints, output.constraints.data());

    // (6) Evaluate gradient of the fitness function
    output.grad_f.resize(output.number_of_variables);
    cppad_nlp.eval_grad_f(output.number_of_variables, x.data(), false, output.grad_f.data());

    // (7) Evaluate gradient of the constraints
    output.jacobian_constraints.resize(output.jacobian_nnz);
    cppad_nlp.eval_jac_g(output.number_of_variables, x.data(), false, output.number_of_constraints, output.jacobian_nnz, nullptr, nullptr, output.jacobian_constraints.data());

    // (8) Evaluate Hessian of the Lagrangian
    output.hessian_lagrangian.resize(output.hessian_nnz);
    cppad_nlp.eval_h(output.number_of_variables, x.data(), false, 1.0, output.number_of_constraints, lambda.data(), true, output.hessian_nnz, nullptr, nullptr, output.hessian_lagrangian.data());

    return output;
}

} // end namespace lioncpp

# endif
