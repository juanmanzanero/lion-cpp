#include "lion/math/sensitivity_analysis.h"
#include "lion/math/ipopt_cppad_handler.hpp"
#include "gtest/gtest.h"


using timeseries = CppAD::AD<double>;

class FG_eval {
public:
    using ADvector = CppAD::vector<timeseries>;

    void operator()(ADvector& fg, const ADvector& x, const ADvector& p)
    {
        throw std::runtime_error("class FG_eval has no parameters");
    }

    void operator()(ADvector& fg, const ADvector& x)
    {   assert( fg.size() == 3 );
        assert( x.size()  == 4 );

        // Fortran style indexing
        timeseries x1 = x[0];
        timeseries x2 = x[1];
        timeseries x3 = x[2];
        timeseries x4 = x[3];
        // f(x)
        fg[0] = x1 * x4 * (x1 + x2 + x3) + x3;
        // g_1 (x)
        fg[1] = x1 * x2 * x3 * x4;
        // g_2 (x)
        fg[2] = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4;
        //
        return;
    }
};


class FG_parameter 
{
 public:
    using ADvector = CppAD::vector<timeseries>;

    FG_parameter(const scalar& a_) : a(a_) {}
    
    void operator()(ADvector& fg, const ADvector& x, const ADvector& p)
    {
        assert(p.size() == 1);
        a = p[0];

        // Call the normal operator()
        (*this)(fg, x);
    }

    void operator()(ADvector& fg, const ADvector& x)
    {
        assert(fg.size() == 2);
        assert(x.size() == 2);

        fg[0] = (x[0] - a)*(x[0] - a) - x[1];
        fg[1] = x[1] - a;
    }

 private:
    CppAD::AD<scalar> a;
};


TEST(Sensitivity_analysis_test, check_consistency_test)
{
    // initial value of the independent variables
    const std::vector<double> x0 = {1.0, 5.0, 5.0, 1.0};

    // lower and upper limits for x
    const std::vector<double> x_lb = {1.0, 1.0, 1.0, 1.0}, x_ub = {5.0, 5.0, 5.0, 5.0};
    const std::vector<double> c_lb = {25.0, 40.0}, c_ub = {1.0e19, 40.0};

    // object that computes objective and constraints
    FG_eval fg_eval;

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "String  sb           yes\n";
    options += "Integer max_iter     10\n";
    options += "Numeric tol          1e-6\n";

    CppAD::ipopt_cppad_result<std::vector<double>> solution;

    // solve the problem
    CppAD::ipopt_cppad_solve<std::vector<double>, FG_eval>(
        options, x0, x_lb, x_ub, c_lb, c_ub, fg_eval, solution
    );

    auto opts = Sensitivity_analysis<FG_eval>::Options{};
    opts.max_error_dual_problem = 1.0e-6;
    EXPECT_TRUE(Sensitivity_analysis(FG_eval{}, solution.x, solution.s, solution.lambda, solution.zl, solution.zu, solution.vl, solution.vu, x_lb, x_ub, c_lb, c_ub, opts).optimality_check.success);
}


TEST(Sensitivity_analysis_test, compute_sensitivity_analysis)
{

    // initial value of the independent variables
    const std::vector<double> x0 = {0.0, 0.0};

    // lower and upper limits for x
    const std::vector<double> x_lb = {-5.0, -5.0};
    const std::vector<double> x_ub = { 5.0,  5.0};
    const std::vector<double> c_lb = {0.0};
    const std::vector<double> c_ub = {1.0};

    // object that computes objective and constraints
    FG_parameter fg(1.0);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "String  sb           yes\n";
    options += "Integer max_iter     10\n";
    options += "Numeric tol          1e-8\n";
    options += "String nlp_scaling_method none\n";
    //options += "Numeric bound_relax_factor 0.0\n";

    CppAD::ipopt_cppad_result<std::vector<double>> solution;

    // solve the problem
    CppAD::ipopt_cppad_solve<std::vector<double>, FG_parameter>(
        options, x0, x_lb, x_ub, c_lb, c_ub, fg, solution
    );

    EXPECT_NEAR(solution.x[0], 1.0, 1.0e-8);
    EXPECT_NEAR(solution.x[1], 2.0, 1.0e-8);
    
    auto opts = Sensitivity_analysis<FG_parameter>::Options{};
    opts.max_error_dual_problem = 1.0e-8;
    auto sensitivity_analysis = Sensitivity_analysis(fg, {1.0}, solution.x, solution.s, solution.lambda, solution.zl, solution.zu, solution.vl, solution.vu, x_lb, x_ub, c_lb, c_ub, opts);

    EXPECT_TRUE(sensitivity_analysis.optimality_check.success);

    EXPECT_NEAR(sensitivity_analysis.dxdparams.front()[0], 1.0, 1.0e-8);
    EXPECT_NEAR(sensitivity_analysis.dxdparams.front()[1], 1.0, 1.0e-8);
}
