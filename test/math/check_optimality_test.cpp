#include "lion/math/check_optimality.h"
#include "lion/math/ipopt_cppad_handler.hpp"
#include "gtest/gtest.h"


using timeseries = CppAD::AD<double>;

class FG_eval {
public:
    using ADvector = CppAD::vector<timeseries>;

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


TEST(Check_optimality_test, aa)
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
    
    auto opts = Check_optimality<FG_eval>::Options{};
    opts.constraint_viol_tolerance = 1.0e-6;
    EXPECT_TRUE(Check_optimality(FG_eval{}, solution.x, solution.s, solution.lambda, solution.zl, solution.zu, solution.vl, solution.vu, x_lb, x_ub, c_lb, c_ub, opts).success);
}
