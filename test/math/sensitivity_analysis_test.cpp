#include "lion/math/sensitivity_analysis.h"
#include "lion/math/ipopt_optimize_nlp_cppad.h"
#include "gtest/gtest.h"


using timeseries = CppAD::AD<double>;

class FG_eval {
public:
    using ADvector = CppAD::vector<timeseries>;

    void operator()(ADvector& fg, const ADvector& x, const ADvector& p)
    {
        if (p.size() > 0) throw lion_exception("class FG_eval has no parameters");

        (*this)(fg,x);
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

//! Auxiliary class that defines the following problem:
//!
//! Minimize (x-a)^2 - y + a + b
//!     s.t. 0 < b.y - a < 1   
//! whose solution is: (for b>0)
//! (x,y) = (a,(1+a)/b)
class FG_parameter 
{
 public:
    using ADvector = CppAD::vector<timeseries>;

    FG_parameter(const scalar& a_, const scalar& b_) : a(a_), b(b_) {}
    
    void operator()(ADvector& fg, const ADvector& x, const ADvector& p)
    {
        assert(p.size() <= 2);
        a = p[0];

        if (p.size() == 2) b = p[1];

        // Call the normal operator()
        (*this)(fg, x);
    }

    void operator()(ADvector& fg, const ADvector& x)
    {
        assert(fg.size() == 2);
        assert(x.size() == 2);

        fg[0] = (x[0] - a)*(x[0] - a) - x[1] + a + b;
        fg[1] = b*x[1] - a;
    }

 private:
    CppAD::AD<scalar> a;
    CppAD::AD<scalar> b;
};


TEST(Sensitivity_analysis_test, check_consistency_test)
{
    using namespace lioncpp;

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

    Optimization_result<std::vector<double>> solution;

    // solve the problem
    ipopt_optimize_nlp_cppad<std::vector<double>, FG_eval>(
        options, x0, x_lb, x_ub, c_lb, c_ub, fg_eval, solution
    );

    auto opts = Sensitivity_analysis<FG_eval>::Options{};
    opts.max_error_dual_problem = 1.0e-6;
    EXPECT_TRUE(Sensitivity_analysis(FG_eval{}, solution.x, solution.s, solution.lambda, solution.zl, solution.zu, solution.vl, solution.vu, x_lb, x_ub, c_lb, c_ub, opts).optimality_check.success);
}


TEST(Sensitivity_analysis_test, compute_sensitivity_analysis_1_parameter)
{
    using namespace lioncpp;

    // initial value of the independent variables
    const std::vector<double> x0 = {0.0, 0.0};

    // lower and upper limits for x
    const std::vector<double> x_lb = {-5.0, -5.0};
    const std::vector<double> x_ub = { 5.0,  5.0};
    const std::vector<double> c_lb = {0.0};
    const std::vector<double> c_ub = {1.0};

    // object that computes objective and constraints
    FG_parameter fg(1.0, 1.0);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "String  sb           yes\n";
    options += "Integer max_iter     10\n";
    options += "Numeric tol          1e-8\n";
    options += "String nlp_scaling_method none\n";
    //options += "Numeric bound_relax_factor 0.0\n";

    Optimization_result<std::vector<double>> solution;

    // solve the problem
    ipopt_optimize_nlp_cppad<std::vector<double>, FG_parameter>(
        options, x0, x_lb, x_ub, c_lb, c_ub, fg, solution
    );

    EXPECT_NEAR(solution.x[0], 1.0, 1.0e-8);
    EXPECT_NEAR(solution.x[1], 2.0, 1.0e-8);
    
    auto opts = Sensitivity_analysis<FG_parameter>::Options{};
    opts.max_error_dual_problem = 1.0e-8;
    auto sensitivity_analysis = Sensitivity_analysis(fg, {1.0}, solution.x, solution.s, solution.lambda, solution.zl, solution.zu, solution.vl, solution.vu, x_lb, x_ub, c_lb, c_ub, opts);

    EXPECT_TRUE(sensitivity_analysis.optimality_check.success);

    EXPECT_EQ(sensitivity_analysis.dxdp.size(), 1);
    EXPECT_NEAR(sensitivity_analysis.dxdp.front()[0], 1.0, 1.0e-8);
    EXPECT_NEAR(sensitivity_analysis.dxdp.front()[1], 1.0, 1.0e-8);
}

TEST(Sensitivity_analysis_test, compute_sensitivity_analysis_2_parameters)
{
    using namespace lioncpp;

    // initial value of the independent variables
    const std::vector<double> x0 = {0.0, 0.0};

    // lower and upper limits for x
    const std::vector<double> x_lb = {-5.0, -5.0};
    const std::vector<double> x_ub = { 5.0,  5.0};
    const std::vector<double> c_lb = {0.0};
    const std::vector<double> c_ub = {1.0};

    // object that computes objective and constraints
    FG_parameter fg(2.0, 2.0);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "String  sb           yes\n";
    options += "Integer max_iter     10\n";
    options += "Numeric tol          1e-8\n";
    options += "String nlp_scaling_method none\n";
    //options += "Numeric bound_relax_factor 0.0\n";

    Optimization_result<std::vector<double>> solution;

    // solve the problem
    ipopt_optimize_nlp_cppad<std::vector<double>, FG_parameter>(
        options, x0, x_lb, x_ub, c_lb, c_ub, fg, solution
    );

    EXPECT_NEAR(solution.x[0], 2.0, 1.0e-8);
    EXPECT_NEAR(solution.x[1], 1.5, 1.0e-8);
    
    auto opts = Sensitivity_analysis<FG_parameter>::Options{};
    opts.max_error_dual_problem = 1.0e-8;
    auto sensitivity_analysis = Sensitivity_analysis(fg, {2.0,2.0}, solution.x, solution.s, solution.lambda, solution.zl, solution.zu, solution.vl, solution.vu, x_lb, x_ub, c_lb, c_ub, opts);

    EXPECT_TRUE(sensitivity_analysis.optimality_check.success);

    EXPECT_EQ(sensitivity_analysis.dxdp.size(), 2);
    EXPECT_NEAR(sensitivity_analysis.dxdp.front()[0], 1.0, 1.0e-8);
    EXPECT_NEAR(sensitivity_analysis.dxdp.front()[1], 0.5, 1.0e-8);
       
    EXPECT_NEAR(sensitivity_analysis.dxdp.back()[0], 0.0, 1.0e-8);
    EXPECT_NEAR(sensitivity_analysis.dxdp.back()[1], -0.75, 1.0e-8);
}
