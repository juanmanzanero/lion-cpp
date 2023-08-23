#include "gtest/gtest.h"
#include "lion/math/ipopt_optimize_nlp_cppad.h"
#include "lion/math/matrix_extensions.h"

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


TEST(ipopt_optimize_nlp_cppad_test,optimize)
{
    using namespace lioncpp;

    size_t i;

    // number of independent variables (domain dimension for f and g)
    size_t nx = 4;
    // number of constraints (range dimension for g)
    size_t ng = 2;
    // initial value of the independent variables
    std::vector<double> xi(nx);
    xi[0] = 1.0;
    xi[1] = 5.0;
    xi[2] = 5.0;
    xi[3] = 1.0;
    // lower and upper limits for x
    std::vector<double> xl(nx), xu(nx);
    for(i = 0; i < nx; i++)
    {   xl[i] = 1.0;
        xu[i] = 5.0;
    }
    // lower and upper limits for g
    std::vector<double> gl(ng), gu(ng);
    gl[0] = 25.0;     gu[0] = 1.0e19;
    gl[1] = 40.0;     gu[1] = 40.0;

    // object that computes objective and constraints
    FG_eval fg_eval;

    // options
    std::string options;
    // turn off any printing
    options += "Integer print_level  0\n";
    options += "String  sb           yes\n";
    // maximum number of iterations
    options += "Integer max_iter     10\n";
    // approximate accuracy in first order necessary conditions;
    // see Mathematical Programming, Volume 106, Number 1,
    // Pages 25-57, Equation (6)
    options += "Numeric tol          1e-6\n";
    // derivative testing
    options += "String  derivative_test            second-order\n";
    // maximum amount of random pertubation; e.g.,
    // when evaluation finite diff
    options += "Numeric point_perturbation_radius  0.\n";

    // place to return solution
    Optimization_result<std::vector<double>> solution;

    // solve the problem
    ipopt_optimize_nlp_cppad<std::vector<double>, FG_eval>(
        options, xi, xl, xu, gl, gu, fg_eval, solution
    );

    //
    // Check some of the solution values
    //
    EXPECT_TRUE(solution.status == Optimization_result<std::vector<double>>::success);

    //
    double check_x[]  = { 1.000000, 4.743000, 3.82115, 1.379408 };
    double check_zl[] = { 1.087871, 0.,       0.,      0.       };
    double check_zu[] = { 0.,       0.,       0.,      0.       };
    for(i = 0; i < nx; i++)
    {   
        EXPECT_NEAR(check_x[i], solution.x[i], 1.0e-6);
        EXPECT_NEAR(check_zl[i], solution.zl[i], 1.0e-6);
        EXPECT_NEAR(check_zu[i], solution.zu[i], 1.0e-6);
    }
}

TEST(ipopt_optimize_nlp_cppad_test,optimize_restart)
{
    using namespace lioncpp;

    size_t i;

    // number of independent variables (domain dimension for f and g)
    size_t nx = 4;
    // number of constraints (range dimension for g)
    size_t ng = 2;
    // lower and upper limits for x
    std::vector<double> xl(nx), xu(nx);
    for(i = 0; i < nx; i++)
    {   xl[i] = 1.0;
        xu[i] = 5.0;
    }
    // lower and upper limits for g
    std::vector<double> gl(ng), gu(ng);
    gl[0] = 25.0;     gu[0] = 1.0e19;
    gl[1] = 40.0;     gu[1] = 40.0;

    std::vector<double> xi = {1.000000073584267, 4.742999579337057, 3.821150081818284, 1.379408168555};
    std::vector<double> lambda_init = {-0.5522936369791295, 0.1614685262869733};
    std::vector<double> zl_init = {1.087871777267788, 2.42876061301246e-08, 3.222458735381987e-08, 2.395768373390596e-07};
    std::vector<double> zu_init = {2.272708130342488e-08, 3.53765983054055e-07, 7.711415533987246e-08, 2.510926239765509e-08};

    // object that computes objective and constraints
    FG_eval fg_eval;

    // options
    std::string options;
    // turn off any printing
    options += "Integer print_level  0\n";
    options += "String  sb           yes\n";
    // maximum number of iterations
    options += "Integer max_iter     10\n";
    // approximate accuracy in first order necessary conditions;
    // see Mathematical Programming, Volume 106, Number 1,
    // Pages 25-57, Equation (6)
    options += "Numeric tol          1e-6\n";
    // derivative testing
    options += "String  derivative_test            second-order\n";
    // maximum amount of random pertubation; e.g.,
    // when evaluation finite diff
    options += "Numeric point_perturbation_radius  0.\n";

    // place to return solution
    Optimization_result<std::vector<double>> solution;

    // solve the problem
    ipopt_optimize_nlp_cppad<std::vector<double>, FG_eval>(
        options, xi, xl, xu, gl, gu, lambda_init, zl_init, zu_init, fg_eval, solution
    );
    
    // The most important check: check that no iterations were performed!
    EXPECT_EQ(solution.iter_count, 0);

    //
    // Check some of the solution values
    //
    EXPECT_TRUE(solution.status == Optimization_result<std::vector<double>>::success);

    //
    double check_x[]  = { 1.000000, 4.743000, 3.82115, 1.379408 };
    double check_zl[] = { 1.087871, 0.,       0.,      0.       };
    double check_zu[] = { 0.,       0.,       0.,      0.       };
    for(i = 0; i < nx; i++)
    {   
        EXPECT_NEAR(check_x[i], solution.x[i], 1.0e-6);
        EXPECT_NEAR(check_zl[i], solution.zl[i], 1.0e-6);
        EXPECT_NEAR(check_zu[i], solution.zu[i], 1.0e-6);
    }
}
