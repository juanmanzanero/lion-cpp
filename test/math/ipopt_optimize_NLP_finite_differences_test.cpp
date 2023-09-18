#include "gtest/gtest.h"
#include "lion/math/ipopt_optimize_nlp_finite_differences.h"
#include "lion/foundation/constants.h"

TEST(Ipopt_optimize_NLP_finite_differences_test, sine_minimum_array)
{
    using namespace lioncpp;
    // This problem is simply to find the minimum of the sine function, but with one intermediate equality constraint

    // Fitness function: simply f(x,y) = y
    class F
    {
     public:
        using argument_type = std::array<double,2>;

        double operator()(const argument_type& x)
        {
            return x[1];
        }
    };
    
    // Equality constraint: y = sin(x)
    class C
    {
     public:
        using argument_type = std::array<double,2>;

        std::array<double,1> operator()(const argument_type& x)
        {
            return { sin(x[0]) - x[1] };
        }
    };

    F f;
    C c;

    // Perform optimisation
    Ipopt_optimize_NLP_finite_differences_options options;
    options.ipopt_findiff_perturbation = 1.0e-7;
    auto result = Ipopt_optimize_NLP_finite_differences<F,C>::optimize(2,1,{0.0,0.0},f,c,{-2.0,-2.0},{2.0,2.00},{0.0},{0.0},options);
    
    // Check results
    EXPECT_EQ(result.solved, true);
    EXPECT_NEAR(result.x[0], -pi*0.5,1.0e-7);
    EXPECT_NEAR(result.x[1], -1, 1.0e-7);
    EXPECT_NEAR(result.cons_err, 0.0, 1.0e-7);
    EXPECT_NEAR(result.f, -1, 1.0e-7);
}

TEST(Ipopt_optimize_NLP_finite_differences_test, sine_minimum_vector)
{
    using namespace lioncpp;
    // This problem is simply to find the minimum of the sine function, but with one intermediate equality constraint

    // Fitness function: simply f(x,y) = y
    class F
    {
     public:
        using argument_type = std::vector<double>;

        double operator()(const argument_type& x)
        {
            return x[1];
        }
    };
    
    // Equality constraint: y = sin(x)
    class C
    {
     public:
        using argument_type = std::vector<double>;

        std::array<double,1> operator()(const argument_type& x)
        {
            return { sin(x[0]) - x[1] };
        }
    };

    F f;
    C c;

    // Perform optimisation
    Ipopt_optimize_NLP_finite_differences_options options;
    options.ipopt_findiff_perturbation = 1.0e-7;
    auto result = Ipopt_optimize_NLP_finite_differences<F,C>::optimize(2,1,{0.0,0.0},f,c,{-2.0,-2.0},{2.0,2.00},{0.0},{0.0},options);
    
    // Check results
    EXPECT_EQ(result.solved, true);
    EXPECT_NEAR(result.x[0], -pi*0.5,1.0e-7);
    EXPECT_NEAR(result.x[1], -1, 1.0e-7);
    EXPECT_NEAR(result.cons_err, 0.0, 1.0e-7);
    EXPECT_NEAR(result.f, -1, 1.0e-7);
}
