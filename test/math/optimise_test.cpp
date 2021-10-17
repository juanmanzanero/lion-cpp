#include "gtest/gtest.h"
#include "lion/math/optimise.h"
#include "lion/foundation/constants.h"

TEST(Optimise_test, sine_minimum_array)
{
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
    Optimise_options options;
    options.findiff_perturbation = 1.0e-7;
    auto result = Optimise<F,C>::optimise(2,1,{0.0,0.0},f,c,{-2.0,-2.0},{2.0,2.00},{0.0},{0.0},options);
    
    // Check results
    EXPECT_EQ(result.solved, true);
    EXPECT_NEAR(result.x[0], -pi*0.5,1.0e-7);
    EXPECT_NEAR(result.x[1], -1, 1.0e-7);
    EXPECT_NEAR(result.cons_err, 0.0, 1.0e-7);
    EXPECT_NEAR(result.f, -1, 1.0e-7);
}

TEST(Optimise_test, sine_minimum_vector)
{
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
    Optimise_options options;
    options.findiff_perturbation = 1.0e-7;
    auto result = Optimise<F,C>::optimise(2,1,{0.0,0.0},f,c,{-2.0,-2.0},{2.0,2.00},{0.0},{0.0},options);
    
    // Check results
    EXPECT_EQ(result.solved, true);
    EXPECT_NEAR(result.x[0], -pi*0.5,1.0e-7);
    EXPECT_NEAR(result.x[1], -1, 1.0e-7);
    EXPECT_NEAR(result.cons_err, 0.0, 1.0e-7);
    EXPECT_NEAR(result.f, -1, 1.0e-7);
}
