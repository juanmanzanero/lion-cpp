#include <vector>

#include "gtest/gtest.h"

#include "lion/math/fsolve_trust_region_dodleg.h"


//
// Defines tests for the "fsolve_trust_region_dodleg" function.
//


template<typename T>
struct Root2d
{
    //
    // The first example from https://www.mathworks.com/help/optim/ug/fsolve.html.
    //

    using ADvector = std::vector<T>;

    void operator()(ADvector &g, const ADvector &x) const
    {
        using std::exp;
        using std::sin;
        using std::cos;

        g[0] = exp(-exp(-(x[0] + x[1]))) - x[1] * (1. + x[0] * x[0]);
        g[1] = x[0] * cos(x[1]) + x[1] * sin(x[0]) - 0.5;
    }
};


TEST(fsolve_trust_region_dodleg_test, solve_Root2d)
{
    const auto solution = CppAD::fsolve_trust_region_dodleg(Root2d<CppAD::AD<double> >{},
        std::vector<double>{ 0., 0. });

    EXPECT_EQ(solution.status, decltype(solution)::success);
    EXPECT_NEAR(solution.g[0], 0., 1e-10);
    EXPECT_NEAR(solution.g[1], 0., 1e-10);
    EXPECT_NEAR(solution.x[0], 0.3532466195, 1e-10);
    EXPECT_NEAR(solution.x[1], 0.6060817366, 1e-10);
    EXPECT_LE(solution.iter_count, 6);
}