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
    const auto max_iters = 100;
    const auto tolf = 1e-10;
    const auto tolx = 1e-10;

    const auto solution = CppAD::fsolve_trust_region_dodleg(Root2d<CppAD::AD<double> >{},
        false, std::vector<double>{ 0., 0. }, max_iters, tolf, tolx);

    EXPECT_EQ(solution.status, decltype(solution)::success);
    EXPECT_NEAR(solution.x[0], 0.35324661960, 1e-10);
    EXPECT_NEAR(solution.x[1], 0.60608173664, 1e-10);
    EXPECT_NEAR(solution.g[0], 0., 1e-10);
    EXPECT_NEAR(solution.g[1], 0., 1e-10);
    EXPECT_NEAR(solution.dg_dx_colmaj[0], -0.16699516067, 1e-10);
    EXPECT_NEAR(solution.dg_dx_colmaj[1], 1.3905452857, 1e-10);
    EXPECT_NEAR(solution.dg_dx_colmaj[2], -0.86358568558, 1e-10);
    EXPECT_NEAR(solution.dg_dx_colmaj[3], 0.14471832326, 1e-10);
    EXPECT_LE(solution.iter_count, 6);
}