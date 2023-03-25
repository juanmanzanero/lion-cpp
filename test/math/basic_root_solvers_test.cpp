#include "gtest/gtest.h"

#include "lion/math/basic_root_solvers.h"


//
// Defines tests for the functions defined
// in the "lion/math/basic_root_solvers.h" header.
//


constexpr auto tolnear = 1e-12;


TEST(basic_root_solvers_test, fzero_seed_test)
{
    const auto x0 = 3.;
    const auto [x, fval, exitflag] = fzero([](auto x){ return std::sin(x); },
        x0, 1000u, 1e-12, 1e-12);
    EXPECT_NEAR(fval, 0., tolnear);
    EXPECT_NEAR(x, pi, tolnear);
    EXPECT_EQ(exitflag, basic_root_solvers_exitflag::success);
}


TEST(basic_root_solvers_test, fzero_interval_test)
{
    const auto xlo = 1.;
    const auto xhi = 2.;
    const auto [x, fval, exitflag] = fzero([](auto x){ return std::cos(x); },
        xlo, xhi, 1000u, 1e-12, 1e-12);
    EXPECT_NEAR(fval, 0., tolnear);
    EXPECT_NEAR(x, 0.5 * pi, tolnear);
    EXPECT_EQ(exitflag, basic_root_solvers_exitflag::success);
}


TEST(basic_root_solvers_test, Newton_method1d_test)
{
    const auto x0 = 1.6;
    const auto[x, fval, exitflag] = Newton_method1d([](auto x) { return std::sin(std::cosh(x)); },
        [](auto x) { return std::cos(std::cosh(x)) * std::sinh(x); },
        x0, 0.9, 1000u, 1e-12, 1e-12);
    EXPECT_NEAR(fval, 0., tolnear);
    EXPECT_NEAR(x, 1.811526272460853, tolnear);
    EXPECT_EQ(exitflag, basic_root_solvers_exitflag::success);
}


TEST(basic_root_solvers_test, fixed_point_iteration_test)
{
    const auto x0 = 0.00079;
    const auto logistic_map = [](auto x) { return 1.5 * x * (1. - x); };
    const auto [x, exitflag] = fixed_point_iteration(logistic_map,
        x0, 1000u, 1e-14);
    EXPECT_NEAR(std::abs(logistic_map(x) - x), 0., tolnear);
    EXPECT_NEAR(x, 0.333333333333333, tolnear);
    EXPECT_EQ(exitflag, basic_root_solvers_exitflag::success);
}


TEST(basic_root_solvers_test, numjac_Newton_method2d_test)
{
    const auto  x0 = std::array<double, 2>{ 0., 1. };
    const auto [x, fval, exitflag] = numjac_Newton_method2d([](const auto &x)
        {
            return std::array<double, 2>{ 2. * x[0] + x[1] - std::exp(-x[0]),
                -x[0] + 2. * x[1] - std::exp(-x[1])};
        },
        x0,
        std::array<double, 2>{ 1e-9, 1e-9 },
        std::array<double, 2>{ 0.9, 0.9 },
        1000, 1e-12, 1e-12);
    EXPECT_NEAR(fval[0], 0., tolnear);
    EXPECT_NEAR(fval[1], 0., tolnear);
    EXPECT_NEAR(x[0], 0.197594329485429, tolnear);
    EXPECT_NEAR(x[1], 0.425514061540109, tolnear);
    EXPECT_EQ(exitflag, basic_root_solvers_exitflag::success);
}


TEST(basic_root_solvers_test, real_root_of_cubic)
{
    const auto x = real_root_of_cubic(1., 0., -2., -5.);
    EXPECT_NEAR(x, 2.094551481542328, tolnear);
}


TEST(basic_root_solvers_test, cubic_roots)
{
    const auto [xs, num_real_roots] = cubic_roots(1., 0., -2., -5.);
    EXPECT_EQ(num_real_roots, 1);
    EXPECT_NEAR(xs[0].real(), 2.094551481542328, tolnear);
    EXPECT_NEAR(xs[0].imag(), 0., tolnear);
    EXPECT_NEAR(xs[1].real(), -1.047275740771163, tolnear);
    EXPECT_NEAR(xs[1].imag(), 1.135939889088928, tolnear);
    EXPECT_NEAR(xs[2].real(), xs[1].real(), tolnear);
    EXPECT_NEAR(xs[2].imag(), -xs[1].imag(), tolnear);
}


TEST(basic_root_solvers_test, quartic_roots)
{
    const auto [xs, num_real_roots] = quartic_roots(1., 2., 1., 4., 2.);
    EXPECT_EQ(num_real_roots, 2);
    EXPECT_NEAR(xs[0].real(), -0.515594615780425, tolnear);
    EXPECT_NEAR(xs[0].imag(), 0., tolnear);
    EXPECT_NEAR(xs[1].real(), -2.187661190106682, tolnear);
    EXPECT_NEAR(xs[1].imag(), 0., tolnear);
    EXPECT_NEAR(xs[2].real(), 0.351627902943554, tolnear);
    EXPECT_NEAR(xs[2].imag(), 1.284325436713318, tolnear);
    EXPECT_NEAR(xs[3].real(), xs[2].real(), tolnear);
    EXPECT_NEAR(xs[3].imag(), -xs[2].imag(), tolnear);
}