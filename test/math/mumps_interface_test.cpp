#include <iostream>
#include "gtest/gtest.h"
#include "lion/math/mumps_interface.h"
#include "lion/math/matrix_extensions.h"


TEST(Mumps_interface_test, solve_3x3_colmaj)
{

    std::vector<double> A1{1,-4,-6,1,-4,-6,1,-4,-8};
    std::vector<double> A2{0,3,6,1,4,7,2,5,8};
    auto A = A1 + A2;
    std::vector<size_t> rows{ (1u - 1u), (2u - 1u), (3u - 1u), (1u - 1u), (2u - 1u), (3u - 1u), (1u - 1u), (2u - 1u), (3u - 1u) };
    std::vector<size_t> cols{ (1u - 1u), (1u - 1u), (1u - 1u), (2u - 1u), (2u - 1u), (2u - 1u), (3u - 1u), (3u - 1u), (3u - 1u) };

    std::vector<std::vector<double>> rhs = { { 1.0, 2.0, 3.0 }, { 1.0, -1.0, 2.0} };

    const auto x = mumps_solve_linear_system<size_t>(3u, 9u, rows, cols, A, rhs, false);

    EXPECT_NEAR(x[0][0], -2.75, 2.0e-15);
    EXPECT_NEAR(x[0][1], 3.0, 2.0e-16);
    EXPECT_NEAR(x[0][2], -0.75, 2.0e-15);

    EXPECT_NEAR(x[1][0], 0.0, 2.0e-15);
    EXPECT_NEAR(x[1][1], 2.0, 2.0e-15);
    EXPECT_NEAR(x[1][2], -1.0, 2.0e-15);

}


TEST(Mumps_interface_test, solve_3x3_repeated_indices)
{
    std::vector<double> A{1,-4,-6,1,-4,-6,1,-4,-8,0,1,2,3,4,5,6,7,8};
    std::vector<size_t> rows{0,1,2,0,1,2,0,1,2,0,0,0,1,1,1,2,2,2};
    std::vector<size_t> cols{0,0,0,1,1,1,2,2,2,0,1,2,0,1,2,0,1,2};

    std::vector<std::vector<double>> rhs = { { 1.0, 2.0, 3.0 }, { 1.0, -1.0, 2.0} };

    const auto x = mumps_solve_linear_system<size_t>(3u, 18u, rows, cols, A, rhs, false);

    EXPECT_NEAR(x[0][0], -2.75, 2.0e-15);
    EXPECT_NEAR(x[0][1], 3.0, 4.0e-15);
    EXPECT_NEAR(x[0][2], -0.75, 2.0e-15);

    EXPECT_NEAR(x[1][0], 0.0, 2.0e-15);
    EXPECT_NEAR(x[1][1], 2.0, 2.0e-15);
    EXPECT_NEAR(x[1][2], -1.0, 2.0e-15);
}
