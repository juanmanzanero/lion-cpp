#include "gtest/gtest.h"
#include "lion/math/legendre_algorithms.h"


TEST(Legendre_algorithms_test, Nodes_and_weights_P4)
{
    const auto xj = gauss_legendre_nodes_and_weights(4);

    std::array<double,5> x = {   -0.906179845938664,
                                 -0.538469310105683,
                                                  0,
                                  0.538469310105683,
                                  0.906179845938664 };

    std::array<double,5> w = {    0.236926885056189,
                                  0.478628670499366,
                                  0.568888888888889,
                                  0.478628670499366,
                                  0.236926885056189 };

    for (size_t j = 0; j < 5; ++j)
    {
        EXPECT_NEAR(xj.first.at(j), x.at(j),7.0e-16);
        EXPECT_NEAR(xj.second.at(j), w.at(j),8.0e-16);
    }
}

TEST(Legendre_algorithms_test, Nodes_and_weights_P7)
{
    const auto xj = gauss_legendre_nodes_and_weights(7);

    std::array<double,8> x = {    -0.960289856497536,
                                  -0.796666477413627,
                                  -0.525532409916329,
                                  -0.183434642495650,
                                   0.183434642495650,
                                   0.525532409916329,
                                   0.796666477413627,
                                   0.960289856497536 };

    std::array<double,8> w = {     0.101228536290377,
                                   0.222381034453374,
                                   0.313706645877887,
                                   0.362683783378362,
                                   0.362683783378362,
                                   0.313706645877887,
                                   0.222381034453374,
                                   0.101228536290377 };

    for (size_t j = 0; j < 8; ++j)
    {
        EXPECT_NEAR(xj.first.at(j), x.at(j),7.0e-16);
        EXPECT_NEAR(xj.second.at(j), w.at(j),7.0e-16);
    }
}


TEST(Legendre_algorithms_test, LGL_Nodes_and_weights_P4)
{
    const size_t N = 4;
    const auto [xj, wj] = gauss_legendre_lobatto_nodes_and_weights(N);

    std::vector<scalar> x =
                  {
                    -1.000000000000000,
                    -0.654653670707977,
                                     0,
                     0.654653670707977,
                     1.000000000000000
                  };

    std::vector<scalar> w =
                  {
                     0.100000000000000,
                     0.544444444444444,
                     0.711111111111111,
                     0.544444444444444,
                     0.100000000000000
                  };

    for (size_t j = 0; j <= N; ++j)
    {
        EXPECT_NEAR(xj.at(j), x.at(j),4.0e-16);
        EXPECT_NEAR(wj.at(j), w.at(j),4.0e-16);
    }
}


TEST(Legendre_algorithms_test, LGL_Nodes_and_weights_P7)
{
    const size_t N = 7;
    const auto [xj, wj] = gauss_legendre_lobatto_nodes_and_weights(N);

    std::vector<scalar> x =
                  {
                    -1.000000000000000,
                    -0.871740148509607,
                    -0.591700181433142,
                    -0.209299217902479,
                     0.209299217902479,
                     0.591700181433142,
                     0.871740148509607,
                     1.000000000000000
                  };


    std::vector<scalar> w =
                  {
                     0.035714285714286,
                     0.210704227143506,
                     0.341122692483504,
                     0.412458794658704,
                     0.412458794658704,
                     0.341122692483504,
                     0.210704227143506,
                     0.035714285714286
                  };

    for (size_t j = 0; j <= N; ++j)
    {
        EXPECT_NEAR(xj.at(j), x.at(j),5.0e-16);
        EXPECT_NEAR(wj.at(j), w.at(j),5.0e-16);
    }
}

