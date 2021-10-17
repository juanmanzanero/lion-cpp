#include "gtest/gtest.h"
#include "lion/math/polynomial.h"
#include "lion/math/matrix3x3.h"
#include <fstream>

template class Polynomial<double>;
template struct Vector3d<double>;
template struct Matrix3x3<double>;

constexpr scalar f(scalar x) { return 1.0 + x + 0.5*x*x + x*x*x/3.0 + x*x*x*x/4.0; }
constexpr scalar df(scalar x) { return 1.0 + x + x*x + x*x*x; }

const inline std::vector<scalar> engine_data = { 1.0000e+4, 1.4251e+1,
                                                 1.0237e+4, 1.7729e+1,
                                                 1.0504e+4, 2.2850e+1,
                                                 1.0760e+4, 2.7198e+1,
                                                 1.1057e+4, 3.1449e+1,
                                                 1.1333e+4, 3.4444e+1,
                                                 1.1659e+4, 3.6860e+1,
                                                 1.1985e+4, 3.9179e+1,
                                                 1.2281e+4, 4.1111e+1,
                                                 1.2598e+4, 4.3816e+1,
                                                 1.2914e+4, 4.6715e+1,
                                                 1.3200e+4, 4.4976e+1,
                                                 1.3407e+4, 4.1884e+1,
                                                 1.3615e+4, 3.7923e+1,
                                                 1.3733e+4, 3.3575e+1,
                                                 1.3812e+4, 3.0290e+1,
                                                 1.3862e+4, 2.6715e+1,
                                                 1.3921e+4, 2.3333e+1,
                                                 1.4000e+4, 1.9952e+1  };

TEST(Polynomial_test, empty_polynomial)
{
    std::vector<timeseries> L;
    std::vector<std::vector<timeseries>> y0;
    sPolynomial p(0.0,L,y0);

    EXPECT_DOUBLE_EQ(p(3.14), 0.0);
}

TEST(Polynomial_test, second_order_polynomial)
{
    std::vector<scalar> x0 = {-1.0,0.0,1.0};
    std::vector<scalar> y0 = {0.0,1.0,2.0};
    Polynomial p(x0,y0,x0.size()-1);

    for (size_t i = 0; i < 3; ++i)
        EXPECT_NEAR(p(x0.at(i)),y0.at(i),2.0e-15);

    EXPECT_NEAR(p.get_coeffs().at(0),1.0,2.0e-15);
    EXPECT_NEAR(p.get_coeffs().at(1),1.0,2.0e-15);
    EXPECT_NEAR(p.get_coeffs().at(2),0.0,2.0e-15);
}


TEST(Polynomial_test, fifth_order_polynomial)
{
    std::vector<scalar> x0 = {-2.0,-1.5,0.0,0.25,0.5,2.0};
    std::vector<scalar> y0 = {0.5,1.0,2.0,0.25,-1.0,0.3};
    Polynomial p(x0,y0,x0.size()-1);

    for (size_t i = 0; i < 6; ++i)
        EXPECT_NEAR(p(x0.at(i)),y0.at(i),5.0e-15);

}


TEST(Polynomial_test, tenth_order_polynomial)
{
    std::vector<scalar> x0 = {-2.0,-1.5,0.0,0.25,0.5,2.0,4.0,4.5,6.0,7.6,8.0};
    std::vector<scalar> y0 = {0.5,1.0,2.0,0.25,-1.0,0.3,-2.0,5.0,9.0,-10.0,-12.0};
    Polynomial p(x0,y0,x0.size()-1);

    for (size_t i = 0; i < 6; ++i)
        EXPECT_NEAR(p(x0.at(i)),y0.at(i),6.0e-15);

}


TEST(Polynomial_test, derivative)
{
    std::vector<scalar> x0 = {-0.5,-0.25,0.0,0.25,0.5};
    std::vector<scalar> y0 = {f(x0[0]), f(x0[1]), f(x0[2]), f(x0[3]), f(x0[4])};

    Polynomial p(x0,y0,x0.size()-1);

    Polynomial dp = p.derivative();

    for (size_t i = 0; i < 5; ++i)
        EXPECT_NEAR(dp(x0.at(i)), df(x0.at(i)),2.0e-14);

}

TEST(Polynomial_test, integral)
{
    std::vector<scalar> x0 = {-0.5,-0.25,0.0,0.25,0.5};
    std::vector<scalar> y0 = {f(x0[0]), f(x0[1]), f(x0[2]), f(x0[3]), f(x0[4])};

    Polynomial p(x0,y0,x0.size()-1);

    Polynomial int_p = p.integral();

    Polynomial d_int_p = int_p.derivative();

    for (size_t i = 0; i < 5; ++i)
        EXPECT_DOUBLE_EQ(d_int_p.get_coeffs().at(i), p.get_coeffs().at(i));
}




TEST(Polynomial_test, smooth_polynomial)
{
    std::vector<scalar> saved_error = { 1.43885e-13, 0.000736931, 0.0159682, -0.0317231, 0.026048, 0.0141357, -0.0149506, 
                                        0.0137441, -0.214174, 0.0111608, 0.0247169, -0.0426472, 0.017571, 0.425088, 
                                        0.412429, -0.704514, -0.14905, 0.0607866, 2.20268e-13};

    std::vector<scalar> x0,y0;

    for (size_t i = 0; i < engine_data.size()/2; ++i)
    {   
        x0.push_back(engine_data[2*i]);
        y0.push_back(engine_data[2*i+1]);
    }

    Polynomial p(x0,y0,20,true);

    for (size_t i = 0; i < x0.size(); ++i)
        EXPECT_NEAR(p(x0[i])-y0[i], saved_error[i],1.0e-6);
}


TEST(Polynomial_test, piecewise_polynomial)
{
    std::vector<scalar> x0,y0;

    for (size_t i = 0; i < engine_data.size()/2; ++i)
    {   
        x0.push_back(engine_data[2*i]);
        y0.push_back(engine_data[2*i+1]);
    }

    Polynomial p(x0,y0,4);

    for (size_t i = 0; i < x0.size(); ++i)
        EXPECT_NEAR(p(x0[i]),y0[i],3.0e-14);
}


TEST(Polynomial_test, piecewise_integral_derivative)
{
    std::vector<scalar> x0,y0;

    for (size_t i = 0; i < engine_data.size()/2; ++i)
    {   
        x0.push_back(engine_data[2*i]);
        y0.push_back(engine_data[2*i+1]);
    }

    Polynomial p(x0,y0,4);

    Polynomial int_p = p.integral();

    Polynomial d_int_p = int_p.derivative();

    for (size_t n = 0; n < p.get_n_blocks(); ++n)
        for (size_t i = 0; i < p.get_coeffs(n).size(); ++i)
            EXPECT_DOUBLE_EQ(d_int_p.get_coeffs(n).at(i), p.get_coeffs(n).at(i));

    // For the integral of the derivative we only check coeffs 1->N (L0 info is lost at derivation)
    Polynomial dp = p.derivative();
    Polynomial int_dp = dp.integral();

    for (size_t n = 0; n < p.get_n_blocks(); ++n)
        for (size_t i = 1; i < p.get_coeffs(n).size(); ++i)
            EXPECT_DOUBLE_EQ(int_dp.get_coeffs(n).at(i), p.get_coeffs(n).at(i));
}


TEST(Polynomial_test, ninety_degrees_corner_test)
{
//    Polynomial(const std::vector<scalar>& L, const std::vector<std::vector<T>>& y0);

    // Part 1: 1m straight
    scalar L1 = 1.0;
    std::vector<sVector3d> y1 = { {0.0,0.0,0.0}, {0.0, 0.5, 0.0}, {0.0,1.0,0.0} };

    // Part 2: 90 degrees bend with 1m radius
    scalar L2 = pi/2.0;
    std::vector<sVector3d> y2;
    const size_t N2 = 10;

    const auto xi2 = std::get<0>(gauss_legendre_lobatto_nodes_and_weights(N2));

    for (size_t i = 0; i <= N2; ++i)
        y2.push_back(sVector3d(1.0-cos(pi/4.0*(xi2[i]+1.0)),1.0+sin(pi/4.0*(xi2[i]+1.0)),0.0));

    // Part 3: 1m horizontal straight
    scalar L3 = 1.0;
    std::vector<sVector3d> y3 = { y2.back(), {1.5, 2.0, 0.0} , {2.0,2.0,0.0} };

    vPolynomial p(0.0, {L1,L2,L3}, {y1,y2,y3});

    Polynomial dp = p.derivative();

    Polynomial d2p = dp.derivative();

    EXPECT_DOUBLE_EQ(dp(0.0).at(X), 0.0);
    EXPECT_DOUBLE_EQ(dp(0.0).at(Y), 1.0);

    EXPECT_DOUBLE_EQ(dp(1.0).at(X), 0.0);
    EXPECT_DOUBLE_EQ(dp(1.0).at(Y), 1.0);

    EXPECT_NEAR(dp(L1+0.5*L2).at(X), cos(pi/4.0), 1.0e-10);
    EXPECT_NEAR(dp(L1+0.5*L2).at(Y), sin(pi/4.0), 1.0e-10);

    EXPECT_NEAR(dp(L1+L2).at(X), 1.0,2.0e-10);
    EXPECT_NEAR(dp(L1+L2).at(Y), 0.0,2.0e-10);

    EXPECT_NEAR(dp(L1+L2+L3).at(X), 1.0,1.0e-15);
    EXPECT_NEAR(dp(L1+L2+L3).at(Y), 0.0,1.0e-15);

    EXPECT_NEAR(d2p(0.0).norm(), 0.0, 1.0e-15);
    EXPECT_NEAR(d2p(L1).norm(), 0.0, 1.0e-15);

    // Curvature of the circular part
    EXPECT_NEAR(cross(dp(L1+0.1*L2),d2p(L1+0.1*L2)).norm()/pow(dp(L1+0.1*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.2*L2),d2p(L1+0.2*L2)).norm()/pow(dp(L1+0.2*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.3*L2),d2p(L1+0.3*L2)).norm()/pow(dp(L1+0.3*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.4*L2),d2p(L1+0.4*L2)).norm()/pow(dp(L1+0.4*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.5*L2),d2p(L1+0.5*L2)).norm()/pow(dp(L1+0.5*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.6*L2),d2p(L1+0.6*L2)).norm()/pow(dp(L1+0.6*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.7*L2),d2p(L1+0.7*L2)).norm()/pow(dp(L1+0.7*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.8*L2),d2p(L1+0.8*L2)).norm()/pow(dp(L1+0.8*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+0.9*L2),d2p(L1+0.9*L2)).norm()/pow(dp(L1+0.9*L2).norm(),3), 1.0, 2.0e-10);
    EXPECT_NEAR(cross(dp(L1+1.0*L2),d2p(L1+1.0*L2)).norm()/pow(dp(L1+1.0*L2).norm(),3), 1.0, 7.0e-09);
}
