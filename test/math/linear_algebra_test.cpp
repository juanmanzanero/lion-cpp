#include <vector>
#include <random>

#include "gtest/gtest.h"

#include "lion/io/Xml_document.h"
#include "lion/math/linear_algebra.h"
#include "lion/math/small_dense_inverse_matrices.h"
#ifdef WITH_EIGEN
#include "cppad/example/atomic_two/eigen_mat_inv.hpp"
#include "cppad/example/atomic_two/eigen_mat_mul.hpp"
#endif


//
// Defines tests for the functions in the
// "lion/math/linear_algebra.h" header.
//


static Xml_document linear_algebra_test_reference_file("data/linear_algebra_test.xml", true);


constexpr auto tolnear_lusolve = 1e3 * std::numeric_limits<double>::epsilon();
constexpr auto tolnear_qrsolve = 1e3 * std::numeric_limits<double>::epsilon();
constexpr auto tolnear_rcond = 1e3 * std::numeric_limits<double>::epsilon();
constexpr auto tolnear_tridiagonalsolve = 1.5e-11; // a bit bigger since we'll compare against reference
                                                   // solutions obtained with Matlab's "linsolve" (i.e., with full matrices)


inline std::vector<double> comma_separated_string2vector_of_doubles(const std::string &str)
{
    //
    // Comma separated string -> vector of doubles.
    // Faster than "Xml_element::get_value(std::vector<double>{})",
    // though less safe.
    //

    std::istringstream ss(str);
    std::string num;
    std::vector<double> vec;
    while (std::getline(ss, num, ',')) {
        vec.emplace_back(std::stod(num));
    }

    return vec;
}


TEST(linear_algebra_test, lusolve_test0)
{
    // 3 x 3 singular problem with 1 column in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    const auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size());
    EXPECT_EQ(n, 3);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    auto A = A_colmaj;
    auto b = b_colmaj;
    const auto ret = lusolve(b.data(), A.data(), n, 1);
    EXPECT_NE(ret, 0);

    auto A_ = A_colmaj;
    std::vector<int> ipiv_A(n);
    const auto ret_ = plufactorize(A_.data(), ipiv_A.data(), n);
    EXPECT_NE(ret_, 0);
    EXPECT_EQ(ret_, ret);
}


TEST(linear_algebra_test, lusolve_test1)
{
    // 4 x 4 problem with 1 column in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    const auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size());
    EXPECT_EQ(n, 4);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    const auto reference_x_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/x")->
        get_value(std::vector<double>{});

    {
        auto A = A_colmaj;
        auto b = b_colmaj;
        EXPECT_EQ(lusolve(b.data(), A.data(), n, 1), 0);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear_lusolve);
        }
    }

    {
        auto A = A_colmaj;
        std::vector<int> ipiv_A(n);
        EXPECT_EQ(plufactorize(A.data(), ipiv_A.data(), n), 0);

        auto b = b_colmaj;
        plusolve(b.data(), A.data(), ipiv_A.data(), n, 1);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear_lusolve);
        }
    }
}


TEST(linear_algebra_test, lusolve_test2)
{
    // 10 x 10 problem with 3 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size()) / 3;
    EXPECT_EQ(n, 10);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    const auto reference_x_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/x")->
        get_value(std::vector<double>{});

    {
        auto A = A_colmaj;
        auto b = b_colmaj;
        EXPECT_EQ(lusolve(b.data(), A.data(), n, 3), 0);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear_lusolve);
        }
    }

    {
        auto A = A_colmaj;
        std::vector<int> ipiv_A(n);
        EXPECT_EQ(plufactorize(A.data(), ipiv_A.data(), n), 0);

        auto b = b_colmaj;
        plusolve(b.data(), A.data(), ipiv_A.data(), n, 3);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear_lusolve);
        }
    }
}


TEST(linear_algebra_test, lusolve_test3)
{
    // 20 x 20 problem with 7 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size()) / 7;
    EXPECT_EQ(n, 20);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    const auto reference_x_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/x")->
        get_value(std::vector<double>{});

    {
        auto A = A_colmaj;
        auto b = b_colmaj;
        EXPECT_EQ(lusolve(b.data(), A.data(), n, 7), 0);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear_lusolve);
        }
    }

    {
        auto A = A_colmaj;
        std::vector<int> ipiv_A(n);
        EXPECT_EQ(plufactorize(A.data(), ipiv_A.data(), n), 0);

        auto b = b_colmaj;
        plusolve(b.data(), A.data(), ipiv_A.data(), n, 7);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear_lusolve);
        }
    }
}


TEST(linear_algebra_test, qrsolve_test0)
{
    // 10 x 4 matrix with 1 column in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto reference_x = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/x")->get_value(std::vector<double>{});
    std::vector<double> x_colmaj(4);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 10, 4, 1);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear_qrsolve);
    }
}


TEST(linear_algebra_test, qrsolve_test1)
{
    // 4 x 10 matrix with 1 column in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto reference_x = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/x")->get_value(std::vector<double>{});
    std::vector<double> x_colmaj(10);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 4, 10, 1);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear_qrsolve);
    }
}


TEST(linear_algebra_test, qrsolve_test2)
{
    // 11 x 6 matrix with 4 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto reference_x = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/x")->get_value(std::vector<double>{});
    std::vector<double> x_colmaj(6 * 4);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 11, 6, 4);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear_qrsolve);
    }
}


TEST(linear_algebra_test, qrsolve_test3)
{
    // 6 x 11 matrix with 7 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value(std::vector<double>{});
    auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/b")->get_value(std::vector<double>{});

    const auto reference_x = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/x")->get_value(std::vector<double>{});
    std::vector<double> x_colmaj(11 * 7);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 6, 11, 7);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear_qrsolve);
    }
}


TEST(linear_algebra_test, rcond_minitest0)
{
    const auto n = 0;
    const auto A_colmaj = std::vector<double>{};
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_TRUE(std::isinf(rcond(A_colmaj.data(), n)));
}


TEST(linear_algebra_test, rcond_minitest1)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ 20 };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 1, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest2)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ 0 };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 0, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest3)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ NaN };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_TRUE(std::isnan(rcond(A_colmaj.data(), n)));
}


TEST(linear_algebra_test, rcond_minitest4)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ -NaN };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_TRUE(std::isnan(rcond(A_colmaj.data(), n)));
}


TEST(linear_algebra_test, rcond_minitest5)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ Inf };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 0, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest6)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ -Inf };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 0, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest7)
{
    const auto n = 3;
    const auto A_colmaj = std::vector<double>{
        Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_TRUE(std::isnan(rcond(A_colmaj.data(), n)));
}


TEST(linear_algebra_test, rcond_minitest8)
{
    const auto n = 5;
    const auto A_colmaj = std::vector<double>{
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
        NaN, NaN, NaN };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_TRUE(std::isnan(rcond(A_colmaj.data(), n)));
}


TEST(linear_algebra_test, rcond_minitest9)
{
    const auto n = 2;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest10)
{
    const auto n = 2;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest11)
{
    const auto n = 3;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest12)
{
    const auto n = 30;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest13)
{
    const auto n = 10;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest14)
{
    const auto n = 5;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest15)
{
    const auto n = 3;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_minitest16)
{
    const auto n = 5;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_randtest0)
{
    const auto n = 263;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_randtest1)
{
    const auto n = 52;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_randtest2)
{
    const auto n = 357;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_randtest3)
{
    const auto n = 474;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_randtest4)
{
    const auto n = 434;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest0)
{
    const auto n = 281;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest1)
{
    const auto n = 396;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest2)
{
    const auto n = 456;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest3)
{
    const auto n = 478;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest4)
{
    const auto n = 334;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest5)
{
    const auto n = 213;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest6)
{
    const auto n = 175;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest7)
{
    const auto n = 444;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest8)
{
    const auto n = 368;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, rcond_gallerytest9)
{
    const auto n = 197;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A")->get_value());
    const auto reference_rcond = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/rcond")->get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear_rcond);
}


TEST(linear_algebra_test, tridiagonalsolve_tests)
{
    auto testcount = 0u;
    while (true) {
        const auto test_name = std::string{ "tridiagonalsolve_test" } + std::to_string(testcount);

        if (!linear_algebra_test_reference_file.get_root_element()->has_child(test_name)) {
            break;
        }
        else {
            ++testcount;
        }

        const auto num_rows_A_rows_B = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/num_rows_A_rows_B")->get_value(int{});
        const auto num_cols_B = linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/num_cols_B")->get_value(int{});
        const auto A_a = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A_a")->get_value());
        const auto A_b = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A_b")->get_value());
        const auto A_c = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/A_c")->get_value());
        const auto B_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/B_colmaj")->get_value());
        const auto reference_X_colmaj = comma_separated_string2vector_of_doubles(linear_algebra_test_reference_file.get_root_element()->get_child(test_name + "/X_colmaj")->get_value());

        ASSERT_EQ(A_a.size(), static_cast<std::size_t>(num_rows_A_rows_B - 1));
        ASSERT_EQ(A_b.size(), static_cast<std::size_t>(num_rows_A_rows_B));
        ASSERT_EQ(A_c.size(), static_cast<std::size_t>(num_rows_A_rows_B - 1));
        ASSERT_EQ(B_colmaj.size(), static_cast<std::size_t>(num_rows_A_rows_B * num_cols_B));
        ASSERT_EQ(reference_X_colmaj.size(), static_cast<std::size_t>(num_rows_A_rows_B * num_cols_B));

        auto X_colmaj = B_colmaj;
        auto A_b_copy = A_b;
        tridiagonalsolve(X_colmaj.data(), A_b_copy.data(), A_a.data(), A_c.data(),
                         num_rows_A_rows_B, num_cols_B);

        EXPECT_EQ(X_colmaj.size(), reference_X_colmaj.size());
        for (auto i = 0u; i < reference_X_colmaj.size(); ++i) {
            EXPECT_NEAR(X_colmaj[i], reference_X_colmaj[i], tolnear_tridiagonalsolve);
        }
    }

    ASSERT_EQ(testcount, 6u);
}


TEST(linear_algebra_test, lufactorize_and_lusolve_Doolittle)
{
    //
    // Tests functions "lufactorize_Doolittle" and
    // "lusolve_Doolittle" for some small matrices.
    // These functions are not safe since they don't
    // do pivoting, but they should suffice for small
    // enough matrices that we trust to be non-singular.
    //

    // the 10 x 10 problem with 3 columns in the rhs from "lusolve_test1"
    auto A_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child("lusolve_test2/A")->get_value(std::vector<double>{});
    auto b_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child("lusolve_test2/b")->get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size()) / 3;
    EXPECT_EQ(n, 10);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    const auto reference_x_colmaj = linear_algebra_test_reference_file.get_root_element()->get_child("lusolve_test2/x")->
        get_value(std::vector<double>{});

    {
        auto LU_A = A_colmaj;
        auto b = b_colmaj;
        lufactorize_Doolittle(LU_A.data(), n);
        lusolve_Doolittle(b.data(), LU_A.data(), n, 3);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear_lusolve);
        }
    }


    // some "n x n" matrices, n <= 7, with variable number "m" of columns, skipped when singular
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<scalar> distr(0., 1.);
    const auto rand_in_limits = [&](scalar lo, scalar hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto reltolnear = 1e-8;
    constexpr auto num_tests = 10000u;
    for (auto t = 0u; t < num_tests; ) {
        const auto n = static_cast<std::size_t>(std::ceil(rand_in_limits(0., 7.)));
        const auto m = static_cast<std::size_t>(std::ceil(rand_in_limits(0., 15.)));

        std::vector<double> A(n * n);
        std::generate(A.begin(), A.end(), [&]() { return rand_in_limits(-100., 100.); });

        std::vector<double> B(n * m);
        std::generate(B.begin(), B.end(), [&]() { return rand_in_limits(-100., 100.); });

        auto A_plu = A;
        auto X_plu = B;
        if (lusolve(X_plu.data(), A_plu.data(), static_cast<int>(n), static_cast<int>(m)) != 0) {
            continue;
        }
        else {
            ++t;
        }

        lufactorize_Doolittle(A.data(), n);
        lusolve_Doolittle(B.data(), A.data(), n, m);
        for (auto k = 0u; k < n * m; ++k) {
            EXPECT_NEAR(B[k], X_plu[k],
                        reltolnear * std::max(1., std::max(std::abs(B[k]), std::abs(X_plu[k]))));
        }
    }
}


TEST(linear_algebra_test, lufactorize_and_lusolve_measure_times)
{
    //
    // Times different ways in which we can invert a small dense matrix.
    //

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<scalar> distr(0., 1.);
    const auto rand_in_limits = [&](scalar lo, scalar hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto num_outer_tests = 1000u;
    constexpr auto num_inner_tests = 2u;
    static constexpr auto n = 6u;

#if defined _OPENMP
    const auto tbegin_plu = omp_get_wtime();
#endif
    for (auto t = 0u; t < num_outer_tests; ++t) {
        std::array<double, n * n> g{};
        std::array<int, n> ipiv_g{};
        std::for_each(g.begin(), g.end(), [&](auto &gi) { gi = rand_in_limits(-100., 100.); });

        plufactorize(&g[0], &ipiv_g[0], n);

        std::array<double, n> b{};
        for (auto t = 0u; t < num_inner_tests; ++t) {
            std::for_each(b.begin(), b.end(), [&](auto &bi) { bi = rand_in_limits(-100., 100.); });
            plusolve(&b[0], &g[0], &ipiv_g[0], n, 1);
        }
    }

#if defined _OPENMP
    std::cout << "PLU done in t = " << omp_get_wtime() - tbegin_plu << std::endl;
#endif


#if defined _OPENMP
    const auto tbegin_lu_Dolittle = omp_get_wtime();
#endif
    for (auto t = 0u; t < num_outer_tests; ++t) {
        std::array<double, n * n> g{};
        std::for_each(g.begin(), g.end(), [&](auto &gi) { gi = rand_in_limits(-100., 100.); });

        lufactorize_Doolittle(&g[0], n);

        std::array<double, n> b{};
        for (auto t = 0u; t < num_inner_tests; ++t) {
            std::for_each(b.begin(), b.end(), [&](auto &bi) { bi = rand_in_limits(-100., 100.); });
            lusolve_Doolittle(&b[0], &g[0], n, 1);
        }
    }

#if defined _OPENMP
    std::cout << "LU-Doolittle done in t = " << omp_get_wtime() - tbegin_lu_Dolittle << std::endl;
#endif


#if defined _OPENMP
    const auto tbegin_inv = omp_get_wtime();
#endif
    for (auto t = 0u; t < num_outer_tests; ++t) {
        std::array<double, n * n> g{};
        std::for_each(g.begin(), g.end(), [&](auto &gi) { gi = rand_in_limits(-100., 100.); });

        const auto ginv = inv6x6(g);
        std::array<double, n> b{};
        for (auto t = 0u; t < num_inner_tests; ++t) {
            std::array<double, n> x{};
            std::for_each(b.begin(), b.end(), [&](auto &bi) { bi = rand_in_limits(-100., 100.); });
            for (auto j = 0u; j < n; ++j) {
                for (auto i = 0u; i < n; ++i) {
                    x[i] += ginv[i + j * n] * b[j];
                }
            }

            // validate...
            //auto gLU = g;
            //lusolve(&b[0], &gLU[0], n, 1);
            //for (auto k = 0u; k < n; ++k) {
            //    EXPECT_NEAR(b[k], x[k], 1e-10 * std::max(1., std::max(std::abs(b[k]), std::abs(x[k]))));
            //}
        }
    }

#if defined _OPENMP
    std::cout << "inv6x6 done in t = " << omp_get_wtime() - tbegin_inv << std::endl;
#endif
}


TEST(linear_algebra_test, lufactorize_and_lusolve_as_cppad_adfun)
{
    //
    // A concept to tape the solution "y(x)" of a linear system
    // "A(x) * y = B(x)" as a CppAD::ADFun.
    //

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<scalar> distr(0., 1.);
    const auto rand_in_limits = [&](scalar lo, scalar hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    const auto n = 6u;
    const auto m = 5u;
    const auto nx = 4u;
    const auto np = 3u;
    const auto maxpow = 5;
    const auto D_A = 4.;

    const auto nn = n * n;
    std::vector<scalar> Acoeffs(nn * np);
    std::vector<int> Apows(nn * np * nx);
    for (auto i = 0u; i < nn; ++i) {
        for (auto j = 0u; j < np; ++j) {
            Acoeffs.at(i + nn * j) = rand_in_limits(-10., 10.);
            for (auto k = 0u; k < nx; ++k) {
                Apows.at(i + nn * (j + np * k)) = static_cast<int>(std::ceil(rand_in_limits(0., maxpow)));
            }
        }
    }

    const auto nm = n * m;
    std::vector<scalar> Bcoeffs(nm * np);
    std::vector<int> Bpows(nm * np * nx);
    for (auto i = 0u; i < nm; ++i) {
        for (auto j = 0u; j < np; ++j) {
            Bcoeffs.at(i + nm * j) = rand_in_limits(-10., 10.);
            for (auto k = 0u; k < nx; ++k) {
                Bpows.at(i + nm * (j + np * k)) = static_cast<int>(std::ceil(rand_in_limits(0., maxpow)));
            }
        }
    }

    const auto F_x = [](const auto &x, auto nF, const auto &coeffs, const auto &pows)
    {
        std::vector<typename std::decay_t<decltype(x)>::value_type> F(nF, scalar{ 0 });
        const auto np = coeffs.size() / nF;
        for (auto i = 0u; i < nF; ++i) {
            for (auto j = 0u; j < np; ++j) {
                typename decltype(F)::value_type monomial = coeffs.at(i + nF * j);
                for (auto k = 0u; k < x.size(); ++k) {
                    monomial *= CppAD::pow(x.at(k), pows.at(i + nF * (j + np * k)));
                }
                F.at(i) += monomial;
            }
        }
        return F;
    };

    const auto AF_x = [F_x, D_A](const auto &x, auto nF, const auto &coeffs, const auto &pows)
    {
        auto A = F_x(x, nF, coeffs, pows);
        const auto n = static_cast<std::size_t>(std::sqrt(static_cast<scalar>(nF)));
        for (auto i = 0u; i < n; ++i) {
            A[i + n * i] += D_A;
        }
        return A;
    };

    std::vector<CppAD::AD<scalar> > a_x(nx, 0.1);

    CppAD::ADFun<scalar> adfun_b;
    {
        CppAD::Independent(a_x);
        auto a_b = F_x(a_x, nm, Bcoeffs, Bpows);
        adfun_b.Dependent(a_x, a_b);
        adfun_b.optimize();
    }

    CppAD::ADFun<scalar> adfun_A;
    {
        CppAD::Independent(a_x);
        auto a_A = AF_x(a_x, nn, Acoeffs, Apows);
        adfun_A.Dependent(a_x, a_A);
        adfun_A.optimize();
    }

    CppAD::ADFun<scalar> adfun_Doolittle;
    {
        CppAD::Independent(a_x);
        auto a_A = AF_x(a_x, nn, Acoeffs, Apows);
        auto a_y = F_x(a_x, nm, Bcoeffs, Bpows);
        lufactorize_Doolittle(&a_A[0], n);
        lusolve_Doolittle(&a_y[0], &a_A[0], n, m);
        adfun_Doolittle.Dependent(a_x, a_y);
        adfun_Doolittle.optimize();
    }

#ifdef WITH_EIGEN
    atomic_eigen_mat_mul<scalar> mat_mul;
    atomic_eigen_mat_inv<scalar> mat_inv;
    CppAD::ADFun<scalar> adfun_eigen;
    {
        using ad_matrix = typename atomic_eigen_mat_inv<scalar>::ad_matrix;

        CppAD::Independent(a_x);
        const auto a_A = AF_x(a_x, nn, Acoeffs, Apows);
        const auto a_b = F_x(a_x, nm, Bcoeffs, Bpows);
        ad_matrix a_eigen_A(n, n);
        ad_matrix a_eigen_b(n, m);
        for (auto i = 0u; i < n; ++i) {
            for (auto j = 0u; j < n; ++j) {
                a_eigen_A(i, j) = a_A[i + j * n];
            }
        }
        for (auto i = 0u; i < n; ++i) {
            for (auto j = 0u; j < m; ++j) {
                a_eigen_b(i, j) = a_b[i + j * n];
            }
        }
        ad_matrix a_eigen_invA = mat_inv.op(a_eigen_A);
        ad_matrix a_eigen_y = mat_mul.op(a_eigen_invA, a_eigen_b);
        std::vector<CppAD::AD<scalar> > a_y(nm);
        for (auto j = 0u; j < m; ++j) {
            for (auto i = 0u; i < n; ++i) {
                a_y[i + j * n] = a_eigen_y(i, j);
            }
        }
        adfun_eigen.Dependent(a_x, a_y);
        // adfun_eigen.optimize(); // CRASHES!!
    }
#endif

    const auto reltolnear = 1e-9;
    const auto num_tests = 1000;
    for (auto t = 0u; t < num_tests; ++t) {
        std::vector<scalar> x(nx);
        std::generate(x.begin(), x.end(), [&]() { return rand_in_limits(-1., 1.); });

        // val
        const auto A = AF_x(x, nn, Acoeffs, Apows);
        const auto b = F_x(x, nm, Bcoeffs, Bpows);
        const auto invA = inv6x6(A);
        std::vector<scalar> val(nm);
        for (auto k = 0u; k < m; ++k) {
            const auto kn = k * n;
            for (auto j = 0u; j < n; ++j) {
                for (auto i = 0u; i < n; ++i) {
                    val[i + kn] += invA[i + j * n] * b[j + kn];
                }
            }
        }

        const auto val_Doolittle = adfun_Doolittle.Forward(0, x);
        ASSERT_EQ(val.size(), val_Doolittle.size());
        for (auto i = 0u; i < val.size(); ++i) {
            EXPECT_NEAR(val[i], val_Doolittle[i],
                        reltolnear * std::max(1., std::max(std::abs(val[i]), std::abs(val_Doolittle[i]))));
        }

#ifdef WITH_EIGEN
        const auto val_eigen = adfun_eigen.Forward(0, x);
        ASSERT_EQ(val.size(), val_eigen.size());
        for (auto i = 0u; i < val.size(); ++i) {
            EXPECT_NEAR(val[i], val_eigen[i],
                reltolnear * std::max(1., std::max(std::abs(val[i]), std::abs(val_eigen[i]))));
        }
#endif


        // jac
        const auto jacA_rowmaj = adfun_A.Jacobian(x);
        const auto jacb_rowmaj = adfun_b.Jacobian(x);
        const auto jac_Doolittle_rowmaj = adfun_Doolittle.Jacobian(x);
#ifdef WITH_EIGEN
        const auto jac_eigen_rowmaj = adfun_eigen.Jacobian(x);
#endif
        for (auto xpos = 0u; xpos < nx; ++xpos) {
            std::vector<scalar> dA_dx(nn);
            for (auto k = 0u; k < nn; ++k) {
                dA_dx[k] = jacA_rowmaj[nx * k + xpos];
            }

            std::vector<scalar> db_dx(nm);
            for (auto k = 0u; k < nm; ++k) {
                db_dx[k] = jacb_rowmaj[nx * k + xpos];
            }

            std::vector<scalar> val_times_dA_dx(nm, 0.);
            for (auto k = 0u; k < m; ++k) {
                const auto kn = k * n;
                for (auto j = 0u; j < n; ++j) {
                    for (auto i = 0u; i < n; ++i) {
                        val_times_dA_dx[i + kn] += dA_dx[i + j * n] * val[j + kn];
                    }
                }
            }

            std::vector<scalar> rhs(nm);
            for (auto k = 0u; k < nm; ++k) {
                rhs[k] = db_dx[k] - val_times_dA_dx[k];
            }

            std::vector<scalar> dval_dx(nm, 0.);
            for (auto k = 0u; k < m; ++k) {
                const auto kn = k * n;
                for (auto j = 0u; j < n; ++j) {
                    for (auto i = 0u; i < n; ++i) {
                        dval_dx[i + kn] += invA[i + j * n] * rhs[j + kn];
                    }
                }
            }

            std::vector<scalar> dval_dx_Doolittle(nm);
            for (auto k = 0u; k < nm; ++k) {
                dval_dx_Doolittle[k] = jac_Doolittle_rowmaj[nx * k + xpos];
            }
            for (auto i = 0u; i < nm; ++i) {
                EXPECT_NEAR(dval_dx[i], dval_dx_Doolittle[i],
                            reltolnear * std::max(1., std::max(std::abs(dval_dx[i]), std::abs(dval_dx_Doolittle[i]))));
            }

#ifdef WITH_EIGEN
            std::vector<scalar> dval_dx_eigen(nm);
            for (auto k = 0u; k < nm; ++k) {
                dval_dx_eigen[k] = jac_eigen_rowmaj[nx * k + xpos];
            }
            for (auto i = 0u; i < nm; ++i) {
                EXPECT_NEAR(dval_dx[i], dval_dx_eigen[i],
                    reltolnear * std::max(1., std::max(std::abs(dval_dx[i]), std::abs(dval_dx_eigen[i]))));
            }
#endif
        }

#ifdef WITH_EIGEN
        // hes
        for (auto k = 0u; k < nm; ++k) {
            const auto hesk_Doolittle = adfun_Doolittle.Hessian(x, k);
            const auto hesk_eigen = adfun_eigen.Hessian(x, k);

            ASSERT_EQ(hesk_Doolittle.size(), nx * nx);
            ASSERT_EQ(hesk_Doolittle.size(), hesk_eigen.size());
            for (auto i = 0u; i < nx * nx; ++i) {
                EXPECT_NEAR(hesk_Doolittle[i], hesk_eigen[i],
                            50. * reltolnear * std::max(1., std::max(std::abs(hesk_Doolittle[i]), std::abs(hesk_eigen[i]))));
            }
        }
#endif
    }
}
