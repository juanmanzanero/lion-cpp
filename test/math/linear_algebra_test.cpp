#include <vector>

#include "gtest/gtest.h"

#include "lion/io/Xml_document.h"
#include "lion/math/linear_algebra.h"


//
// Defines tests for the functions in the
// "lion/math/linear_algebra.h" header.
//


constexpr auto tolnear = 1e3 * std::numeric_limits<double>::epsilon();


static Xml_document reference_file("data/linear_algebra_test.xml", true);


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
    const auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    const auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

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
    const auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    const auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size());
    EXPECT_EQ(n, 4);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    const auto reference_x_colmaj = reference_file.get_root_element().get_child(test_name + "/x").
        get_value(std::vector<double>{});

    {
        auto A = A_colmaj;
        auto b = b_colmaj;
        EXPECT_EQ(lusolve(b.data(), A.data(), n, 1), 0);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear);
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
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear);
        }
    }
}


TEST(linear_algebra_test, lusolve_test2)
{
    // 10 x 10 problem with 3 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size()) / 3u;
    EXPECT_EQ(n, 10);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    const auto reference_x_colmaj = reference_file.get_root_element().get_child(test_name + "/x").
        get_value(std::vector<double>{});

    {
        auto A = A_colmaj;
        auto b = b_colmaj;
        EXPECT_EQ(lusolve(b.data(), A.data(), n, 3), 0);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear);
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
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear);
        }
    }
}


TEST(linear_algebra_test, lusolve_test3)
{
    // 20 x 20 problem with 7 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

    const auto n = static_cast<int>(b_colmaj.size()) / 7u;
    EXPECT_EQ(n, 20);
    EXPECT_EQ(std::sqrt(A_colmaj.size()), n);

    const auto reference_x_colmaj = reference_file.get_root_element().get_child(test_name + "/x").
        get_value(std::vector<double>{});

    {
        auto A = A_colmaj;
        auto b = b_colmaj;
        EXPECT_EQ(lusolve(b.data(), A.data(), n, 7), 0);
        const auto &x = b;

        EXPECT_EQ(reference_x_colmaj.size(), x.size());
        for (auto i = 0u; i < reference_x_colmaj.size(); ++i) {
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear);
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
            EXPECT_NEAR(x[i], reference_x_colmaj[i], tolnear);
        }
    }
}


TEST(linear_algebra_test, qrsolve_test0)
{
    // 10 x 4 matrix with 1 column in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

    const auto reference_x = reference_file.get_root_element().get_child(test_name + "/x").get_value(std::vector<double>{});
    std::vector<double> x_colmaj(4);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 10, 4, 1);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear);
    }
}


TEST(linear_algebra_test, qrsolve_test1)
{
    // 4 x 10 matrix with 1 column in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

    const auto reference_x = reference_file.get_root_element().get_child(test_name + "/x").get_value(std::vector<double>{});
    std::vector<double> x_colmaj(10);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 4, 10, 1);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear);
    }
}


TEST(linear_algebra_test, qrsolve_test2)
{
    // 11 x 6 matrix with 4 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

    const auto reference_x = reference_file.get_root_element().get_child(test_name + "/x").get_value(std::vector<double>{});
    std::vector<double> x_colmaj(6 * 4);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 11, 6, 4);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear);
    }
}


TEST(linear_algebra_test, qrsolve_test3)
{
    // 6 x 11 matrix with 7 columns in the rhs
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    auto A_colmaj = reference_file.get_root_element().get_child(test_name + "/A").get_value(std::vector<double>{});
    auto b_colmaj = reference_file.get_root_element().get_child(test_name + "/b").get_value(std::vector<double>{});

    const auto reference_x = reference_file.get_root_element().get_child(test_name + "/x").get_value(std::vector<double>{});
    std::vector<double> x_colmaj(11 * 7);
    EXPECT_EQ(reference_x.size(), x_colmaj.size());
    qrsolve(x_colmaj.data(), A_colmaj.data(), b_colmaj.data(), 6, 11, 7);
    for (auto i = 0u; i < reference_x.size(); ++i) {
        EXPECT_NEAR(x_colmaj[i], reference_x[i], tolnear);
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

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 1, tolnear);
}


TEST(linear_algebra_test, rcond_minitest2)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ 0 };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 0, tolnear);
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

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 0, tolnear);
}


TEST(linear_algebra_test, rcond_minitest6)
{
    const auto n = 1;
    const auto A_colmaj = std::vector<double>{ -Inf };
    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));

    EXPECT_NEAR(rcond(A_colmaj.data(), n), 0, tolnear);
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
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_minitest10)
{
    const auto n = 2;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_minitest11)
{
    const auto n = 3;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_minitest12)
{
    const auto n = 30;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_minitest13)
{
    const auto n = 10;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_minitest14)
{
    const auto n = 5;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_minitest15)
{
    const auto n = 3;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_minitest16)
{
    const auto n = 5;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_randtest0)
{
    const auto n = 263;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_randtest1)
{
    const auto n = 52;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_randtest2)
{
    const auto n = 357;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_randtest3)
{
    const auto n = 474;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_randtest4)
{
    const auto n = 434;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest0)
{
    const auto n = 281;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest1)
{
    const auto n = 396;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest2)
{
    const auto n = 456;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest3)
{
    const auto n = 478;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest4)
{
    const auto n = 334;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest5)
{
    const auto n = 213;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest6)
{
    const auto n = 175;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest7)
{
    const auto n = 444;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest8)
{
    const auto n = 368;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}


TEST(linear_algebra_test, rcond_gallerytest9)
{
    const auto n = 197;
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const auto A_colmaj = comma_separated_string2vector_of_doubles(reference_file.get_root_element().get_child(test_name + "/A").get_value());
    const auto reference_rcond = reference_file.get_root_element().get_child(test_name + "/rcond").get_value(double{});

    EXPECT_EQ(n, static_cast<int>(std::sqrt(A_colmaj.size())));
    EXPECT_NEAR(rcond(A_colmaj.data(), n), reference_rcond, tolnear);
}
