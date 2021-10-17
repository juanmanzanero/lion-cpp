#include "gtest/gtest.h" 
#include "lion/foundation/type_traits.h"
#include <string>


TEST(Type_traits_test, is_vector_test)
{
    EXPECT_EQ((is_vector<double>::value), false);
    EXPECT_EQ((is_vector<double*>::value), false);
    EXPECT_EQ((is_vector<double**>::value), false);
    EXPECT_EQ((is_vector<std::vector<double>>::value), true);
    EXPECT_EQ((is_vector<std::array<double,5>>::value), true);
    EXPECT_EQ((is_vector<std::string>::value), false);
    EXPECT_EQ((is_vector<std::vector<std::vector<double>>>::value), false);
    EXPECT_EQ((is_vector<std::array<std::array<double,5>,6>>::value), false);
    EXPECT_EQ((is_matrix<std::vector<std::vector<std::vector<double>>>>::value), false);
    EXPECT_EQ((is_matrix<std::array<std::array<std::array<double,2>,5>,6>>::value), false);
}

TEST(Type_traits_test, is_matrix_test)
{
    EXPECT_EQ((is_matrix<double>::value), false);
    EXPECT_EQ((is_matrix<double*>::value), false);
    EXPECT_EQ((is_matrix<double**>::value), false);
    EXPECT_EQ((is_matrix<std::vector<double>>::value), false);
    EXPECT_EQ((is_matrix<std::array<double,5>>::value), false);
    EXPECT_EQ((is_matrix<std::string>::value), false);
    EXPECT_EQ((is_matrix<std::vector<std::vector<double>>>::value), true);
    EXPECT_EQ((is_matrix<std::array<std::array<double,5>,6>>::value), true);
    EXPECT_EQ((is_matrix<std::vector<std::vector<std::vector<double>>>>::value), false);
    EXPECT_EQ((is_matrix<std::array<std::array<std::array<double,2>,5>,6>>::value), false);
}
