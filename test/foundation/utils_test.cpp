#include "gtest/gtest.h"
#include "lion/foundation/utils.h"


TEST(utils_test, string_to_double_vector)
{
    const std::string s = "1.02324355, 5.0624692 \n -7.354634, \n 5.23035234392, 243.242532425";

    std::vector<double> v = string_to_double_vector(s);
  
    std::vector<double> v_expected = { 1.02324355,5.0624692,-7.354634,5.23035234392,243.242532425 };

    EXPECT_EQ(v.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v[i], v_expected[i]);

}
