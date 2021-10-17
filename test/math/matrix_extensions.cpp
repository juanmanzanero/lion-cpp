#include "gtest/gtest.h"
#include "lion/math/matrix_extensions.h"



class Matrix_operations_test : public ::testing::Test
{
 protected:
    const std::array<std::array<double,3>,3> a = 
            {{{ 0.7537,   0.0759,   0.7792 },
              { 0.3804,   0.0540,   0.9340 },
              { 0.5678,   0.5308,   0.1299 }}};


    const std::array<std::array<double,3>,3> b = 
            {{{ 0.5688,   0.3371,   0.3112 },
              { 0.4694,   0.1622,   0.5285 },
              { 0.0119,   0.7943,   0.1656 }}};

    const std::array<std::array<double,3>,3> c =
            {{{ 0.6020,   0.6892,   0.0838 },
              { 0.2630,   0.7482,   0.2290 },
              { 0.6541,   0.4505,   0.9133 }}}; 

    const std::array<std::array<double,3>,3> a_plus_b = 
       {{{ 1.3225,   0.4130,   1.0904 },
         { 0.8498,   0.2162,   1.4625 },
         { 0.5797,   1.3251,   0.2955 }}};

    const std::array<std::array<double,3>,3> a_plus_c = 
       {{{ 1.3557,   0.7651,   0.8630 },
         { 0.6434,   0.8022,   1.1630 },
         { 1.2219,   0.9813,   1.0432 }}};

    const std::array<std::array<double,3>,3> b_plus_c = 
       {{{ 1.1708,   1.0263,   0.3950 },
         { 0.7324,   0.9104,   0.7575 },
         { 0.6660,   1.2448,   1.0789 }}};

    const std::array<std::array<double,4>,2> u =
           {{{ 0.1524,   0.5383,   0.0782,   0.1067 },
             { 0.8258,   0.9961,   0.4427,   0.9619 }}};

    const std::array<std::array<double,3>,4> v =
       {{{ 0.0046,   0.0844,   0.4314 },
         { 0.7749,   0.3998,   0.9106 },
         { 0.8173,   0.2599,   0.1818 },
         { 0.8687,   0.8001,   0.2638 }}};

    const std::array<double,4> vec1 = {{0.0046, 0.7749, 0.8173, 0.8687}};

    const std::array<double,2> vec2 = {{ 0.8147, 0.9058 }};
};


TEST_F(Matrix_operations_test, matrix_addition)
{
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a+b)[i][j], a_plus_b[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a+c)[i][j], a_plus_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (b+c)[i][j], b_plus_c[i][j]);
}


TEST_F(Matrix_operations_test, matrix_subtraction)
{
    constexpr const std::array<std::array<double,3>,3> a_minus_b = 
      {{{  0.1849,  -0.2612,   0.4680 },
        { -0.0890,  -0.1082,   0.4055 },
        {  0.5559,  -0.2635,  -0.0357 }}};

    constexpr const std::array<std::array<double,3>,3> a_minus_c = 
      {{{  0.1517,  -0.6133,   0.6954 },
        {  0.1174,  -0.6942,   0.7050 },
        { -0.0863,   0.0803,  -0.7834 }}};

    constexpr const std::array<std::array<double,3>,3> b_minus_c = 
      {{{ -0.0332,  -0.3521,   0.2274 },
        {  0.2064,  -0.5860,   0.2995 },
        { -0.6422,   0.3438,  -0.7477 }}};


    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a-(-b))[i][j], a_plus_b[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a-(-c))[i][j], a_plus_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (b-(-c))[i][j], b_plus_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a-b)[i][j], a_minus_b[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a-c)[i][j], a_minus_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (b-c)[i][j], b_minus_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (b-a)[i][j], -a_minus_b[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (c-a)[i][j], -a_minus_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (c-b)[i][j], -b_minus_c[i][j]);
}

TEST_F(Matrix_operations_test, Matrix_product)
{
    constexpr const std::array<std::array<double,3>,3> a_times_b = 
        {{{  0.473604500,  0.885301810,  0.403700110 },
          {  0.252833720,  0.878867840,  0.301589880 },
          {  0.573667970,  0.380680710,  0.478738600 }}};

    constexpr const std::array<std::array<double,3>,3> a_times_c = 
      {{{ 0.983363820,  0.927268020,  0.792184520 },
        { 0.854132200,  0.723341480,  0.897265720 },
        { 0.566383590,  0.846992270,  0.287772510 }}};

    constexpr const std::array<std::array<double,3>,3> b_times_c = 
      {{{ 0.634630820,  0.784430780,  0.409080300 },
        { 0.670929250,  0.682957770,  0.559158570 },
        { 0.324383660,  0.677099540,  0.334134400 }}};

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a*b)[i][j], a_times_b[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (a*c)[i][j], a_times_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (b*c)[i][j], b_times_c[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (transpose(b)*transpose(a))[i][j], a_times_b[j][i]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (transpose(c)*transpose(a))[i][j], a_times_c[j][i]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (transpose(c)*transpose(b))[i][j], b_times_c[j][i]);

    constexpr const std::array<std::array<double,3>,2> u_times_v = 
       {{{ 0.574432860,  0.333769750,  0.598285560 },
         { 1.973097810,  1.352612220,  1.597530860 }}};

    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( (u*v)[i][j], u_times_v[i][j]);

    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ( ((-u)*v)[i][j], -u_times_v[i][j]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 2; ++j)
            EXPECT_DOUBLE_EQ( (transpose(v)*transpose(u))[i][j], u_times_v[j][i]);
}


TEST_F(Matrix_operations_test, Matrix_vector_products)
{
    constexpr const std::array<double,2> u_times_vec1 = {{ 0.57443286, 1.97309781 }};

    constexpr const std::array<double,4> vec2_times_u = {{0.87216992, 1.34082039, 0.46470720, 0.95821751}};

    for (size_t i = 0; i < 2; ++i)
        EXPECT_DOUBLE_EQ( (u*vec1).at(i), u_times_vec1.at(i) );

    for (size_t i = 0; i < 2; ++i)
        EXPECT_DOUBLE_EQ( (vec1*transpose(u)).at(i), u_times_vec1.at(i) );

    for (size_t i = 0; i < 4; ++i)
        EXPECT_DOUBLE_EQ( (vec2*u).at(i), vec2_times_u.at(i) );

    for (size_t i = 0; i < 4; ++i)
        EXPECT_DOUBLE_EQ( (transpose(u)*vec2).at(i), vec2_times_u.at(i) );

}




TEST_F(Matrix_operations_test, linsolve_10x10)
{
    std::array<std::array<double,10>,10> A = {{
                                              {0.06698,2.74249,6.50832,3.67703,2.71744,8.81633,2.94511,3.12523,6.15616,4.47479},
                                              {2.71338,3.80140,2.46664,1.60183,5.91685,4.49776,2.10519,1.61385,6.03287,7.13949},
                                              {8.01550,8.16740,0.50531,3.97456,8.66729,7.29586,1.63632,7.27620,2.38327,8.21052},
                                              {6.17009,7.83230,0.02025,1.64468,5.17865,1.32359,1.22712,2.63549,2.40153,7.03108},
                                              {8.28342,7.89850,5.23441,1.92203,4.14574,7.77246,7.47021,0.04083,1.93083,0.41663},
                                              {8.79967,5.75105,7.35315,2.23842,5.50054,3.48890,8.34791,5.49344,0.82640,8.99430},
                                              {8.66048,1.66353,5.53240,8.56250,7.13401,7.70530,7.81393,5.40602,2.02684,1.12470},
                                              {4.28197,2.32865,4.31847,7.54178,6.93049,6.23996,5.62008,8.69732,4.89048,2.17583},
                                              {6.58578,4.26765,5.14290,8.34002,4.81665,8.95123,4.54251,8.12440,5.11957,0.00861},
                                              {0.86879,1.27197,1.48969,6.06280,5.17434,5.33063,1.69299,7.32451,4.42106,0.50510}
                                             }};

    std::array<double,10> b = {6.67186,1.31601,0.68879,0.62785,4.68025,4.72920,8.54279,7.51548,1.09594,0.94729};

    auto x = linsolve(10,A,b);
    auto b1 = A*x;

    for (size_t i = 0; i < 10; ++i)
        EXPECT_NEAR(b1.at(i), b.at(i),5.0e-14);
}

TEST_F(Matrix_operations_test, Max_array_and_double)
{
    const std::array<double, 10> a = {1.0,2.0,3.0,4.0,-1.0,-6.0,5.0,-10.0,40.0,0.0};
    const double k = 0.1;

    const auto result = ext::max(a,k);

    const std::array<double,10> expected = {1.0,2.0,3.0,4.0,0.1,0.1,5.0,0.1,40.0,0.1};

    for (size_t i = 0; i < 10; ++i)
        EXPECT_DOUBLE_EQ(result[i], expected[i]);
}
