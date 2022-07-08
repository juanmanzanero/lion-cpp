#include "lion/io/database_parameters.h"
#include "gtest/gtest.h"

class Database_parameter_test : public ::testing::Test
{
 protected:
    Database_parameter_test() { doc.load(); }
    Xml_document doc = {"data/example.xml"};
};


TEST_F(Database_parameter_test, double_parameter)
{
    class Sample_class
    {
     public:
        double d1;
        double d2;
        DECLARE_PARAMS({"double_parameter",d1},{"double_parameter_2",d2});
    } s;

    EXPECT_EQ(s.__used_parameters.size(), 2);
    EXPECT_FALSE(s.__used_parameters[0]);
    EXPECT_FALSE(s.__used_parameters[1]);

    read_parameters(doc, "xml_doc/parameters/", s.get_parameters(), s.__used_parameters);

    EXPECT_TRUE(s.__used_parameters[0]);
    EXPECT_TRUE(s.__used_parameters[1]);

    EXPECT_DOUBLE_EQ(s.d1, 3.14);
    EXPECT_DOUBLE_EQ(s.d2, 6.28);
}


TEST_F(Database_parameter_test, int_parameter)
{

    class Sample_class
    {
     public:
        int i;
        DECLARE_PARAMS({"int_parameter",i});
    } s;

    EXPECT_EQ(s.__used_parameters.size(), 1);
    EXPECT_FALSE(s.__used_parameters[0]);

    read_parameters(doc, "xml_doc/parameters/", s.get_parameters(), s.__used_parameters );

    EXPECT_TRUE(s.__used_parameters[0]);
    EXPECT_DOUBLE_EQ(s.i, 100.0);
}

TEST_F(Database_parameter_test, vector_parameter)
{
    class Sample_class
    {
     public:
        std::vector<double> v;
        DECLARE_PARAMS({"vector_parameter",v});
    } s;

    std::vector<double> v_expected = { 1.0, 3.0, 5.0, 5.0, 6.0, -6 };

    EXPECT_EQ(s.__used_parameters.size(), 1);
    EXPECT_FALSE(s.__used_parameters[0]);

    read_parameters(doc, "xml_doc/parameters/", s.get_parameters(), s.__used_parameters );

    EXPECT_TRUE(s.__used_parameters[0]);
    EXPECT_EQ(s.v.size(), v_expected.size());

    for (size_t i = 0; i < s.v.size(); ++i)
        EXPECT_DOUBLE_EQ(s.v[i], v_expected[i]);
}


TEST_F(Database_parameter_test, vector3_parameter)
{
    class Sample_class
    {
     public:
        sVector3d v;
        DECLARE_PARAMS({ "vector3_parameter", v });
    } s;

    sVector3d v_expected = { 0.6,0.8,-1.0 };

    EXPECT_EQ(s.__used_parameters.size(), 1);
    EXPECT_FALSE(s.__used_parameters[0]);

    read_parameters(doc, "xml_doc/parameters/", s.get_parameters(), s.__used_parameters);

    EXPECT_TRUE(s.__used_parameters[0]);

    for (size_t i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(s.v[i], v_expected[i]);
}


TEST_F(Database_parameter_test, matrix3x3_parameter)
{
    class Sample_class
    {
     public:
        sMatrix3x3 m;
        DECLARE_PARAMS({ "matrix3_parameter", m });
    } s;


    sMatrix3x3 m_expected = { 1,2,3,4,5,6,7,8,9 };
    EXPECT_EQ(s.__used_parameters.size(), 1);
    EXPECT_FALSE(s.__used_parameters[0]);

    read_parameters(doc, "xml_doc/parameters/", s.get_parameters(), s.__used_parameters);

    EXPECT_TRUE(s.__used_parameters[0]);
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(s.m(i,j), m_expected(i,j));
}
