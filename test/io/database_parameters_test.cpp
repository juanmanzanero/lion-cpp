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
    double d;
    Database_parameter_mutable p = { "double_parameter", d };

    read_parameters(doc, "xml_doc/parameters/", {p} );

    EXPECT_DOUBLE_EQ(d, 3.14);
}


TEST_F(Database_parameter_test, int_parameter)
{
    int i;
    Database_parameter_mutable p = { "int_parameter", i };

    read_parameters(doc, "xml_doc/parameters/", {p} );

    EXPECT_DOUBLE_EQ(i, 100.0);
}

TEST_F(Database_parameter_test, vector_parameter)
{
    std::vector<double> v;
    std::vector<double> v_expected = { 1.0, 3.0, 5.0, 5.0, 6.0, -6 };
    Database_parameter_mutable p = { "vector_parameter", v };

    read_parameters(doc, "xml_doc/parameters/", {p} );

    EXPECT_EQ(v.size(), v_expected.size());

    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v[i], v_expected[i]);
}


TEST_F(Database_parameter_test, vector3_parameter)
{
    sVector3d v;
    sVector3d v_expected = { 0.6,0.8,-1.0 };
    Database_parameter_mutable p = { "vector3_parameter", v };

    read_parameters(doc, "xml_doc/parameters/", {p} );

    for (size_t i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(v[i], v_expected[i]);
}



TEST_F(Database_parameter_test, matrix3x3_parameter)
{
    sMatrix3x3 m;
    sMatrix3x3 m_expected = { 1,2,3,4,5,6,7,8,9 };
    Database_parameter_mutable p = { "matrix3_parameter", m };

    read_parameters(doc, "xml_doc/parameters/", {p} );

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(m(i,j), m_expected(i,j));
}


