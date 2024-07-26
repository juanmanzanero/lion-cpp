#include "lion/io/database_parameters.h"
#include "gtest/gtest.h"

class Database_parameter_test : public ::testing::Test
{
 protected:
    Database_parameter_test()
    {
        doc.load();
        doc2.load();
    }
    Xml_document doc = {"data/example.xml"};
    Xml_document doc2 = {"data/example2.xml"};
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

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters);

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

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters );

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

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters );

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

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters);

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

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters);

    EXPECT_TRUE(s.__used_parameters[0]);
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(s.m(i,j), m_expected(i,j));
}

TEST_F(Database_parameter_test, string_parameter)
{
    class Sample_class
    {
     public:
        std::string s;
        DECLARE_PARAMS({ "string_parameter", s });
    } s;


    std::string s_expected = " String parameter test ";
    EXPECT_EQ(s.__used_parameters.size(), 1);
    EXPECT_FALSE(s.__used_parameters[0]);

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters);

    EXPECT_TRUE(s.__used_parameters[0]);
    EXPECT_EQ(s.s, s_expected);
}


TEST_F(Database_parameter_test, unused_parameter)
{
    class Sample_class
    {
     public:
        double d1;
        double d2;
        DECLARE_PARAMS({"double_parameter",d1},{"double_parameter_2",d2});
    } s;

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters);

    Xml_element element = doc.get_element("xml_doc/parameters").cast<Xml_element>();
    std::string attribute, expected("false");

    for (size_t i = 0; i < s.__get_parameters().size(); ++i)
    {
        element = doc.get_element("xml_doc/parameters/"+ s.__get_parameters()[i].name).cast<Xml_element>();
        // Check until you get
        while(element.has_parent())
        {
            attribute = element.get_attribute("__unused__");
            EXPECT_EQ(attribute, expected);
            element = element.get_parent();
        }
    }
}

TEST_F(Database_parameter_test, unusued_parameter_2)
{
    class Sample_class
    {
     public:
        double d1;
        double d2;
        DECLARE_PARAMS({"double_parameter",d1});
    } s;

    read_parameters(doc, "xml_doc/parameters/", s.__get_parameters(), s.__used_parameters);

    Xml_element element = doc.get_element("xml_doc/parameters/double_parameter_2").cast<Xml_element>();
    bool attribute_existance = element.has_attribute("__unused__");
    EXPECT_FALSE(attribute_existance);

    element = doc.get_element("xml_doc/parameters/double_parameter").cast<Xml_element>();
    attribute_existance = element.has_attribute("__unused__");
    EXPECT_TRUE(attribute_existance);
}

TEST_F(Database_parameter_test, parameters_all_used_check)
{
    class Sample_class
    {
     public:
        double d1, d2, d3, d4;
        int i, i1;
        std::vector<double> v;
        sVector3d v2;
        sMatrix3x3 m, m1;
        DECLARE_PARAMS({"double_parameter",d1}, {"double_parameter_2",d2}, {"int_parameter",i},
            {"vector_parameter",v}, {"vector3_parameter",v2}, {"matrix3_parameter",m},
            {"child1/double_parameter", d3}, {"child1/child2/double_parameter", d4},
            {"child1/child2/child3/matrix3_parameter", m1}, {"child1/child2/child3/child4/int_parameter", i1}
        );
    } s;

    read_parameters(doc2, "xml_doc/", s.__get_parameters(), s.__used_parameters);
    auto root = doc2.get_root_element().cast<Xml_element>();
    bool parameters_all_used = database_parameters_all_used(root);
    EXPECT_TRUE(parameters_all_used);
}

TEST_F(Database_parameter_test, parameters_all_used_check_2)
{
    class Sample_class
    {
     public:
        double d1;
        DECLARE_PARAMS({"double_parameter",d1});
    } s;

    read_parameters(doc2, "xml_doc/", s.__get_parameters(), s.__used_parameters);
    auto root = doc2.get_root_element().cast<Xml_element>();
    bool parameters_all_used = database_parameters_all_used(root);
    EXPECT_FALSE(parameters_all_used);
}
