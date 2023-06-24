#include "gtest/gtest.h"
#include "lion/io/Xml_document.h"

class Xml_test : public ::testing::Test
{
 protected:
    Xml_document doc = {"data/example.xml"};
};


TEST_F(Xml_test, open_document)
{
    try
    {
        doc.load();
        SUCCEED();
    }
    catch(...)
    {
        FAIL();
    }
}


TEST_F(Xml_test, root_node)
{
    doc.load();
    Xml_element root = doc.get_root_element();

    EXPECT_EQ(root.get_name(), "xml_doc");

    EXPECT_EQ(root.get_value(std::string()), "");
}


TEST_F(Xml_test, root_node_children)
{
    doc.load();
    auto root_children = doc.get_root_element().get_children();

    EXPECT_EQ(root_children.size(), 3);
    EXPECT_EQ(root_children[0].get_name(), "child1");
    EXPECT_EQ(root_children[1].get_name(), "child2");
    EXPECT_EQ(root_children[2].get_name(), "parameters");
}


TEST_F(Xml_test, read_vector1)
{
    doc.load();

    // Find node using full path
    auto vector_child = doc.get_root_element().get_child("child2/vector1");

    std::vector<double> v = vector_child.get_value(std::vector<double>());
    std::vector<double> v_expected = { 1, 3, 5, 7 };

    EXPECT_EQ(v.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v[i], v_expected[i]);

    // Find node by node
    auto child2 = doc.get_root_element().get_child("child2");
    auto vector1 = child2.get_child("vector1");

    std::vector<double> v2 = vector1.get_value(std::vector<double>());

    EXPECT_EQ(v2.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v2[i], v_expected[i]);
    
}

TEST_F(Xml_test, read_vector2)
{
    doc.load();

    // Find node using full path
    auto vector_child = doc.get_root_element().get_child("child2/vector2");

    std::vector<double> v = vector_child.get_value(std::vector<double>());
    std::vector<double> v_expected = { 1, 2, 3, 4 };

    EXPECT_EQ(v.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v[i], v_expected[i]);

    // Find node by node
    auto child2 = doc.get_root_element().get_child("child2");
    auto vector1 = child2.get_child("vector2");

    std::vector<double> v2 = vector1.get_value(std::vector<double>());

    EXPECT_EQ(v2.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v2[i], v_expected[i]);
    
}


TEST_F(Xml_test, Matrix3x3)
{
    doc.load();

    sMatrix3x3 mat_expected = 
    {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0 
    };

    sMatrix3x3 mat = doc.get_root_element().get_child("child2/matrix3x3").get_value(sMatrix3x3());

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(mat_expected(i,j), mat(i,j));
}


TEST_F(Xml_test, print_doc)
{
    doc.load();
    std::ostringstream sOut, sOutRef;

    sOutRef << "<xml_doc>"                                                        << std::endl;
    sOutRef << "    <child1> Text for child 1 </child1>"                          << std::endl;
    sOutRef << "    <child2>"                                                     << std::endl;
    sOutRef << "        <vector1>1.0 3.0 5.0 7.0</vector1>"                       << std::endl;
    sOutRef << "        <vector2>"                                                << std::endl;
    sOutRef << "            1.0"                                                  << std::endl;
    sOutRef << "            2.0"                                                  << std::endl;
    sOutRef << "            3.0"                                                  << std::endl;
    sOutRef << "            4.0"                                                  << std::endl;
    sOutRef << "        </vector2>"                                               << std::endl;
    sOutRef << "        <matrix3x3>"                                              << std::endl;
    sOutRef << "            1.0 2.0 3.0"                                          << std::endl;
    sOutRef << "            4.0 5.0 6.0"                                          << std::endl;
    sOutRef << "            7.0 8.0 9.0"                                          << std::endl;
    sOutRef << "        </matrix3x3>"                                             << std::endl;
    sOutRef << "    </child2>"                                                    << std::endl;
    sOutRef << "    <!-- Parameters for the parameters test -->"                  << std::endl;
    sOutRef << "    <parameters>"                                                 << std::endl;
    sOutRef << "        <double_parameter>3.14</double_parameter>"                << std::endl;
    sOutRef << "        <double_parameter_2>6.28</double_parameter_2>"            << std::endl;
    sOutRef << "        <int_parameter>100</int_parameter>"                       << std::endl;
    sOutRef << "        <vector_parameter>1.0 3.0 5.0 "                           << std::endl;
    sOutRef << "                          5.0 6.0 -6 </vector_parameter>"         << std::endl;
    sOutRef << "        <vector3_parameter>0.6 0.8 -1.0</vector3_parameter>"      << std::endl;
    sOutRef << "        <matrix3_parameter>1 2 3 4 5 6 7 8 9</matrix3_parameter>" << std::endl;
    sOutRef << "        <string_parameter> String parameter test </string_parameter>" << std::endl;
    sOutRef << "    </parameters>"                                                << std::endl;
    sOutRef << "</xml_doc>"                                                       << std::endl;
    sOut << doc ;

    EXPECT_EQ(sOut.str(), sOutRef.str());

}



TEST_F(Xml_test, print_element)
{
    doc.load();
    std::ostringstream sOut;
    sOut << "<child2>" << std::endl;
    sOut << "    <vector1>1.0 3.0 5.0 7.0</vector1>" << std::endl;
    sOut << "    <vector2>" << std::endl;
    sOut << "            1.0" << std::endl;
    sOut << "            2.0" << std::endl;
    sOut << "            3.0" << std::endl;
    sOut << "            4.0" << std::endl;
    sOut << "        </vector2>" << std::endl;
    sOut << "    <matrix3x3>" << std::endl;
    sOut << "            1.0 2.0 3.0" << std::endl;
    sOut << "            4.0 5.0 6.0" << std::endl;
    sOut << "            7.0 8.0 9.0" << std::endl;
    sOut << "        </matrix3x3>" << std::endl;
    sOut << "</child2>" << std::endl;

    std::ostringstream sOutXml;
    sOutXml << doc.get_root_element().get_child("child2");
    EXPECT_EQ(sOutXml.str(), sOut.str());


}

TEST_F(Xml_test, create_root)
{
    Xml_document doc_out("dummy.xml");
    doc_out.create_root_element("root");
    std::ostringstream sOut;
    sOut << doc_out;
    EXPECT_EQ(sOut.str(), "<root/>\n");
}


TEST_F(Xml_test, create_child)
{
    Xml_document doc_out("dummy.xml");
    doc_out.create_root_element("root");
    Xml_element child = doc_out.get_root_element().add_child("child");
    child.set_value("value for child");
    
    std::ostringstream sOut;
    sOut << doc_out;
    EXPECT_EQ(sOut.str(), "<root>\n    <child>value for child</child>\n</root>\n");
}

TEST_F(Xml_test, get_value_as_bool)
{
    Xml_document doc;
    doc.parse("<doc> <a> true </a> <b> false </b> </doc>");

    EXPECT_TRUE(doc.get_element("doc/a").get_value(bool()));
    EXPECT_FALSE(doc.get_element("doc/b").get_value(bool()));
}


TEST_F(Xml_test, add_child)
{
    Xml_document doc;
    doc.create_root_element("root");

    doc.add_element("root/a/b/../b/./c/d").set_value("my value");

    EXPECT_EQ(doc.get_element("root/a/b/c/d").get_value(), "my value");
}


TEST_F(Xml_test, deep_copy)
{
    Xml_document doc1;
    doc1.parse("<doc> <a> <a1> true </a1> <empty_node/> </a> </doc>");

    Xml_document doc2;
    doc2.parse("<doc2> <content1> <content2> content </content2> </content1> </doc2>");

    doc1.get_element("doc/a/empty_node").copy_contents(doc2.get_element("doc2"));

    std::ostringstream s_out_ref;
    s_out_ref << "<doc>" << std::endl;
    s_out_ref << "    <a>" << std::endl;
    s_out_ref << "        <a1> true </a1>" << std::endl;
    s_out_ref << "        <empty_node>" << std::endl;
    s_out_ref << "            <content1>" << std::endl;
    s_out_ref << "                <content2> content </content2>" << std::endl;
    s_out_ref << "            </content1>" << std::endl;
    s_out_ref << "        </empty_node>" << std::endl;
    s_out_ref << "    </a>" << std::endl;
    s_out_ref << "</doc>" << std::endl;

    std::ostringstream s_out;
    doc1.print(s_out);
    EXPECT_EQ(s_out.str(), s_out_ref.str());
    
}
