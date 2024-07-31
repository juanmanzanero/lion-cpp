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
    auto root = doc.get_root_element();

    EXPECT_EQ(root->get_name(), "xml_doc");

    EXPECT_EQ(root->get_value(std::string()), "");
}


TEST_F(Xml_test, root_node_children)
{
    doc.load();
    auto root_children = doc.get_root_element()->get_children();

    EXPECT_EQ(root_children.size(), 3u);
    EXPECT_EQ(root_children[0]->get_name(), "child1");
    EXPECT_EQ(root_children[1]->get_name(), "child2");
    EXPECT_EQ(root_children[2]->get_name(), "parameters");
}


TEST_F(Xml_test, read_vector1)
{
    doc.load();

    // Find node using full path
    auto vector_child = doc.get_root_element()->get_child("child2/vector1");

    std::vector<double> v = vector_child->get_value(std::vector<double>());
    std::vector<double> v_expected = { 1, 3, 5, 7 };

    EXPECT_EQ(v.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v[i], v_expected[i]);

    // Find node by node
    auto child2 = doc.get_root_element()->get_child("child2");
    auto vector1 = child2->get_child("vector1");

    std::vector<double> v2 = vector1->get_value(std::vector<double>());

    EXPECT_EQ(v2.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v2[i], v_expected[i]);
    
}

TEST_F(Xml_test, read_vector2)
{
    doc.load();

    // Find node using full path
    auto vector_child = doc.get_root_element()->get_child("child2/vector2");

    std::vector<double> v = vector_child->get_value(std::vector<double>());
    std::vector<double> v_expected = { 1, 2, 3, 4 };

    EXPECT_EQ(v.size(), v_expected.size());
    
    for (size_t i = 0; i < v.size(); ++i)
        EXPECT_DOUBLE_EQ(v[i], v_expected[i]);

    // Find node by node
    auto child2 = doc.get_root_element()->get_child("child2");
    auto vector1 = child2->get_child("vector2");

    std::vector<double> v2 = vector1->get_value(std::vector<double>());

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

    sMatrix3x3 mat = doc.get_root_element()->get_child("child2/matrix3x3")->get_value(sMatrix3x3());

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(mat_expected(i,j), mat(i,j));
}


TEST_F(Xml_test, setget_stdarray_Matrix3x3_Vector3d)
{
    constexpr auto ref_array = std::array<double, 10>{ -0.747873558438421,
                                                       0.152910391091676,
                                                       0.020299178810689,
                                                       -0.462828559182863,
                                                       0.348012509296819,
                                                       0.963524857198698,
                                                      -0.848449653427868,
                                                       0.106551410006897,
                                                       -0.351271659559704,
                                                       0.214090352577798 };

    constexpr auto ref_matrix3x3 = sMatrix3x3{ -0.230420379098272, -0.729266329894491, -0.822025202317546,
                                               -0.794278216759767, 0.653874348334355, -0.840119645921260,
                                               -0.699372948125601, 0.179282943946838, 0.858659940414170 };

    constexpr auto ref_vector3d = sVector3d{ -0.383548029744673,
                                              0.278827481911275,
                                             -0.984115682466366 };

    Xml_document doc;
    doc.create_root_element("doc");
    doc.get_root_element()->add_child("array")->set_value(ref_array);
    doc.get_root_element()->add_child("Matrix3x3")->set_value(ref_matrix3x3);
    doc.get_root_element()->add_child("Vector3d")->set_value(ref_vector3d);

    const auto array_ = doc.get_root_element()->get_child("array")->get_value(std::remove_const_t<decltype(ref_array)>{});
    ASSERT_EQ(array_.size(), ref_array.size());
    for (auto i = 0u; i < array_.size(); ++i) {
        EXPECT_DOUBLE_EQ(array_[i], ref_array[i]);
    }

    const auto array_try = doc.get_root_element()->try_get_child_value("array", std::remove_const_t<decltype(ref_array)>{});
    for (auto i = 0u; i < array_try.size(); ++i) {
        EXPECT_NE(array_try[i], scalar{ 0 });
    }

    const auto matrix3x3 = doc.get_element("doc/Matrix3x3")->get_value(std::remove_const_t<decltype(ref_matrix3x3)>{});
    ASSERT_EQ(matrix3x3.base_type::size(), ref_matrix3x3.base_type::size());
    for (auto i = 0u; i < matrix3x3.base_type::size(); ++i) {
        EXPECT_DOUBLE_EQ(matrix3x3[i], ref_matrix3x3[i]);
    }

    const auto matrix3x3_try = doc.get_root_element()->try_get_child_value("Matrix3x3", std::remove_const_t<decltype(ref_matrix3x3)>{});
    for (auto i = 0u; i < matrix3x3_try.base_type::size(); ++i) {
        EXPECT_NE(matrix3x3_try[i], scalar{ 0 });
    }

    const auto vector3d = doc.get_root_element()->get_child("Vector3d")->get_value(std::remove_const_t<decltype(ref_vector3d)>{});
    ASSERT_EQ(vector3d.base_type::size(), ref_vector3d.base_type::size());
    for (auto i = 0u; i < vector3d.base_type::size(); ++i) {
        EXPECT_DOUBLE_EQ(vector3d[i], ref_vector3d[i]);
    }

    const auto vector3d_try = doc.get_root_element()->try_get_child_value("Vector3d", std::remove_const_t<decltype(ref_vector3d)>{});
    for (auto i = 0u; i < vector3d_try.base_type::size(); ++i) {
        EXPECT_NE(vector3d_try[i], scalar{ 0 });
    }
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
    sOutXml << doc.get_root_element()->get_child("child2").cast<Xml_element>();
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
    auto child = doc_out.get_root_element()->add_child("child");
    child->set_value("value for child");
    
    std::ostringstream sOut;
    sOut << doc_out;
    EXPECT_EQ(sOut.str(), "<root>\n    <child>value for child</child>\n</root>\n");
}


TEST_F(Xml_test, add_child)
{
    Xml_document doc;
    doc.create_root_element("root");

    doc.add_element("root/a/b/../b/./c/d").set_value("my value");

    EXPECT_EQ(doc.get_element("root/a/b/c/d")->get_value(), "my value");
}


TEST_F(Xml_test, deep_copy)
{
    Xml_document doc1;
    doc1.parse("<doc> <a> <a1> true </a1> <empty_node/> </a> </doc>");

    Xml_document doc2;
    doc2.parse("<doc2> <content1> <content2> content </content2> </content1> </doc2>");

    auto element_to_copy = doc2.get_element("doc2");
    doc1.get_element("doc/a/empty_node")->copy_contents(element_to_copy);

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


TEST_F(Xml_test, this_and_child_have_attributes)
{
    auto xmldef = "<test>"
        "<abc>"
            "<dfg>"
                "<qwer aaaa=\"\"/>"
            "</dfg>"
        "</abc>"
        "<qwe abc=\"auuu\">"
        "</qwe>"
    "</test>";

    Xml_document doc;
    doc.parse(xmldef);

    ASSERT_TRUE(doc.get_root_element()->this_and_childs_have_attributes());

    doc.get_element("test/abc/dfg/qwer").cast<Xml_element>().delete_attribute("aaaa");
    ASSERT_TRUE(doc.get_root_element()->this_and_childs_have_attributes());

    doc.get_element("test/qwe").cast<Xml_element>().delete_attribute("abc");
    ASSERT_FALSE(doc.get_root_element()->this_and_childs_have_attributes());
}


TEST_F(Xml_test, this_and_childs_have_both_value_and_children)
{

    auto xmldef = "<test>"
            "<abc>"
                "<auu>"
                "  valueee "
                "  <childd/>"
                "</auu>"
            "</abc>"
        "</test>";

    Xml_document doc;
    doc.parse(xmldef);

    ASSERT_TRUE(doc.get_root_element()->this_and_childs_have_both_value_and_children());

    doc.get_element("test/abc/auu")->set_value("");
    ASSERT_FALSE(doc.get_root_element()->this_and_childs_have_both_value_and_children());

}
