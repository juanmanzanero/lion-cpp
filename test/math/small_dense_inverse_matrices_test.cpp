#include "lion/math/small_dense_inverse_matrices.h"

#include "gtest/gtest.h"

#include "lion/io/Xml_document.h"


//
// Defines tests for the functions in header
// "src/core/math/small_dense_inverse_matrices.h".
//


static Xml_document small_dense_inverse_matrices_test_reference_file("data/small_dense_inverse_matrices.xml", true);


template<typename DetAndInverseMatrixFun>
void test_matrices(Document_element_ptr reference_xml_element, DetAndInverseMatrixFun &&det_and_inverse_matrix_fun)
{
    constexpr auto reltolnear = 1e-12;
    constexpr auto abstolnear = 1e-10;

    for (auto &test : reference_xml_element->get_children()) {
        const auto reference_matrix = test->get_child("matrix_colmaj")->get_value(std::vector<double>{});
        const auto reference_det = test->get_child("det")->get_value(double{});
        const auto reference_inverse_matrix = test->get_child("inverse_matrix_colmaj")->get_value(std::vector<double>{});

        const auto ret = det_and_inverse_matrix_fun(reference_matrix);
        const auto &det = std::get<0>(ret);
        const auto &inverse_matrix = std::get<1>(ret);

        EXPECT_NEAR(det, reference_det, reltolnear * std::max(1., std::abs(reference_det)));

        EXPECT_EQ(inverse_matrix.size(), reference_inverse_matrix.size());
        for (auto i = 0u; i < reference_inverse_matrix.size(); ++i) {
            EXPECT_NEAR(inverse_matrix[i], reference_inverse_matrix[i], abstolnear);
        }
    }
}


TEST(small_dense_inverse_matrices, inv2x2)
{
    test_matrices(small_dense_inverse_matrices_test_reference_file.get_root_element().get_child("inv2x2"),
        [](const auto &mat)
        {
            EXPECT_EQ(mat.size(), 4u);
            return std::make_pair(det2x2(mat), inv2x2(mat));
        });
}


TEST(small_dense_inverse_matrices, inv3x3)
{
    test_matrices(small_dense_inverse_matrices_test_reference_file.get_root_element().get_child("inv3x3"),
        [](const auto &mat)
        {
            EXPECT_EQ(mat.size(), 9u);
            return std::make_pair(det3x3(mat), inv3x3(mat));
        });
}


TEST(small_dense_inverse_matrices, inv4x4)
{
    test_matrices(small_dense_inverse_matrices_test_reference_file.get_root_element().get_child("inv4x4"),
        [](const auto &mat)
        {
            EXPECT_EQ(mat.size(), 16u);
            return std::make_pair(det4x4(mat), inv4x4(mat));
        });
}


TEST(small_dense_inverse_matrices, inv5x5)
{
    test_matrices(small_dense_inverse_matrices_test_reference_file.get_root_element().get_child("inv5x5"),
        [](const auto &mat)
        {
            EXPECT_EQ(mat.size(), 25u);
            return std::make_pair(det5x5(mat), inv5x5(mat));
        });
}


TEST(small_dense_inverse_matrices, inv6x6)
{
    test_matrices(small_dense_inverse_matrices_test_reference_file.get_root_element().get_child("inv6x6"),
        [](const auto &mat)
        {
            EXPECT_EQ(mat.size(), 36u);
            return std::make_pair(det6x6(mat), inv6x6(mat));
        });
}
