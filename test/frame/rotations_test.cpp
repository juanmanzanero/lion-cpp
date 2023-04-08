#include "gtest/gtest.h"
#include "lion/frame/rotations.h"

TEST(Rotations_test, omega_vector_to_matrix)
{
    Vector3d<scalar> omega_vector = {1.0, 2.0, 3.0};

    Matrix3x3<scalar> omega_matrix_expected = {0.0, -3.0, 2.0,
                                               3.0,  0.0,-1.0,
                                              -2.0,  1.0, 0.0};


    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).xx(), omega_matrix_expected.xx());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).xy(), omega_matrix_expected.xy());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).xz(), omega_matrix_expected.xz());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).yx(), omega_matrix_expected.yx());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).yy(), omega_matrix_expected.yy());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).yz(), omega_matrix_expected.yz());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).zx(), omega_matrix_expected.zx());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).zy(), omega_matrix_expected.zy());
    EXPECT_DOUBLE_EQ(omega_tensor_from_vector(omega_vector).zz(), omega_matrix_expected.zz());

    EXPECT_DOUBLE_EQ(omega_vector_from_matrix(omega_tensor_from_vector(omega_vector)).x(), omega_vector.x());
    EXPECT_DOUBLE_EQ(omega_vector_from_matrix(omega_tensor_from_vector(omega_vector)).y(), omega_vector.y());
    EXPECT_DOUBLE_EQ(omega_vector_from_matrix(omega_tensor_from_vector(omega_vector)).z(), omega_vector.z());
}
