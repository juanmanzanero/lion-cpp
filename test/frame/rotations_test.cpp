#include "gtest/gtest.h"

#include "lion/frame/rotations.h"


//
// Defines tests for the functions in header
// "lion/frame/rotations.h".
//


TEST(rotations_test, omega_vector_to_matrix)
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


TEST(rotations_test, ea2rotmat_and_rotmat2ea)
{
    constexpr auto tolnear = 10. * std::numeric_limits<double>::epsilon();

    // we'll take values in [-pi, pi) to isolate a bit the test
    // from the "wrap_to_pi" function
    for (auto yaw_rad : linspace(-pi, pi - tolnear, 100)) {
        for (auto pitch_rad : linspace(-0.5 * pi, 0.5 * pi, 100)) {
            for (auto roll_rad : linspace(-pi, pi - tolnear, 100)) {

                const auto ea = Euler_angles{ yaw_rad, pitch_rad, roll_rad };
                const auto rotmat = ea2rotmat(ea.yaw(), ea.pitch(), ea.roll());
                const auto identity_ = rotmat * rotmat.t();
                EXPECT_NEAR(identity_[0], 1., tolnear);
                EXPECT_NEAR(identity_[1], 0., tolnear);
                EXPECT_NEAR(identity_[2], 0., tolnear);
                EXPECT_NEAR(identity_[3], 0., tolnear);
                EXPECT_NEAR(identity_[4], 1., tolnear);
                EXPECT_NEAR(identity_[5], 0., tolnear);
                EXPECT_NEAR(identity_[6], 0., tolnear);
                EXPECT_NEAR(identity_[7], 0., tolnear);
                EXPECT_NEAR(identity_[8], 1., tolnear);

                const auto ea_ = Euler_angles{ rotmat2ea(rotmat) };
                EXPECT_NEAR(ea_.yaw(), ea.yaw(), tolnear);
                EXPECT_NEAR(ea_.pitch(), ea.pitch(), tolnear);
                EXPECT_NEAR(ea_.roll(), ea.roll(), tolnear);
            }
        }
    }
}


TEST(rotations_test, scs2rotmat)
{
    constexpr auto tolnear = 10. * std::numeric_limits<double>::epsilon();

    constexpr auto T_ea2scs = Matrix3x3<double>{ 0., -1., 0.,
        1., 0., 0.,
        0., 0., 1. };

    for (auto yaw_rad : linspace(-pi, pi, 100)) {
        for (auto pitch_rad : linspace(-0.5 * pi, 0.5 * pi, 100)) {
            for (auto roll_rad : linspace(-pi, pi, 100)) {

                const auto T_by_scs = T_ea2scs *
                    ea2rotmat(yaw_rad, pitch_rad, roll_rad) * T_ea2scs.t();

                const auto steer_rad = yaw_rad;
                const auto camber_rad = -pitch_rad;
                const auto spin_rad = roll_rad;
                const auto T_by_scs_ = scs2rotmat(steer_rad, camber_rad, spin_rad);

                for (auto i = 0u; i < 9u; ++i) {
                    EXPECT_NEAR(T_by_scs[i], T_by_scs_[i], tolnear);
                }
            }
        }
    }
}


TEST(rotations_test, angular_kinematic_relationships_and_inverse)
{
    constexpr auto tolnear = 1e-12;

    for (auto yaw_rad : linspace(-pi, pi, 20)) {
        for (auto pitch_rad : linspace(-0.5 * pi + 0.02, 0.5 * pi - 0.02, 20)) {
            for (auto roll_rad : linspace(-pi, pi, 20)) {

                const auto ea = Euler_angles{ yaw_rad, pitch_rad, roll_rad };

                for (auto thing0 : linspace(-100., 100., 20)) {
                    for (auto thing1 : linspace(-100., 100., 20)) {
                        for (auto thing2 : linspace(-100., 100., 20)) {

                            {
                                const auto pqr = angular_kinematic_relationships(
                                    thing0, thing1, thing2,
                                    ea.yaw(), ea.pitch(), ea.roll());

                                const auto yprdot = inverse_angular_kinematic_relationships(
                                    pqr[0], pqr[1], pqr[2],
                                    ea.yaw(), ea.pitch(), ea.roll());

                                EXPECT_NEAR(yprdot[0], thing0, tolnear);
                                EXPECT_NEAR(yprdot[1], thing1, tolnear);
                                EXPECT_NEAR(yprdot[2], thing2, tolnear);
                            }

                            {
                                const auto yprdot = inverse_angular_kinematic_relationships(
                                    thing0, thing1, thing2,
                                    ea.yaw(), ea.pitch(), ea.roll());

                                const auto pqr = angular_kinematic_relationships(
                                    yprdot[0], yprdot[1], yprdot[2],
                                    ea.yaw(), ea.pitch(), ea.roll());

                                EXPECT_NEAR(pqr[0], thing0, tolnear);
                                EXPECT_NEAR(pqr[1], thing1, tolnear);
                                EXPECT_NEAR(pqr[2], thing2, tolnear);
                            }
                        }
                    }
                }
            }
        }
    }
}


