#include "lion/frame/rotations.h"

#include <random>

#include "gtest/gtest.h"

#include "cppad/cppad.hpp"


//
// Defines tests for the functions in header
// "lion/frame/rotations.h".
//


TEST(rotations_test, omega_vector_to_matrix)
{
    Vector3d<double> omega_vector = {1.0, 2.0, 3.0};

    Matrix3x3<double> omega_matrix_expected = {0.0, -3.0, 2.0,
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


TEST(rotations_test, sis2rotmat_and_rotmat2sis)
{
    //
    // Think of this test as of rotating a body with a UNIQUE body-axes
    // frame, but on which we change the "labels" of the axes (so
    // that we rotate it using either Euler angles or "sis" angles).
    //

    constexpr auto tolnear = 10. * std::numeric_limits<double>::epsilon();

    constexpr auto T_ea2sis = Matrix3x3<double>{ 0., -1., 0.,
        1., 0., 0.,
        0., 0., 1. };

    for (auto yaw_rad : linspace(-pi, pi - tolnear, 100)) {
        for (auto pitch_rad : linspace(-0.5 * pi, 0.5 * pi - tolnear, 100)) {
            for (auto roll_rad : linspace(-pi, pi - tolnear, 100)) {

                const auto T_by_sis = T_ea2sis *
                    ea2rotmat(yaw_rad, pitch_rad, roll_rad) * T_ea2sis.t();

                const auto steer_rad = yaw_rad;
                const auto inclination_rad = -pitch_rad;
                const auto spin_rad = roll_rad;
                const auto T_by_sis_ = sis2rotmat(steer_rad, inclination_rad, spin_rad);

                for (auto i = 0u; i < 9u; ++i) {
                    EXPECT_NEAR(T_by_sis[i], T_by_sis_[i], tolnear);
                }

                const auto [steer_rad_, inclination_rad_, spin_rad_] =
                    rotmat2sis(T_by_sis_);
                EXPECT_NEAR(steer_rad_, steer_rad, tolnear);
                EXPECT_NEAR(inclination_rad_, inclination_rad, tolnear);
                EXPECT_NEAR(spin_rad_, spin_rad, tolnear);
            }
        }
    }
}


TEST(rotations_test, angular_kinematic_relationships)
{
    // tape the "ea2rotmat" function to get the derivative
    // of the rotation matrix w.r.t. the Euler angles
    using ADvector = std::vector<CppAD::AD<double> >;

    ADvector a_ea_rad(3u, 0.);
    CppAD::Independent(a_ea_rad);
    const auto a_T_Matrix3x3 = ea2rotmat(a_ea_rad[0], a_ea_rad[1], a_ea_rad[2]);
    ADvector a_T(a_T_Matrix3x3.cbegin(), a_T_Matrix3x3.cend());
    CppAD::ADFun<double> adfun;
    adfun.Dependent(a_ea_rad, a_T);
    adfun.optimize();


    // compare the formula "Omega_crossmatrix = transpose(T) * Tdot"
    // with the angular kinematic relationships
    constexpr auto tolnear = 1e-12;
    for (auto yaw_rad : linspace(-pi, pi, 20)) {
        for (auto pitch_rad : linspace(-0.5 * pi + 0.02, 0.5 * pi - 0.02, 20)) {
            for (auto roll_rad : linspace(-pi, pi, 20)) {

                const auto ea_rad = std::vector<double>{ yaw_rad, pitch_rad, roll_rad };
                const auto T_from_cppad = adfun.Forward(0, ea_rad);
                const auto T = ea2rotmat(yaw_rad, pitch_rad, roll_rad);
                for (auto i = 0u; i < 9u; ++i) {
                    EXPECT_DOUBLE_EQ(T[i], T_from_cppad[i]);
                }

                const auto dT_dypr = adfun.Jacobian(ea_rad);
                EXPECT_EQ(dT_dypr.size(), 9u * 3u);

                for (auto yawdot : linspace(-100., 100., 20)) {
                    for (auto pitchdot : linspace(-100., 100., 20)) {
                        for (auto rolldot : linspace(-100., 100., 20)) {

                            sMatrix3x3 Tdot;
                            for (auto i = 0u; i < 9u; ++i) {
                                Tdot[i] = dT_dypr[3 * i] * yawdot +
                                          dT_dypr[3 * i + 1] * pitchdot +
                                          dT_dypr[3 * i + 2] * rolldot;
                            }

                            const auto Omega_crossmatrix = transpose(T) * Tdot;
                            EXPECT_NEAR(Omega_crossmatrix[0], 0., tolnear);
                            EXPECT_NEAR(Omega_crossmatrix[4], 0., tolnear);
                            EXPECT_NEAR(Omega_crossmatrix[8], 0., tolnear);
                            EXPECT_NEAR(Omega_crossmatrix[1], -Omega_crossmatrix[3], tolnear);
                            EXPECT_NEAR(Omega_crossmatrix[2], -Omega_crossmatrix[6], tolnear);
                            EXPECT_NEAR(Omega_crossmatrix[5], -Omega_crossmatrix[7], tolnear);

                            const auto pqr = angular_kinematic_relationships(
                                yawdot, pitchdot, rolldot,
                                yaw_rad, pitch_rad, roll_rad);

                            EXPECT_NEAR(pqr[0], Omega_crossmatrix[5], tolnear);
                            EXPECT_NEAR(pqr[1], Omega_crossmatrix[6], tolnear);
                            EXPECT_NEAR(pqr[2], Omega_crossmatrix[1], tolnear);
                        }
                    }
                }
            }
        }
    }
}


TEST(rotations_test, inverse_angular_kinematic_relationships)
{
    constexpr auto tolnear = 1e-12;

    for (auto yaw_rad : linspace(-pi, pi, 20)) {
        for (auto pitch_rad : linspace(-0.5 * pi + 0.02, 0.5 * pi - 0.02, 20)) {
            for (auto roll_rad : linspace(-pi, pi, 20)) {

                const auto ea = Euler_angles{ yaw_rad, pitch_rad, roll_rad };

                for (auto thingdot0 : linspace(-100., 100., 20)) {
                    for (auto thingdot1 : linspace(-100., 100., 20)) {
                        for (auto thingdot2 : linspace(-100., 100., 20)) {

                            {
                                const auto pqr = angular_kinematic_relationships(
                                    thingdot0, thingdot1, thingdot2,
                                    ea.yaw(), ea.pitch(), ea.roll());

                                const auto yprdot = inverse_angular_kinematic_relationships(
                                    pqr[0], pqr[1], pqr[2],
                                    ea.yaw(), ea.pitch(), ea.roll());

                                EXPECT_NEAR(yprdot[0], thingdot0, tolnear);
                                EXPECT_NEAR(yprdot[1], thingdot1, tolnear);
                                EXPECT_NEAR(yprdot[2], thingdot2, tolnear);
                            }

                            {
                                const auto yprdot = inverse_angular_kinematic_relationships(
                                    thingdot0, thingdot1, thingdot2,
                                    ea.yaw(), ea.pitch(), ea.roll());

                                const auto pqr = angular_kinematic_relationships(
                                    yprdot[0], yprdot[1], yprdot[2],
                                    ea.yaw(), ea.pitch(), ea.roll());

                                EXPECT_NEAR(pqr[0], thingdot0, tolnear);
                                EXPECT_NEAR(pqr[1], thingdot1, tolnear);
                                EXPECT_NEAR(pqr[2], thingdot2, tolnear);
                            }
                        }
                    }
                }
            }
        }
    }
}


TEST(rotations_test, sis_kinematic_relationships_and_inverse)
{
    //
    // Think of this test as of rotating a body with a UNIQUE body-axes
    // frame, but on which we change the "labels" of the axes (so
    // that we rotate it using either Euler angles or "sis" angles).
    //

    constexpr auto T_ea2sis = Matrix3x3<double>{ 0., -1., 0.,
        1., 0., 0.,
        0., 0., 1. };

    constexpr auto tolnear = 1e-10;
    for (auto yaw_rad : linspace(-pi, pi, 20)) {
        for (auto pitch_rad : linspace(-0.5 * pi + 0.02, 0.5 * pi - 0.02, 20)) {
            for (auto roll_rad : linspace(-pi, pi, 20)) {

                const auto ea = Euler_angles{ yaw_rad, pitch_rad, roll_rad };

                for (auto thingdot0 : linspace(-100., 100., 20)) {
                    for (auto thingdot1 : linspace(-100., 100., 20)) {
                        for (auto thingdot2 : linspace(-100., 100., 20)) {

                            {
                                const auto pqr = angular_kinematic_relationships(
                                    thingdot0, thingdot1, thingdot2,
                                    ea.yaw(), ea.pitch(), ea.roll());

                                const auto pqr_sisaxes = T_ea2sis * sVector3d(pqr);

                                const auto steer_rad = yaw_rad;
                                const auto inclination_rad = -pitch_rad;
                                const auto spin_rad = roll_rad;

                                const auto steerdot = thingdot0;
                                const auto inclinationdot = -thingdot1;
                                const auto spindot = thingdot2;

                                const auto pqr_by_sis = sis_kinematic_relationships(
                                    steerdot, inclinationdot, spindot,
                                    steer_rad, inclination_rad, spin_rad);

                                EXPECT_NEAR(pqr_by_sis[0], pqr_sisaxes[0], tolnear);
                                EXPECT_NEAR(pqr_by_sis[1], pqr_sisaxes[1], tolnear);
                                EXPECT_NEAR(pqr_by_sis[2], pqr_sisaxes[2], tolnear);

                                const auto sisdot_alt = inverse_sis_kinematic_relationships(
                                    pqr_by_sis[0], pqr_by_sis[1], pqr_by_sis[2],
                                    steer_rad, inclination_rad, spin_rad);

                                EXPECT_NEAR(steerdot, sisdot_alt[0], tolnear);
                                EXPECT_NEAR(inclinationdot, sisdot_alt[1], tolnear);
                                EXPECT_NEAR(spindot, sisdot_alt[2], tolnear);
                            }

                            {
                                const auto yprdot = inverse_angular_kinematic_relationships(
                                    thingdot0, thingdot1, thingdot2,
                                    ea.yaw(), ea.pitch(), ea.roll());

                                const auto pqr_sisaxes = T_ea2sis * sVector3d{ thingdot0, thingdot1, thingdot2 };

                                const auto steer_rad = yaw_rad;
                                const auto inclination_rad = -pitch_rad;
                                const auto spin_rad = roll_rad;

                                const auto sisdot = inverse_sis_kinematic_relationships(
                                    pqr_sisaxes[0], pqr_sisaxes[1], pqr_sisaxes[2],
                                    steer_rad, inclination_rad, spin_rad);

                                EXPECT_NEAR(sisdot[0], yprdot[0], tolnear);
                                EXPECT_NEAR(sisdot[1], -yprdot[1], tolnear);
                                EXPECT_NEAR(sisdot[2], yprdot[2], tolnear);

                                const auto pqr_sisaxes_alt = sis_kinematic_relationships(
                                    sisdot[0], sisdot[1], sisdot[2],
                                    steer_rad, inclination_rad, spin_rad);

                                EXPECT_NEAR(pqr_sisaxes_alt[0], pqr_sisaxes[0], tolnear);
                                EXPECT_NEAR(pqr_sisaxes_alt[1], pqr_sisaxes[1], tolnear);
                                EXPECT_NEAR(pqr_sisaxes_alt[2], pqr_sisaxes[2], tolnear);
                            }
                        }
                    }
                }
            }
        }
    }
}


TEST(rotations_test, angular_kinematic_relationships_derivatives)
{
    // tape the "angular_kinematic_relationships" function to get
    // the derivative of the angular velocity w.r.t. the Euler angles
    // and their derivatives
    using ADvector = std::vector<CppAD::AD<double> >;

    ADvector a_eadot_and_ea_rad(6u, 0.);
    CppAD::Independent(a_eadot_and_ea_rad);
    const auto a_pqr_array = angular_kinematic_relationships(a_eadot_and_ea_rad[0],
                                                             a_eadot_and_ea_rad[1],
                                                             a_eadot_and_ea_rad[2],
                                                             a_eadot_and_ea_rad[3],
                                                             a_eadot_and_ea_rad[4],
                                                             a_eadot_and_ea_rad[5]);
    ADvector a_pqr(a_pqr_array.cbegin(), a_pqr_array.cend());
    CppAD::ADFun<double> adfun;
    adfun.Dependent(a_eadot_and_ea_rad, a_pqr);
    adfun.optimize();


    // compare the adfun's results with the angular kinematic relationships
    // derivatives
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distr(0., 1.);
    const auto rand_in_limits = [&](double lo, double hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto tolnear = 1e-12;
    constexpr auto num_tests = 5e6;
    for (auto t = 0u; t < num_tests; ++t) {

        const auto yaw_rad = rand_in_limits(-pi, pi);
        const auto pitch_rad = rand_in_limits(-0.5 * pi + 0.02, 0.5 * pi - 0.02);
        const auto roll_rad = rand_in_limits(-pi, pi);

        const auto yawdot = rand_in_limits(-5., 5.);
        const auto pitchdot = rand_in_limits(-5., 5.);
        const auto rolldot = rand_in_limits(-5., 5.);

        const auto pqr_cppad = adfun.Forward(0,
            std::vector<double>{ yawdot, pitchdot, rolldot, yaw_rad, pitch_rad, roll_rad });

        const auto pqr = angular_kinematic_relationships(
            yawdot, pitchdot, rolldot,
            yaw_rad, pitch_rad, roll_rad);

        EXPECT_NEAR(pqr[0], pqr_cppad[0], tolnear);
        EXPECT_NEAR(pqr[1], pqr_cppad[1], tolnear);
        EXPECT_NEAR(pqr[2], pqr_cppad[2], tolnear);

        const auto dpqr_dypryprdot = adfun.Jacobian(
            std::vector<double>{ yawdot, pitchdot, rolldot, yaw_rad, pitch_rad, roll_rad });

        const auto yawdotdot = rand_in_limits(-20., 20.);
        const auto pitchdotdot = rand_in_limits(-20., 20.);
        const auto rolldotdot = rand_in_limits(-20., 20.);

        sVector3d pqrdot_cppad;
        for (auto i = 0u; i < 3u; ++i) {
            pqrdot_cppad[i] = dpqr_dypryprdot[6 * i] * yawdotdot +
                dpqr_dypryprdot[6 * i + 1] * pitchdotdot +
                dpqr_dypryprdot[6 * i + 2] * rolldotdot +
                dpqr_dypryprdot[6 * i + 3] * yawdot +
                dpqr_dypryprdot[6 * i + 4] * pitchdot +
                dpqr_dypryprdot[6 * i + 5] * rolldot;
        }

        const auto pqrdot = angular_kinematic_relationships_derivatives(
            yawdotdot, pitchdotdot, rolldotdot,
            yawdot, pitchdot, rolldot,
            yaw_rad, pitch_rad, roll_rad);

        EXPECT_NEAR(pqrdot[0], pqrdot_cppad[0], tolnear);
        EXPECT_NEAR(pqrdot[1], pqrdot_cppad[1], tolnear);
        EXPECT_NEAR(pqrdot[2], pqrdot_cppad[2], tolnear);
    }
}


TEST(rotations_test, inverse_angular_kinematic_relationships_derivatives)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distr(0., 1.);
    const auto rand_in_limits = [&](double lo, double hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto tolnear = 1e-9;
    constexpr auto num_tests = 5e6;
    for (auto t = 0u; t < num_tests; ++t) {
        const auto yaw_rad = rand_in_limits(-pi, pi);
        const auto pitch_rad = rand_in_limits(-0.5 * pi + 0.02, 0.5 * pi - 0.02);
        const auto roll_rad = rand_in_limits(-pi, pi);

        const auto yawdot = rand_in_limits(-100., 100.);
        const auto pitchdot = rand_in_limits(-100., 100.);
        const auto rolldot = rand_in_limits(-100., 100.);

        const auto thingdotdot0 = rand_in_limits(-100., 100.);
        const auto thingdotdot1 = rand_in_limits(-100., 100.);
        const auto thingdotdot2 = rand_in_limits(-100., 100.);

        {
            const auto pqrdot = angular_kinematic_relationships_derivatives(thingdotdot0, thingdotdot1, thingdotdot2,
                                                                            yawdot, pitchdot, rolldot,
                                                                            yaw_rad, pitch_rad, roll_rad);

            const auto yprdotdot = inverse_angular_kinematic_relationships_derivatives(pqrdot[0], pqrdot[1], pqrdot[2],
                                                                                       yawdot, pitchdot, rolldot,
                                                                                       yaw_rad, pitch_rad, roll_rad);

            EXPECT_NEAR(yprdotdot[0], thingdotdot0, tolnear);
            EXPECT_NEAR(yprdotdot[1], thingdotdot1, tolnear);
            EXPECT_NEAR(yprdotdot[2], thingdotdot2, tolnear);
        }

        {
            const auto yprdotdot = inverse_angular_kinematic_relationships_derivatives(thingdotdot0, thingdotdot1, thingdotdot2,
                                                                                       yawdot, pitchdot, rolldot,
                                                                                       yaw_rad, pitch_rad, roll_rad);

            const auto pqrdot = angular_kinematic_relationships_derivatives(yprdotdot[0], yprdotdot[1], yprdotdot[2],
                                                                            yawdot, pitchdot, rolldot,
                                                                            yaw_rad, pitch_rad, roll_rad);

            EXPECT_NEAR(pqrdot[0], thingdotdot0, tolnear);
            EXPECT_NEAR(pqrdot[1], thingdotdot1, tolnear);
            EXPECT_NEAR(pqrdot[2], thingdotdot2, tolnear);
        }
    }
}


TEST(rotations_test, sis_kinematic_relationships_derivatives_and_inverse)
{
    //
    // Think of this test as of rotating a body with a UNIQUE body-axes
    // frame, but on which we change the "labels" of the axes (so
    // that we rotate it using either Euler angles or "sis" angles).
    //

    constexpr auto T_ea2sis = Matrix3x3<double>{ 0., -1., 0.,
        1., 0., 0.,
        0., 0., 1. };

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distr(0., 1.);
    const auto rand_in_limits = [&](double lo, double hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto tolnear = 1e-10;
    constexpr auto num_tests = 5e6;
    for (auto t = 0u; t < num_tests; ++t) {

        const auto yaw_rad = rand_in_limits(-pi, pi);
        const auto pitch_rad = rand_in_limits(-0.5 * pi + 0.02, 0.5 * pi - 0.02);
        const auto roll_rad = rand_in_limits(-pi, pi);

        const auto yawdot = rand_in_limits(-5., 5.);
        const auto pitchdot = rand_in_limits(-5., 5.);
        const auto rolldot = rand_in_limits(-5., 5.);

        const auto thingdotdot0 = rand_in_limits(-20., 20.);
        const auto thingdotdot1 = rand_in_limits(-20., 20.);
        const auto thingdotdot2 = rand_in_limits(-20., 20.);

        {
            const auto pqrdot = angular_kinematic_relationships_derivatives(thingdotdot0, thingdotdot1, thingdotdot2,
                                                                            yawdot, pitchdot, rolldot,
                                                                            yaw_rad, pitch_rad, roll_rad);

            const auto pqrdot_sisaxes = T_ea2sis * sVector3d(pqrdot);

            const auto steer_rad = yaw_rad;
            const auto inclination_rad = -pitch_rad;
            const auto spin_rad = roll_rad;

            const auto steerdot = yawdot;
            const auto inclinationdot = -pitchdot;
            const auto spindot = rolldot;

            const auto steerdotdot = thingdotdot0;
            const auto inclinationdotdot = -thingdotdot1;
            const auto spindotdot = thingdotdot2;

            const auto pqrdot_by_sis = sis_kinematic_relationships_derivatives(
                steerdotdot, inclinationdotdot, spindotdot,
                steerdot, inclinationdot, spindot,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(pqrdot_by_sis[0], pqrdot_sisaxes[0], tolnear);
            EXPECT_NEAR(pqrdot_by_sis[1], pqrdot_sisaxes[1], tolnear);
            EXPECT_NEAR(pqrdot_by_sis[2], pqrdot_sisaxes[2], tolnear);

            const auto sisdotdot_alt = inverse_sis_kinematic_relationships_derivatives(
                pqrdot_by_sis[0], pqrdot_by_sis[1], pqrdot_by_sis[2],
                steerdot, inclinationdot, spindot,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(sisdotdot_alt[0], steerdotdot, tolnear);
            EXPECT_NEAR(sisdotdot_alt[1], inclinationdotdot, tolnear);
            EXPECT_NEAR(sisdotdot_alt[2], spindotdot, tolnear);
        }

        {
            const auto yprdotdot = inverse_angular_kinematic_relationships_derivatives(
                thingdotdot0, thingdotdot1, thingdotdot2,
                yawdot, pitchdot, rolldot,
                yaw_rad, pitch_rad, roll_rad);

            const auto pqrdot_sisaxes = T_ea2sis * sVector3d{ thingdotdot0, thingdotdot1, thingdotdot2 };

            const auto steer_rad = yaw_rad;
            const auto inclination_rad = -pitch_rad;
            const auto spin_rad = roll_rad;

            const auto steerdot = yawdot;
            const auto inclinationdot = -pitchdot;
            const auto spindot = rolldot;

            const auto sisdotdot = inverse_sis_kinematic_relationships_derivatives(
                pqrdot_sisaxes[0], pqrdot_sisaxes[1], pqrdot_sisaxes[2],
                steerdot, inclinationdot, spindot,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(sisdotdot[0], yprdotdot[0], tolnear);
            EXPECT_NEAR(sisdotdot[1], -yprdotdot[1], tolnear);
            EXPECT_NEAR(sisdotdot[2], yprdotdot[2], tolnear);

            const auto pqrdot_sisaxes_alt = sis_kinematic_relationships_derivatives(
                sisdotdot[0], sisdotdot[1], sisdotdot[2],
                steerdot, inclinationdot, spindot,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(pqrdot_sisaxes_alt[0], pqrdot_sisaxes[0], tolnear);
            EXPECT_NEAR(pqrdot_sisaxes_alt[1], pqrdot_sisaxes[1], tolnear);
            EXPECT_NEAR(pqrdot_sisaxes_alt[2], pqrdot_sisaxes[2], tolnear);
        }
    }
}


TEST(rotations_test, angular_kinematic_relationships_crossderivatives)
{
    // tape the "angular_kinematic_relationships" function to get
    // the derivative of the angular velocity w.r.t. the Euler angles
    // and their derivatives
    using ADvector = std::vector<CppAD::AD<double> >;

    ADvector a_eadot_and_ea_rad(6u, 0.);
    CppAD::Independent(a_eadot_and_ea_rad);
    const auto a_pqr_array = angular_kinematic_relationships(a_eadot_and_ea_rad[0],
                                                             a_eadot_and_ea_rad[1],
                                                             a_eadot_and_ea_rad[2],
                                                             a_eadot_and_ea_rad[3],
                                                             a_eadot_and_ea_rad[4],
                                                             a_eadot_and_ea_rad[5]);
    ADvector a_pqr(a_pqr_array.cbegin(), a_pqr_array.cend());
    CppAD::ADFun<double> adfun;
    adfun.Dependent(a_eadot_and_ea_rad, a_pqr);
    adfun.optimize();


    // compare the adfun's results with the angular kinematic relationships
    // derivatives
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distr(0., 1.);
    const auto rand_in_limits = [&](double lo, double hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto tolnear = 1e-12;
    constexpr auto num_tests = 5e6;
    for (auto t = 0u; t < num_tests; ++t) {

        const auto yaw_rad = rand_in_limits(-pi, pi);
        const auto pitch_rad = rand_in_limits(-0.5 * pi + 0.02, 0.5 * pi - 0.02);
        const auto roll_rad = rand_in_limits(-pi, pi);

        const auto yawdot0 = rand_in_limits(-5., 5.);
        const auto pitchdot0 = rand_in_limits(-5., 5.);
        const auto rolldot0 = rand_in_limits(-5., 5.);

        const auto yawdot1 = rand_in_limits(-5., 5.);
        const auto pitchdot1 = rand_in_limits(-5., 5.);
        const auto rolldot1 = rand_in_limits(-5., 5.);

        const auto Dpqr01 = adfun.Jacobian(
            std::vector<double>{ yawdot0, pitchdot0, rolldot0, yaw_rad, pitch_rad, roll_rad });

        const auto yawdot0dot1 = rand_in_limits(-20., 20.);
        const auto pitchdot0dot1 = rand_in_limits(-20., 20.);
        const auto rolldot0dot1 = rand_in_limits(-20., 20.);

        sVector3d pqrdot01_cppad;
        for (auto i = 0u; i < 3u; ++i) {
            pqrdot01_cppad[i] = Dpqr01[6 * i] * yawdot0dot1 +
                Dpqr01[6 * i + 1] * pitchdot0dot1 +
                Dpqr01[6 * i + 2] * rolldot0dot1 +
                Dpqr01[6 * i + 3] * yawdot1 +
                Dpqr01[6 * i + 4] * pitchdot1 +
                Dpqr01[6 * i + 5] * rolldot1;
        }

        const auto pqrdot01 = angular_kinematic_relationships_crossderivatives(
            yawdot0dot1, pitchdot0dot1, rolldot0dot1,
            yawdot0, pitchdot0, rolldot0,
            yawdot1, pitchdot1, rolldot1,
            yaw_rad, pitch_rad, roll_rad);

        EXPECT_NEAR(pqrdot01[0], pqrdot01_cppad[0], tolnear);
        EXPECT_NEAR(pqrdot01[1], pqrdot01_cppad[1], tolnear);
        EXPECT_NEAR(pqrdot01[2], pqrdot01_cppad[2], tolnear);
    }
}


TEST(rotations_test, inverse_angular_kinematic_relationships_crossderivatives)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distr(0., 1.);
    const auto rand_in_limits = [&](double lo, double hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto tolnear = 1e-9;
    constexpr auto num_tests = 5e6;
    for (auto t = 0u; t < num_tests; ++t) {
        const auto yaw_rad = rand_in_limits(-pi, pi);
        const auto pitch_rad = rand_in_limits(-0.5 * pi + 0.02, 0.5 * pi - 0.02);
        const auto roll_rad = rand_in_limits(-pi, pi);

        const auto yawdot0 = rand_in_limits(-100., 100.);
        const auto pitchdot0 = rand_in_limits(-100., 100.);
        const auto rolldot0 = rand_in_limits(-100., 100.);

        const auto yawdot1 = rand_in_limits(-100., 100.);
        const auto pitchdot1 = rand_in_limits(-100., 100.);
        const auto rolldot1 = rand_in_limits(-100., 100.);

        const auto thingdotdotA = rand_in_limits(-100., 100.);
        const auto thingdotdotB = rand_in_limits(-100., 100.);
        const auto thingdotdotC = rand_in_limits(-100., 100.);

        {
            const auto pqrdot01 = angular_kinematic_relationships_crossderivatives(
                thingdotdotA, thingdotdotB, thingdotdotC,
                yawdot0, pitchdot0, rolldot0,
                yawdot1, pitchdot1, rolldot1,
                yaw_rad, pitch_rad, roll_rad);

            const auto yprdot0dot1 = inverse_angular_kinematic_relationships_crossderivatives(
                pqrdot01[0], pqrdot01[1], pqrdot01[2],
                yawdot0, pitchdot0, rolldot0,
                yawdot1, pitchdot1, rolldot1,
                yaw_rad, pitch_rad, roll_rad);

            EXPECT_NEAR(yprdot0dot1[0], thingdotdotA, tolnear);
            EXPECT_NEAR(yprdot0dot1[1], thingdotdotB, tolnear);
            EXPECT_NEAR(yprdot0dot1[2], thingdotdotC, tolnear);
        }

        {
            const auto yprdot0dot1 = inverse_angular_kinematic_relationships_crossderivatives(
                thingdotdotA, thingdotdotB, thingdotdotC,
                yawdot0, pitchdot0, rolldot0,
                yawdot1, pitchdot1, rolldot1,
                yaw_rad, pitch_rad, roll_rad);

            const auto pqrdot01 = angular_kinematic_relationships_crossderivatives(
                yprdot0dot1[0], yprdot0dot1[1], yprdot0dot1[2],
                yawdot0, pitchdot0, rolldot0,
                yawdot1, pitchdot1, rolldot1,
                yaw_rad, pitch_rad, roll_rad);

            EXPECT_NEAR(pqrdot01[0], thingdotdotA, tolnear);
            EXPECT_NEAR(pqrdot01[1], thingdotdotB, tolnear);
            EXPECT_NEAR(pqrdot01[2], thingdotdotC, tolnear);
        }
    }
}


TEST(rotations_test, sis_kinematic_relationships_crossderivatives_and_inverse)
{
    //
    // Think of this test as of rotating a body with a UNIQUE body-axes
    // frame, but on which we change the "labels" of the axes (so
    // that we rotate it using either Euler angles or "sis" angles).
    //

    constexpr auto T_ea2sis = Matrix3x3<double>{ 0., -1., 0.,
        1., 0., 0.,
        0., 0., 1. };

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distr(0., 1.);
    const auto rand_in_limits = [&](double lo, double hi)
    {
        return (hi - lo) * distr(generator) + lo;
    };

    constexpr auto tolnear = 1e-10;
    constexpr auto num_tests = 5e6;
    for (auto t = 0u; t < num_tests; ++t) {

        const auto yaw_rad = rand_in_limits(-pi, pi);
        const auto pitch_rad = rand_in_limits(-0.5 * pi + 0.02, 0.5 * pi - 0.02);
        const auto roll_rad = rand_in_limits(-pi, pi);

        const auto yawdot0 = rand_in_limits(-5., 5.);
        const auto pitchdot0 = rand_in_limits(-5., 5.);
        const auto rolldot0 = rand_in_limits(-5., 5.);

        const auto yawdot1 = rand_in_limits(-5., 5.);
        const auto pitchdot1 = rand_in_limits(-5., 5.);
        const auto rolldot1 = rand_in_limits(-5., 5.);

        const auto thingdotdotA = rand_in_limits(-20., 20.);
        const auto thingdotdotB = rand_in_limits(-20., 20.);
        const auto thingdotdotC = rand_in_limits(-20., 20.);

        {
            const auto pqrdot01 = angular_kinematic_relationships_crossderivatives(
                thingdotdotA, thingdotdotB, thingdotdotC,
                yawdot0, pitchdot0, rolldot0,
                yawdot1, pitchdot1, rolldot1,
                yaw_rad, pitch_rad, roll_rad);

            const auto pqrdot01_sisaxes = T_ea2sis * sVector3d(pqrdot01);

            const auto steer_rad = yaw_rad;
            const auto inclination_rad = -pitch_rad;
            const auto spin_rad = roll_rad;

            const auto steerdot0 = yawdot0;
            const auto inclinationdot0 = -pitchdot0;
            const auto spindot0 = rolldot0;

            const auto steerdot1 = yawdot1;
            const auto inclinationdot1 = -pitchdot1;
            const auto spindot1 = rolldot1;

            const auto steerdot0dot1 = thingdotdotA;
            const auto inclinationdot0dot1 = -thingdotdotB;
            const auto spindot0dot1 = thingdotdotC;

            const auto pqrdot01_by_sis = sis_kinematic_relationships_crossderivatives(
                steerdot0dot1, inclinationdot0dot1, spindot0dot1,
                steerdot0, inclinationdot0, spindot0,
                steerdot1, inclinationdot1, spindot1,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(pqrdot01_by_sis[0], pqrdot01_sisaxes[0], tolnear);
            EXPECT_NEAR(pqrdot01_by_sis[1], pqrdot01_sisaxes[1], tolnear);
            EXPECT_NEAR(pqrdot01_by_sis[2], pqrdot01_sisaxes[2], tolnear);

            const auto sisdot0dot1_alt = inverse_sis_kinematic_relationships_crossderivatives(
                pqrdot01_by_sis[0], pqrdot01_by_sis[1], pqrdot01_by_sis[2],
                steerdot0, inclinationdot0, spindot0,
                steerdot1, inclinationdot1, spindot1,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(sisdot0dot1_alt[0], steerdot0dot1, tolnear);
            EXPECT_NEAR(sisdot0dot1_alt[1], inclinationdot0dot1, tolnear);
            EXPECT_NEAR(sisdot0dot1_alt[2], spindot0dot1, tolnear);
        }

        {
            const auto yprdot0dot1 = inverse_angular_kinematic_relationships_crossderivatives(
                thingdotdotA, thingdotdotB, thingdotdotC,
                yawdot0, pitchdot0, rolldot0,
                yawdot1, pitchdot1, rolldot1,
                yaw_rad, pitch_rad, roll_rad);

            const auto pqrdot01_sisaxes = T_ea2sis * sVector3d{ thingdotdotA, thingdotdotB, thingdotdotC };

            const auto steer_rad = yaw_rad;
            const auto inclination_rad = -pitch_rad;
            const auto spin_rad = roll_rad;

            const auto steerdot0 = yawdot0;
            const auto inclinationdot0 = -pitchdot0;
            const auto spindot0 = rolldot0;

            const auto steerdot1 = yawdot1;
            const auto inclinationdot1 = -pitchdot1;
            const auto spindot1 = rolldot1;

            const auto sisdot0dot1 = inverse_sis_kinematic_relationships_crossderivatives(
                pqrdot01_sisaxes[0], pqrdot01_sisaxes[1], pqrdot01_sisaxes[2],
                steerdot0, inclinationdot0, spindot0,
                steerdot1, inclinationdot1, spindot1,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(sisdot0dot1[0], yprdot0dot1[0], tolnear);
            EXPECT_NEAR(sisdot0dot1[1], -yprdot0dot1[1], tolnear);
            EXPECT_NEAR(sisdot0dot1[2], yprdot0dot1[2], tolnear);

            const auto pqrdot01_sisaxes_alt = sis_kinematic_relationships_crossderivatives(
                sisdot0dot1[0], sisdot0dot1[1], sisdot0dot1[2],
                steerdot0, inclinationdot0, spindot0,
                steerdot1, inclinationdot1, spindot1,
                steer_rad, inclination_rad, spin_rad);

            EXPECT_NEAR(pqrdot01_sisaxes_alt[0], pqrdot01_sisaxes[0], tolnear);
            EXPECT_NEAR(pqrdot01_sisaxes_alt[1], pqrdot01_sisaxes[1], tolnear);
            EXPECT_NEAR(pqrdot01_sisaxes_alt[2], pqrdot01_sisaxes[2], tolnear);
        }
    }
}
