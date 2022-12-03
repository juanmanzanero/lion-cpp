#include "gtest/gtest.h"
#include "lion/frame/frame.h"
#include "lion/foundation/constants.h"
#include "lion/foundation/types.h"
#include <cmath>
#include "lion/thirdparty/include/logger.hpp"

  
void generic_inertial_frame_tests(Frame<scalar>& inertial_frame);

class Rotated_frames_test : public ::testing::Test
{
 protected:
    Rotated_frames_test(): inertial_frame(),
                           frame_rot_x(sVector3d(0.0), sVector3d(0.0),{theta_x},{w_x},{X},inertial_frame),
                           frame_rot_y(sVector3d(0.0), sVector3d(0.0),{theta_y},{w_y},{Y},inertial_frame),
                           frame_rot_z(sVector3d(0.0), sVector3d(0.0),{theta_z},{w_z},{Z},inertial_frame),
                           frame_rot_xy(sVector3d(0.0), sVector3d(0.0),{theta_x,theta_y},{w_x,w_y},{X,Y},inertial_frame),
                           frame_rot_yz(sVector3d(0.0), sVector3d(0.0),{theta_y,theta_z},{w_y,w_z},{Y,Z},inertial_frame),
                           frame_rot_zx(sVector3d(0.0), sVector3d(0.0),{theta_z,theta_x},{w_z,w_x},{Z,X},inertial_frame),
                           frame_rot_xyz(sVector3d(0.0), sVector3d(0.0),{theta_x,theta_y,theta_z},{w_x,w_y,w_z},{X,Y,Z},inertial_frame),
                           euler_frame(sVector3d(0.0), sVector3d(0.0), {phi,theta,psi},{dphi,dtheta,dpsi},{Z,X,Z},inertial_frame)
                  {} ;

    // Data used
    const double theta_x = 45.0*DEG;
    const double theta_y = 135.0*DEG;
    const double theta_z = -90.0*DEG;
    const double w_x = 8.3;
    const double w_y = -5.6;
    const double w_z = 0.12;

    const scalar phi = 10.0*DEG;
    const scalar theta = 50.0*DEG;
    const scalar psi = -300.0*DEG;

    const scalar dphi = pi;
    const scalar dtheta = pi/3.0;
    const scalar dpsi = 2.0*pi;


    Frame<scalar> inertial_frame;
    const Frame<scalar,decltype(inertial_frame),1> frame_rot_x;
    const Frame<scalar,decltype(inertial_frame),1> frame_rot_y;
    const Frame<scalar,decltype(inertial_frame),1> frame_rot_z;
    const Frame<scalar,decltype(inertial_frame),2> frame_rot_xy;
    const Frame<scalar,decltype(inertial_frame),2> frame_rot_yz;
    const Frame<scalar,decltype(inertial_frame),2> frame_rot_zx;
    const Frame<scalar,decltype(inertial_frame),3> frame_rot_xyz;
    const Frame<scalar,decltype(inertial_frame),3> euler_frame;
};

class Rotation_wrt_target_test : public :: testing::Test
/*
*       This test defines a frame tree such that:
*                       1           11
*                       |
*                       2
*                     /   \
*                    3     4
*                    |     |      
*                    5     6
*                   / \   / \
*                  7   8 9  10
*/

{
 protected:
    Rotation_wrt_target_test() {}


    // Data used
    const std::array<scalar,9> theta = {1.0, -0.4, 2.0, 1.2, -0.333, 0.89, 0.5, -0.6667, 0.707};
    const std::array<scalar,9> dtheta = {0.1, -0.25, pi/12.0, -pi/12.0, pi, -pi/4.0, pi/8.0, -pi/8.0, 0.668};
    const std::array<Axis,9> axis = {X,Y,Z,Z,Y,X,X,Z,X};

    const Frame<scalar> frame1;
    const Frame<scalar, decltype(frame1), 1> frame2 = { sVector3d(0.0), sVector3d(0.0), {theta[0]}, {dtheta[0]}, {axis[0]}, frame1};
    const Frame<scalar, decltype(frame2), 1> frame3 = { sVector3d(0.0), sVector3d(0.0), {theta[1]}, {dtheta[1]}, {axis[1]}, frame2};
    const Frame<scalar, decltype(frame2), 1> frame4 = { sVector3d(0.0), sVector3d(0.0), {theta[2]}, {dtheta[2]}, {axis[2]}, frame2};
    const Frame<scalar, decltype(frame3), 1> frame5 = { sVector3d(0.0), sVector3d(0.0), {theta[3]}, {dtheta[3]}, {axis[3]}, frame3};
    const Frame<scalar, decltype(frame4), 1> frame6 = { sVector3d(0.0), sVector3d(0.0), {theta[4]}, {dtheta[4]}, {axis[4]}, frame4};
    const Frame<scalar, decltype(frame5), 1> frame7 = { sVector3d(0.0), sVector3d(0.0), {theta[5]}, {dtheta[5]}, {axis[5]}, frame5};
    const Frame<scalar, decltype(frame5), 1> frame8 = { sVector3d(0.0), sVector3d(0.0), {theta[6]}, {dtheta[6]}, {axis[6]}, frame5};
    const Frame<scalar, decltype(frame6), 1> frame9 = { sVector3d(0.0), sVector3d(0.0), {theta[7]}, {dtheta[7]}, {axis[7]}, frame6};
    const Frame<scalar, decltype(frame6), 1> frame10 = { sVector3d(0.0), sVector3d(0.0), {theta[8]}, {dtheta[8]}, {axis[8]}, frame6};

    const Frame<scalar> frame11; // A disconnected frame

    // Rotation matrices
    const sMatrix3x3 Q12 = rotation_matrix_x(theta[0]);
    const sMatrix3x3 Q23 = rotation_matrix_y(theta[1]);
    const sMatrix3x3 Q35 = rotation_matrix_z(theta[3]);
    const sMatrix3x3 Q57 = rotation_matrix_x(theta[5]);
    const sMatrix3x3 Q58 = rotation_matrix_x(theta[6]);
    const sMatrix3x3 Q24 = rotation_matrix_z(theta[2]);
    const sMatrix3x3 Q46 = rotation_matrix_y(theta[4]);
    const sMatrix3x3 Q69 = rotation_matrix_z(theta[7]);
    const sMatrix3x3 Q610 = rotation_matrix_x(theta[8]);

    // Rotation matrices derivatives
    const sMatrix3x3 dQ12 = dtheta[0] * drotation_matrix_x(theta[0]);
    const sMatrix3x3 dQ23 = dtheta[1] * drotation_matrix_y(theta[1]);
    const sMatrix3x3 dQ35 = dtheta[3] * drotation_matrix_z(theta[3]);
    const sMatrix3x3 dQ57 = dtheta[5] * drotation_matrix_x(theta[5]);
    const sMatrix3x3 dQ58 = dtheta[6] * drotation_matrix_x(theta[6]);
    const sMatrix3x3 dQ24 = dtheta[2] * drotation_matrix_z(theta[2]);
    const sMatrix3x3 dQ46 = dtheta[4] * drotation_matrix_y(theta[4]);
    const sMatrix3x3 dQ69 = dtheta[7] * drotation_matrix_z(theta[7]);
    const sMatrix3x3 dQ610 = dtheta[8] * drotation_matrix_x(theta[8]);

    // Full path rotation matrices
    const sMatrix3x3 Q47  = transpose(Q24)*Q23*Q35*Q57 ;
    const sMatrix3x3 Q29  = Q24*Q46*Q69;
    const sMatrix3x3 Q19  = Q12*Q29;
    const sMatrix3x3 Q28  = Q23*Q35*Q58;
    const sMatrix3x3 Q210 = Q24*Q46*Q610;
    const sMatrix3x3 Q810 = transpose(Q28)*Q210;
    const sMatrix3x3 Q78  = transpose(Q57)*Q58 ;

    // Full path rotation matrices derivatives
    const sMatrix3x3 dQ47 = transpose(dQ24)*Q23*Q35*Q57 + transpose(Q24)*(dQ23*Q35*Q57+Q23*dQ35*Q57+Q23*Q35*dQ57);

    const sMatrix3x3 dQ78 = transpose(dQ57)*Q58 + transpose(Q57)*dQ58;

    const sMatrix3x3 dQ28 = dQ23*Q35*Q58 + Q23*dQ35*Q58 + Q23*Q35*dQ58;
    const sMatrix3x3 dQ210 = dQ24*Q46*Q610 + Q24*dQ46*Q610 + Q24*Q46*dQ610;
    const sMatrix3x3 dQ810 = transpose(dQ28)*Q210 + transpose(Q28)*dQ210;

    const sMatrix3x3 dQ29 = dQ24*Q46*Q69 + Q24*dQ46*Q69 + Q24*Q46*dQ69;

    const sMatrix3x3 dQ19 = dQ12*Q29 + Q12*dQ29;

    // Full path omega matrices
    const sMatrix3x3 M_omega_47_7 = transpose(Q47)*dQ47;
    const sMatrix3x3 M_omega_47_4 = dQ47*transpose(Q47);

    const sMatrix3x3 M_omega_78_8 = transpose(Q78)*dQ78;
    const sMatrix3x3 M_omega_78_7 = dQ78*transpose(Q78);
   
    const sMatrix3x3 M_omega_810_10 = transpose(Q810)*dQ810;
    const sMatrix3x3 M_omega_810_8  = dQ810*transpose(Q810);

    const sMatrix3x3 M_omega_29_9 = transpose(Q29)*dQ29;
    const sMatrix3x3 M_omega_29_2 = dQ29*transpose(Q29);

    const sMatrix3x3 M_omega_19_9 = transpose(Q19)*dQ19;
    const sMatrix3x3 M_omega_19_1 = dQ19*transpose(Q19);

    // Full path omega vectors
    const sVector3d omega_47_7 = omega_vector_from_matrix(M_omega_47_7);
    const sVector3d omega_47_4 = omega_vector_from_matrix(M_omega_47_4);

    const sVector3d omega_78_8 = omega_vector_from_matrix(M_omega_78_8);
    const sVector3d omega_78_7 = omega_vector_from_matrix(M_omega_78_7);

    const sVector3d omega_810_10 = omega_vector_from_matrix(M_omega_810_10);
    const sVector3d omega_810_8  = omega_vector_from_matrix(M_omega_810_8);

    const sVector3d omega_29_9 = omega_vector_from_matrix(M_omega_29_9);
    const sVector3d omega_29_2  = omega_vector_from_matrix(M_omega_29_2);

    const sVector3d omega_19_9 = omega_vector_from_matrix(M_omega_19_9);
    const sVector3d omega_19_1  = omega_vector_from_matrix(M_omega_19_1);
};


inline sVector3d omega_euler(const double phi, const double theta, const double psi, 
                                       const double dphi, const double dtheta, const double dpsi)
/* https://mathworld.wolfram.com/EulerAngles.html
 * Purpose: to check that after three rotations (phi,theta,psi) around (Z,X,Z), the 
 * omega in the body frame is   (s(psi)s(theta)dphi+c(psi)dtheta)
 *              ----------  w = (c(psi)s(theta)dphi-s(psi)dtheta)
 *                              (c(theta)dphi+dpsi              )
*/

{
    return { sin(psi)*sin(theta)*dphi + cos(psi)*dtheta , 
             cos(psi)*sin(theta)*dphi - sin(psi)*dtheta ,
             cos(theta)*dphi + dpsi                       };
}

inline sMatrix3x3 Qcp_euler(const double phi, const double theta, const double psi)
{
    return { cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi),
             cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi),
             sin(psi)*sin(theta),
            -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi),
            -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi),
             cos(psi)*sin(theta),
             sin(theta)*sin(phi),
            -sin(theta)*cos(phi),
             cos(theta) } ;
}


inline sMatrix3x3 Qpc_euler(const double phi, const double theta, const double psi)
{
    return transpose(Qcp_euler(phi,theta,psi));
}


template<typename Frame_type, typename Parent_frame_type>
void test_frame_1_rotation(const Frame_type& frame, const Parent_frame_type& parent, Axis axis, double theta, double w)
{
    EXPECT_EQ(&frame.get_parent(), &parent);
    EXPECT_EQ(frame.get_parent_ptr(), &parent);
    EXPECT_EQ(frame.get_rotation_angles().size(), 1);
    EXPECT_EQ(frame.get_rotation_angles_derivative().size(), 1);
    EXPECT_EQ(frame.get_rotation_axis().size(), 1);

    EXPECT_DOUBLE_EQ(frame.get_rotation_angles().at(0), theta);
    EXPECT_DOUBLE_EQ(frame.get_rotation_angles_derivative().at(0), w);
    EXPECT_DOUBLE_EQ(frame.get_rotation_axis().at(0), axis);

    EXPECT_EQ(frame.is_inertial, false);
    EXPECT_EQ(frame.is_updated(), true);


    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(frame.get_origin().at(i), 0.0);
        EXPECT_EQ(frame.get_relative_velocity().at(i), 0.0);

        // Check that omega=(w_x,0,0) in both axis
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_body().at(i),w*(double)(i==axis));
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_parent().at(i),w*(double)(i==axis));
    }

    // Qpc and Qcp must be inverse matrices
    sMatrix3x3 QQinv = frame.get_rotation_matrix()*frame.get_back_rotation_matrix();

    for ( int i = 0; i < 3; ++i)
        for ( int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(QQinv(i,j), (double)(i==j));

    // Check Qpc
    sMatrix3x3 Qpc;
    sMatrix3x3 dQpc;
    switch(axis)
    {
     case X:
        Qpc = sMatrix3x3( 1.0, 0.0, 0.0,
                          0.0, std::cos(theta), -std::sin(theta),
                          0.0, std::sin(theta),  std::cos(theta));
        dQpc = w*drotation_matrix_x(theta);
        break;

     case Y:
        Qpc = sMatrix3x3(  std::cos(theta), 0.0, std::sin(theta),
                           0.0,             1.0, 0.0,
                          -std::sin(theta), 0.0, std::cos(theta));
        dQpc = w*drotation_matrix_y(theta);
        break;

     case Z:
        Qpc = sMatrix3x3( std::cos(theta), -std::sin(theta), 0.0,
                          std::sin(theta),  std::cos(theta), 0.0,
                          0.0,             0.0,              1.0 );
        dQpc = w*drotation_matrix_z(theta);
        break;
    }

    for ( int i = 0; i < 3; ++i)
        for ( int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(frame.get_rotation_matrix()(i,j), Qpc(i,j));

    // Check that omega is the result of a skew-symmetric tensor
    sVector3d omega_c = omega_vector_from_matrix(transpose(Qpc)*dQpc);
    sVector3d omega_p = omega_vector_from_matrix(dQpc*transpose(Qpc));
    
    for ( int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_body().at(i),omega_c[i]);
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_parent().at(i),omega_p[i]);
    }

    // Check that absolute rotations are the same than relative rotations (i.e. parent is inertial)
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_omega_absolute_in_body().at(i),omega_c[i]);
        EXPECT_DOUBLE_EQ(frame.get_omega_absolute_in_parent().at(i),omega_p[i]);
        EXPECT_DOUBLE_EQ(frame.get_omega_absolute_in_inertial().at(i),omega_p[i]);

        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(frame.get_absolute_rotation_matrix()(i,j),Qpc(i,j));
    }


        
}


template<typename Frame_type, typename Parent_frame_type>
void test_frame_2_rotations(const Frame_type& frame, const Parent_frame_type& parent, const std::array<Axis,2>& axis, 
                              const std::array<double,2>& theta, const std::array<double,2>& w)
{
    EXPECT_EQ(&frame.get_parent(), &parent);
    EXPECT_EQ(frame.get_parent_ptr(), &parent);
    EXPECT_EQ(frame.get_rotation_angles().size(), 2);
    EXPECT_EQ(frame.get_rotation_angles_derivative().size(), 2);
    EXPECT_EQ(frame.get_rotation_axis().size(), 2);

    for ( int i = 0; i < 2; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_rotation_angles().at(i), theta[i]);
        EXPECT_DOUBLE_EQ(frame.get_rotation_angles_derivative().at(i), w[i]);
        EXPECT_DOUBLE_EQ(frame.get_rotation_axis().at(i), axis[i]);
    }

    EXPECT_EQ(frame.is_inertial, false);
    EXPECT_EQ(frame.is_updated(), true);


    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(frame.get_origin().at(i), 0.0);
        EXPECT_EQ(frame.get_relative_velocity().at(i), 0.0);

    }


    // Qpc and Qcp must be inverse matrices
    sMatrix3x3 QQinv = frame.get_rotation_matrix()*frame.get_back_rotation_matrix();

    for ( int i = 0; i < 3; ++i)
        for ( int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(QQinv(i,j)+1.0, (double)(i==j)+1.0);

    // Check Qpc
    std::array<sMatrix3x3,2> Qij;
    std::array<sMatrix3x3,2> dQij;
    for ( int i = 0; i < 2; ++i) 
    {
        switch(axis[i])
        {
         case X:
            Qij[i] = sMatrix3x3( 1.0, 0.0               , 0.0,
                                0.0, std::cos(theta[i]), -std::sin(theta[i]),
                                0.0, std::sin(theta[i]),  std::cos(theta[i]));
            dQij[i] = w[i]*drotation_matrix_x(theta[i]);
            break;

         case Y:
            Qij[i] = sMatrix3x3(  std::cos(theta[i]), 0.0, std::sin(theta[i]),
                                  0.0,                1.0, 0.0,
                                 -std::sin(theta[i]), 0.0, std::cos(theta[i]));
            dQij[i] = w[i]*drotation_matrix_y(theta[i]);
            break;

         case Z:
            Qij[i] = sMatrix3x3( std::cos(theta[i]), -std::sin(theta[i]), 0.0,
                                 std::sin(theta[i]),  std::cos(theta[i]), 0.0,
                                 0.0,                 0.0,                1.0 );
            dQij[i] = w[i]*drotation_matrix_z(theta[i]);
            break;
        }
    }

    const sMatrix3x3 Qpc = Qij[0]*Qij[1];

    
    for ( int i = 0; i < 3; ++i)
        for ( int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(frame.get_rotation_matrix()(i,j), (Qij[0]*Qij[1])(i,j));

    // Test omega: omega|parent = w0 e_axis0 + Q01 w1 e_axis1
    //             omega|child  = Qcp omega|parent
    constexpr std::array<sVector3d,3> unitary{UX,UY,UZ};
    const sVector3d omega_p = w[0]*unitary[axis[0]] + Qij[0]*w[1]*unitary[axis[1]];
    const sVector3d omega_c = frame.get_back_rotation_matrix()*omega_p;

    for (int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_parent().at(i),omega_p[i]);
        EXPECT_NEAR(frame.get_omega_wrt_parent_in_body().at(i),omega_c[i], omega_c.norm()*3.0e-16);
    }

    // Check that omega is the result of a skew-symmetric tensor: tr(Q).dQ / dQ.tr(Q)
    sVector3d omega_c_from_tens = omega_vector_from_matrix(transpose(Qij[0]*Qij[1])*(dQij[0]*Qij[1]+Qij[0]*dQij[1]));
    sVector3d omega_p_from_tens = omega_vector_from_matrix((dQij[0]*Qij[1]+Qij[0]*dQij[1])*transpose(Qij[0]*Qij[1]));

    for ( int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_body().at(i),omega_c_from_tens[i]);
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_parent().at(i),omega_p_from_tens[i]);
    }

    // Check that absolute rotations are the same than relative rotations (i.e. parent is inertial)
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(frame.get_omega_absolute_in_body().at(i),omega_c[i], omega_c.norm()*3.0e-16);
        EXPECT_DOUBLE_EQ(frame.get_omega_absolute_in_parent().at(i),omega_p[i]);
        EXPECT_NEAR(frame.get_omega_absolute_in_inertial().at(i),omega_p[i], omega_p.norm()*3.0e-16);

        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(frame.get_absolute_rotation_matrix()(i,j),Qpc(i,j));
    }

}


template<typename Frame_type, typename Parent_frame_type>
void test_frame_3_rotations(const Frame_type& frame, const Parent_frame_type& parent, const std::array<Axis,3>& axis, 
                              const std::array<double,3>& theta, const std::array<double,3>& w)
{
    EXPECT_EQ(&frame.get_parent(), &parent);
    EXPECT_EQ(frame.get_parent_ptr(), &parent);
    EXPECT_EQ(frame.get_rotation_angles().size(), 3);
    EXPECT_EQ(frame.get_rotation_angles_derivative().size(), 3);
    EXPECT_EQ(frame.get_rotation_axis().size(), 3);

    for ( int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_rotation_angles().at(i), theta[i]);
        EXPECT_DOUBLE_EQ(frame.get_rotation_angles_derivative().at(i), w[i]);
        EXPECT_DOUBLE_EQ(frame.get_rotation_axis().at(i), axis[i]);
    }

    EXPECT_EQ(frame.is_inertial, false);
    EXPECT_EQ(frame.is_updated(), true);


    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(frame.get_origin().at(i), 0.0);
        EXPECT_EQ(frame.get_relative_velocity().at(i), 0.0);

    }


    // Qpc and Qcp must be inverse matrices
    sMatrix3x3 QQinv = frame.get_rotation_matrix()*frame.get_back_rotation_matrix();

    for ( int i = 0; i < 3; ++i)
        for ( int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(QQinv(i,j)+1.0, (double)(i==j)+1.0);

    // Check Qpc
    std::array<sMatrix3x3,3> Qij;
    std::array<sMatrix3x3,3> dQij;
    for ( int i = 0; i < 3; ++i) 
    {
        switch(axis[i])
        {
         case X:
            Qij[i] = sMatrix3x3( 1.0, 0.0               , 0.0,
                                0.0, std::cos(theta[i]), -std::sin(theta[i]),
                                0.0, std::sin(theta[i]),  std::cos(theta[i]));
            dQij[i] = w[i]*drotation_matrix_x(theta[i]);
            break;

         case Y:
            Qij[i] = sMatrix3x3(  std::cos(theta[i]), 0.0, std::sin(theta[i]),
                                  0.0,                1.0, 0.0,
                                 -std::sin(theta[i]), 0.0, std::cos(theta[i]));
            dQij[i] = w[i]*drotation_matrix_y(theta[i]);
            break;

         case Z:
            Qij[i] = sMatrix3x3( std::cos(theta[i]), -std::sin(theta[i]), 0.0,
                                 std::sin(theta[i]),  std::cos(theta[i]), 0.0,
                                 0.0,                 0.0,                1.0 );
            dQij[i] = w[i]*drotation_matrix_z(theta[i]);
            break;
        }
    }

    
    for ( int i = 0; i < 3; ++i)
        for ( int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(frame.get_rotation_matrix()(i,j), (Qij[0]*Qij[1]*Qij[2])(i,j));

    // Test omega: omega|parent = w0 e_axis0 + Q01 w1 e_axis1 + Q02 w2 e_axis2
    //             omega|child  = Qcp omega|parent
    constexpr std::array<sVector3d,3> unitary{UX,UY,UZ};
    const sVector3d omega_p = w[0]*unitary[axis[0]] + Qij[0]*(w[1]*unitary[axis[1]]+Qij[1]*w[2]*unitary[axis[2]]);
    const sVector3d omega_c = frame.get_back_rotation_matrix()*omega_p;

    for (int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_parent().at(i),omega_p[i]);
        EXPECT_NEAR(frame.get_omega_wrt_parent_in_body().at(i),omega_c[i], omega_c.norm()*3.0e-16);
    }

    // Check that omega is the result of a skew-symmetric tensor: tr(Q).dQ / dQ.tr(Q)
    sVector3d omega_c_from_tens = omega_vector_from_matrix(transpose(Qij[0]*Qij[1]*Qij[2])*(dQij[0]*Qij[1]*Qij[2]+Qij[0]*dQij[1]*Qij[2]+Qij[0]*Qij[1]*dQij[2]));
    sVector3d omega_p_from_tens = omega_vector_from_matrix((dQij[0]*Qij[1]*Qij[2]+Qij[0]*dQij[1]*Qij[2]+Qij[0]*Qij[1]*dQij[2])*transpose(Qij[0]*Qij[1]*Qij[2]));

    for ( int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_body().at(i),omega_c_from_tens[i]);
        EXPECT_DOUBLE_EQ(frame.get_omega_wrt_parent_in_parent().at(i),omega_p_from_tens[i]);
    }

    // Check that absolute rotations are the same than relative rotations (i.e. parent is inertial)
    sMatrix3x3 Qpc = Qij[0]*Qij[1]*Qij[2];
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(frame.get_omega_absolute_in_body().at(i),omega_c[i], omega_p.norm()*3.0e-16);
        EXPECT_DOUBLE_EQ(frame.get_omega_absolute_in_parent().at(i),omega_p[i]);
        EXPECT_DOUBLE_EQ(frame.get_omega_absolute_in_inertial().at(i),omega_p[i]);

        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(frame.get_absolute_rotation_matrix()(i,j),Qpc(i,j));
    }


}


TEST_F(Rotated_frames_test, Bad_construction_test)
{
    try
    {
        Frame<scalar,decltype(inertial_frame),2> error_frame(sVector3d(0.0), sVector3d(0.0), {0.0,1.0}, {0.0,1.0}, {X,static_cast<Axis>(5)}, inertial_frame);
        FAIL();
    }
    catch(const lion_exception& error)
    {
        SUCCEED();
    }
}


TEST_F(Rotated_frames_test, Frame_rotation_wrt_x) 
{ 
    test_frame_1_rotation(frame_rot_x, inertial_frame, X, theta_x, w_x); 

    auto local_frame(frame_rot_x);

    local_frame.set_rotation_angle(0, 60.0*DEG);
    test_frame_1_rotation(local_frame, inertial_frame, X, 60.0*DEG, w_x); 

    local_frame.set_rotation_angle(0, 120.0*DEG, 3.14);
    test_frame_1_rotation(local_frame, inertial_frame, X, 120.0*DEG, 3.14); 

    local_frame.set_angular_speed(0, 1.7);
    test_frame_1_rotation(local_frame, inertial_frame, X, 120.0*DEG, 1.7); 

    // For good measure
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Rotated_frames_test, Frame_rotation_wrt_y) 
{ 
    test_frame_1_rotation(frame_rot_y, inertial_frame, Y, theta_y, w_y); 

    auto local_frame(frame_rot_y);

    local_frame.set_rotation_angle(0, 40.0*DEG);
    test_frame_1_rotation(local_frame, inertial_frame, Y, 40.0*DEG, w_y); 

    local_frame.set_rotation_angle(0, -120.0*DEG, -3.14);
    test_frame_1_rotation(local_frame, inertial_frame, Y, -120.0*DEG, -3.14); 

    // For good measure
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Rotated_frames_test, Frame_rotation_wrt_z) 
{ 
    test_frame_1_rotation(frame_rot_z, inertial_frame, Z, theta_z, w_z); 

    auto local_frame(frame_rot_z);

    local_frame.set_rotation_angle(0, 150.0*DEG);
    test_frame_1_rotation(local_frame, inertial_frame, Z, 150.0*DEG, w_z); 

    local_frame.set_rotation_angle(0, 225.0*DEG, pi);
    test_frame_1_rotation(local_frame, inertial_frame, Z, 225.0*DEG, pi); 

    // For good measure
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Rotated_frames_test, Frame_rotation_wrt_xy) 
{ 
    test_frame_2_rotations(frame_rot_xy, inertial_frame, {X,Y}, {theta_x,theta_y}, {w_x,w_y} ); 

    auto local_frame(frame_rot_xy);

    local_frame.set_rotation_angle(1, 36.0*DEG, 40.0, false);
    EXPECT_EQ(local_frame.is_updated(), false);

    local_frame.set_rotation_angle(0, 15.0*DEG, -10.0);
    EXPECT_EQ(local_frame.is_updated(), true);

    test_frame_2_rotations(local_frame, inertial_frame, {X,Y},{15.0*DEG,36.0*DEG}, {-10.0,40.0} ); 

    try
    { 
        local_frame.set_rotation_angle(2, 36.0*DEG, 40.0);
        FAIL();
    }
    catch( const std::out_of_range& ex) 
    {
        SUCCEED();
    }

    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Rotated_frames_test, Frame_rotation_wrt_yz) 
{ 
    test_frame_2_rotations(frame_rot_yz, inertial_frame, {Y,Z}, {theta_y,theta_z}, {w_y,w_z} ); 

    // For good measure
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Rotated_frames_test, Frame_rotation_wrt_zx) 
{ 
    test_frame_2_rotations(frame_rot_zx, inertial_frame, {Z,X}, {theta_z,theta_x}, {w_z,w_x} ); 

    // For good measure
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Rotated_frames_test, Frame_rotation_wrt_xyz) 
{ 
    test_frame_3_rotations(frame_rot_xyz, inertial_frame, {X,Y,Z}, {theta_x,theta_y,theta_z}, {w_x,w_y,w_z} ); 

    // For good measure
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Rotated_frames_test, Euler_Angles_test)
{
    using namespace std;

    test_frame_3_rotations(euler_frame, inertial_frame, {Z,X,Z}, {phi,theta,psi},{dphi,dtheta,dpsi});

    for (size_t i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(euler_frame.get_omega_wrt_parent_in_body().at(i), omega_euler(phi,theta,psi,dphi,dtheta,dpsi).at(i));


    for ( size_t i = 0; i < 3; ++i )
        for ( size_t j = 0; j < 3; ++j ) 
        {
            EXPECT_DOUBLE_EQ(euler_frame.get_rotation_matrix()(i,j), Qpc_euler(phi,theta,psi)(i,j));
            EXPECT_DOUBLE_EQ(euler_frame.get_back_rotation_matrix()(i,j), Qcp_euler(phi,theta,psi)(i,j));
        }
}


TEST_F(Rotated_frames_test, Nested_rotation_frames)
{
    using namespace std;

    const Frame<scalar,decltype(inertial_frame),1> euler_1st_rot(sVector3d(0.0), sVector3d(0.0), {phi}, {dphi}, {Z}, inertial_frame);

    const Frame<scalar,decltype(euler_1st_rot),1> euler_2nd_rot(sVector3d(0.0), sVector3d(0.0), {theta}, {dtheta}, {X}, euler_1st_rot);

    const Frame<scalar,decltype(euler_2nd_rot),1> euler_3rd_rot(sVector3d(0.0), sVector3d(0.0), {psi}, {dpsi}, {Z}, euler_2nd_rot);

    // Check the three rotations
    for ( size_t i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(euler_3rd_rot.get_omega_absolute_in_body().at(i), omega_euler(phi,theta,psi,dphi,dtheta,dpsi).at(i));
        EXPECT_DOUBLE_EQ(euler_3rd_rot.get_omega_absolute_in_inertial().at(i),(Qpc_euler(phi,theta,psi) * omega_euler(phi,theta,psi,dphi,dtheta,dpsi)).at(i));
        EXPECT_DOUBLE_EQ(euler_3rd_rot.get_omega_absolute_in_parent().at(i),(rotation_matrix_z(psi) * omega_euler(phi,theta,psi,dphi,dtheta,dpsi)).at(i));

        EXPECT_DOUBLE_EQ(euler_3rd_rot.get_omega_in_body(inertial_frame).at(i), omega_euler(phi,theta,psi,dphi,dtheta,dpsi).at(i));
        EXPECT_DOUBLE_EQ(euler_3rd_rot.get_omega_in_target(inertial_frame).at(i),(Qpc_euler(phi,theta,psi) * omega_euler(phi,theta,psi,dphi,dtheta,dpsi)).at(i));
        EXPECT_DOUBLE_EQ(euler_3rd_rot.get_omega_in_parent(inertial_frame).at(i),(rotation_matrix_z(psi) * omega_euler(phi,theta,psi,dphi,dtheta,dpsi)).at(i));
    }

    for ( size_t i = 0; i < 3; ++i )
        for ( size_t j = 0; j < 3; ++j ) 
        {
            EXPECT_DOUBLE_EQ(euler_3rd_rot.get_absolute_rotation_matrix()(i,j), Qpc_euler(phi,theta,psi)(i,j));
            EXPECT_DOUBLE_EQ(euler_3rd_rot.get_rotation_matrix(inertial_frame)(i,j), Qpc_euler(phi,theta,psi)(i,j));
        }

    // Check the two first rotations
    for ( size_t i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(euler_2nd_rot.get_omega_absolute_in_body().at(i), omega_euler(phi,theta,0.0,dphi,dtheta,0.0).at(i));
        EXPECT_DOUBLE_EQ(euler_2nd_rot.get_omega_absolute_in_inertial().at(i),(Qpc_euler(phi,theta,0.0) * omega_euler(phi,theta,0.0,dphi,dtheta,0.0)).at(i));
        EXPECT_NEAR( euler_2nd_rot.get_omega_absolute_in_parent().at(i),
                    (rotation_matrix_x(theta) * omega_euler(phi,theta,0.0,dphi,dtheta,0.0)).at(i),
                     omega_euler(phi,theta,0.0,dphi,dtheta,0.0).norm()*3.0e-16);


        EXPECT_DOUBLE_EQ(euler_2nd_rot.get_omega_in_body(inertial_frame).at(i), omega_euler(phi,theta,0.0,dphi,dtheta,0.0).at(i));
        EXPECT_DOUBLE_EQ(euler_2nd_rot.get_omega_in_target(inertial_frame).at(i),(Qpc_euler(phi,theta,0.0) * omega_euler(phi,theta,0.0,dphi,dtheta,0.0)).at(i));
        EXPECT_NEAR( euler_2nd_rot.get_omega_in_parent(inertial_frame).at(i),
                    (rotation_matrix_x(theta) * omega_euler(phi,theta,0.0,dphi,dtheta,0.0)).at(i),
                     omega_euler(phi,theta,0.0,dphi,dtheta,0.0).norm()*3.0e-16);
    }

    for ( size_t i = 0; i < 3; ++i )
        for ( size_t j = 0; j < 3; ++j ) 
        {
            EXPECT_DOUBLE_EQ(euler_2nd_rot.get_absolute_rotation_matrix()(i,j), Qpc_euler(phi,theta,0.0)(i,j));
            EXPECT_DOUBLE_EQ(euler_2nd_rot.get_rotation_matrix(inertial_frame)(i,j), Qpc_euler(phi,theta,0.0)(i,j));
        }
}


TEST_F(Rotation_wrt_target_test, rotation_matrix_wrt_target)
{
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(frame7.get_rotation_matrix(frame4)(i,j), Q47(i,j));
            EXPECT_DOUBLE_EQ(frame4.get_rotation_matrix(frame7)(i,j), Q47(j,i));
            EXPECT_EQ(lioncpp::detail::get_crossing_generation(frame7,frame4),1);

            // Test the self rotation
            EXPECT_DOUBLE_EQ(frame4.get_rotation_matrix(frame4)(i,j), (double)(i==j));
            EXPECT_EQ(lioncpp::detail::get_crossing_generation(frame4,frame4),2);
        }


    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(frame8.get_rotation_matrix(frame7)(i,j)+1.0, Q78(i,j)+1.0);
            EXPECT_DOUBLE_EQ(frame7.get_rotation_matrix(frame8)(i,j)+1.0, Q78(j,i)+1.0);
            EXPECT_EQ(lioncpp::detail::get_crossing_generation(frame7,frame8),3);
        }


    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(frame10.get_rotation_matrix(frame8)(i,j)+1.0, Q810(i,j)+1.0);
            EXPECT_DOUBLE_EQ(frame8.get_rotation_matrix(frame10)(i,j)+1.0, Q810(j,i)+1.0);
            EXPECT_EQ(lioncpp::detail::get_crossing_generation(frame8,frame10),1);
        }


    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(frame9.get_rotation_matrix(frame2)(i,j)+1.0, Q29(i,j)+1.0);
            EXPECT_DOUBLE_EQ(frame2.get_rotation_matrix(frame9)(i,j)+1.0, Q29(j,i)+1.0);
            EXPECT_EQ(lioncpp::detail::get_crossing_generation(frame2,frame9),1);

            EXPECT_DOUBLE_EQ(frame9.get_rotation_matrix(frame1)(i,j)+1.0, Q19(i,j)+1.0);
            EXPECT_DOUBLE_EQ(frame1.get_rotation_matrix(frame9)(i,j)+1.0, Q19(j,i)+1.0);
            EXPECT_EQ(lioncpp::detail::get_crossing_generation(frame1,frame9),0);

        }

    try
    {
        frame4.get_rotation_matrix(frame11);
        FAIL();
    }
    catch (const lion_exception& error)
    {
        SUCCEED();
    }

    try
    {
        const size_t dummy(lioncpp::detail::get_crossing_generation(frame4, frame11));
        out(2) << dummy << std::endl;
        FAIL();
    }
    catch (const lion_exception& error)
    {
        SUCCEED();
    }
}


TEST_F(Rotation_wrt_target_test, omega_wrt_target)
{
    // Omega of 7 wrt 4 and viceversa
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame7.get_omega_in_body(frame4).at(i), omega_47_7.at(i));
        EXPECT_DOUBLE_EQ(frame7.get_omega_in_target(frame4).at(i), omega_47_4.at(i));

        EXPECT_DOUBLE_EQ(frame4.get_omega_in_body(frame7).at(i), -omega_47_4.at(i));
        EXPECT_DOUBLE_EQ(frame4.get_omega_in_target(frame7).at(i), -omega_47_7.at(i));
        EXPECT_DOUBLE_EQ(frame4.get_omega_in_parent(frame7).at(i), -(Q24*omega_47_4).at(i));

        // Test the self rotation
        EXPECT_DOUBLE_EQ(frame4.get_omega_in_body(frame4).at(i), 0.0);
        EXPECT_DOUBLE_EQ(frame4.get_omega_in_target(frame4).at(i), 0.0);
    }

    // Omega of 8 wrt 7 and viceversa
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame8.get_omega_in_body(frame7).at(i), omega_78_8.at(i));
        EXPECT_DOUBLE_EQ(frame8.get_omega_in_target(frame7).at(i), omega_78_7.at(i));
        EXPECT_DOUBLE_EQ(frame8.get_omega_in_parent(frame7).at(i),(Q58*omega_78_8).at(i));

        EXPECT_DOUBLE_EQ(frame7.get_omega_in_body(frame8).at(i), -omega_78_7.at(i));
        EXPECT_DOUBLE_EQ(frame7.get_omega_in_target(frame8).at(i), -omega_78_8.at(i));
        EXPECT_DOUBLE_EQ(frame7.get_omega_in_parent(frame8).at(i), -(Q57*omega_78_8).at(i));
    }

    // Omega of 10 wrt 8 and viceversa
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame10.get_omega_in_body(frame8).at(i), omega_810_10.at(i));
        EXPECT_DOUBLE_EQ(frame10.get_omega_in_target(frame8).at(i), omega_810_8.at(i));

        EXPECT_DOUBLE_EQ(frame8.get_omega_in_body(frame10).at(i), -omega_810_8.at(i));
        EXPECT_DOUBLE_EQ(frame8.get_omega_in_target(frame10).at(i), -omega_810_10.at(i));
    }

    // Omega of 9 wrt 2 and 1 and viceversa
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(frame9.get_omega_in_body(frame2).at(i), omega_29_9.at(i));
        EXPECT_NEAR(frame9.get_omega_in_target(frame2).at(i), omega_29_2.at(i),5.0e-16);

        EXPECT_NEAR(frame2.get_omega_in_body(frame9).at(i), -omega_29_2.at(i),5.0e-16);
        EXPECT_DOUBLE_EQ(frame2.get_omega_in_target(frame9).at(i), -omega_29_9.at(i));

        EXPECT_DOUBLE_EQ(frame9.get_omega_in_body(frame1).at(i), omega_19_9.at(i));
        EXPECT_DOUBLE_EQ(frame9.get_omega_in_target(frame1).at(i), omega_19_1.at(i));

        EXPECT_DOUBLE_EQ(frame1.get_omega_in_body(frame9).at(i), -omega_19_1.at(i));
        EXPECT_DOUBLE_EQ(frame1.get_omega_in_target(frame9).at(i), -omega_19_9.at(i));
    }

    try
    {
        frame4.get_omega_in_body(frame11);
        FAIL();
    }
    catch (const lion_exception& error)
    {
        SUCCEED();
    }

    try
    {
        frame4.get_omega_in_target(frame11);
        FAIL();
    }
    catch (const lion_exception& error)
    {
        SUCCEED();
    }

    try
    {
        frame4.get_omega_in_parent(frame11);
        FAIL();
    }
    catch (const lion_exception& error)
    {
        SUCCEED();
    }



}
