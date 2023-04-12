#include "gtest/gtest.h"
#include "lion/frame/frame.h"
#include "lion/math/matrix_extensions.h"

using Matrix4x4 = std::array<std::array<double,4>,4>;

constexpr Matrix4x4 total_transformation_matrix(const sMatrix3x3& Q, const sVector3d& x)
{
    return {{{ Q(0,0), Q(0,1), Q(0,2), x[0] },
             { Q(1,0), Q(1,1), Q(1,2), x[1] },
             { Q(2,0), Q(2,1), Q(2,2), x[2] },
             {    0.0,    0.0,    0.0,  1.0 }}};

}

constexpr Matrix4x4 dtotal_transformation_matrix(const sMatrix3x3& Q, const sVector3d& x)
{
    return {{{ Q(0,0), Q(0,1), Q(0,2), x[0] },
             { Q(1,0), Q(1,1), Q(1,2), x[1] },
             { Q(2,0), Q(2,1), Q(2,2), x[2] },
             {    0.0,    0.0,    0.0,  0.0 }}};

}


class Generic_position_velocity_test : public ::testing::Test
{
/*
*       This test defines a frame tree such that:
*                       0           10
*                       |
*                       1
*                     /   \
*                    2     3
*                    |     |      
*                    4     5
*                   / \   / \
*                  6   7 8  9
*/
 protected:
    Generic_position_velocity_test() {};

    const std::array<sVector3d,10> positions = {{{ 0.8147,0.1576,0.6557},
                                                 { 0.9058,0.9706,0.0357},
                                                 { 0.1270,0.9572,0.8491},
                                                 { 0.9134,0.4854,0.9340},
                                                 { 0.6324,0.8003,0.6787},
                                                 { 0.0975,0.1419,0.7577},
                                                 { 0.2785,0.4218,0.7431},
                                                 { 0.5469,0.9157,0.3922},
                                                 { 0.9575,0.7922,0.6555},
                                                 { 0.9649,0.9595,0.1712}}};


    const std::array<sVector3d,10> velocities = {{{ 0.7060,0.4387,0.2760},
                                                  { 0.0318,0.3816,0.6797},
                                                  { 0.2769,0.7655,0.6551},
                                                  { 0.0462,0.7952,0.1626},
                                                  { 0.0971,0.1869,0.1190},
                                                  { 0.8235,0.4898,0.4984},
                                                  { 0.6948,0.4456,0.9597},
                                                  { 0.3171,0.6463,0.3404},
                                                  { 0.9502,0.7094,0.5853},
                                                  { 0.0344,0.7547,0.2238}}};

    const std::array<double,10> angles = {{ 0.7513,
                                            0.2551,
                                            0.5060,
                                            0.6991,
                                            0.8909,
                                            0.9593,
                                            0.5472,
                                            0.1386,
                                            0.1493,
                                            0.2575 }};

    const std::array<double,10> dangles = {{ 0.8407,
                                             0.2543,
                                             0.8143,
                                             0.2435,
                                             0.9293,
                                             0.3500,
                                             0.1966,
                                             0.2511,
                                             0.6160,
                                             0.4733  }};

    const std::array<Axis,10> axis = {{ Y,Z,X,X,Y,Z,Y,Y,Z,X }};
   
    // Frames to be tested    
    const sFrame frame0 = sFrame(); // Inertial frame
    const sFrame frame1 = sFrame(positions[1],velocities[1], {angles[1]}, {dangles[1]}, {axis[1]}, frame0, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame2 = sFrame(positions[2],velocities[2], {angles[2]}, {dangles[2]}, {axis[2]}, frame1, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame3 = sFrame(positions[3],velocities[3], {angles[3]}, {dangles[3]}, {axis[3]}, frame1, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame4 = sFrame(positions[4],velocities[4], {angles[4]}, {dangles[4]}, {axis[4]}, frame2, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame5 = sFrame(positions[5],velocities[5], {angles[5]}, {dangles[5]}, {axis[5]}, frame3, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame6 = sFrame(positions[6],velocities[6], {angles[6]}, {dangles[6]}, {axis[6]}, frame4, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame7 = sFrame(positions[7],velocities[7], {angles[7]}, {dangles[7]}, {axis[7]}, frame4, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame8 = sFrame(positions[8],velocities[8], {angles[8]}, {dangles[8]}, {axis[8]}, frame5, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame9 = sFrame(positions[9],velocities[9], {angles[9]}, {dangles[9]}, {axis[9]}, frame5, sFrame::Frame_velocity_types::parent_frame);
    const sFrame frame10 = sFrame(); // Different inertial frame

    // Rotation matrices
    const sMatrix3x3 Q01 = rotation_matrix_z(angles[1]);
    const sMatrix3x3 Q12 = rotation_matrix_x(angles[2]);
    const sMatrix3x3 Q13 = rotation_matrix_x(angles[3]);
    const sMatrix3x3 Q24 = rotation_matrix_y(angles[4]);
    const sMatrix3x3 Q35 = rotation_matrix_z(angles[5]);
    const sMatrix3x3 Q46 = rotation_matrix_y(angles[6]);
    const sMatrix3x3 Q47 = rotation_matrix_y(angles[7]);
    const sMatrix3x3 Q58 = rotation_matrix_z(angles[8]);
    const sMatrix3x3 Q59 = rotation_matrix_x(angles[9]);

    // Rotation matrices derivatives
    const sMatrix3x3 dQ01 = dangles[1] * drotation_matrix_z(angles[1]);
    const sMatrix3x3 dQ12 = dangles[2] * drotation_matrix_x(angles[2]);
    const sMatrix3x3 dQ13 = dangles[3] * drotation_matrix_x(angles[3]);
    const sMatrix3x3 dQ24 = dangles[4] * drotation_matrix_y(angles[4]);
    const sMatrix3x3 dQ35 = dangles[5] * drotation_matrix_z(angles[5]);
    const sMatrix3x3 dQ46 = dangles[6] * drotation_matrix_y(angles[6]);
    const sMatrix3x3 dQ47 = dangles[7] * drotation_matrix_y(angles[7]);
    const sMatrix3x3 dQ58 = dangles[8] * drotation_matrix_z(angles[8]);
    const sMatrix3x3 dQ59 = dangles[9] * drotation_matrix_x(angles[9]);

    // Total transformation matrices
    const Matrix4x4 T01 = total_transformation_matrix(Q01, positions[1]);
    const Matrix4x4 T12 = total_transformation_matrix(Q12, positions[2]);
    const Matrix4x4 T13 = total_transformation_matrix(Q13, positions[3]);
    const Matrix4x4 T24 = total_transformation_matrix(Q24, positions[4]);
    const Matrix4x4 T35 = total_transformation_matrix(Q35, positions[5]);
    const Matrix4x4 T46 = total_transformation_matrix(Q46, positions[6]);
    const Matrix4x4 T47 = total_transformation_matrix(Q47, positions[7]);
    const Matrix4x4 T58 = total_transformation_matrix(Q58, positions[8]);
    const Matrix4x4 T59 = total_transformation_matrix(Q59, positions[9]);
    
    // Total transformation matrices inverses
    const Matrix4x4 T10 = total_transformation_matrix(transpose(Q01), -transpose(Q01)*positions[1]);
    const Matrix4x4 T21 = total_transformation_matrix(transpose(Q12), -transpose(Q12)*positions[2]);
    const Matrix4x4 T31 = total_transformation_matrix(transpose(Q13), -transpose(Q13)*positions[3]);
    const Matrix4x4 T42 = total_transformation_matrix(transpose(Q24), -transpose(Q24)*positions[4]);
    const Matrix4x4 T53 = total_transformation_matrix(transpose(Q35), -transpose(Q35)*positions[5]);
    const Matrix4x4 T64 = total_transformation_matrix(transpose(Q46), -transpose(Q46)*positions[6]);
    const Matrix4x4 T74 = total_transformation_matrix(transpose(Q47), -transpose(Q47)*positions[7]);
    const Matrix4x4 T85 = total_transformation_matrix(transpose(Q58), -transpose(Q58)*positions[8]);
    const Matrix4x4 T95 = total_transformation_matrix(transpose(Q59), -transpose(Q59)*positions[9]);

    // Total transformation matrices derivatives
    const Matrix4x4 dT01 = dtotal_transformation_matrix(dQ01, velocities[1]);
    const Matrix4x4 dT12 = dtotal_transformation_matrix(dQ12, velocities[2]);
    const Matrix4x4 dT13 = dtotal_transformation_matrix(dQ13, velocities[3]);
    const Matrix4x4 dT24 = dtotal_transformation_matrix(dQ24, velocities[4]);
    const Matrix4x4 dT35 = dtotal_transformation_matrix(dQ35, velocities[5]);
    const Matrix4x4 dT46 = dtotal_transformation_matrix(dQ46, velocities[6]);
    const Matrix4x4 dT47 = dtotal_transformation_matrix(dQ47, velocities[7]);
    const Matrix4x4 dT58 = dtotal_transformation_matrix(dQ58, velocities[8]);
    const Matrix4x4 dT59 = dtotal_transformation_matrix(dQ59, velocities[9]);

    // Total transformation matrices inverses derivatives
    const Matrix4x4 dT10 = dtotal_transformation_matrix(transpose(dQ01), -transpose(dQ01)*positions[1]-transpose(Q01)*velocities[1]);
    const Matrix4x4 dT21 = dtotal_transformation_matrix(transpose(dQ12), -transpose(dQ12)*positions[2]-transpose(Q12)*velocities[2]);
    const Matrix4x4 dT31 = dtotal_transformation_matrix(transpose(dQ13), -transpose(dQ13)*positions[3]-transpose(Q13)*velocities[3]);
    const Matrix4x4 dT42 = dtotal_transformation_matrix(transpose(dQ24), -transpose(dQ24)*positions[4]-transpose(Q24)*velocities[4]);
    const Matrix4x4 dT53 = dtotal_transformation_matrix(transpose(dQ35), -transpose(dQ35)*positions[5]-transpose(Q35)*velocities[5]);
    const Matrix4x4 dT64 = dtotal_transformation_matrix(transpose(dQ46), -transpose(dQ46)*positions[6]-transpose(Q46)*velocities[6]);
    const Matrix4x4 dT74 = dtotal_transformation_matrix(transpose(dQ47), -transpose(dQ47)*positions[7]-transpose(Q47)*velocities[7]);
    const Matrix4x4 dT85 = dtotal_transformation_matrix(transpose(dQ58), -transpose(dQ58)*positions[8]-transpose(Q58)*velocities[8]);
    const Matrix4x4 dT95 = dtotal_transformation_matrix(transpose(dQ59), -transpose(dQ59)*positions[9]-transpose(Q59)*velocities[9]);
};


TEST_F(Generic_position_velocity_test, frames_1_0)
{
    constexpr const sVector3d x(0.3517, 0.8308, 0.5853);
    constexpr const sVector3d dx(0.5497, 0.9172, 0.2858);

    constexpr const std::array<double,4> xT{x[0],x[1],x[2],1.0};
    constexpr const std::array<double,4> vT{dx[0],dx[1],dx[2],0.0};

    sVector3d position, velocity;

    // Frame 1 wrt frame 0
    std::tie(position,velocity) = frame1.get_position_and_velocity_in_target(frame0, x, dx);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position.at(i), (T01*xT).at(i),1.0e-16 );
        EXPECT_NEAR(velocity.at(i), (dT01*xT + T01*vT).at(i), 1.0e-16);
    }

    // Frame 0 wrt frame 1
    std::tie(position,velocity) = frame0.get_position_and_velocity_in_target(frame1, x, dx);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position.at(i), (T10*xT).at(i),5.0e-16 );
        EXPECT_NEAR(velocity.at(i), (dT10*xT + T10*vT).at(i), 5.0e-16);
    }

}


TEST_F(Generic_position_velocity_test, frames_6_1)
{
    constexpr const sVector3d x(0.3517, 0.8308, 0.5853);
    constexpr const sVector3d dx(0.5497, 0.9172, 0.2858);

    std::array<double,4> xT{x[0],x[1],x[2],1.0};
    std::array<double,4> vT{dx[0],dx[1],dx[2],0.0};

    sVector3d position, velocity;
    sVector3d position2, velocity2;


    // Frame 6 wrt frame 1
    std::tie(position,velocity) = frame6.get_position_and_velocity_in_target(frame1, x, dx);

    const Matrix4x4 T16(T12*T24*T46);
    const Matrix4x4 dT16(dT12*T24*T46 + T12*dT24*T46 + T12*T24*dT46);
  
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position.at(i), (T16*xT).at(i),5.0e-16);
        EXPECT_NEAR(velocity.at(i), (dT16*xT + T16*vT).at(i), 5.0e-16);
    }

    // Frame 1 wrt frame 6
    std::tie(position2,velocity2) = frame1.get_position_and_velocity_in_target(frame6, position, velocity);

    const Matrix4x4 T61(T64*T42*T21);
    const Matrix4x4 dT61(dT64*T42*T21 + T64*dT42*T21 + T64*T42*dT21);

    xT = {position[0], position[1], position[2], 1.0};
    vT = {velocity[0], velocity[1], velocity[2], 0.0};

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position2.at(i), x.at(i), 5.0e-16);
        EXPECT_NEAR(velocity2.at(i), dx.at(i), 5.0e-16);

        EXPECT_NEAR(position2.at(i), (T61*xT).at(i), 5.0e-16);
        EXPECT_NEAR(velocity2.at(i), (dT61*xT + T61*vT).at(i), 5.0e-16);
    }

}


TEST_F(Generic_position_velocity_test, frames_7_8)
{
    constexpr const sVector3d x(0.3517, 0.8308, 0.5853);
    constexpr const sVector3d dx(0.5497, 0.9172, 0.2858);

    std::array<double,4> xT{x[0],x[1],x[2],1.0};
    std::array<double,4> vT{dx[0],dx[1],dx[2],0.0};

    sVector3d position, velocity;
    sVector3d position2, velocity2;


    // Frame 7 wrt frame 8
    std::tie(position,velocity) = frame7.get_position_and_velocity_in_target(frame8, x, dx);

    const Matrix4x4 T87(T85*T53*T31*T12*T24*T47);
    const Matrix4x4 dT87(dT85*T53*T31*T12*T24*T47 + T85*dT53*T31*T12*T24*T47 + T85*T53*dT31*T12*T24*T47
                       + T85*T53*T31*dT12*T24*T47 + T85*T53*T31*T12*dT24*T47 + T85*T53*T31*T12*T24*dT47);
  
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position.at(i), (T87*xT).at(i),2.0e-15);
        EXPECT_NEAR(velocity.at(i), (dT87*xT + T87*vT).at(i), 2.0e-15);
    }

    // Frame 8 wrt frame 7
    std::tie(position2,velocity2) = frame8.get_position_and_velocity_in_target(frame7, position, velocity);

    const Matrix4x4 T78(T74*T42*T21*T13*T35*T58);
    const Matrix4x4 dT78(dT74*T42*T21*T13*T35*T58 + T74*dT42*T21*T13*T35*T58 + T74*T42*dT21*T13*T35*T58
                        + T74*T42*T21*dT13*T35*T58 + T74*T42*T21*T13*dT35*T58 + T74*T42*T21*T13*T35*dT58);

    xT = {position[0], position[1], position[2], 1.0};
    vT = {velocity[0], velocity[1], velocity[2], 0.0};

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position2.at(i), x.at(i), 2.0e-15);
        EXPECT_NEAR(velocity2.at(i), dx.at(i), 2.0e-15);

        EXPECT_NEAR(position2.at(i), (T78*xT).at(i), 2.0e-15);
        EXPECT_NEAR(velocity2.at(i), (dT78*xT + T78*vT).at(i), 2.0e-15);
    }
}


TEST_F(Generic_position_velocity_test, frames_9_0)
{
    constexpr const sVector3d x(0.3517, 0.8308, 0.5853);
    constexpr const sVector3d dx(0.5497, 0.9172, 0.2858);

    std::array<double,4> xT{x[0],x[1],x[2],1.0};
    std::array<double,4> vT{dx[0],dx[1],dx[2],0.0};

    sVector3d position, velocity;
    sVector3d position2, velocity2;


    // Frame 9 wrt frame 0
    std::tie(position,velocity) = frame9.get_position_and_velocity_in_target(frame0, x, dx);

    const Matrix4x4 T09(T01*T13*T35*T59);
    const Matrix4x4 dT09(dT01*T13*T35*T59 + T01*(dT13*T35*T59 + T13*(dT35*T59 + T35*dT59)));
  
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position.at(i), (T09*xT).at(i),2.0e-15);
        EXPECT_NEAR(velocity.at(i), (dT09*xT + T09*vT).at(i), 2.0e-15);
    }

    // Frame 0 wrt frame 9
    std::tie(position2,velocity2) = frame0.get_position_and_velocity_in_target(frame9, position, velocity);

    const Matrix4x4 T90(T95*T53*T31*T10);
    const Matrix4x4 dT90(dT95*T53*T31*T10 + T95*(dT53*T31*T10 + T53*(dT31*T10 + T31*dT10)));

    xT = {position[0], position[1], position[2], 1.0};
    vT = {velocity[0], velocity[1], velocity[2], 0.0};

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(position2.at(i), x.at(i), 2.0e-15);
        EXPECT_NEAR(velocity2.at(i), dx.at(i), 2.0e-15);

        EXPECT_NEAR(position2.at(i), (T90*xT).at(i), 2.0e-15);
        EXPECT_NEAR(velocity2.at(i), (dT90*xT + T90*vT).at(i), 2.0e-15);
    }

    // Using the absolute functions
    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(position.at(i), frame9.get_absolute_position(x).at(i));
        EXPECT_DOUBLE_EQ(velocity.at(i), frame9.get_absolute_velocity_in_inertial(x,dx).at(i));
        
        EXPECT_DOUBLE_EQ((transpose(frame9.get_absolute_rotation_matrix())*velocity).at(i),
                         frame9.get_absolute_velocity_in_body(x,dx).at(i));
    }
}
