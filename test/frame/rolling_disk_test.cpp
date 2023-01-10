#include "gtest/gtest.h"
#include "lion/frame/frame.h"
#include "lion/foundation/constants.h"
#include "lion/foundation/types.h"
#include <cmath>
#include <tuple>

void generic_inertial_frame_tests(sFrame& inertial_frame);

class Rolling_disk_test : public ::testing::TestWithParam<std::tuple<double,double>>
{
 protected:

    const scalar R = 0.25;
    const scalar x0 = 3.43;
    const scalar omega = pi;
    const scalar theta = pi/6.0;
    const scalar alpha = std::get<0>(GetParam());
    const scalar dalpha = std::get<1>(GetParam());

    sFrame inertial_frame = {};
    sFrame inclined_frame = {sVector3d(0.0), sVector3d(0.0), {alpha}, {dalpha}, {Z}, inertial_frame, sFrame::Frame_velocity_types::parent_frame};
    sFrame frame1         = {{x0,R,0.0}, {omega*R,0.0,0.0},{}, {}, {}, inclined_frame, sFrame::Frame_velocity_types::parent_frame};
    sFrame frame2         = {sVector3d(0.0),sVector3d(0.0),{theta},{-omega},{Z}, frame1, sFrame::Frame_velocity_types::parent_frame};
    sFrame frame3         = {sVector3d(0.0), sVector3d(0.0), {-theta}, {0.0}, {Z}, frame2, sFrame::Frame_velocity_types::parent_frame};
};


INSTANTIATE_TEST_SUITE_P(Rolling_disk_no_centrifugal, Rolling_disk_test, ::testing::Values(
                                                                            std::make_tuple(0.0,0.0),
                                                                            std::make_tuple(-pi/4.0,0.0),
                                                                            std::make_tuple(pi/6.0,0.0)));

INSTANTIATE_TEST_SUITE_P(Rolling_disk_centrifugal, Rolling_disk_test, ::testing::Values(
                                                                            std::make_tuple(pi/4.0,0.25),
                                                                            std::make_tuple(-pi/6.0,-0.10)));

TEST_P(Rolling_disk_test, get_origin) 
{
    EXPECT_DOUBLE_EQ(frame1.get_origin().at(0), x0);
    EXPECT_DOUBLE_EQ(frame1.get_origin().at(1),  R);
    EXPECT_DOUBLE_EQ(frame1.get_origin().at(2),0.0);

    for (size_t i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(frame2.get_origin().at(i),0.0);

    for (size_t i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(frame3.get_origin().at(i),0.0);
}



TEST_P(Rolling_disk_test, get_absolute_position)
{
    EXPECT_DOUBLE_EQ(frame1.get_absolute_position().at(0),x0*cos(alpha)-R*sin(alpha));
    EXPECT_DOUBLE_EQ(frame1.get_absolute_position().at(1),x0*sin(alpha)+R*cos(alpha));
    EXPECT_DOUBLE_EQ(frame1.get_absolute_position().at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame2.get_absolute_position().at(0), x0*cos(alpha) - R*sin(alpha));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_position().at(1), x0*sin(alpha) + R*cos(alpha));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_position().at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame2.get_absolute_position({-R*sin(theta),-R*cos(theta),0.0}).at(0), x0*cos(alpha));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_position({-R*sin(theta),-R*cos(theta),0.0}).at(1), x0*sin(alpha));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_position({-R*sin(theta),-R*cos(theta),0.0}).at(2),           0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_position().at(0),  x0*cos(alpha) - R*sin(alpha));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_position().at(1),  x0*sin(alpha) + R*cos(alpha));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_position().at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_position({0.0,-R,0.0}).at(0), x0*cos(alpha));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_position({0.0,-R,0.0}).at(1), x0*sin(alpha));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_position({0.0,-R,0.0}).at(2),           0.0);
}

TEST_P(Rolling_disk_test, get_relative_velocity_in_parent)
{
    EXPECT_DOUBLE_EQ(frame1.get_relative_velocity_in_parent().at(0),  omega*R);
    EXPECT_DOUBLE_EQ(frame1.get_relative_velocity_in_parent().at(1),      0.0);
    EXPECT_DOUBLE_EQ(frame1.get_relative_velocity_in_parent().at(2),      0.0);

    for (size_t i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(frame2.get_relative_velocity_in_parent().at(i), 0.0);

    for (size_t i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(frame3.get_relative_velocity_in_parent().at(i), 0.0);
}


TEST_P(Rolling_disk_test, get_absolute_velocity_in_inertial)
{
    EXPECT_DOUBLE_EQ(frame1.get_absolute_velocity_in_inertial().at(0),  omega*R*cos(alpha) - dalpha*(x0*sin(alpha)+R*cos(alpha)));
    EXPECT_DOUBLE_EQ(frame1.get_absolute_velocity_in_inertial().at(1),  omega*R*sin(alpha) + dalpha*(x0*cos(alpha)-R*sin(alpha)));
    EXPECT_DOUBLE_EQ(frame1.get_absolute_velocity_in_inertial().at(2),                 0.0);

    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_inertial().at(0),  omega*R*cos(alpha) - dalpha*(x0*sin(alpha)+R*cos(alpha)));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_inertial().at(1),  omega*R*sin(alpha) + dalpha*(x0*cos(alpha)-R*sin(alpha)));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_inertial().at(2),                 0.0);

    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_inertial({-R*sin(theta),-R*cos(theta),0.0}).at(0), -dalpha*x0*sin(alpha));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_inertial({-R*sin(theta),-R*cos(theta),0.0}).at(1),  dalpha*x0*cos(alpha));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_inertial({-R*sin(theta),-R*cos(theta),0.0}).at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial().at(0),  omega*R*cos(alpha) - dalpha*(x0*sin(alpha)+R*cos(alpha)));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial().at(1),  omega*R*sin(alpha) + dalpha*(x0*cos(alpha)-R*sin(alpha)));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial().at(2),                 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial({0.0,-R,0.0}).at(0), -dalpha*x0*sin(alpha));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial({0.0,-R,0.0}).at(1),  dalpha*x0*cos(alpha));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial({0.0,-R,0.0}).at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial({0.0,R,0.0}).at(0),  2.0*omega*R*cos(alpha) - dalpha*(x0*sin(alpha)+2.0*R*cos(alpha)));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial({0.0,R,0.0}).at(1),  2.0*omega*R*sin(alpha) + dalpha*(x0*cos(alpha)-2.0*R*sin(alpha)));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_inertial({0.0,R,0.0}).at(2),                 0.0);
}


TEST_P(Rolling_disk_test, get_absolute_velocity_in_body)
{
    EXPECT_DOUBLE_EQ(frame1.get_absolute_velocity_in_body().at(0),  omega*R-dalpha*R);
    EXPECT_NEAR     (frame1.get_absolute_velocity_in_body().at(1), dalpha*x0, 2.0e-16);
    EXPECT_DOUBLE_EQ(frame1.get_absolute_velocity_in_body().at(2),      0.0);

    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_body().at(0),  omega*R*cos(theta) + dalpha*x0*sin(theta)-dalpha*R*cos(theta));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_body().at(1), -omega*R*sin(theta) + dalpha*x0*cos(theta)+dalpha*R*sin(theta));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_body().at(2),                 0.0);

    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_body({-R*sin(theta),-R*cos(theta),0.0}).at(0), dalpha*x0*sin(theta));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_body({-R*sin(theta),-R*cos(theta),0.0}).at(1), dalpha*x0*cos(theta));
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_body({-R*sin(theta),-R*cos(theta),0.0}).at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_body().at(0),  omega*R - dalpha*R);
    EXPECT_NEAR     (frame3.get_absolute_velocity_in_body().at(1),  dalpha*x0, 2.0e-16 );
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_body().at(2),  0.0 );

    EXPECT_NEAR     (frame3.get_absolute_velocity_in_body({0.0,-R,0.0}).at(0), 0.0, 2.0e-16);
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_body({0.0,-R,0.0}).at(1), dalpha*x0);
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_body({0.0,-R,0.0}).at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_body({0.0,R,0.0}).at(0), 2.0*omega*R - 2.0*dalpha*R);
    EXPECT_NEAR     (frame3.get_absolute_velocity_in_body({0.0,R,0.0}).at(1), dalpha*x0, 4.0e-16);
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_body({0.0,R,0.0}).at(2), 0.0);
}


TEST_P(Rolling_disk_test, get_absolute_velocity_in_parent)
{
    EXPECT_DOUBLE_EQ(frame1.get_absolute_velocity_in_parent().at(0),  omega*R -dalpha*R);
    EXPECT_NEAR     (frame1.get_absolute_velocity_in_parent().at(1), dalpha*x0, 2.0e-16);
    EXPECT_DOUBLE_EQ(frame1.get_absolute_velocity_in_parent().at(2),      0.0);

    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_parent().at(0), omega*R-dalpha*R);
    EXPECT_NEAR     (frame2.get_absolute_velocity_in_parent().at(1), dalpha*x0, 2.0e-16);
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_parent().at(2), 0.0);

    EXPECT_NEAR     (frame2.get_absolute_velocity_in_parent({-R*sin(theta),-R*cos(theta),0.0}).at(0), 0.0, 2.0e-16);
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_parent({-R*sin(theta),-R*cos(theta),0.0}).at(1), dalpha*x0);
    EXPECT_DOUBLE_EQ(frame2.get_absolute_velocity_in_parent({-R*sin(theta),-R*cos(theta),0.0}).at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent().at(0),  omega*R*cos(theta) + dalpha*x0*sin(theta)-dalpha*R*cos(theta));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent().at(1), -omega*R*sin(theta) + dalpha*x0*cos(theta)+dalpha*R*sin(theta));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_body().at(2),                 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent({0.0,-R,0.0}).at(0), dalpha*x0*sin(theta));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent({0.0,-R,0.0}).at(1), dalpha*x0*cos(theta));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent({0.0,-R,0.0}).at(2), 0.0);

    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent({0.0,R,0.0}).at(0),   2.0*omega*R*cos(theta) 
                                                                                + dalpha*x0*sin(theta)-2.0*dalpha*R*cos(theta));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent({0.0,R,0.0}).at(1), - 2.0*omega*R*sin(theta) 
                                                                                + dalpha*x0*cos(theta)+2.0*dalpha*R*sin(theta));
    EXPECT_DOUBLE_EQ(frame3.get_absolute_velocity_in_parent({0.0,R,0.0}).at(2), 0.0);
}
