#include "gtest/gtest.h"
#include "lion/frame/frame.h"
#include "lion/foundation/constants.h"
#include <cmath>
#include "lion/thirdparty/include/logger.hpp"

class Inertial_frame_test : public ::testing::Test
{
 protected:
    Inertial_frame_test(): inertial_frame()
                  {} ;


    Frame inertial_frame;
};


void generic_inertial_frame_tests(Frame& inertial_frame)
{
    try
    { 
        const Frame& parent(inertial_frame.get_parent());
        out(2) << &parent << std::endl;
        FAIL();
    } 
    catch( const std::runtime_error& error )
    {
        SUCCEED();
    } 

    inertial_frame.update();

    EXPECT_EQ(inertial_frame.get_parent_ptr(), nullptr);

    for ( size_t i = 0; i < 3; i++ )
    {
        EXPECT_EQ(inertial_frame.get_origin().at(i), 0.0 );
        EXPECT_EQ(inertial_frame.get_relative_velocity().at(i), 0.0);
    } 

    EXPECT_EQ(inertial_frame.get_rotation_angles().size(), 0u);
    EXPECT_EQ(inertial_frame.get_rotation_axis().size(), 0u);

    EXPECT_EQ(inertial_frame.is_inertial(), true);
    EXPECT_EQ(inertial_frame.is_updated(), true);

    for ( size_t i = 0; i < 3; i++ ) {
        for ( size_t j = 0; j < 3; j++) 
        {
            EXPECT_EQ(inertial_frame.get_rotation_matrix()(i,j), (int)(i==j)); 
            EXPECT_EQ(inertial_frame.get_absolute_rotation_matrix()(i,j), (int)(i==j)); 
        }
    }

    for ( size_t i = 0; i < 3; i++ )
    { 
        EXPECT_EQ(inertial_frame.get_omega_wrt_parent_in_body().at(i), 0.0);
        EXPECT_EQ(inertial_frame.get_omega_wrt_parent_in_parent().at(i), 0.0);
        EXPECT_EQ(inertial_frame.get_omega_absolute_in_body().at(i), 0.0);
        EXPECT_EQ(inertial_frame.get_omega_absolute_in_parent().at(i), 0.0);
        EXPECT_EQ(inertial_frame.get_omega_absolute_in_inertial().at(i), 0.0);
    }
}


TEST_F(Inertial_frame_test, Inertial_frame_test_parent) 
{ 
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Inertial_frame_test, Inertial_frame_test_set_origin) 
{ 
    inertial_frame.set_origin({ 1.0, 2.0, 3.0} );
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Inertial_frame_test, Inertial_frame_test_set_origin_and_velocity) 
{ 
    inertial_frame.set_origin({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0} );
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Inertial_frame_test, Inertial_frame_test_set_velocity) 
{ 
    inertial_frame.set_velocity({9.0, 10.0, 11.0} );
    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Inertial_frame_test, Inertial_frame_test_set_rot_angle) 
{ 
    try
    {
        inertial_frame.set_rotation_angle(0, 40.0);
        FAIL();
    } 
    catch( const std::out_of_range& ex)
    {
        SUCCEED();
    } 

    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Inertial_frame_test, Inertial_frame_test_set_rot_angle_and_speed) 
{ 
    try
    {
        inertial_frame.set_rotation_angle(0, 40.0, 100.0);
        FAIL();
    } 
    catch( const std::out_of_range& ex)
    {
        SUCCEED();
    } 

    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Inertial_frame_test, Inertial_frame_test_set_rot_angular_speed) 
{ 
    try
    {
        inertial_frame.set_angular_speed(0, 80.0);
        FAIL();
    } 
    catch( const std::out_of_range& ex)
    {
        SUCCEED();
    } 

    generic_inertial_frame_tests(inertial_frame);
}


TEST_F(Inertial_frame_test, Inertial_frame_test_add_rotation) 
{ 
    inertial_frame.add_rotation(4.0, 5.0, X);

    generic_inertial_frame_tests(inertial_frame);
}
