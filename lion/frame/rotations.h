#ifndef ROTATIONS_H
#define ROTATIONS_H

#include "lion/math/vector3d.h"
#include "lion/math/matrix3x3.h"

//! Get a rotation matrix around the X axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> rotation_matrix_x(const T phi)
{
    return Matrix3x3<T>( 1.0,      0.0,       0.0,
                         0.0, cos(phi), -sin(phi),
                         0.0, sin(phi), cos(phi) );
}


//! Get a rotation matrix around the Y axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> rotation_matrix_y(const T phi)
{
    return Matrix3x3<T>(  cos(phi), 0.0, sin(phi),
		                 0.0, 1.0,           0.0,
                      -sin(phi), 0.0, cos(phi) ); 
}


//! Get a rotation matrix around the Z axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> rotation_matrix_z(const T phi)
{
    return Matrix3x3<T>(  cos(phi), -sin(phi), 0.0,
                       sin(phi),  cos(phi), 0.0,
                                 0.0,            0.0, 1.0 ); 
}


//! Get the derivative of the rotation matrix around the X axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> drotation_matrix_x(const T phi)
{
    return Matrix3x3<T>( 0.0,            0.0,            0.0,
                      0.0, -sin(phi), -cos(phi),
                      0.0,  cos(phi), -sin(phi) );
}


//! Get the derivative of the rotation matrix around the Y axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> drotation_matrix_y(const T phi)
{
    return Matrix3x3<T>( -sin(phi), 0.0,  cos(phi),
		                 0.0        , 0.0,            0.0,
                      -cos(phi), 0.0, -sin(phi) ); 
}


//! Get the derivative of the rotation matrix around the Z axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> drotation_matrix_z(const T phi)
{
    return Matrix3x3<T>( -sin(phi), -cos(phi), 0.0,
                       cos(phi), -sin(phi), 0.0,
                                 0.0,            0.0, 0.0 ); 
}


//! Get the omega vector from the omega matrix
//! @param[in] matrix: omega matrix tensor
template<typename T>
constexpr Vector3d<T> omega_vector_from_matrix(const Matrix3x3<T>& matrix)
{
    return Vector3d<T>(matrix(2,1),matrix(0,2),matrix(1,0));
}


//! Get the omega tensor from the omega vector
//! @param[in] vector: omega vector
template<typename T>
constexpr Matrix3x3<T> omega_tensor_from_vector(const Vector3d<T>& vector)
{
    return Matrix3x3<T>(0.0,      -vector[2], vector[1],
                        vector[2],       0.0, -vector[0],
                       -vector[1], vector[0],       0.0);
}

#endif
