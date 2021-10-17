#ifndef __ROTATIONS_H__
#define __ROTATIONS_H__

//! Get a rotation matrix around the X axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> rotation_matrix_x(const T phi);

//! Get a rotation matrix around the Y axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> rotation_matrix_y(const T phi);

//! Get a rotation matrix around the Z axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> rotation_matrix_z(const T phi);

//! Get the derivative of the rotation matrix around the X axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> drotation_matrix_x(const T phi);

//! Get the derivative of the rotation matrix around the Y axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> drotation_matrix_y(const T phi);

//! Get the derivative of the rotation matrix around the Z axis
//! @param[in] phi: rotation angle [rad]
template<typename T>
constexpr Matrix3x3<T> drotation_matrix_z(const T phi);

//! Get the omega vector from the omega matrix
//! @param[in] matrix: omega matrix tensor
template<typename T>
constexpr Vector3d<T> omega_vector_from_matrix(const Matrix3x3<T>& matrix);

#include "rotations.hpp"

#endif
