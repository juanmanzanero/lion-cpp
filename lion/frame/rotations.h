#ifndef LION_FRAME_ROTATIONS_H
#define LION_FRAME_ROTATIONS_H
#pragma once


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


template<typename T,
    typename Matrix3x3Type = Matrix3x3<T> >
constexpr Matrix3x3Type ea2rotmat(T yaw_rad, T pitch_rad, T roll_rad)
{
    //
    // Converts the input Euler angles (in radians and in standard
    // "yaw-pitch-roll" "Z-Y-X" sequence) to a rotation matrix. This
    // rotation matrix transforms from the rotated frame onto the
    // original one, i.e.,
    // "x_original = ea2rotmat(yaw_rad, pitch_rad, roll_rad) * x_rotated".
    //

    using std::cos;
    using std::sin;

    const auto cpsi = cos(yaw_rad);
    const auto spsi = sin(yaw_rad);
    const auto ctheta = cos(pitch_rad);
    const auto stheta = sin(pitch_rad);
    const auto cphi = cos(roll_rad);
    const auto sphi = sin(roll_rad);

    return Matrix3x3Type{ ctheta * cpsi, sphi * stheta * cpsi - spsi * cphi, cphi * cpsi * stheta + sphi * spsi,
        ctheta * spsi, sphi * stheta * spsi + cpsi * cphi, cphi * spsi * stheta - sphi * cpsi,
        -stheta, sphi * ctheta, cphi * ctheta };
}


template<typename Matrix3x3Type,
    typename Array3Type = std::array<typename Matrix3x3Type::value_type, 3u> >
constexpr Array3Type rotmat2ea(const Matrix3x3Type &M)
{
    //
    // Returns an std::array holding the Euler angles in 'ZYX' sequence
    // "[yaw, pitch, roll]", in rad, from the input rotation matrix,
    // given as an std::array in column-major format. This input rotation
    // matrix should transform from the rotated frame onto the original
    // one, i.e., "x_original = M * x_rotated".
    //

    using T = typename Array3Type::value_type;
    using std::atan2;
    using std::abs;
    using std::asin;

    // gymbal lock tolerance
    constexpr auto gymbal_lock_tol_deg = T{ 0.001 };
    constexpr auto gymbal_lock_zone = std::clamp(
        std::cos(gymbal_lock_tol_deg * DEG), T{ 0 }, T{ 1 });

    const auto sinp = -M[2];
    T pitch_rad;
    if (abs(sinp) <= gymbal_lock_zone) {
        pitch_rad = asin(sinp);
    }
    else if (sinp >= T{ 0 }) {
        pitch_rad = T{ 0.5 } * pi_T<T>;
    }
    else {
        pitch_rad = -T{ 0.5 } * pi_T<T>;
    }

    const auto yaw_rad = atan2(M[1], M[0]);
    const auto roll_rad = atan2(M[5], M[8]);

    return Array3Type{ wrap_to_pi(yaw_rad),
        wrap_to_pi(pitch_rad),
        wrap_to_pi(roll_rad) };
}


template<typename T,
    typename Matrix3x3Type = Matrix3x3<T> >
constexpr Matrix3x3Type tcs2rotmat(T toe_rad, T camber_rad, T spin_rad)
{
    //
    // Converts the input Euler angles (in radians and in "toe-camber-spin"
    // "Z-X-Y" sequence) to a rotation matrix. This rotation matrix
    // transforms from the rotated frame onto the original one, i.e.,
    // "x_original = tcs2rotmat(toe_rad, camber_rad, spin_rad) * x_rotated".
    //

    using std::cos;
    using std::sin;

    const auto ctoe = cos(toe_rad);
    const auto stoe = sin(toe_rad);
    const auto ccamber = cos(camber_rad);
    const auto scamber = sin(camber_rad);
    const auto cspin = cos(spin_rad);
    const auto sspin = sin(spin_rad);

    return Matrix3x3Type{ cspin * ctoe - scamber * sspin * stoe, -ccamber * stoe, ctoe * sspin + cspin * scamber * stoe,
        cspin * stoe + ctoe * scamber * sspin, ccamber * ctoe, sspin * stoe - cspin * ctoe * scamber,
        -ccamber * sspin, scamber, ccamber * cspin };
}

#endif