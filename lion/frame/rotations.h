#ifndef LION_FRAME_ROTATIONS_H
#define LION_FRAME_ROTATIONS_H
#pragma once


#include "lion/math/vector3d.h"
#include "lion/math/matrix3x3.h"


template<typename T>
struct Euler_angles : std::array<T, 3>
{
    using base_type = std::array<T, 3>;

    //! Default constructor
    constexpr Euler_angles() : base_type{} {}

    //! Constructor from coordinates
    constexpr Euler_angles(T yaw, T pitch, T roll) : base_type{ yaw, pitch, roll } {};

    //! Constructor from size 3 array
    constexpr explicit Euler_angles(const base_type& arr) : base_type{ arr } {};

    //! Get yaw
    constexpr const T& yaw() const { return (*this)[0]; }
    constexpr       T& yaw()       { return (*this)[0]; }

    //! Get pitch
    constexpr const T& pitch() const { return (*this)[1]; }
    constexpr       T& pitch()       { return (*this)[1]; }

    //! Get roll
    constexpr const T& roll() const { return (*this)[2]; }
    constexpr       T& roll()       { return (*this)[2]; }
};


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
    // "yaw-pitch-roll" "ZYX" sequence) to a rotation matrix. This
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

    return Matrix3x3Type{
        ctheta * cpsi, sphi * stheta * cpsi - spsi * cphi, cphi * cpsi * stheta + sphi * spsi,
        ctheta * spsi, sphi * stheta * spsi + cpsi * cphi, cphi * spsi * stheta - sphi * cpsi,
        -stheta, sphi * ctheta, cphi * ctheta };
}


template<typename Matrix3x3Type,
    typename Array3Type = std::array<typename Matrix3x3Type::value_type, 3u> >
constexpr Array3Type rotmat2ea(const Matrix3x3Type &M)
{
    //
    // Returns an std::array holding the Euler angles in "ZYX" sequence
    // "[yaw, pitch, roll]", in rad, from the input rotation matrix
    // (specified as an array in column-major order). This input rotation
    // matrix should transform from the rotated frame onto the original
    // one, i.e., "x_original = M * x_rotated".
    //

    using T = typename Array3Type::value_type;
    using std::abs;
    using std::atan2;
    using std::asin;

    // gymbal lock tolerance
    //constexpr auto gymbal_lock_tol_deg = T{ 1e-3 };
    //constexpr auto gymbal_lock_zone = std::clamp(
    //    std::cos(gymbal_lock_tol_deg * DEG), T{ 0 }, T{ 1 });
    constexpr auto gymbal_lock_zone = T{ 0.99999999984769128 };

    const auto spitch = -M[2];
    T pitch_rad;
    if (abs(spitch) <= gymbal_lock_zone) {
        pitch_rad = asin(spitch);
    }
    else if (spitch >= T{ 0 }) {
        pitch_rad = T{ 0.5 } * pi_T<T>;
    }
    else {
        pitch_rad = -T{ 0.5 } * pi_T<T>;
    }

    const auto yaw_rad = atan2(M[1], M[0]);
    const auto roll_rad = atan2(M[5], M[8]);

    return Array3Type{
        wrap_to_pi(yaw_rad),
        wrap_to_pi(pitch_rad),
        wrap_to_pi(roll_rad) };
}


template<typename T, typename Vector3dType,
    typename Matrix3x3Type = Matrix3x3<T> >
constexpr Matrix3x3Type angax2rotmat(T angle_rad, const Vector3dType &axis)
{
    //
    // Converts the input angle-axis pair (with the angle given in
    // radians) to a rotation matrix. This rotation matrix transforms
    // from the rotated frame onto the original one, i.e.,
    // "x_original = angax2rotmat(angle_rad, axis) * x_rotated".
    //

    using std::sqrt;
    using std::cos;
    using std::sin;

    const auto invnorm = sqrt(axis[0] * axis[0] +
        axis[1] * axis[1] + axis[2] * axis[2]);

    const auto x = invnorm * axis[0];
    const auto y = invnorm * axis[1];
    const auto z = invnorm * axis[2];

    const auto c = cos(angle_rad);
    const auto s = sin(angle_rad);
    const auto C = T{ 1 } - c;

    const auto xyC = x * y * C;
    const auto xzC = x * z * C;
    const auto yzC = y * z * C;

    const auto xs = x * s;
    const auto ys = y * s;
    const auto zs = z * s;

    return Matrix3x3Type{
        x * x * C + c, xyC - zs, xzC + ys,
        xyC + zs, y * y * C + c, yzC - xs,
        xzC - ys, yzC + xs, z * z * C + c };
}


template<typename T,
    typename Matrix3x3Type = Matrix3x3<T> >
constexpr Matrix3x3Type sis2rotmat(T steer_rad, T inclination_rad, T spin_rad)
{
    //
    // Converts the input Euler angles (in radians and in "steer-inclination-spin"
    // "ZXY" sequence) to a rotation matrix. This rotation matrix
    // transforms from the rotated frame onto the original one, i.e.,
    // "x_original = sis2rotmat(steer_rad, inclination_rad, spin_rad) * x_rotated".
    //

    using std::cos;
    using std::sin;

    const auto csteer = cos(steer_rad);
    const auto ssteer = sin(steer_rad);
    const auto cinclination = cos(inclination_rad);
    const auto sinclination = sin(inclination_rad);
    const auto cspin = cos(spin_rad);
    const auto sspin = sin(spin_rad);

    return Matrix3x3Type{
        cspin * csteer - sinclination * sspin * ssteer, -cinclination * ssteer, csteer * sspin + cspin * sinclination * ssteer,
        cspin * ssteer + csteer * sinclination * sspin, cinclination * csteer, sspin * ssteer - cspin * csteer * sinclination,
        -cinclination * sspin, sinclination, cinclination * cspin };
}


template<typename Matrix3x3Type,
    typename Array3Type = std::array<typename Matrix3x3Type::value_type, 3u> >
constexpr Array3Type rotmat2sis(const Matrix3x3Type &M)
{
    //
    // Returns an std::array holding the Euler angles in "ZXY" sequence
    // "[steer, inclination, spin]", in rad, from the input rotation matrix
    // (specified as an array in column-major order). This input rotation
    // matrix should transform from the rotated frame onto the original
    // one, i.e., "x_original = M * x_rotated".
    //

    using T = typename Array3Type::value_type;
    using std::abs;
    using std::atan2;
    using std::asin;

    // gymbal lock tolerance
    //constexpr auto gymbal_lock_tol_deg = T{ 1e-3 };
    //constexpr auto gymbal_lock_zone = std::clamp(
    //    std::cos(gymbal_lock_tol_deg * DEG), T{ 0 }, T{ 1 });
    constexpr auto gymbal_lock_zone = T{ 0.99999999984769128 };

    const auto sinclination = M[5];
    T inclination_rad;
    if (abs(sinclination) <= gymbal_lock_zone) {
        inclination_rad = asin(sinclination);
    }
    else if (sinclination >= T{ 0 }) {
        inclination_rad = T{ 0.5 } * pi_T<T>;
    }
    else {
        inclination_rad = -T{ 0.5 } * pi_T<T>;
    }

    const auto steer_rad = atan2(-M[3], M[4]);
    const auto spin_rad = atan2(-M[2], M[8]);

    return Array3Type{
        wrap_to_pi(steer_rad),
        wrap_to_pi(inclination_rad),
        wrap_to_pi(spin_rad) };
}


template<typename T,
    typename Array3Type = std::array<T, 3u> >
constexpr Array3Type angular_kinematic_relationships(T yawdot, T pitchdot, T rolldot,
    T yaw_rad, T pitch_rad, T roll_rad)
{
    //
    // Returns an array holding the angular velocity components
    // "[p, q, r]" of a rotating frame (projected onto this same
    // frame), from the time derivatives of its Euler angles in
    // "ZYX" sequence "[yawdot, pitchdot, rolldot]" and the values
    // of these angles "[yaw_rad, pitch_rad, roll_rad]" (given in
    // radians). The resulting angular velocity components will
    // share the units of the input Euler angle derivatives.
    //

    using std::cos;
    using std::sin;

    const auto cphi = cos(roll_rad);
    const auto sphi = sin(roll_rad);
    const auto yawdot_ctheta = yawdot * cos(pitch_rad);

    return Array3Type{
        rolldot - yawdot * sin(pitch_rad),
        yawdot_ctheta * sphi + pitchdot * cphi,
        yawdot_ctheta * cphi - pitchdot * sphi };
}



template<typename T,
    typename Array3Type = std::array<T, 3u> >
constexpr Array3Type angular_kinematic_relationships_derivative(const T& yawdotdot, const T& pitchdotdot, const T& rolldotdot, 
    const T& yawdot, const T& pitchdot, const T& rolldot,
    const T& yaw_rad, const T& pitch_rad, const T& roll_rad)
{
    //
    // Returns an array holding the derivative of the angular velocity components
    // "[p, q, r]" of a rotating frame (projected onto this same
    // frame), from the time derivatives of its Euler angles in
    // "ZYX" sequence "[yawdot, pitchdot, rolldot]" and the values
    // of these angles "[yaw_rad, pitch_rad, roll_rad]" (given in
    // radians). The resulting angular velocity components will
    // share the units of the input Euler angle derivatives.
    //

    using std::cos;
    using std::sin;

    const auto spitch = sin(pitch_rad);
    const auto cpitch = cos(pitch_rad);

    const auto sroll = sin(roll_rad);
    const auto croll = cos(roll_rad);

    return Array3Type{
        rolldotdot - yawdotdot * spitch - pitchdot * yawdot * cpitch,
        pitchdotdot * croll - rolldot * (pitchdot * sroll - yawdot * cpitch * croll) + yawdotdot * cpitch * sroll - pitchdot * yawdot * spitch * sroll,
        yawdotdot * cpitch * croll - pitchdotdot * sroll - rolldot * (pitchdot * croll + yawdot * cpitch * sroll) - pitchdot * yawdot * croll * spitch
    };
}




template<typename T,
    typename Array3Type = std::array<T, 3u> >
constexpr Array3Type inverse_angular_kinematic_relationships(T p, T q, T r,
    T yaw_rad, T pitch_rad, T roll_rad)
{
    //
    // Returns an array holding the time derivatives of a rotating
    // frame's Euler angles in "ZYX" sequence "[yawdot, pitchdot,
    // rolldot]" from its angular velocity components "[p, q, r]"
    // (projected onto the rotating frame) and its Euler angles
    // "[yaw_rad, pitch_rad, roll_rad]" (given in radians). These
    // output derivatives will share the units of the input angular
    // velocity components.
    //

    using std::cos;
    using std::sin;

    const auto cphi = cos(roll_rad);
    const auto sphi = sin(roll_rad);
    const auto yawdot = (q * sphi + r * cphi) / cos(pitch_rad);

    return Array3Type{
        yawdot,
        q * cphi - r * sphi,
        p + yawdot * sin(pitch_rad) };
}

#endif
