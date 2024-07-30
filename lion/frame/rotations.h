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
constexpr Matrix3x3Type ea2rotmat(const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
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

    const auto cyaw = cos(yaw_rad);
    const auto syaw = sin(yaw_rad);
    const auto cpitch = cos(pitch_rad);
    const auto spitch = sin(pitch_rad);
    const auto croll = cos(roll_rad);
    const auto sroll = sin(roll_rad);

    return Matrix3x3Type{
        cpitch * cyaw, sroll * spitch * cyaw - syaw * croll, croll * cyaw * spitch + sroll * syaw,
        cpitch * syaw, sroll * spitch * syaw + cyaw * croll, croll * syaw * spitch - sroll * cyaw,
        -spitch, sroll * cpitch, croll * cpitch };
}


template<typename T,
         typename Matrix3x3Type = Matrix3x3<T> >
constexpr Matrix3x3Type sis2rotmat(const T &steer_rad, const T &inclination_rad, const T &spin_rad)
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
constexpr Array3Type rotmat2ea(const Matrix3x3Type &T_colmaj)
{
    //
    // Returns an std::array holding the Euler angles in "ZYX" sequence
    // "[yaw, pitch, roll]", in rad, from the input rotation matrix
    // (specified as an array in column-major order). This input rotation
    // matrix should transform from the rotated frame onto the original
    // one, i.e., "x_original = T_colmaj * x_rotated".
    //

    using value_type = typename Array3Type::value_type;

    using std::abs;
    using std::atan2;
    using std::asin;

    // gymbal lock tolerance
    //constexpr auto gymbal_lock_tol_deg = value_type{ 1e-3 };
    //constexpr auto gymbal_lock_zone = std::clamp(
    //    std::cos(gymbal_lock_tol_deg * DEG), value_type{ 0 }, value_type{ 1 });
    constexpr auto gymbal_lock_zone = value_type{ 0.99999999984769128 };

    const auto spitch = -T_colmaj[2];
    value_type pitch_rad;
    if (abs(spitch) <= gymbal_lock_zone) {
        pitch_rad = asin(spitch);
    }
    else if (spitch >= value_type{ 0 }) {
        pitch_rad = value_type{ 0.5 } * pi_T<value_type>;
    }
    else {
        pitch_rad = -value_type{ 0.5 } * pi_T<value_type>;
    }

    const auto yaw_rad = atan2(T_colmaj[1], T_colmaj[0]);
    const auto roll_rad = atan2(T_colmaj[5], T_colmaj[8]);

    return Array3Type{
        wrap_to_pi(yaw_rad),
        wrap_to_pi(pitch_rad),
        wrap_to_pi(roll_rad) };
}


template<typename Matrix3x3Type,
         typename Array3Type = std::array<typename Matrix3x3Type::value_type, 3u> >
constexpr Array3Type rotmat2sis(const Matrix3x3Type &T_colmaj)
{
    //
    // Returns an std::array holding the Euler angles in "ZXY" sequence
    // "[steer, inclination, spin]", in rad, from the input rotation matrix
    // (specified as an array in column-major order). This input rotation
    // matrix should transform from the rotated frame onto the original
    // one, i.e., "x_original = T_colmaj * x_rotated".
    //

    using value_type = typename Array3Type::value_type;

    using std::abs;
    using std::atan2;
    using std::asin;

    // gymbal lock tolerance
    //constexpr auto gymbal_lock_tol_deg = value_type{ 1e-3 };
    //constexpr auto gymbal_lock_zone = std::clamp(
    //    std::cos(gymbal_lock_tol_deg * DEG), value_type{ 0 }, value_type{ 1 });
    constexpr auto gymbal_lock_zone = value_type{ 0.99999999984769128 };

    const auto sinclination = T_colmaj[5];
    value_type inclination_rad;
    if (abs(sinclination) <= gymbal_lock_zone) {
        inclination_rad = asin(sinclination);
    }
    else if (sinclination >= value_type{ 0 }) {
        inclination_rad = value_type{ 0.5 } * pi_T<value_type>;
    }
    else {
        inclination_rad = -value_type{ 0.5 } * pi_T<value_type>;
    }

    const auto steer_rad = atan2(-T_colmaj[3], T_colmaj[4]);
    const auto spin_rad = atan2(-T_colmaj[2], T_colmaj[8]);

    return Array3Type{
        wrap_to_pi(steer_rad),
        wrap_to_pi(inclination_rad),
        wrap_to_pi(spin_rad) };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type angular_kinematic_relationships(const T &yawdot, const T &pitchdot, const T &rolldot,
                                                     const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
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

    (void)yaw_rad;

    using std::cos;
    using std::sin;

    const auto croll = cos(roll_rad);
    const auto sroll = sin(roll_rad);
    const auto yawdot_ctheta = yawdot * cos(pitch_rad);

    return Array3Type{
        rolldot - yawdot * sin(pitch_rad),
        yawdot_ctheta * sroll + pitchdot * croll,
        yawdot_ctheta * croll - pitchdot * sroll };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type sis_kinematic_relationships(const T &steerdot, const T &inclinationdot, const T &spindot,
                                                 const T &steer_rad, const T &inclination_rad, const T &spin_rad)
{
    //
    // Returns an array holding the angular velocity components
    // "[p, q, r]" of a rotating frame (projected onto this same
    // frame), from the time derivatives of its Euler angles in
    // "ZXY" sequence "[steerdot, inclinationdot, spindot]" and the
    // values of these angles "[steer_rad, inclination_rad, spin_rad]"
    // (given in radians). The resulting angular velocity components
    // will share the units of the input Euler angle derivatives.
    //

    (void)steer_rad;

    using std::cos;
    using std::sin;

    const auto cspin = cos(spin_rad);
    const auto sspin = sin(spin_rad);
    const auto steerdot_cinclination = steerdot * cos(inclination_rad);

    return Array3Type{
        inclinationdot * cspin - steerdot_cinclination * sspin,
        spindot + steerdot * sin(inclination_rad),
        inclinationdot * sspin + steerdot_cinclination * cspin };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type inverse_angular_kinematic_relationships(const T &p, const T &q, const T &r,
                                                             const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
{
    //
    // Returns an array holding the time derivatives of a rotating
    // frame's Euler angles in "ZYX" sequence "[yawdot, pitchdot,
    // rolldot]" from its angular velocity components "[p, q, r]"
    // (projected onto the rotating frame) and its Euler angles
    // "[yaw_rad, pitch_rad, roll_rad]" (given in radians). These
    // output angle derivatives will share the units of the input
    // angular velocity components.
    //

    (void)yaw_rad;

    using std::cos;
    using std::sin;

    const auto croll = cos(roll_rad);
    const auto sroll = sin(roll_rad);
    const auto yawdot = (q * sroll + r * croll) / cos(pitch_rad);

    return Array3Type{
        yawdot,
        q * croll - r * sroll,
        p + yawdot * sin(pitch_rad) };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type inverse_sis_kinematic_relationships(const T &p, const T &q, const T &r,
                                                         const T &steer_rad, const T &inclination_rad, const T &spin_rad)
{
    //
    // Returns an array holding the time derivatives of a rotating
    // frame's Euler angles in "ZXY" sequence "[steerdot, inclinationdot,
    // spindot]" from its angular velocity components "[p, q, r]"
    // (projected onto the rotating frame) and its "[steer_rad,
    // inclination_rad, spin_rad]" angles (given in radians). These
    // output angle derivatives will share the units of the input
    // angular velocity components.
    //

    (void)steer_rad;

    using std::cos;
    using std::sin;

    const auto cspin = cos(spin_rad);
    const auto sspin = sin(spin_rad);
    const auto steerdot = (r * cspin - p * sspin) / cos(inclination_rad);

    return Array3Type{
        steerdot,
        p * cspin + r * sspin,
        q - steerdot * sin(inclination_rad) };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type angular_kinematic_relationships_derivatives(const T &yawdotdot, const T &pitchdotdot, const T &rolldotdot,
                                                                 const T &yawdot, const T &pitchdot, const T &rolldot,
                                                                 const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
{
    //
    // Returns an array holding the derivatives of the angular velocity
    // components "[pdot, qdot, rdot]" of a rotating frame (projected
    // onto this same frame), from the second & first time derivatives
    // of its Euler angles in "ZYX" sequence "[yawdotdot, pitchdotdot, rolldotdot]"
    // & "[yawdot, pitchdot, rolldot]" (whose units MUST be consistent) and
    // the values of these angles "[yaw_rad, pitch_rad, roll_rad]" (given in
    // radians). The resulting angular acceleration components will
    // share the units of the input Euler angle second derivatives.
    //

    (void)yaw_rad;

    using std::cos;
    using std::sin;

    const auto spitch = sin(pitch_rad);
    const auto cpitch = cos(pitch_rad);

    const auto sroll = sin(roll_rad);
    const auto croll = cos(roll_rad);

    const auto yawdot_spitch = yawdot * spitch;
    const auto yawdot_cpitch = yawdot * cpitch;

    const auto pitchdot_sroll = pitchdot * sroll;
    const auto pitchdot_croll = pitchdot * croll;

    const auto yawdotdot_cpitch = yawdotdot * cpitch;

    return Array3Type{
        rolldotdot - yawdotdot * spitch - pitchdot * yawdot_cpitch,
        pitchdotdot * croll - rolldot * (pitchdot_sroll - yawdot_cpitch * croll) + yawdotdot_cpitch * sroll - pitchdot_sroll * yawdot_spitch,
        yawdotdot_cpitch * croll - pitchdotdot * sroll - rolldot * (pitchdot_croll + yawdot_cpitch * sroll) - pitchdot_croll * yawdot_spitch };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type sis_kinematic_relationships_derivatives(const T &steerdotdot, const T &inclinationdotdot, const T &spindotdot,
                                                             const T &steerdot, const T &inclinationdot, const T &spindot,
                                                             const T &steer_rad, const T &inclination_rad, const T &spin_rad)
{
    //
    // Returns an array holding the derivatives of the angular velocity
    // components "[pdot, qdot, rdot]" of a rotating frame (projected
    // onto this same frame), from the second & first time derivatives
    // of its Euler angles in "ZXY" sequence "[steerdotdot, inclinationdotdot,
    // spindotdot]" & "[steerdot, inclinationdot, spindot]" (whose units MUST
    // be consistent) and the values of these angles "[steer_rad, inclination_rad,
    // spin_rad]" (given in radians). The resulting angular acceleration
    // components will share the units of the input Euler angle second derivatives.
    //

    (void)steer_rad;

    using std::cos;
    using std::sin;

    const auto sinclination = sin(inclination_rad);
    const auto cinclination = cos(inclination_rad);

    const auto sspin = sin(spin_rad);
    const auto cspin = cos(spin_rad);

    const auto steerdot_sinclination = steerdot * sinclination;
    const auto steerdot_cinclination = steerdot * cinclination;

    const auto inclinationdot_sspin = inclinationdot * sspin;
    const auto inclinationdot_cspin = inclinationdot * cspin;

    const auto steerdotdot_cinclination = steerdotdot * cinclination;

    return Array3Type{
        inclinationdotdot * cspin - spindot * (inclinationdot_sspin + steerdot_cinclination * cspin) - steerdotdot_cinclination * sspin + inclinationdot_sspin * steerdot_sinclination,
        spindotdot + steerdotdot * sinclination + inclinationdot * steerdot_cinclination,
        steerdotdot_cinclination * cspin + inclinationdotdot * sspin + spindot * (inclinationdot_cspin - steerdot_cinclination * sspin) - inclinationdot_cspin * steerdot_sinclination };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type inverse_angular_kinematic_relationships_derivatives(const T &pdot, const T &qdot, const T &rdot,
                                                                         const T &yawdot, const T &pitchdot, const T &rolldot,
                                                                         const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
{
    //
    // Returns an array holding the second time derivatives of a rotating
    // frame's Euler angles in "ZYX" sequence "[yawdotdot, pitchdotdot,
    // rolldotdot]" from its angular acceleration components "[pdot, qdot, rdot]"
    // (projected onto the rotating frame), the derivatives of these angles
    // "[yawdot, pitchdot, rolldot]" (whose units MUST be consistent with those
    // of the angular acceleration components), and the values of these angles
    // "[yaw_rad, pitch_rad, roll_rad]" (given in radians). These output angle
    // second derivatives will share the units of the input angular acceleration
    // components.
    //

    (void)yaw_rad;

    using std::cos;
    using std::sin;

    const auto croll = cos(roll_rad);
    const auto sroll = sin(roll_rad);

    const auto cpitch = cos(pitch_rad);
    const auto inv_cpitch = T{ 1 } / cpitch;
    const auto spitch = sin(pitch_rad);

    const auto yawdotdotA = (rolldot * pitchdot + rdot * croll + qdot * sroll) * inv_cpitch;
    const auto yawdotdotB = yawdot * pitchdot * inv_cpitch;

    return Array3Type{
        yawdotdotA + spitch * yawdotdotB,
        qdot * croll - rdot * sroll - yawdot * rolldot * cpitch,
        pdot + spitch * yawdotdotA + yawdotdotB };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type inverse_sis_kinematic_relationships_derivatives(const T &pdot, const T &qdot, const T &rdot,
                                                                     const T &steerdot, const T &inclinationdot, const T &spindot,
                                                                     const T &steer_rad, const T &inclination_rad, const T &spin_rad)
{
    //
    // Returns an array holding the second time derivatives of a rotating
    // frame's Euler angles in "ZXY" sequence "[steerdotdot, inclinationdotdot,
    // spindotdot]" from its angular acceleration components "[pdot, qdot, rdot]"
    // (projected onto the rotating frame), the derivatives of these angles
    // "[steerdot, inclinationdot, spindot]" (whose units MUST be consistent with
    // those of the angular acceleration components), and the values of these angles
    // "[steer_rad, inclination_rad, spin_rad]" (given in radians). These output angle
    // second derivatives will share the units of the input angular acceleration
    // components.
    //

    (void)steer_rad;

    using std::cos;
    using std::sin;

    const auto cspin = cos(spin_rad);
    const auto sspin = sin(spin_rad);

    const auto cinclination = cos(inclination_rad);
    const auto inv_cinclination = T{ 1 } / cinclination;
    const auto sinclination = sin(inclination_rad);

    const auto steerdotdotA = (inclinationdot * spindot - rdot * cspin + pdot * sspin) * inv_cinclination;
    const auto steerdotdotB =  steerdot * inclinationdot * inv_cinclination;

    return Array3Type{
        sinclination * steerdotdotB - steerdotdotA,
        pdot * cspin + rdot * sspin + spindot * steerdot * cinclination,
        qdot - steerdotdotB + sinclination * steerdotdotA };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type angular_kinematic_relationships_crossderivatives(const T &yawdot0dot1, const T &pitchdot0dot1, const T &rolldot0dot1,
                                                                      const T &yawdot0, const T &pitchdot0, const T &rolldot0,
                                                                      const T &yawdot1, const T &pitchdot1, const T &rolldot1,
                                                                      const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
{
    //
    // Like function "angular_kinematic_relationships_derivatives", but with
    // two distinct time increments ("0" and "1") in the Euler angles ("ZYX"
    // sequence) to form the second derivatives. Returns the angular
    // acceleration "[pdot01, qdot01, rdot01]" that these two increments yield.
    //

    (void)rolldot0;
    (void)yawdot1;
    (void)yaw_rad;

    using std::cos;
    using std::sin;

    const auto spitch = sin(pitch_rad);
    const auto cpitch = cos(pitch_rad);

    const auto sroll = sin(roll_rad);
    const auto croll = cos(roll_rad);

    const auto yawdot0_spitch = yawdot0 * spitch;
    const auto yawdot0_cpitch = yawdot0 * cpitch;

    const auto yawdot0dot1_cpitch = yawdot0dot1 * cpitch;

    return Array3Type{
        rolldot0dot1 - yawdot0dot1 * spitch - pitchdot1 * yawdot0_cpitch,
        pitchdot0dot1 * croll - rolldot1 * (pitchdot0 * sroll - yawdot0_cpitch * croll) + yawdot0dot1_cpitch * sroll - pitchdot1 * sroll * yawdot0_spitch,
        yawdot0dot1_cpitch * croll - pitchdot0dot1 * sroll - rolldot1 * (pitchdot0 * croll + yawdot0_cpitch * sroll) - pitchdot1 * croll * yawdot0_spitch };
}

template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type sis_kinematic_relationships_crossderivatives(const T &steerdot0dot1, const T &inclinationdot0dot1, const T &spindot0dot1,
                                                                  const T &steerdot0, const T &inclinationdot0, const T &spindot0,
                                                                  const T &steerdot1, const T &inclinationdot1, const T &spindot1,
                                                                  const T &steer_rad, const T &inclination_rad, const T &spin_rad)
{
    //
    // Like function "sis_kinematic_relationships_derivatives", but
    // with two distinct time increments ("0" and "1") in the
    // "steer-inclination-spin" angles ("ZXY" sequence) to form the
    // second derivatives. Returns the angular acceleration
    // "[pdot01, qdot01, rdot01]" that these two increments yield.
    //

    (void)spindot0;
    (void)steerdot1;
    (void)steer_rad;

    using std::cos;
    using std::sin;

    const auto sinclination = sin(inclination_rad);
    const auto cinclination = cos(inclination_rad);

    const auto sspin = sin(spin_rad);
    const auto cspin = cos(spin_rad);

    const auto steerdot0_sinclination = steerdot0 * sinclination;
    const auto steerdot0_cinclination = steerdot0 * cinclination;

    const auto steerdot0dot1_cinclination = steerdot0dot1 * cinclination;

    return Array3Type{
        inclinationdot0dot1 * cspin - spindot1 * (inclinationdot0 * sspin + steerdot0_cinclination * cspin) - steerdot0dot1_cinclination * sspin + inclinationdot1 * sspin * steerdot0_sinclination,
        spindot0dot1 + steerdot0dot1 * sinclination + inclinationdot1 * steerdot0_cinclination,
        steerdot0dot1_cinclination * cspin + inclinationdot0dot1 * sspin + spindot1 * (inclinationdot0 * cspin - steerdot0_cinclination * sspin) - inclinationdot1 * cspin * steerdot0_sinclination };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type inverse_angular_kinematic_relationships_crossderivatives(const T &pdot01, const T &qdot01, const T &rdot01,
                                                                              const T &yawdot0, const T &pitchdot0, const T &rolldot0,
                                                                              const T &yawdot1, const T &pitchdot1, const T &rolldot1,
                                                                              const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
{
    //
    // Like function "inverse_angular_kinematic_relationships_derivatives", but with
    // two distinct time increments ("0" and "1") in the Euler angles ("ZYX"
    // sequence) to form the angular acceleration. Returns the angular second
    // derivatives "[yawdot0dot1, pitchdot0dot1, rolldot0dot1]" that these two increments
    // yield.
    //

    (void)rolldot0;
    (void)yawdot1;
    (void)yaw_rad;

    using std::cos;
    using std::sin;

    const auto croll = cos(roll_rad);
    const auto sroll = sin(roll_rad);

    const auto cpitch = cos(pitch_rad);
    const auto inv_cpitch = T{ 1 } / cpitch;
    const auto spitch = sin(pitch_rad);

    const auto yawdotdotA = (rolldot1 * pitchdot0 + rdot01 * croll + qdot01 * sroll) * inv_cpitch;
    const auto yawdotdotB = yawdot0 * pitchdot1 * inv_cpitch;

    return Array3Type{
        yawdotdotA + spitch * yawdotdotB,
        qdot01 * croll - rdot01 * sroll - yawdot0 * rolldot1 * cpitch,
        pdot01 + spitch * yawdotdotA + yawdotdotB };
}


template<typename T,
         typename Array3Type = std::array<T, 3u> >
constexpr Array3Type inverse_sis_kinematic_relationships_crossderivatives(const T &pdot01, const T &qdot01, const T &rdot01,
                                                                          const T &steerdot0, const T &inclinationdot0, const T &spindot0,
                                                                          const T &steerdot1, const T &inclinationdot1, const T &spindot1,
                                                                          const T &steer_rad, const T &inclination_rad, const T &spin_rad)
{
    //
    // Like function "inverse_sis_kinematic_relationships_derivatives", but with
    // two distinct time increments ("0" and "1") in the "steer-inclination-spin"
    // angles ("ZXY" sequence) to form the angular acceleration. Returns the
    // angular second derivatives "[steerdot0dot1, inclinationdot0dot1, spindot0dot1]"
    // that these two increments yield.
    //

    (void)spindot0;
    (void)steerdot1;
    (void)steer_rad;

    using std::cos;
    using std::sin;

    const auto cspin = cos(spin_rad);
    const auto sspin = sin(spin_rad);

    const auto cinclination = cos(inclination_rad);
    const auto inv_cinclination = T{ 1 } / cinclination;
    const auto sinclination = sin(inclination_rad);

    const auto steerdotdotA = (inclinationdot0 * spindot1 - rdot01 * cspin + pdot01 * sspin) * inv_cinclination;
    const auto steerdotdotB =  steerdot0 * inclinationdot1 * inv_cinclination;

    return Array3Type{
        sinclination * steerdotdotB - steerdotdotA,
        pdot01 * cspin + rdot01 * sspin + spindot1 * steerdot0 * cinclination,
        qdot01 - steerdotdotB + sinclination * steerdotdotA };
}

#endif
