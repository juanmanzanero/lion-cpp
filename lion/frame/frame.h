#ifndef FRAME_H
#define FRAME_H

#include "lion/math/vector3d.hpp"
#include "lion/math/matrix3x3.h"
#include "lion/foundation/type_traits.h"
#include "rotations.h"

namespace lioncpp {
namespace detail {

    template<typename T>
    class Inertial_parent_frame 
    {
        using basic_type = T;
        using aggregated_type = T;
    };

    template<typename P>
    constexpr const size_t get_frame_generation()
    {
     if constexpr (!std::is_same_v<P,Inertial_parent_frame>)
     {
         return P::generation+1;
     }
     else
     {
         return 0;
     }
    }
    
    /** 
     *  This function takes two frames \p f1 and \p f2, and computes the generation
     *  of the earliest common parent they share. 
     *  Throws an exception if they do not share any (belong to different families)
     *  <pre>
     *    gen 0.............0......0*
     *                      |           For 4 and 8 would be gen=1 (frame 1)
     *    gen 1.............1           For 6 and 7 would be gen=3 (frame 4)
     *                    /   \         For 8 and 9 would be gen=3 (frame 5)
     *    gen 2..........2.....3        For 2 and 5 would be gen=5 (frame 1)
     *                   |     |        For 6 and 1 would be gen=1 (frame 1)
     *    gen 3..........4.....5        For 4 and 0* it will throw an exception
     *                  / \   / \
     *    gen 4........6...7.8...9
     *  </pre>
     * 
     *  @param[in] f1: first frame
     *  @param[in] f2: second frame
     *  @param[out] number of the first shared generation.
     */
    template<typename Input_frame_1, typename Input_frame_2>
    constexpr size_t get_crossing_generation(const Input_frame_1& f1, const Input_frame_2& f2);
}
}


//!  Class defining moving frames
template<typename T, typename Parent_frame_type = lioncpp::detail::Inertial_parent_frame<T>, size_t Number_of_rotations = 0>
class Frame
{
 public:
    using basic_type                         = T;
    using parent_frame_type                  = Parent_frame_type;
    using number_of_rotations                = Number_of_rotations;
    constexpr const static size_t generation = lioncpp::detail::get_frame_generation<Parent_frame_type>();
    using aggregated_type                    = combine_types<T, Parent_frame_type::aggregated_type>::type;
    constexpr const static bool is_inertial  = (generation == 0);

    //!  Default constructor. Returns an inertial frame
    Frame() = default;

    //!  Constructor for non-inertial frames
    //!
    //! @param[in] x: Origin coordinates in \p parent Frame [m]
    //! @param[in] dx: Origin velocity in \p parent Frame [m/s]
    //! @param[in] angles: Vector of rotation angles [rad]
    //! @param[in] dangles: Vector of rotation angles derivatives [rad/s]
    //! @param[in] axis: Vector of axis of rotations: X, Y, or Z
    //! @param[in] parent: Frame for which all the quantities above are referred to
    Frame(const Vector3d<T>& x, 
          const Vector3d<T>& dx, 
          const std::array<T,Number_of_rotations>& angles, 
          const std::array<T,Number_of_rotations>& dangles,
          const std::array<Axis,Number_of_rotations>& axis, 
          const Parent_frame_type& parent
         );

    //!  Updates all internal data members: rotation matrices and angular velocity vectors
    void update();

    // Getters
    
    //!  Returns a const reference to the parent Frame
    std::enable_if_t<!is_inertial, const Parent_frame_type&> get_parent() const { return *_parent; }

    //!  Returns a const pointer to the parent Frame
    constexpr const Parent_frame_type* get_parent_ptr() const { return _parent; }

    //! Sets the parent pointer to a new parent
    //! @param[in] parent: the new parent
    constexpr void set_parent(const Parent_frame_type& parent) { _parent = &parent; } 

    //!  Returns a const reference to the origin coordinates (in \p parent Frame)
    constexpr const Vector3d<T>& get_origin() const { return _x; }

    //!  Returns a const reference to the origin velocity (in \p parent Frame)
    constexpr const Vector3d<T>& get_relative_velocity() const { return _dx; } 

    //!  Returns a const reference to the rotation angles vector
    constexpr const std::array<T,number_of_rotations>& get_rotation_angles() const { return _angles; } 

    //!  Returns a const reference to the rotation angles derivatives
    constexpr const std::array<T,number_of_rotations>& get_rotation_angles_derivative() const { return _dangles; } 

    //!  Returns a const reference to the rotation axis  
    constexpr const std::array<Axis,number_of_rotations>& get_rotation_axis() const { return _axis; } 

    //!  Whether a frame has been updated or not. When a frame is modified, it is updated
    //! automatically, unless is is told otherwise.
    //!  @param[out] is_updated: true if updated, false otherwise
    constexpr const bool& is_updated() const { return _updated; }

    //!  Returns the absolute position of a vector \p x of this frame
    //! Applies x = xO + Qpc*x' over all the chain of parent frames until the inertial frame
    //! is reached
    //! @param[in] x: the vector of this frame to be computed. Defaults to the origin (0,0,0)
    template<typename U>
    auto get_absolute_position(const Vector3d<U>& x = Vector3d<U>(0.0)) const
        -> Vector3d<typename combine_types<aggregated_type,U>::type>;

    //!  Returns the absolute velocity of a point \p x with velocity \p dx of this frame
    //! AND projects it back to this frame
    //! Applies v = v0 + Qpc(omega^x' + v') over all the chain of parent frames until the 
    //! inertial frame is reached. Then it is projected back to this frame (aka "body")
    //! xBody = Q(Body|Inertial)xInertial
    //! @param[in] x: position of the point of this frame. Defaults to the origin (0,0,0)
    //! @param[in] dx: local velocity of the point of this frame. Defaults to fixed point (0,0,0)
    template<typename U>
    auto get_absolute_velocity_in_body(const Vector3d<U>& x = Vector3d<U>(0.0),
                                       const Vector3d<U>& dx = Vector3d<U>(0.0)
                                      ) const
        -> Vector3d<typename combine_types<aggregated_type, U>::type>;

    //!  Returns the absolute velocity of a point \p x with velocity \p dx of this frame
    //! AND projects to this frame's parent.
    //! Applies v = v0 + Qpc(omega^x' + v') over all the chain of parent frames until the 
    //! inertial frame is reached. Then it is projected back to this frame's parent
    //! xParent = Q(Parent|Inertial)xInertial
    //! @param[in] x: position of the point of this frame. Defaults to the origin (0,0,0)
    //! @param[in] dx: local velocity of the point of this frame. Defaults to fixed point (0,0,0)
    template<typename U>
    auto get_absolute_velocity_in_parent(const Vector3d<U>& x = Vector3d<U>(0.0),
                                         const Vector3d<U>& dx = Vector3d<U>(0.0)
                                        ) const
        -> Vector3d<typename combine_types<aggregated_type, U>::type>;
    
    //!  Returns the absolute velocity of a point \p x with velocity \p dx of this frame
    //! projected in inertial frame.
    //! Applies v = v0 + Qpc(omega^x' + v') over all the chain of parent frames until the 
    //! inertial frame is reached.
    //! @param[in] x: position of the point of this frame. Defaults to the origin (0,0,0)
    //! @param[in] dx: local velocity of the point of this frame. Defaults to fixed point (0,0,0)
    template<typename U>
    auto get_absolute_velocity_in_inertial(const Vector3d<U>& x = Vector3d<U>(0.0),
                                      const Vector3d<U>& dx = Vector3d<U>(0.0)
                                     ) const
        -> Vector3d<typename combine_types<aggregated_type, U>::type>;

    //!  Returns the position and velocity of a point \p x with velocity \p dx relative to another
    //! moving frame \p target.
    //! Applies x = x0 + Qpc x' and v = v0 + Qpc(omega^x' + v') over all the chain of parent frames
    //! until the target is reached.
    //! If the target belongs to a generation younger than the body (i.e. the target is a children
    //! of the body), the formulas are:
    //! x = Qcp(-x0 + x')
    //! v = Qcp(-v0 + v' + omega^(-x0+x'))
    //! The origin needs to be transformed to the child coordinates because it is written on the parent
    //! coordinates, and it is introduced in the coriolis term for the same reason.
    //! @param[in] x: position of the point of this frame. Defaults to the origin (0,0,0)
    //! @param[in] dx: local velocity of the point of this frame. Defaults to fixed point (0,0,0)
    //! @param[out] returns the tuple (x,v). One can use std::tie(x,v) = f() to unwrap the variables.
    //!                 or std::get<0|1>f() to only get one of the two.
    template<typename Input_frame>
    auto get_position_and_velocity_in_target(const Input_frame& target, 
                                        const Vector3d<T>& x = Vector3d<T>(0.0), 
                                        const Vector3d<T>& dx = Vector3d<T>(0.0)
                                       ) const
        -> std::pair<Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>, 
                     Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>>;

    //! Returns the rotation matrix from the children to the parent
    //! This matrix is such that x|Parent = Q(Parent|Body)x|Body
    constexpr const Matrix3x3<T>& get_rotation_matrix() const { return _Qpc; }

    //! Returns the rotation matrix from the parent to the children
    //! This matrix is such that x|Body = Q(Body|Parent)x|Parent
    constexpr const Matrix3x3<T>& get_back_rotation_matrix() const { return _Qcp; }

    //! Returns the angular velocity of the frame wrt its parent written in this frame
    constexpr const Vector3d<T>& get_omega_wrt_parent_in_body() const { return _omega_pc_self; }

    //! Returns the angular velocity of the frame wrt its parent written in the parent frame
    constexpr const Vector3d<T>& get_omega_wrt_parent_in_parent() const { return _omega_pc_parent; }

    //! Returns the rotation matrix of the frame wrt a target frame
    //! x|Target = Q(target|body)x|body
    //! @param[in] target: target Frame
    template<typename Input_frame>
    Matrix3x3<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type> get_rotation_matrix(const Input_frame& target) const;

    //! Returns the angular velocity of the frame wrt a target frame 
    //! written in this frame (aka body)
    //! @param[in] target: target Frame
    template<typename Input_frame>
    auto get_omega_in_body(const Input_frame& target) const
        ->Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>;

    //! Returns the angular velocity of the frame wrt a target frame
    //! written in this frame's parent 
    //! @param[in] target: target Frame
    template<typename Input_frame>
    auto get_omega_in_parent(const Input_frame& target) const
        ->Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>
    { return _Qpc * get_omega_in_body(target); }

    //! Returns the angular velocity of the frame wrt a target frame
    //! written in the target frame
    //! @param[in] target: target Frame
    template<typename Input_frame>
    Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type> get_omega_in_target(const Input_frame& target) const;

    //! Computes the rotation matrix from this frame to the inertial frame
    //! x|Inertial = Q(Inertial|Body)x|Body
    Matrix3x3<aggregated_type> get_absolute_rotation_matrix() const;

    //! Computes the total (i.e. wrt the inertial frame) angular velocity of 
    //! the frame, written in this frame
    Vector3d<aggregated_type> get_omega_absolute_in_body() const;

    //! Computes the total (i.e. wrt the inertial frame) angular velocity
    //! of the frame, written in this frame's parent
    Vector3d<aggregated_type> get_omega_absolute_in_parent() const;

    //! Computes the total (i.e. wrt the inertial frame) angular velocity
    //! of the frame, written in the inertial frame
    Vector3d<aggregated_type> get_omega_absolute_in_inertial() const;

    //! Set the origin of the frame
    //! @param[in] x: new origin [m]
    //! @param[in] bUpdate: if true, will call update()
    constexpr void set_origin(const Vector3d<T>& x, const bool bUpdate = true);

    //! Set the origin and velocity of the frame
    //! @param[in] x: new origin [m]
    //! @param[in] dx: new velocity [m/s]
    //! @param[in] bUpdate: if true, will call update()
    constexpr void set_origin(const Vector3d<T>& x, const Vector3d<T>& dx, const bool bUpdate = true);

    //! Set the velocity of the frame
    //! @param[in] dx: new velocity [m/s]
    //! @param[in] bUpdate: if true, will call update()
    constexpr void set_velocity(const Vector3d<T>& dx, const bool bUpdate = true);
  
    //! Set a rotation angle
    //! @param[in] which: index of the angle
    //! @param[in] angle: new angle [rad]
    //! @param[in] bUpdate: if true, will call update()
    void set_rotation_angle(const size_t which, const T angle, const bool bUpdate = true);

    //! Set a rotation angle and speed
    //! @param[in] which: index of the angle
    //! @param[in] angle: new angle [rad]
    //! @param[in] dangle: new angular velocity [rad/s]
    //! @param[in] bUpdate: if true, will call update()
    void set_rotation_angle(const size_t which, const T angle, const T dangle, const bool bUpdate = true);

    //! Set a rotation angular speed
    //! @param[in] which: index of the angle
    //! @param[in] dangle: new angular velocity [rad/s]
    //! @param[in] bUpdate: if true, will call update()
    void set_angular_speed(const size_t which, const T dangle, const bool bUpdate = true);

 private:

     template<typename U>
     struct Frame_kinematics
     {
         Vector3d<U>  x;
         Vector3d<U>  dx;
         Matrix3x3<U> Qpc;
         Matrix3x3<U> Qcp;
         Vector3d<U>  omega_pc_self;
         Vector3d<U>  omega_pc_parent;
     };


    template<typename return_type = aggregated_type>
    auto get_frame_kinematics_at_generation(size_t requested_generation) const->Frame_kinematics<return_type>;
   
    //! Settable parameters
    const Parent_frame_type* _parent = nullptr; //! Parent frame
    Vector3d<T> _x;                   //! Origin coordinates in parent frame
    Vector3d<T> _dx;                  //! Origin velocity in parent frame

    std::array<T,Number_of_rotations>    _angles;  //! Rotation angles
    std::array<T,Number_of_rotations>    _dangles; //! Rotation angles derivatives
    std::array<Axis,Number_of_rotations> _axis;    //! Rotation axes

    //! Updated parameters
    bool _updated = true;                 //! If the frame is updated
    Matrix3x3<T> _Qpc = Matrix3x3<T>::eye();  //! Rotation matrix (from parent Xp=QpcXc)
    Matrix3x3<T> _Qcp = Matrix3x3<T>::eye();  //! Rotation matrix (from children Xc=QcpXp)
    Vector3d<T>  _omega_pc_self = Vector3d<T>::zeros();    //! Omega from parent (in this frame)
    Vector3d<T>  _omega_pc_parent = Vector3d<T>::zeros();  //! Omega from parent (in parent)
};


#include "frame.hpp"
#endif
