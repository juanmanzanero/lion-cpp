#ifndef __FRAME_H__
#define __FRAME_H__

#include "lion/math/vector3d.h"
#include "lion/math/matrix3x3.h"
#include "rotations.h"

//!  Class defining moving frames
class Frame
{
 public:
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
    Frame(const tVector3d& x, 
          const tVector3d& dx, 
          const std::vector<timeseries> angles, 
          const std::vector<timeseries> dangles,
          const std::vector<Axis> axis, 
          const Frame& parent
         );

    //!  Updates all internal data members: rotation matrices and angular velocity vectors
    void update();

    // Getters
    
    //!  Returns a const reference to the parent Frame
    const Frame& get_parent() const;

    //!  Returns a const pointer to the parent Frame
    constexpr const Frame* get_parent_ptr() const { return _parent; }

    //!  Computes the generation of the frame
    constexpr size_t generation() const ;

    //! Sets the parent pointer to a new parent
    //! @param[in] parent: the new parent
    constexpr void set_parent(const Frame& parent) { _parent = &parent; } 

    //!  Returns a const reference to the origin coordinates (in \p parent Frame)
    constexpr const tVector3d& get_origin() const { return _x; }

    //!  Returns a const reference to the origin velocity (in \p parent Frame)
    constexpr const tVector3d& get_relative_velocity() const { return _dx; } 

    //!  Returns a const reference to the rotation angles vector
    constexpr const std::vector<timeseries>& get_rotation_angles() const { return _angles; } 

    //!  Returns a const reference to the rotation angles derivatives
    constexpr const std::vector<timeseries>& get_rotation_angles_derivative() const { return _dangles; } 

    //!  Returns a const reference to the rotation axis  
    constexpr const std::vector<Axis>& get_rotation_axis() const { return _axis; } 

    //!  Whether a frame is inertial or not
    //!  @param[out] is_inertial: true if inertial, false otherwise. 
    constexpr bool is_inertial() const { return generation() == 0; } 
   
    //!  Whether a frame has been updated or not. When a frame is modified, it is updated
    //! automatically, unless is is told otherwise.
    //!  @param[out] is_updated: true if updated, false otherwise
    constexpr const bool& is_updated() const { return _updated; }

    //!  Returns the absolute position of a vector \p x of this frame
    //! Applies x = xO + Qpc*x' over all the chain of parent frames until the inertial frame
    //! is reached
    //! @param[in] x: the vector of this frame to be computed. Defaults to the origin (0,0,0)
    constexpr tVector3d get_absolute_position(const tVector3d& x = tVector3d(0.0)) const;

    //!  Returns the absolute velocity of a point \p x with velocity \p dx of this frame
    //! AND projects it back to this frame
    //! Applies v = v0 + Qpc(omega^x' + v') over all the chain of parent frames until the 
    //! inertial frame is reached. Then it is projected back to this frame (aka "body")
    //! xBody = Q(Body|Inertial)xInertial
    //! @param[in] x: position of the point of this frame. Defaults to the origin (0,0,0)
    //! @param[in] dx: local velocity of the point of this frame. Defaults to fixed point (0,0,0)
    constexpr tVector3d get_absolute_velocity_in_body(const tVector3d& x = tVector3d(0.0),
                                                      const tVector3d& dx = tVector3d(0.0) 
                                                     ) const;

    //!  Returns the absolute velocity of a point \p x with velocity \p dx of this frame
    //! AND projects to this frame's parent.
    //! Applies v = v0 + Qpc(omega^x' + v') over all the chain of parent frames until the 
    //! inertial frame is reached. Then it is projected back to this frame's parent
    //! xParent = Q(Parent|Inertial)xInertial
    //! @param[in] x: position of the point of this frame. Defaults to the origin (0,0,0)
    //! @param[in] dx: local velocity of the point of this frame. Defaults to fixed point (0,0,0)
    constexpr tVector3d get_absolute_velocity_in_parent(const tVector3d& x = tVector3d(0.0),
                                                        const tVector3d& dx = tVector3d(0.0)
                                                       ) const;
    
    //!  Returns the absolute velocity of a point \p x with velocity \p dx of this frame
    //! projected in inertial frame.
    //! Applies v = v0 + Qpc(omega^x' + v') over all the chain of parent frames until the 
    //! inertial frame is reached.
    //! @param[in] x: position of the point of this frame. Defaults to the origin (0,0,0)
    //! @param[in] dx: local velocity of the point of this frame. Defaults to fixed point (0,0,0)
    constexpr tVector3d get_absolute_velocity_in_inertial(const tVector3d& x = tVector3d(0.0),
                                                          const tVector3d& dx = tVector3d(0.0)
                                                         ) const;

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
    std::pair<tVector3d,tVector3d> get_position_and_velocity_in_target(const Frame& target, 
                                                                       const tVector3d& x = tVector3d(0.0), 
                                                                       const tVector3d& dx = tVector3d(0.0)
                                                                      ) const;

    //! Returns the rotation matrix from the children to the parent
    //! This matrix is such that x|Parent = Q(Parent|Body)x|Body
    constexpr const tMatrix3x3& get_rotation_matrix() const { return _Qpc; }

    //! Returns the rotation matrix from the parent to the children
    //! This matrix is such that x|Body = Q(Body|Parent)x|Parent
    constexpr const tMatrix3x3& get_back_rotation_matrix() const { return _Qcp; }

    //! Returns the angular velocity of the frame wrt its parent written in this frame
    constexpr const tVector3d& get_omega_wrt_parent_in_body() const { return _omega_pc_self; }

    //! Returns the angular velocity of the frame wrt its parent written in the parent frame
    constexpr const tVector3d& get_omega_wrt_parent_in_parent() const { return _omega_pc_parent; }

    //! Returns the rotation matrix of the frame wrt a target frame
    //! x|Target = Q(target|body)x|body
    //! @param[in] target: target Frame
    tMatrix3x3 get_rotation_matrix(const Frame& target) const;

    //! Returns the angular velocity of the frame wrt a target frame 
    //! written in this frame (aka body)
    //! @param[in] target: target Frame
    tVector3d get_omega_in_body(const Frame& target) const;

    //! Returns the angular velocity of the frame wrt a target frame
    //! written in this frame's parent 
    //! @param[in] target: target Frame
    tVector3d get_omega_in_parent(const Frame& target) const;

    //! Returns the angular velocity of the frame wrt a target frame
    //! written in the target frame
    //! @param[in] target: target Frame
    tVector3d get_omega_in_target(const Frame& target) const;

    //! Computes the rotation matrix from this frame to the inertial frame
    //! x|Inertial = Q(Inertial|Body)x|Body
    constexpr tMatrix3x3 get_absolute_rotation_matrix() const;

    //! Computes the total (i.e. wrt the inertial frame) angular velocity of 
    //! the frame, written in this frame
    constexpr tVector3d get_omega_absolute_in_body() const;

    //! Computes the total (i.e. wrt the inertial frame) angular velocity
    //! of the frame, written in this frame's parent
    constexpr tVector3d get_omega_absolute_in_parent() const;

    //! Computes the total (i.e. wrt the inertial frame) angular velocity
    //! of the frame, written in the inertial frame
    constexpr tVector3d get_omega_absolute_in_inertial() const;

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
    static size_t get_crossing_generation(const Frame& f1, const Frame& f2);

    //! Set the origin of the frame
    //! @param[in] x: new origin [m]
    //! @param[in] bUpdate: if true, will call update()
    constexpr void set_origin(const tVector3d& x, const bool bUpdate = true);

    //! Set the origin and velocity of the frame
    //! @param[in] x: new origin [m]
    //! @param[in] dx: new velocity [m/s]
    //! @param[in] bUpdate: if true, will call update()
    constexpr void set_origin(const tVector3d& x, const tVector3d& dx, const bool bUpdate = true);

    //! Set the velocity of the frame
    //! @param[in] dx: new velocity [m/s]
    //! @param[in] bUpdate: if true, will call update()
    constexpr void set_velocity(const tVector3d& dx, const bool bUpdate = true);
  
    //! Set a rotation angle
    //! @param[in] which: index of the angle
    //! @param[in] angle: new angle [rad]
    //! @param[in] bUpdate: if true, will call update()
    void set_rotation_angle(const size_t which, const timeseries angle, const bool bUpdate = true);

    //! Set a rotation angle and speed
    //! @param[in] which: index of the angle
    //! @param[in] angle: new angle [rad]
    //! @param[in] dangle: new angular velocity [rad/s]
    //! @param[in] bUpdate: if true, will call update()
    void set_rotation_angle(const size_t which, const timeseries angle, const timeseries dangle, const bool bUpdate = true);

    //! Set a rotation angular speed
    //! @param[in] which: index of the angle
    //! @param[in] dangle: new angular velocity [rad/s]
    //! @param[in] bUpdate: if true, will call update()
    void set_angular_speed(const size_t which, const timeseries dangle, const bool bUpdate = true);

    //! Append a new rotation angular
    //! @param[in] angle: new angle [rad]
    //! @param[in] dangle: new angular velocity [rad/s]
    //! @param[in] axis: which axis (X/Y/Z)
    //! @param[in] bUpdate: if true, will call update()
    void add_rotation(const timeseries angle, const timeseries dangle, const Axis axis, const bool bUpdate = true);

 private:
    
    //! Settable parameters
    const Frame* _parent = nullptr; //! Parent frame
    tVector3d _x;                   //! Origin coordinates in parent frame
    tVector3d _dx;                  //! Origin velocity in parent frame

    std::vector<timeseries> _angles;  //! Rotation angles
    std::vector<timeseries> _dangles; //! Rotation angles derivatives
    std::vector<Axis> _axis;    //! Rotation axes

    //! Updated parameters
    bool _updated = true;                 //! If the frame is updated
    tMatrix3x3 _Qpc = tMatrix3x3::eye();  //! Rotation matrix (from parent Xp=QpcXc)
    tMatrix3x3 _Qcp = tMatrix3x3::eye();  //! Rotation matrix (from children Xc=QcpXp)
    tVector3d  _omega_pc_self = tVector3d::zeros();    //! Omega from parent (in this frame)
    tVector3d  _omega_pc_parent = tVector3d::zeros();  //! Omega from parent (in parent)
};


#include "frame.hpp"
#endif
