#ifndef FRAME_HPP
#define FRAME_HPP

#include <tuple>

template<typename Input_frame_1, typename Input_frame_2>
constexpr size_t lioncpp::detail::get_crossing_generation(const Input_frame_1& f1, const Input_frame_2& f2)
{
    constexpr const size_t first_common_generation = std::min(f1.generation, f2.generation);  
    
    if constexpr (f1.generation == 0 && f2.generation == 0)
    {
        if constexpr (std::is_same_v<Input_frame_1, Input_frame_2>)
        {
            if (&f1 == &f2)
            {
                return 0;
            }
            else
            {
                throw lion_exception("[ERROR] get_crossing_generation -> f1 and f2 do not share a common ancestor");  
            }
        }
        else
        {
            throw lion_exception("[ERROR] get_crossing_generation -> f1 and f2 do not share a common ancestor");  
        }
    }
    else if constexpr (f1.generation > first_common_generation)
    {
        return get_crossing_generation(f1.get_parent(), f2);
    }
    else if constexpr (f2.generation > first_common_generation)
    {
        return get_crossing_generation(f1, f2.get_parent());
    }
    else
    {
        if constexpr (std::is_same_v<Input_frame_1,Input_frame_2>)
        {
            if (&f1 == &f2)
                return first_common_generation;
            else
                return get_crossing_generation(f1.get_parent(), f2.get_parent());
        }
        else
            return get_crossing_generation(f1.get_parent(), f2.get_parent());
    }
}

template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline Frame<T,Parent_frame_type,Number_of_rotations>::Frame(const Vector3d<T>& x, 
    const Vector3d<T>& dx, const std::array<T,Number_of_rotations>& angles, 
    const std::array<T,Number_of_rotations>& dangles, const std::array<Axis,Number_of_rotations>& axis, 
    const Parent_frame_type& parent) 
: _parent(&parent),
  _x(x), 
  _dx(dx), 
  _angles(angles), 
  _dangles(dangles), 
  _axis(axis),
  _updated(false)
{
    for (size_t i = 0; i < _axis.size(); ++i)
    {
        switch(_axis[i])
        {
         case X: break;
         case Y: break;
         case Z: break;

         default:
            throw lion_exception("Invalid rotation axis: "+std::to_string(_axis[i])+", in position "+std::to_string(i)); 
         }
    }

    update();
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline void Frame<T,Parent_frame_type,Number_of_rotations>::update()
{
    if constexpr (is_inertial)
        return;

    if ( _angles.size() == 0 )
    {
        _Qpc = Matrix3x3<T>::eye();
        _Qcp = Matrix3x3<T>::eye();
        _omega_pc_self = Vector3d<T>::zeros();
        _omega_pc_parent = Vector3d<T>::zeros();

        _updated = true;
        return;
    }

//  Compute the rotation matrix
    std::array<Matrix3x3<T>,Number_of_rotations> fwd_rot_matrices;
    std::array<Matrix3x3<T>,Number_of_rotations> bwd_rot_matrices;

//  Accumulated forward rotation matrices: X_parent = Q_{12}Q_{23}...Q_{j-1,j}Xj
    std::array<Matrix3x3<T>,Number_of_rotations> accumulated_fwd_rot_matrices;

//  Accumulated backward rotation matrices: X_child = Q_{N,N-1}Q_{N-1,N-2}...Q_{j+1,j}Xj
    std::array<Matrix3x3<T>,Number_of_rotations> accumulated_bwd_rot_matrices;

    for (size_t i = 0; i < Number_of_rotations; ++i)
    {
        switch(_axis[i])
        {
         case(X):
            fwd_rot_matrices[i] = rotation_matrix_x(_angles[i]);
            break;

         case(Y):
            fwd_rot_matrices[i] = rotation_matrix_y(_angles[i]);
            break;

         case(Z):
            fwd_rot_matrices[i] = rotation_matrix_z(_angles[i]);
            break;
        } 

        bwd_rot_matrices[i] = transpose(fwd_rot_matrices.at(i));
    }

    accumulated_fwd_rot_matrices.front() = fwd_rot_matrices.front();

    for (size_t i=1; i< Number_of_rotations; ++i)
        accumulated_fwd_rot_matrices[i] = accumulated_fwd_rot_matrices[i-1]*fwd_rot_matrices[i];

    accumulated_bwd_rot_matrices.back() = bwd_rot_matrices.back();
    for (int i = Number_of_rotations-2; i >= 0; --i)
        accumulated_bwd_rot_matrices[i] = accumulated_bwd_rot_matrices[i+1]*bwd_rot_matrices[i];

    _Qpc = accumulated_fwd_rot_matrices.back();
    _Qcp = accumulated_bwd_rot_matrices.front();

    // Compute omega
    _omega_pc_self = Vector3d<T>::zeros();
    _omega_pc_parent = Vector3d<T>::zeros();

    for (size_t i = 0; i < Number_of_rotations; ++i)
    {
        Vector3d<T> new_rot(Vector3d<T>::zeros());
        new_rot[_axis[i]] = _dangles[i];
        _omega_pc_self += accumulated_bwd_rot_matrices[i]*new_rot;
        _omega_pc_parent += accumulated_fwd_rot_matrices[i]*new_rot;
    }

    _updated = true;

    return;
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename U>
inline auto Frame<T, Parent_frame_type, Number_of_rotations>::get_absolute_position(const Vector3d<U>& x) const
    -> Vector3d<typename combine_types<aggregated_type,U>::type>
{
    if constexpr (is_inertial)
        return x;
    else
        return get_parent().get_absolute_position(_x + _Qpc * x);
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename U>
inline auto Frame<T, Parent_frame_type, Number_of_rotations>::get_absolute_velocity_in_body(const Vector3d<U>& x, const Vector3d<U>& dx) const
    -> Vector3d<typename combine_types<aggregated_type,U>::type>
{
    if constexpr (is_inertial)
        return dx;
    else
        return _Qcp*get_parent().get_absolute_velocity_in_body(_x + _Qpc*x, _dx + _Qpc*(dx+cross(_omega_pc_self,x)));
}

    
template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename U>
inline auto Frame<T, Parent_frame_type, Number_of_rotations>::get_absolute_velocity_in_parent
    (const Vector3d<U>& x, const Vector3d<U>& dx) const
    -> Vector3d<typename combine_types<aggregated_type, U>::type>
{
    if constexpr (is_inertial)
        return dx;
    else
        return _Qpc*get_absolute_velocity_in_body(x, dx);
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename U>
inline auto Frame<T, Parent_frame_type, Number_of_rotations>::get_absolute_velocity_in_inertial
    (const Vector3d<U>& x, const Vector3d<U>& dx) const
    -> Vector3d<typename combine_types<aggregated_type, U>::type>
{
    if constexpr (is_inertial)
        return dx;
    else
        return get_parent().get_absolute_velocity_in_inertial(_x + _Qpc*x, _dx + _Qpc*(dx+cross(_omega_pc_self,x)));
}
 
template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename Input_frame>
inline auto Frame<T, Parent_frame_type, Number_of_rotations>::get_position_and_velocity_in_target
    (const Input_frame& target, const Vector3d<T>& x, const Vector3d<T>& dx) const
    -> std::pair<Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>, 
                 Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>> 
{
    using return_type = typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type;

    // Target is self
    if constexpr (std::is_same_v<Frame<T, Parent_frame_type, Number_of_rotations>, Input_frame>)
        if ( this == &target ) 
            return {x,dx};

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation == common_generation ) 
    {
        Vector3d<return_type> current_velocity = dx;
        Vector3d<return_type> current_position = x;

        for (size_t i_generation = common_generation + 1; i_generation <= target.generation; ++i_generation)
        {
            auto frame_kinematics = target.template get_frame_kinematics_at_generation<return_type>(i_generation);

            current_velocity = frame_kinematics.Qcp*(-frame_kinematics.dx-cross(frame_kinematics.omega_pc_parent,-frame_kinematics.x+current_position)+current_velocity);
            current_position = frame_kinematics.Qcp*(-frame_kinematics.x + current_position); 
        }

        return {current_position, current_velocity};
    }
    else if ( target.generation == common_generation )
    {
        // Only go from body to target
        Vector3d<return_type> current_position = _x + _Qpc*x;
        Vector3d<return_type> current_velocity = _dx + _Qpc*(cross(_omega_pc_self,x) + dx);

        for (size_t i_generation = generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = get_frame_kinematics_at_generation<return_type>(i_generation);
            
            current_velocity = frame_kinematics.dx + frame_kinematics.Qpc*(cross(frame_kinematics.omega_pc_self,current_position) + current_velocity);
            current_position = frame_kinematics.x + frame_kinematics.Qpc*current_position;
        }
        
        return {current_position, current_velocity};
    }
    else
    {
        // Go from body to common generation
        Vector3d<return_type> current_position = _x + _Qpc*x;
        Vector3d<return_type> current_velocity = _dx + _Qpc*(cross(_omega_pc_self,x) + dx);

        for (size_t i_generation = generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = get_frame_kinematics_at_generation<return_type>(i_generation);
            
            current_velocity = frame_kinematics.dx + frame_kinematics.Qpc*(cross(frame_kinematics.omega_pc_self,current_position) + current_velocity);
            current_position = frame_kinematics.x + frame_kinematics.Qpc*current_position;
        }

        // Go from common frame to target
        for (size_t i_generation = common_generation + 1; i_generation <= target.generation; ++i_generation)
        {
            auto frame_kinematics = target.template get_frame_kinematics_at_generation<return_type>(i_generation);

            current_velocity = frame_kinematics.Qcp*(-frame_kinematics.dx-cross(frame_kinematics.omega_pc_parent,-frame_kinematics.x+current_position)+current_velocity);
            current_position = frame_kinematics.Qcp*(-frame_kinematics.x + current_position); 
        }

        return {current_position, current_velocity};
    }
}

template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline auto Frame<T, Parent_frame_type, Number_of_rotations>::get_absolute_rotation_matrix() const -> Matrix3x3<aggregated_type>
{
    if constexpr (is_inertial) 
        return Matrix3x3<T>::eye();

    else if constexpr (!Parent_frame_type::is_inertial)
        return get_parent().get_absolute_rotation_matrix() * get_rotation_matrix();

    else
        return get_rotation_matrix();
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline auto Frame<T,Parent_frame_type,Number_of_rotations>::get_omega_absolute_in_body() const -> Vector3d<aggregated_type>
{
    if constexpr (is_inertial)
        return Vector3d<T>::zeros();
    else
    {
        const auto& omega_wrt_parent = get_omega_wrt_parent_in_body();
        const auto omega_parent = get_back_rotation_matrix() * get_parent().get_omega_absolute_in_body();

        return omega_wrt_parent + omega_parent;
    }
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline auto Frame<T,Parent_frame_type,Number_of_rotations>::get_omega_absolute_in_parent() const -> Vector3d<aggregated_type>
{
    if constexpr (is_inertial)
        return Vector3d<aggregated_type>::zeros();
    else
    {
        const auto& omega_wrt_parent = get_omega_wrt_parent_in_parent();
        const auto  omega_parent = get_parent().get_omega_absolute_in_body();

        return omega_wrt_parent + omega_parent;
    }
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename Input_frame>
inline auto Frame<T,Parent_frame_type,Number_of_rotations>::get_rotation_matrix(const Input_frame& target) const 
    -> Matrix3x3<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>
/*
*       frame_i.get_rotation_matrix(frame_j) = Qji
*/
{
    using return_type = typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type;

    // Self rotation...
    if constexpr (std::is_same_v<Frame<T,Parent_frame_type,Number_of_rotations>, Input_frame>)
        if ( this == &target ) 
            return Matrix3x3<return_type>::eye();

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation == common_generation ) 
    {
        // Go through the path only for the target
        Matrix3x3<return_type> Qtarget = target.get_rotation_matrix();

        for (size_t i_generation = target.generation - 1; i_generation > common_generation; --i_generation)
        {
            Qtarget = target.template get_frame_kinematics_at_generation<return_type>(i_generation).Qpc * Qtarget;
        }

        return transpose(Qtarget);
    }
    else if ( target.generation == common_generation )
    {
        // Go through the path only for the body
        auto Qbody = get_rotation_matrix();

        for ( size_t i_generation = generation-1; i_generation > common_generation; --i_generation)
        {
            Qbody = get_frame_kinematics_at_generation<return_type>(i_generation).Qpc * Qbody;
        }

        return Qbody;
    }
    else
    {
        Matrix3x3<return_type> Qtarget = target.get_rotation_matrix();

        for (size_t i_generation = target.generation - 1; i_generation > common_generation; --i_generation)
        {
            Qtarget = target.template get_frame_kinematics_at_generation<return_type>(i_generation).Qpc * Qtarget;
        }

        Matrix3x3<return_type> Qbody = get_rotation_matrix();

        for (size_t i_generation = generation - 1; i_generation > common_generation; --i_generation)
        {
            Qbody = get_frame_kinematics_at_generation<return_type>(i_generation).Qpc * Qbody;
        }

        return transpose(Qtarget)*Qbody;
    }
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename Input_frame>
inline auto Frame<T,Parent_frame_type,Number_of_rotations>::get_omega_in_body(const Input_frame& target) const 
    -> Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>
{
    using return_type = typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type;

    // Self rotation...
    if constexpr (std::is_same_v<Frame<T,Parent_frame_type,Number_of_rotations>, Input_frame>)
        if ( this == &target ) 
            return Vector3d<return_type>::zeros();

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation == common_generation ) 
    {
        // Go through the path only for the target
        auto omega_target = -target.get_omega_wrt_parent_in_parent();

        for (size_t i_generation = target.generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = target.template get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_target = -frame_kinematics.omega_pc_parent + frame_kinematics.Qpc * omega_target;
        }

        return omega_target;
    }
    else if ( target.generation == common_generation )
    {
        // Go through the path only for the body
        auto omega_body = get_omega_wrt_parent_in_body();
        auto Qbody_currentparent = _Qcp;

        for (size_t i_generation = generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_body = Qbody_currentparent * frame_kinematics.omega_pc_self + omega_body;
            Qbody_currentparent = Qbody_currentparent * frame_kinematics.Qcp;
        }

        return omega_body;
    }
    else
    {
        // Go through the path only for the target
        auto omega_target = -target.get_omega_wrt_parent_in_parent();

        for (size_t i_generation = target.generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = target.template get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_target = -frame_kinematics.omega_pc_parent + frame_kinematics.Qpc * omega_target;
        }

        // Go through the path only for the body
        auto omega_body = get_omega_wrt_parent_in_body();
        auto Qbody_currentparent = _Qcp;

        for (size_t i_generation = generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_body = Qbody_currentparent * frame_kinematics.omega_pc_self + omega_body;
            Qbody_currentparent = Qbody_currentparent * frame_kinematics.Qcp;
        }

        return Qbody_currentparent*omega_target + omega_body;
    }
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename Input_frame>
inline auto Frame<T,Parent_frame_type,Number_of_rotations>::get_omega_in_target(const Input_frame& target) const
    -> Vector3d<typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type>
{
    using return_type = typename combine_types<aggregated_type, typename Input_frame::aggregated_type>::type;

    // Self rotation...
    if constexpr (std::is_same_v<Frame<T,Parent_frame_type,Number_of_rotations>, Input_frame>)
        if ( this == &target ) 
            return Vector3d<return_type>::zeros();

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation == common_generation ) 
    {
        // Go through the path only for the target
        auto omega_target = -target.get_omega_wrt_parent_in_body();
        auto Qcp = target.get_back_rotation_matrix();

        for (size_t i_generation = target.generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = target.template get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_target = -Qcp * frame_kinematics.omega_pc_self + omega_target;
            Qcp *= frame_kinematics.Qcp;
        }
    
        return omega_target;
    }
    else if ( target.generation == common_generation )
    {
        // Go through the path only for the body
        auto omega_body = get_omega_wrt_parent_in_parent();

        for (size_t i_generation = generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_body = frame_kinematics.omega_pc_parent + frame_kinematics.Qpc * omega_body;
        }

        return omega_body;
    }
    else
    {
        // Go through the path only for the target
        auto omega_target = -target.get_omega_wrt_parent_in_body();
        auto Qcp = target.get_back_rotation_matrix();

        for (size_t i_generation = target.generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = target.template get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_target = -Qcp * frame_kinematics.omega_pc_self + omega_target;
            Qcp *= frame_kinematics.Qcp;
        }

        // Go through the path only for the body
        auto omega_body = get_omega_wrt_parent_in_parent();

        for (size_t i_generation = generation - 1; i_generation > common_generation; --i_generation)
        {
            auto frame_kinematics = get_frame_kinematics_at_generation<return_type>(i_generation);
            omega_body = frame_kinematics.omega_pc_parent + frame_kinematics.Qpc * omega_body;
        }

        return omega_target + Qcp*omega_body;
    }
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline auto Frame<T,Parent_frame_type,Number_of_rotations>::get_omega_absolute_in_inertial() const -> Vector3d<aggregated_type>
{
    if constexpr (is_inertial)
        return Vector3d<aggregated_type>::zeros();

    return get_absolute_rotation_matrix()*get_omega_absolute_in_body();
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
constexpr inline void Frame<T,Parent_frame_type,Number_of_rotations>::set_origin(const Vector3d<T>& x, const bool bUpdate) 
{ 
    if constexpr (is_inertial)
        throw lion_exception("[ERROR] Frame::set_origin -> origin cannot be set for an inertial frame");

    _x = x; 
    _updated = false; 

    if ( bUpdate ) 
        update(); 

    return;
} 


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
constexpr inline void Frame<T,Parent_frame_type,Number_of_rotations>::set_origin(const Vector3d<T>& x, const Vector3d<T>& dx, const bool bUpdate) 
{ 
    if constexpr (is_inertial)
        throw lion_exception("[ERROR] Frame::set_origin -> origin cannot be set for an inertial frame");

    _x = x;
    _dx = dx; 
    _updated = false; 

    if ( bUpdate )
        update(); 

    return;
}  


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
constexpr inline void Frame<T,Parent_frame_type,Number_of_rotations>::set_velocity(const Vector3d<T>& dx, const bool bUpdate)
{ 
    if constexpr (is_inertial)
        throw lion_exception("[ERROR] Frame::set_velocity -> origin cannot be set for an inertial frame");
    
    _dx = dx; 
    _updated = false; 
        
    if ( bUpdate )
        update();

    return;
}


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline void Frame<T,Parent_frame_type,Number_of_rotations>::set_rotation_angle(const size_t which, const T angle, const bool bUpdate) 
{ 
    if constexpr (is_inertial)
        throw lion_exception("[ERROR] Frame::set_rotation_angle -> rotation angles cannot be set for an inertial frame");

    _angles.at(which) = angle; 
    _updated = false; 

    if ( bUpdate ) 
        update(); 

    return;
} 


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline void Frame<T,Parent_frame_type,Number_of_rotations>::set_rotation_angle(const size_t which, const T angle, const T dangle, const bool bUpdate) 
{
    if constexpr (is_inertial)
        throw lion_exception("[ERROR] Frame::set_rotation_angle -> rotation angles cannot be set for an inertial frame");

    _angles.at(which) = angle; 
    _dangles.at(which) = dangle; 
    _updated = false;

    if ( bUpdate )
        update();
    
    return;
} 


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
inline void Frame<T,Parent_frame_type,Number_of_rotations>::set_angular_speed(const size_t which, const T dangle, const bool bUpdate) 
{
    if constexpr (is_inertial)
        throw lion_exception("[ERROR] Frame::set_angular_speed -> angular speed cannot be set for an inertial frame");

    _dangles.at(which) = dangle; 
    _updated = false; 

    if ( bUpdate )
        update();

    return;
}  


template<typename T, typename Parent_frame_type, size_t Number_of_rotations>
template<typename return_type>
inline auto Frame<T, Parent_frame_type, Number_of_rotations>::get_frame_kinematics_at_generation(size_t requested_generation) const -> lioncpp::detail::Frame_kinematics<return_type>
{
    if (requested_generation > generation)
    {
        throw lion_exception("[ERROR] Frame::get_frame_kinematics_at_generation -> input generation was lower than actual generation");
    }

    if constexpr (generation == 0)
    {
        if (requested_generation == 0)
        {
            return lioncpp::detail::Frame_kinematics<return_type>
            {
                .x = _x,
                .dx = _dx,
                .Qpc = _Qpc,
                .Qcp = _Qcp,
                .omega_pc_self = _omega_pc_self,
                .omega_pc_parent = _omega_pc_parent
            };
        }
    }
    else
    {
        if (requested_generation == generation)
        {
            return lioncpp::detail::Frame_kinematics<return_type>
            {
                .x = _x,
                .dx = _dx,
                .Qpc = _Qpc,
                .Qcp = _Qcp,
                .omega_pc_self = _omega_pc_self,
                .omega_pc_parent = _omega_pc_parent
            };
        }
        else
        {
            return get_parent().template get_frame_kinematics_at_generation<return_type>(requested_generation);
        }
    }
}

#endif
