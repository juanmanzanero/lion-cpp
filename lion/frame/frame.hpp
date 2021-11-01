#include <tuple>

template<typename T>
inline Frame<T>::Frame(const typename Frame<T>::tVector3d& x, const typename Frame<T>::tVector3d& dx, const std::vector<T> angles, 
                    const std::vector<T> dangles, const std::vector<Axis> axis, const Frame<T>& parent) 
: _parent(&parent),
  _x(x), 
  _dx(dx), 
  _angles(angles), 
  _dangles(dangles), 
  _axis(axis),
  _updated(false)
{
 try
 { 
    if ( _angles.size() != _dangles.size() ) 
    {
        throw std::runtime_error("angles and angles derivatives vectors do not have the same size"); 
    } 
    if ( _angles.size() != _axis.size() ) 
    {
        throw std::runtime_error("angles and rotation axis vectors do not have the same size"); 
    } 

    for (size_t i = 0; i < _axis.size(); ++i)
    {
        switch(_axis[i])
        {
         case X: break;
         case Y: break;
         case Z: break;

         default:
            throw std::runtime_error("Invalid rotation axis: "+std::to_string(_axis[i])+", in position "+std::to_string(i)); 
         }
    }

    update();
 } 
 catch( const std::runtime_error& error )
 {
    throw;
 }
}


template<typename T>
inline void Frame<T>::update()
{
    if ( is_inertial() )
        return;

    if ( _angles.size() == 0 )
    {
        _Qpc = tMatrix3x3::eye();
        _Qcp = tMatrix3x3::eye();
        _omega_pc_self = tVector3d::zeros();
        _omega_pc_parent = tVector3d::zeros();

        _updated = true;
        return;
    }

//  Compute the rotation matrix
    std::vector<tMatrix3x3> fwd_rot_matrices;
    std::vector<tMatrix3x3> bwd_rot_matrices;

//  Accumulated forward rotation matrices: X_parent = Q_{12}Q_{23}...Q_{j-1,j}Xj
    std::vector<tMatrix3x3> accumulated_fwd_rot_matrices(_angles.size());

//  Accumulated backward rotation matrices: X_child = Q_{N,N-1}Q_{N-1,N-2}...Q_{j+1,j}Xj
    std::vector<tMatrix3x3> accumulated_bwd_rot_matrices(_angles.size());

    for (size_t i = 0; i < _angles.size(); ++i)
    {
        switch(_axis[i])
        {
         case(X):
            fwd_rot_matrices.push_back(rotation_matrix_x(_angles[i]));
            break;

         case(Y):
            fwd_rot_matrices.push_back(rotation_matrix_y(_angles[i]));
            break;

         case(Z):
            fwd_rot_matrices.push_back(rotation_matrix_z(_angles[i]));
            break;
        } 

        bwd_rot_matrices.push_back(transpose(fwd_rot_matrices.at(i)));
    }

    accumulated_fwd_rot_matrices.front() = fwd_rot_matrices.front();

    for (size_t i=1; i< _angles.size(); ++i)
        accumulated_fwd_rot_matrices[i] = accumulated_fwd_rot_matrices[i-1]*fwd_rot_matrices[i];


    accumulated_bwd_rot_matrices.back() = bwd_rot_matrices.back();
    for (int i = _angles.size()-2; i >= 0; --i)
        accumulated_bwd_rot_matrices[i] = accumulated_bwd_rot_matrices[i+1]*bwd_rot_matrices[i];

    _Qpc = accumulated_fwd_rot_matrices.back();
    _Qcp = accumulated_bwd_rot_matrices.front();

    // Compute omega
    _omega_pc_self = tVector3d::zeros();
    _omega_pc_parent = tVector3d::zeros();

    for (size_t i = 0; i < _angles.size(); ++i)
    {
        tVector3d new_rot(tVector3d::zeros());
        new_rot[_axis[i]] = _dangles[i];
        _omega_pc_self += accumulated_bwd_rot_matrices[i]*new_rot;
        _omega_pc_parent += accumulated_fwd_rot_matrices[i]*new_rot;
    }

    _updated = true;

    return;
}


template<typename T>
inline const Frame<T>& Frame<T>::get_parent() const 
{
 try
 { 
    if ( !is_inertial() ) 
        return *_parent; 
    else
        throw std::runtime_error("Inertial frames have no parents");
    
 }
 catch( const std::runtime_error& error )
 {
    throw;
 }
}


template<typename T>
constexpr inline size_t Frame<T>::generation() const
{
    size_t gen = 0;
    const Frame* current = this;
    while ( current->get_parent_ptr() != nullptr )
    {
        ++gen;
        current = current->get_parent_ptr();
    }

    return gen;
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_absolute_position(const typename Frame<T>::tVector3d& x) const
{
    if (is_inertial())
        return x;

    return get_parent_ptr()->get_absolute_position(_x + _Qpc*x);
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_absolute_velocity_in_body(const typename Frame<T>::tVector3d& x, const typename Frame<T>::tVector3d& dx) const
{
    if ( is_inertial() )
        return dx;

    return _Qcp*get_parent_ptr()->get_absolute_velocity_in_body(_x + _Qpc*x, _dx + _Qpc*(dx+cross(_omega_pc_self,x)));
}

    
template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_absolute_velocity_in_parent(const typename Frame<T>::tVector3d& x, const typename Frame<T>::tVector3d& dx) const
{
    if ( is_inertial() )
        return dx;

    return _Qpc*get_absolute_velocity_in_body(x, dx);
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_absolute_velocity_in_inertial(const typename Frame<T>::tVector3d& x, const typename Frame<T>::tVector3d& dx) const
{
    if ( is_inertial() )
        return typename Frame<T>::tVector3d(dx);

    return get_parent_ptr()->get_absolute_velocity_in_inertial(_x + _Qpc*x, _dx + _Qpc*(dx+cross(_omega_pc_self,x)));
}
 
template<typename T>
inline std::pair<typename Frame<T>::tVector3d,typename Frame<T>::tVector3d> Frame<T>::get_position_and_velocity_in_target(const Frame<T>& target, const typename Frame<T>::tVector3d& x, const typename Frame<T>::tVector3d& dx) const
{
 try
 {
    // Target is self
    if ( this == &target ) 
        return {x,dx};

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation() == common_generation ) 
    {
        // Get the list of frames from common_generation to the target frame: start=common_gen_frame+1 -> end=target
        std::vector<const Frame*> frames(target.generation() - common_generation);

        const Frame* cf = &target;
        for (size_t gen = target.generation(); gen > common_generation; --gen)
        {
            frames.at(gen-common_generation-1) = cf;
            cf = cf->get_parent_ptr();
        }

        tVector3d current_velocity = dx;
        tVector3d current_position = x;

        for (size_t f = 0; f < target.generation() - common_generation; ++f)
        {
            current_velocity = frames.at(f)->_Qcp*(-frames.at(f)->_dx-cross(frames.at(f)->_omega_pc_parent,-frames.at(f)->_x+current_position)+current_velocity);
            current_position = frames.at(f)->_Qcp*(-frames.at(f)->_x + current_position); 
        }

        return {current_position, current_velocity};
    }
    else if ( target.generation() == common_generation )
    {
        // Only go from body to target
        tVector3d current_position = _x + _Qpc*x;
        tVector3d current_velocity = _dx + _Qpc*(cross(_omega_pc_self,x) + dx);
        const Frame* cf = this;
        
        for (size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cf = cf->get_parent_ptr();

            current_velocity = cf->_dx + cf->_Qpc*(cross(cf->_omega_pc_self,current_position) + current_velocity);
            current_position = cf->_x + cf->_Qpc*current_position;
        }

        return {current_position, current_velocity};
    }
    else
    {
        // Go from body to common generation
        tVector3d current_position = _x + _Qpc*x;
        tVector3d current_velocity = _dx + _Qpc*(cross(_omega_pc_self,x) + dx);
        const Frame* cf = this;
        
        for (size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cf = cf->get_parent_ptr();

            current_velocity = cf->_dx + cf->_Qpc*(cross(cf->_omega_pc_self,current_position) + current_velocity);
            current_position = cf->_x + cf->_Qpc*current_position;
        }

        // Go from common frame to target

        // Get the list of frames from common_generation to the target frame: start=common_gen_frame+1 -> end=target
        std::vector<const Frame*> frames(target.generation() - common_generation);

        cf = &target;
        for (size_t gen = target.generation(); gen > common_generation; --gen)
        {
            frames.at(gen-common_generation-1) = cf;
            cf = cf->get_parent_ptr();
        }

        for (size_t f = 0; f < target.generation() - common_generation; ++f)
        {
            current_velocity = frames.at(f)->_Qcp*(-frames.at(f)->_dx-cross(frames.at(f)->_omega_pc_parent,-frames.at(f)->_x+current_position)+current_velocity);
            current_position = frames.at(f)->_Qcp*(-frames.at(f)->_x + current_position); 
        }

        return {current_position, current_velocity};
    }
 }
 catch (const std::runtime_error& error)
 {
    throw;
 }
}

template<typename T>
inline typename Frame<T>::tMatrix3x3 Frame<T>::get_absolute_rotation_matrix() const
{
    if ( is_inertial() ) 
        return tMatrix3x3::eye();

    else if ( !(get_parent_ptr()->is_inertial()) )
        return (get_parent_ptr()->get_absolute_rotation_matrix()) * get_rotation_matrix();

    else
        return get_rotation_matrix();

}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_omega_absolute_in_body() const
{
    if ( is_inertial() )
        return tVector3d::zeros();

    const tVector3d omega_wrt_parent(get_omega_wrt_parent_in_body());
    const tVector3d omega_parent(get_back_rotation_matrix()*get_parent_ptr()->get_omega_absolute_in_body());

    return omega_wrt_parent + omega_parent;
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_omega_absolute_in_parent() const
{
    if ( is_inertial() )
        return tVector3d::zeros();

    const tVector3d omega_wrt_parent = get_omega_wrt_parent_in_parent();
    const tVector3d omega_parent = get_parent_ptr()->get_omega_absolute_in_body();

    return omega_wrt_parent + omega_parent;
}


template<typename T>
inline typename Frame<T>::tMatrix3x3 Frame<T>::get_rotation_matrix(const Frame<T>& target) const
/*
*       frame_i.get_rotation_matrix(frame_j) = Qji
*/
{
 try
 {
    // Self rotation...
    if ( this == &target ) 
        return tMatrix3x3::eye();

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation() == common_generation ) 
    {
        // Go through the path only for the target
        tMatrix3x3 Qtarget = target.get_rotation_matrix();
        const Frame* cft = &target;

        for ( size_t gen = target.generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            Qtarget = cft->get_rotation_matrix()*Qtarget;
        }

        return transpose(Qtarget);
    }
    else if ( target.generation() == common_generation )
    {
        // Go through the path only for the body
        tMatrix3x3 Qbody = get_rotation_matrix();
        const Frame* cfb = this;

        for ( size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cfb = cfb->get_parent_ptr();
            Qbody = cfb->get_rotation_matrix()*Qbody;
        }

        return Qbody;
    }
    else
    {
        tMatrix3x3 Qtarget = target.get_rotation_matrix();
        const Frame* cft = &target;

        for ( size_t gen = target.generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            Qtarget = cft->get_rotation_matrix()*Qtarget;
        }

        tMatrix3x3 Qbody = get_rotation_matrix();
        const Frame* cfb = this;

        for ( size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cfb = cfb->get_parent_ptr();
            Qbody = cfb->get_rotation_matrix()*Qbody;
        }

        return transpose(Qtarget)*Qbody;
    }
 }
 catch( const std::runtime_error& error )
 {
    throw;
 }
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_omega_in_body(const Frame<T>& target) const
{
 try
 {
    // Self rotation...
    if ( this == &target ) 
        return tVector3d::zeros();

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation() == common_generation ) 
    {
        // Go through the path only for the target
        tVector3d omega_target = -target.get_omega_wrt_parent_in_parent();
        const Frame* cft = &target;

        for ( size_t gen = target.generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_target = -cft->get_omega_wrt_parent_in_parent() + cft->_Qpc*omega_target;
        }

        return omega_target;
    }
    else if ( target.generation() == common_generation )
    {
        // Go through the path only for the body
        tVector3d omega_body = get_omega_wrt_parent_in_body();
        const Frame* cft = this;
        tMatrix3x3 Qbody_currentparent = _Qcp;

        for ( size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_body = Qbody_currentparent * cft->get_omega_wrt_parent_in_body() + omega_body;
            Qbody_currentparent = Qbody_currentparent*cft->_Qcp;
        }

        return omega_body;
    }
    else
    {
        // Go through the path for the target
        tVector3d omega_target = -target.get_omega_wrt_parent_in_parent();
        const Frame* cft = &target;

        for ( size_t gen = target.generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_target = -cft->get_omega_wrt_parent_in_parent() + cft->_Qpc*omega_target;
        }

        // Go through the path for the body
        tVector3d omega_body = get_omega_wrt_parent_in_body();
        cft = this;
        tMatrix3x3 Qbody_currentparent = _Qcp;

        for ( size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_body = Qbody_currentparent * cft->get_omega_wrt_parent_in_body() + omega_body;
            Qbody_currentparent = Qbody_currentparent*cft->_Qcp;
        }

        return Qbody_currentparent*omega_target + omega_body;
    }
 }
 catch( const std::runtime_error& error )
 {
    throw;
 }
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_omega_in_parent(const Frame<T>& target) const
{
 try
 {
    return _Qpc*get_omega_in_body(target);
 }
 catch(const std::runtime_error& error)
 { 
    throw;
 }
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_omega_in_target(const Frame<T>& target) const
{
 try
 {
    // Self rotation...
    if ( this == &target ) 
        return tVector3d::zeros();

    const size_t common_generation(get_crossing_generation(*this, target));

    if ( generation() == common_generation ) 
    {
        // Go through the path only for the target
        tVector3d omega_target = -target.get_omega_wrt_parent_in_body();
        const Frame* cft = &target;
        tMatrix3x3 Qcp = target._Qcp;
    
        for ( size_t gen = target.generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_target = -Qcp*cft->get_omega_wrt_parent_in_body() + omega_target;
            Qcp *= cft->_Qcp;
        }

        return omega_target;
    }
    else if ( target.generation() == common_generation )
    {
        // Go through the path only for the body
        tVector3d omega_body = get_omega_wrt_parent_in_parent();
        const Frame* cft = this;

        for ( size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_body = cft->get_omega_wrt_parent_in_parent() + cft->_Qpc*omega_body;
        }

        return omega_body;
    }
    else
    {
        // Go through the path for the target
        tVector3d omega_target = -target.get_omega_wrt_parent_in_body();
        const Frame* cft = &target;
        tMatrix3x3 Qcp = target._Qcp;
    
        for ( size_t gen = target.generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_target = -Qcp*cft->get_omega_wrt_parent_in_body() + omega_target;
            Qcp *= cft->_Qcp;
        }

        // Go through the path for the body
        tVector3d omega_body = get_omega_wrt_parent_in_parent();
        cft = this;

        for ( size_t gen = generation()-1; gen > common_generation; --gen)
        {
            cft = cft->get_parent_ptr();
            omega_body = cft->get_omega_wrt_parent_in_parent() + cft->_Qpc*omega_body;
        }

        return omega_target + Qcp*omega_body;
    }
 }
 catch( const std::runtime_error& error )
 {
    throw;
 }
}


template<typename T>
inline typename Frame<T>::tVector3d Frame<T>::get_omega_absolute_in_inertial() const
{
    if ( is_inertial() )
        return tVector3d::zeros();

    return get_absolute_rotation_matrix()*get_omega_absolute_in_body();
}


template<typename T>
inline size_t Frame<T>::get_crossing_generation(const Frame<T>& f1, const Frame<T>& f2)
{
 try
 {
    const size_t common_generation = std::min(f1.generation(), f2.generation());

    // Complete path up to generation
    const Frame* cf1 = &f1;
    for (size_t gen = f1.generation(); gen > common_generation; --gen)
        cf1 = cf1->get_parent_ptr();

    const Frame* cf2 = &f2;
    for (size_t gen = f2.generation(); gen > common_generation; --gen)
        cf2 = cf2->get_parent_ptr();

    // Get all the way to the end until two frames coincide, or the inertial frame is reached
    if ( cf1 == cf2 ) 
        return common_generation;

    for (size_t gen = common_generation; gen>0; --gen)
    {
        cf1 = cf1->get_parent_ptr();
        cf2 = cf2->get_parent_ptr();

        if ( cf1 == cf2 ) 
            return gen-1;
    }

    // If the code reaches here, means that the inertial frames are not the same. Throw an exception
    throw std::runtime_error("The two frames do not share an ancestor"); 
 }
 catch(const std::runtime_error& error)
 {
    throw;
 }
}


template<typename T>
constexpr inline void Frame<T>::set_origin(const typename Frame<T>::tVector3d& x, const bool bUpdate) 
{ 
    if (is_inertial())
        return;

    _x = x; 
    _updated = false; 

    if ( bUpdate ) 
        update(); 

    return;
} 


template<typename T>
constexpr inline void Frame<T>::set_origin(const typename Frame<T>::tVector3d& x, const typename Frame<T>::tVector3d& dx, const bool bUpdate) 
{ 
    if (is_inertial())
        return;

    _x = x;
    _dx = dx; 
    _updated = false; 

    if ( bUpdate )
        update(); 

    return;
}  


template<typename T>
constexpr inline void Frame<T>::set_velocity(const typename Frame<T>::tVector3d& dx, const bool bUpdate)
{ 
    if (is_inertial())
        return;
    
    _dx = dx; 
    _updated = false; 
        
    if ( bUpdate )
        update();

    return;
}


template<typename T>
inline void Frame<T>::set_rotation_angle(const size_t which, const T angle, const bool bUpdate) 
{ 
    _angles.at(which) = angle; 
    _updated = false; 

    if ( bUpdate ) 
        update(); 

    return;
} 


template<typename T>
inline void Frame<T>::set_rotation_angle(const size_t which, const T angle, const T dangle, const bool bUpdate) 
{
    _angles.at(which) = angle; 
    _dangles.at(which) = dangle; 
    _updated = false;

    if ( bUpdate )
        update();
    
    return;
} 


template<typename T>
inline void Frame<T>::set_angular_speed(const size_t which, const T dangle, const bool bUpdate) 
{
    _dangles.at(which) = dangle; 
    _updated = false; 

    if ( bUpdate )
        update();

    return;
}  


template<typename T>
inline void Frame<T>::add_rotation(const T angle, const T dangle, const Axis axis, const bool bUpdate) 
{ 
    if (is_inertial())
        return;

    _angles.push_back(angle) ; 
    _dangles.push_back(dangle); 
    _axis.push_back(axis); 
    _updated = false;

    if ( bUpdate )
        update();
}
