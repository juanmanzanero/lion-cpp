#ifndef __ROTATIONS_HPP__
#define __ROTATIONS_HPP__

template<typename T>
constexpr Matrix3x3<T> rotation_matrix_x(const T phi)
{
    return Matrix3x3( 1.0,           0.0,            0.0,
                      0.0, std::cos(phi), -std::sin(phi),
                      0.0, std::sin(phi), std::cos(phi) );
}


template<typename T>
constexpr Matrix3x3<T> rotation_matrix_y(const T phi)
{
    return Matrix3x3(  std::cos(phi), 0.0, std::sin(phi),
		                 0.0, 1.0,           0.0,
                      -std::sin(phi), 0.0, std::cos(phi) ); 
}


template<typename T>
constexpr Matrix3x3<T> rotation_matrix_z(const T phi)
{
    return Matrix3x3(  std::cos(phi), -std::sin(phi), 0.0,
                       std::sin(phi),  std::cos(phi), 0.0,
                                 0.0,            0.0, 1.0 ); 
}


template<typename T>
constexpr Matrix3x3<T> drotation_matrix_x(const T phi)
{
    return Matrix3x3( 0.0,            0.0,            0.0,
                      0.0, -std::sin(phi), -std::cos(phi),
                      0.0,  std::cos(phi), -std::sin(phi) );
}


template<typename T>
constexpr Matrix3x3<T> drotation_matrix_y(const T phi)
{
    return Matrix3x3( -std::sin(phi), 0.0,  std::cos(phi),
		                 0.0        , 0.0,            0.0,
                      -std::cos(phi), 0.0, -std::sin(phi) ); 
}


template<typename T>
constexpr Matrix3x3<T> drotation_matrix_z(const T phi)
{
    return Matrix3x3( -std::sin(phi), -std::cos(phi), 0.0,
                       std::cos(phi), -std::sin(phi), 0.0,
                                 0.0,            0.0, 0.0 ); 
}


template<typename T>
constexpr Vector3d<T> omega_vector_from_matrix(const Matrix3x3<T>& matrix)
{
    return Vector3d<T>(matrix(2,1),matrix(0,2),matrix(1,0));
}


#endif
