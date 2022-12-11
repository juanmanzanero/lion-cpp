#ifndef EULER_ANGLES_H
#define EULER_ANGLES_H

#include "lion/math/euler_angles.h"

template<typename T>
struct Euler_angles : public std::array<T, 3>
{
public:
    using base_type = std::array<T, 3>;

    //! Default constructor
    constexpr Euler_angles() : base_type{} {}

    //! Constructor from coordinates
    constexpr Euler_angles(T yaw, T pitch, T roll) : base_type{ yaw,pitch,roll } {};

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

#endif
