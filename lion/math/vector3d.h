#ifndef _TR7_BASE_VECTOR3D_H_
#define _TR7_BASE_VECTOR3D_H_


#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>
#include <cmath>
#include "lion/foundation/utils.h"
#include "lion/foundation/constants.h"
#include "lion/foundation/types.h"

//
// Defines the typical 3 x 1 array of Ts.
//
template <typename T>
struct Vector3d : std::array<T, 3>
{
    //
    // 3D Vector [x; y; z], and its algebra
    //
    using base_type = std::array<T, 3>;

    static constexpr std::size_t size = std::tuple_size<base_type>::value;


    constexpr Vector3d() : base_type{} {}
    constexpr Vector3d(T x, T y, T z);
    constexpr explicit Vector3d(T val);
    constexpr explicit Vector3d(const base_type &arr);
    explicit Vector3d(const std::vector<T> &v) noexcept(false);

    constexpr const T& x() const { return (*this)[0]; }
    constexpr const T& y() const { return (*this)[1]; }
    constexpr const T& z() const { return (*this)[2]; }

    constexpr T& x() { return (*this)[0]; }
    constexpr T& y() { return (*this)[1]; }
    constexpr T& z() { return (*this)[2]; }

    constexpr Vector3d& normalize();

    constexpr T norm() const;
    constexpr T safe_norm(T tol = eps) const;

    constexpr static Vector3d zeros();
    constexpr static Vector3d ones();

    constexpr Vector3d& operator+=(const Vector3d &other);
    constexpr Vector3d& operator-=(const Vector3d &other);
    constexpr Vector3d& operator*=(const Vector3d &other);
    constexpr Vector3d& operator/=(const Vector3d &other);

    constexpr Vector3d& operator=(T val);
    constexpr Vector3d& operator+=(T val);
    constexpr Vector3d& operator-=(T val);
    constexpr Vector3d& operator*=(T val);
    constexpr Vector3d& operator/=(T val);

    constexpr Vector3d& operator=(const base_type &arr);
};


template<typename T>
constexpr Vector3d<T> normalize(Vector3d<T> arg);

template<typename T>
constexpr T norm(const Vector3d<T> &arg);

template<typename T>
constexpr T angle(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr T dot(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> cross(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr T distance(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> sqrt(const Vector3d<T> &arg);

template<typename T>
constexpr Vector3d<T> abs(const Vector3d<T> &arg);

template<typename T>
constexpr T sum(const Vector3d<T> &arg);

template<typename T>
constexpr Vector3d<T> operator+(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator-(const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator-(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator*(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator/(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator+(const Vector3d<T> &lhs, T val);

template<typename T>
constexpr Vector3d<T> operator+(T val, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator-(const Vector3d<T> &lhs, T val);

template<typename T>
constexpr Vector3d<T> operator-(T val, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator*(const Vector3d<T> &lhs, T val);

template<typename T>
constexpr Vector3d<T> operator*(T val, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> operator/(const Vector3d<T> &lhs, T val);

template<typename T>
constexpr Vector3d<T> operator/(T val, const Vector3d<T> &rhs);

template<typename T>
std::ostream& operator<<(std::ostream &os, const Vector3d<T> &v);

template<typename T>
std::istream& operator>>(std::istream &is, Vector3d<T> &v);

using sVector3d = Vector3d<scalar>;
using tVector3d = Vector3d<timeseries>;

#include "vector3d.hpp"

constexpr sVector3d UX{1.0, 0.0, 0.0};
constexpr sVector3d UY{0.0, 1.0, 0.0};
constexpr sVector3d UZ{0.0, 0.0, 1.0};

#endif
