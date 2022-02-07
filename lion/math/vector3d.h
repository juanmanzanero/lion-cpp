#ifndef _TR7_BASE_VECTOR3D_H_
#define _TR7_BASE_VECTOR3D_H_


#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>
#include <cmath>
#include "lion/foundation/constants.h"
#include "lion/foundation/type_traits.h"
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

    template<typename U>
    explicit constexpr Vector3d(const Vector3d<U>& other);

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

template<typename U, typename V, typename W = typename combine_types<U,V>::type>
constexpr W dot(const Vector3d<U> &lhs, const Vector3d<V> &rhs);

template<typename U, typename V, typename W = typename combine_types<U,V>::type>
constexpr Vector3d<W> cross(const Vector3d<U> &lhs, const Vector3d<V> &rhs);

template<typename T>
constexpr T distance(const Vector3d<T> &lhs, const Vector3d<T> &rhs);

template<typename T>
constexpr Vector3d<T> sqrt(const Vector3d<T> &arg);

template<typename T>
constexpr Vector3d<T> abs(const Vector3d<T> &arg);

template<typename T>
constexpr T sum(const Vector3d<T> &arg);

template<typename U, typename V, typename W = typename combine_types<U,V>::type>
constexpr Vector3d<W> operator+(const Vector3d<U> &lhs, const Vector3d<V> &rhs);

template<typename T>
constexpr Vector3d<T> operator-(const Vector3d<T> &rhs);

template<typename U, typename V, typename W = typename combine_types<U,V>::type>
constexpr Vector3d<W> operator-(const Vector3d<U> &lhs, const Vector3d<V> &rhs);

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

template<typename U, typename V, typename W = typename combine_types<U,V>::type>
constexpr Vector3d<W> operator*(const Vector3d<U> &lhs, V val);

template<typename U, typename V, typename W = typename combine_types<U,V>::type>
constexpr Vector3d<W> operator*(U val, const Vector3d<V> &rhs);

template<typename T>
constexpr Vector3d<T> operator/(const Vector3d<T> &lhs, T val);

template<typename T>
constexpr Vector3d<T> operator/(T val, const Vector3d<T> &rhs);

template<typename T>
std::ostream& operator<<(std::ostream &os, const Vector3d<T> &v);

template<typename T>
std::istream& operator>>(std::istream &is, Vector3d<T> &v);

using sVector3d = Vector3d<scalar>;

//! Add the capability to add vectors of different type to scalar vectors
template<typename T>
typename std::enable_if<!std::is_same<T,scalar>::value,Vector3d<T>>::type operator+(Vector3d<T> lhs, const sVector3d& rhs);

template<typename T>
typename std::enable_if<!std::is_same<T,scalar>::value,Vector3d<T>>::type operator+(const sVector3d& lhs, const Vector3d<T>& rhs) { return rhs+lhs; }

template<typename T>
typename std::enable_if<!std::is_same<T,scalar>::value,Vector3d<T>>::type operator/(const Vector3d<T>& lhs, const scalar rhs) { return {lhs[0]/rhs, lhs[1]/rhs, lhs[2]/rhs}; }

#endif
