#ifndef _TR7_BASE_MATRIX3X3_H_
#define _TR7_BASE_MATRIX3X3_H_


#include <initializer_list>


#include "vector3d.h"


//
// Defines the typical 3 x 3 matrix of Ts, stored in column-major format.
//
template <typename T>
struct Matrix3x3 : std::array<T, 9>
{
    //
    // 3 x 3 Matrix and its algebra, [xx, xy, xz; ...
    //                                yx, yy, yz; ...
    //                                zx, zy, zz], 
    // stored in column-major format.
    //

    using base_type = std::array<T, 9>;

    static constexpr std::size_t size = std::tuple_size<base_type>::value;


    constexpr Matrix3x3() : base_type{} {}
    constexpr Matrix3x3(T xx, T xy, T xz,
                        T yx, T yy, T yz,
                        T zx, T zy, T zz);
    constexpr explicit Matrix3x3(T val);
    constexpr Matrix3x3(const Vector3d<T> &row_x, const Vector3d<T> &row_y, const Vector3d<T> &row_z);
    constexpr explicit Matrix3x3(const base_type &arr);
    explicit Matrix3x3(const std::vector<T> &m_colmaj) noexcept(false);
    explicit Matrix3x3(std::initializer_list<std::initializer_list<T>> ilist_rows) noexcept(false);

    constexpr const T& xx() const { return (*this)[0]; }
    constexpr const T& yx() const { return (*this)[1]; }
    constexpr const T& zx() const { return (*this)[2]; }

    constexpr const T& xy() const { return (*this)[3]; }
    constexpr const T& yy() const { return (*this)[4]; }
    constexpr const T& zy() const { return (*this)[5]; }

    constexpr const T& xz() const { return (*this)[6]; }
    constexpr const T& yz() const { return (*this)[7]; }
    constexpr const T& zz() const { return (*this)[8]; }

    constexpr T& xx() { return (*this)[0]; }
    constexpr T& yx() { return (*this)[1]; }
    constexpr T& zx() { return (*this)[2]; }

    constexpr T& xy() { return (*this)[3]; }
    constexpr T& yy() { return (*this)[4]; }
    constexpr T& zy() { return (*this)[5]; }

    constexpr T& xz() { return (*this)[6]; }
    constexpr T& yz() { return (*this)[7]; }
    constexpr T& zz() { return (*this)[8]; }

    constexpr Matrix3x3 t() const;
    constexpr T trace() const;
    constexpr T det() const;

    constexpr bool is_orthogonal() const;
    constexpr bool is_rotmat() const;
    constexpr bool is_eye() const;

    constexpr static Matrix3x3 eye();
    constexpr static Matrix3x3 zeros();
    constexpr static Matrix3x3 ones();
    constexpr static Matrix3x3 columns(const Vector3d<T> &col_x,
                                       const Vector3d<T> &col_y,
                                       const Vector3d<T> &col_z);

    constexpr const T& operator()(std::size_t i, std::size_t j) const;
    constexpr T& operator()(std::size_t i, std::size_t j);

    constexpr Matrix3x3& operator+=(const Matrix3x3 &other);
    constexpr Matrix3x3& operator-=(const Matrix3x3 &other);
    constexpr Matrix3x3& operator*=(const Matrix3x3 &other);

    constexpr Matrix3x3& operator=(T val);
    constexpr Matrix3x3& operator+=(T val);
    constexpr Matrix3x3& operator-=(T val);
    constexpr Matrix3x3& operator*=(T val);
    constexpr Matrix3x3& operator/=(T val);

    constexpr Matrix3x3& operator=(const base_type &arr);


};

template <typename T>
constexpr Matrix3x3<T> transpose(const Matrix3x3<T> &arg);

template <typename T>
constexpr T trace(const Matrix3x3<T> &arg);

template <typename T>
constexpr T det(const Matrix3x3<T> &arg);

template <typename T>
constexpr bool is_orthogonal(const Matrix3x3<T> &arg);

template <typename T>
constexpr bool is_rotmat(const Matrix3x3<T> &arg);

template <typename T>
constexpr bool is_eye(const Matrix3x3<T> &arg);


template <typename T>
constexpr Matrix3x3<T> inv(const Matrix3x3<T> &arg);

template <typename T>
constexpr T dot(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs);

template <typename T>
constexpr T sumsqr(const Matrix3x3<T> &arg);

template <typename T>
constexpr T sumabs(const Matrix3x3<T> &arg);

template <typename T>
constexpr Matrix3x3<T> sqrt(const Matrix3x3<T> &arg);

template <typename T>
constexpr Matrix3x3<T> abs(const Matrix3x3<T> &arg);

template <typename T>
constexpr Matrix3x3<T> operator+(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs);

template <typename T>
constexpr Matrix3x3<T> operator-(const Matrix3x3<T> &rhs);

template <typename T>
constexpr Matrix3x3<T> operator-(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs);

template <typename T>
constexpr Matrix3x3<T> operator*(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs);


template <typename T>
constexpr Matrix3x3<T> operator+(const Matrix3x3<T> &lhs, T val);

template <typename T>
constexpr Matrix3x3<T> operator+(T lhs, const Matrix3x3<T> &rhs);

template <typename T>
constexpr Matrix3x3<T> operator-(const Matrix3x3<T> &lhs, T val);

template <typename T>
constexpr Matrix3x3<T> operator-(T lhs, const Matrix3x3<T> &rhs);

template <typename T>
constexpr Matrix3x3<T> operator*(const Matrix3x3<T> &lhs, T val);

template <typename T>
constexpr Matrix3x3<T> operator*(T lhs, const Matrix3x3<T> &rhs);

template <typename T>
constexpr Matrix3x3<T> operator/(const Matrix3x3<T> &lhs, T val);

template <typename T>
std::ostream& operator<<(std::ostream &os, const Matrix3x3<T> &m);

template <typename T>
std::istream& operator>>(std::istream &is, Matrix3x3<T> &m);


// outer operations
template <typename T>
constexpr Vector3d<T> operator*(const Matrix3x3<T> &m, const Vector3d<T> &v);

template <typename T>
constexpr Vector3d<T> operator*(const Vector3d<T> &v, const Matrix3x3<T> &m);

template <typename T>
Vector3d<T> linsolve(const Matrix3x3<T> &m_lhs, const Vector3d<T> &v);

template <typename T>
Vector3d<T> linsolve(const Vector3d<T> &v, const Matrix3x3<T> &m_rhs);

template <typename T>
constexpr Matrix3x3<T> outer(const Vector3d<T> &v_lhs, const Vector3d<T> &v_rhs);

template <typename T>
constexpr Matrix3x3<T> diag(const Vector3d<T> &v);

template <typename T>
constexpr Matrix3x3<T> crossmat(const Vector3d<T> &v);

using sMatrix3x3 = Matrix3x3<scalar>;

// implementation header
#include "matrix3x3.hpp"
#endif
