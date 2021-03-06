#pragma once

#include "lion/foundation/utils.hpp"

//
// tr7::Matrix3x3 implementation header, #included by "tr7/base/matrix3x3.h"
//
template<typename T>
constexpr Matrix3x3<T>::Matrix3x3(T xx, T xy, T xz,
                                  T yx, T yy, T yz,
                                  T zx, T zy, T zz) :
    base_type{ xx, yx, zx,
               xy, yy, zy,
               xz, yz, zz } {}

template<typename T>
constexpr Matrix3x3<T>::Matrix3x3(T val) :
    base_type{ val, val, val,
               val, val, val,
               val, val, val } {}

template<typename T>
constexpr Matrix3x3<T>::Matrix3x3(const Vector3d<T> &row_x, const Vector3d<T> &row_y, const Vector3d<T> &row_z) :
    Matrix3x3<T>{ row_x.x(), row_x.y(), row_x.z(),
               row_y.x(), row_y.y(), row_y.z(),
               row_z.x(), row_z.y(), row_z.z() } {}

template<typename T>
constexpr Matrix3x3<T>::Matrix3x3(const base_type &arr) : base_type(arr) {}

template<typename T>
inline Matrix3x3<T>::Matrix3x3(const std::vector<T> &m_colmaj) noexcept(false)
{
    const bool dims_ok{ m_colmaj.size() == Matrix3x3<T>::size };
    if (dims_ok) {
        *this =
            Matrix3x3<T>{ m_colmaj[0], m_colmaj[3], m_colmaj[6],
            m_colmaj[1], m_colmaj[4], m_colmaj[7],
            m_colmaj[2], m_colmaj[5], m_colmaj[8] };
    }
    else {
        lion_exception excp("tr7::Matrix3x3::Matrix3x3: invalid argument, the input std::vector is not of size 9.");
        std::cerr << excp.what() << std::endl;
        throw excp;
    }
}

template<typename T>
inline Matrix3x3<T>::Matrix3x3(std::initializer_list<std::initializer_list<T>> ilist_rows) noexcept(false)
{
    bool cols_ok{ ilist_rows.size() == 3 };
    if (cols_ok) {
        auto count{ 0 };
        for (const auto &row : ilist_rows) {
            auto rows_ok{ row.size() == 3 };
            if (rows_ok) {
                for (auto col : row) {
                    this->operator[](count) = col;
                    ++count;
                }
            }
            else {
                lion_exception excp("tr7::Matrix3x3::Matrix3x3: row dimensions of the input std::initializer_list must agree.");
                std::cerr << excp.what() << std::endl;
                throw excp;
            }
        }

        *this = this->t();

    }
    else {
        lion_exception excp("tr7::Matrix3x3::Matrix3x3: column dimensions of the input std::initializer_list must agree.");
        std::cerr << excp.what() << std::endl;
        throw excp;
    }

}

template<typename T>
constexpr Matrix3x3<T> Matrix3x3<T>::t() const
{
    return Matrix3x3<T>{ xx(), yx(), zx(),
        xy(), yy(), zy(),
        xz(), yz(), zz() };
}

template<typename T>
constexpr T Matrix3x3<T>::trace() const
{
    return xx() + yy() + zz();
}

template<typename T>
constexpr T Matrix3x3<T>::det() const
{
    return xx() * yy() * zz() +
        yx() * zy() * xz() +
        xy() * yz() * zx() -
        zx() * xz() * yy() -
        zy() * yz() * xx() -
        yx() * xy() * zz();
}

template<typename T>
constexpr bool Matrix3x3<T>::is_orthogonal() const
{
    // A * transpose(A) = eye(3)
    constexpr auto tol{ 1e3 };

    const auto temp = *this * this->t();
    return eq_fp<T>(temp.xx(), 1., tol) && eq_fp<double>(temp.xy(), 0., tol) && eq_fp<double>(temp.xz(), 0., tol) &&
        eq_fp<T>(temp.yx(), 0., tol) && eq_fp<double>(temp.yy(), 1., tol) && eq_fp<double>(temp.yz(), 0., tol) &&
        eq_fp<T>(temp.zx(), 0., tol) && eq_fp<double>(temp.zy(), 0., tol) && eq_fp<double>(temp.zz(), 1., tol);

}

template<typename T>
constexpr bool Matrix3x3<T>::is_rotmat() const
{
    // rotation Matrix3x3<T> == "special orthogonal Matrix3x3<T>" == [orthogonal Matrix3x3<T> with det == +1.0]
    return is_orthogonal() && (det() > eps);
}

template<typename T>
constexpr bool Matrix3x3<T>::is_eye() const
{
    constexpr auto tol{ 1e3 };
    return eq_fp<T>(xy(), 0., tol) &&
        eq_fp<T>(yx(), 0., tol) &&
        eq_fp<T>(xz(), 0., tol) &&
        eq_fp<T>(zx(), 0., tol) &&
        eq_fp<T>(yz(), 0., tol) &&
        eq_fp<T>(zy(), 0., tol) &&
        eq_fp<T>(xx(), 1., tol) &&
        eq_fp<T>(yy(), 1., tol) &&
        eq_fp<T>(zz(), 1., tol);
}

template<typename T>
constexpr Matrix3x3<T> Matrix3x3<T>::eye()
{
    return Matrix3x3<T>{ 1., 0., 0.,
                      0., 1., 0.,
                      0., 0., 1. };
}

template<typename T>
constexpr Matrix3x3<T> Matrix3x3<T>::zeros()
{
    return Matrix3x3<T>{};
}

template<typename T>
constexpr Matrix3x3<T> Matrix3x3<T>::ones()
{
    return Matrix3x3<T>{ 1. };
}

template<typename T>
constexpr Matrix3x3<T> Matrix3x3<T>::columns(const Vector3d<T> &col_x,
                                              const Vector3d<T> &col_y,
                                              const Vector3d<T> &col_z)
{
    return Matrix3x3<T>{ col_x.x(), col_y.x(), col_z.x(),
        col_x.y(), col_y.y(), col_z.y(),
        col_x.z(), col_y.z(), col_z.z(), };
}

template<typename T>
constexpr const T& Matrix3x3<T>::operator()(std::size_t i, std::size_t j) const //noexcept(false)
{
    //if ((i < 3u) && (j < 3u)) {
    return (*this)[i + 3u * j];
    //}
    //else {
    //    lion_exception excp("naiveDG::Matrix3x3::operator(): out of bounds.");
    //    std::cerr << excp.what() << std::endl;
    //    throw excp;
    //}
}

template<typename T>
constexpr T& Matrix3x3<T>::operator()(std::size_t i, std::size_t j) //noexcept(false)
{
    return (*this)[i + 3u * j];
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator+=(const Matrix3x3<T> &other)
{
    auto o_ij = other.cbegin();
    for (auto &&m_ij : *this) {
        m_ij += *o_ij++;
    }
    return *this;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator-=(const Matrix3x3<T> &other)
{
    auto o_ij = other.cbegin();
    for (auto &&m_ij : *this) {
        m_ij -= *o_ij++;
    }
    return *this;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator*=(const Matrix3x3<T> &other)
{
    return *this = (*this) * other;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator=(T val)
{
    for (auto &&m_ij : *this) {
        m_ij = val;
    }
    return *this;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator+=(T val)
{
    for (auto &&m_ij : *this) {
        m_ij += val;
    }
    return *this;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator-=(T val)
{
    for (auto &&m_ij : *this) {
        m_ij -= val;
    }
    return *this;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator*=(T val)
{
    for (auto &&m_ij : *this) {
        m_ij *= val;
    }
    return *this;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator/=(T val)
{
    const auto inv_val{ 1. / val };
    for (auto &&m_ij : *this) {
        m_ij *= inv_val;
    }
    return *this;
}

template<typename T>
constexpr Matrix3x3<T>& Matrix3x3<T>::operator=(const base_type &arr)
{
    base_type::operator=(arr);
    return *this;
}

template<typename T>
constexpr Matrix3x3<T> transpose(const Matrix3x3<T> &arg)
{
    return arg.t();
}

template<typename T>
constexpr T trace(const Matrix3x3<T> &arg)
{
    return arg.trace();
}

template<typename T>
constexpr T det(const Matrix3x3<T> &arg)
{
    return arg.det();
}

template<typename T>
constexpr bool is_orthogonal(const Matrix3x3<T> &arg)
{
    return arg.is_orthogonal();
}

template<typename T>
constexpr bool is_rotmat(const Matrix3x3<T> &arg)
{
    return arg.is_rotmat();
}

template<typename T>
constexpr bool is_eye(const Matrix3x3<T> &arg)
{
    return arg.is_eye();
}

template<typename T>
constexpr Matrix3x3<T> inv(const Matrix3x3<T> &arg)
{
    const auto d = arg.det();

    if (fabs(d) <= eps) {
        std::cerr << "tr7::inv: warning, input Matrix3x3<T> is singular to working precision" << std::endl;
        return Matrix3x3<T>{ inf };
    }
    else {
        return Matrix3x3<T>{ arg.yy() * arg.zz() - arg.zy() * arg.yz(),
            arg.xz() * arg.zy() - arg.xy() * arg.zz(),
            arg.xy() * arg.yz() - arg.xz() * arg.yy(),

            arg.yz() * arg.zx() - arg.yx() * arg.zz(),
            arg.xx() * arg.zz() - arg.xz() * arg.zx(),
            arg.xz() * arg.yx() - arg.xx() * arg.yz(),

            arg.yx() * arg.zy() - arg.yy() * arg.zx(),
            arg.xy() * arg.zx() - arg.xx() * arg.zy(),
            arg.xx() * arg.yy() - arg.xy() * arg.yx() } / d;
    }

}

template<typename T>
constexpr T dot(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs)
{
    return lhs.xx() * rhs.xx() + lhs.xy() * rhs.xy() + lhs.xz() * rhs.xz() +
        lhs.yx() * rhs.yx() + lhs.yy() * rhs.yy() + lhs.yz() * rhs.yz() +
        lhs.zx() * rhs.zx() + lhs.zy() * rhs.zy() + lhs.zz() * rhs.zz();
}

template<typename T>
constexpr T sumsqr(const Matrix3x3<T> &arg)
{
    return arg.xx() * arg.xx() + arg.xy() * arg.xy() + arg.xz() * arg.xz() +
        arg.yx() * arg.yx() + arg.yy() * arg.yy() + arg.yz() * arg.yz() +
        arg.zx() * arg.zx() + arg.zy() * arg.zy() + arg.zz() * arg.zz();
}

template<typename T>
constexpr T sumabs(const Matrix3x3<T> &arg)
{
    return fabs(arg.xx()) + fabs(arg.xy()) + fabs(arg.xz()) +
        fabs(arg.yx()) + fabs(arg.yy()) + fabs(arg.yz()) +
        fabs(arg.zx()) + fabs(arg.zy()) + fabs(arg.zz());
}

template<typename T>
constexpr Matrix3x3<T> sqrt(const Matrix3x3<T> &arg)
{
    return Matrix3x3<T>{ std::sqrt(arg.xx()), std::sqrt(arg.xy()), std::sqrt(arg.xz()),
        std::sqrt(arg.yx()), std::sqrt(arg.yy()), std::sqrt(arg.yz()),
        std::sqrt(arg.zx()), std::sqrt(arg.zy()), std::sqrt(arg.zz()) };
}

template<typename T>
constexpr Matrix3x3<T> abs(const Matrix3x3<T> &arg)
{
    return Matrix3x3<T>{ fabs(arg.xx()), fabs(arg.xy()), fabs(arg.xz()),
        fabs(arg.yx()), fabs(arg.yy()), fabs(arg.yz()),
        fabs(arg.zx()), fabs(arg.zy()), fabs(arg.zz()) };
}

template<typename T>
constexpr Matrix3x3<T> operator+(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs)
{
    return Matrix3x3<T>{ lhs.xx() + rhs.xx(), lhs.xy() + rhs.xy(), lhs.xz() + rhs.xz(),
        lhs.yx() + rhs.yx(), lhs.yy() + rhs.yy(), lhs.yz() + rhs.yz(),
        lhs.zx() + rhs.zx(), lhs.zy() + rhs.zy(), lhs.zz() + rhs.zz() };
}

template<typename T>
constexpr Matrix3x3<T> operator-(const Matrix3x3<T> &rhs)
{
    return Matrix3x3<T>{ -rhs.xx(), -rhs.xy(), -rhs.xz(),
        -rhs.yx(), -rhs.yy(), -rhs.yz(),
        -rhs.zx(), -rhs.zy(), -rhs.zz() };
}

template<typename T>
constexpr Matrix3x3<T> operator-(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs)
{
    return Matrix3x3<T>{ lhs.xx() - rhs.xx(), lhs.xy() - rhs.xy(), lhs.xz() - rhs.xz(),
        lhs.yx() - rhs.yx(), lhs.yy() - rhs.yy(), lhs.yz() - rhs.yz(),
        lhs.zx() - rhs.zx(), lhs.zy() - rhs.zy(), lhs.zz() - rhs.zz() };
}

template<typename T>
constexpr Matrix3x3<T> operator*(const Matrix3x3<T> &lhs, const Matrix3x3<T> &rhs)
{
    return Matrix3x3<T>{ lhs.xx() * rhs.xx() + lhs.xy() * rhs.yx() + lhs.xz() * rhs.zx(),
        lhs.xx() * rhs.xy() + lhs.xy() * rhs.yy() + lhs.xz() * rhs.zy(),
        lhs.xx() * rhs.xz() + lhs.xy() * rhs.yz() + lhs.xz() * rhs.zz(),

        lhs.yx() * rhs.xx() + lhs.yy() * rhs.yx() + lhs.yz() * rhs.zx(),
        lhs.yx() * rhs.xy() + lhs.yy() * rhs.yy() + lhs.yz() * rhs.zy(),
        lhs.yx() * rhs.xz() + lhs.yy() * rhs.yz() + lhs.yz() * rhs.zz(),

        lhs.zx() * rhs.xx() + lhs.zy() * rhs.yx() + lhs.zz() * rhs.zx(),
        lhs.zx() * rhs.xy() + lhs.zy() * rhs.yy() + lhs.zz() * rhs.zy(),
        lhs.zx() * rhs.xz() + lhs.zy() * rhs.yz() + lhs.zz() * rhs.zz() };
}

template<typename T>
constexpr Matrix3x3<T> operator+(const Matrix3x3<T> &lhs, T val)
{
    return Matrix3x3<T>{ lhs.xx() + val, lhs.xy() + val, lhs.xz() + val,
        lhs.yx() + val, lhs.yy() + val, lhs.yz() + val,
        lhs.zx() + val, lhs.zy() + val, lhs.zz() + val };
}

template<typename T>
constexpr Matrix3x3<T> operator+(T val, const Matrix3x3<T> &rhs)
{
    return Matrix3x3<T>{ val + rhs.xx(), val + rhs.xy(), val + rhs.xz(),
        val + rhs.yx(), val + rhs.yy(), val + rhs.yz(),
        val + rhs.zx(), val + rhs.zy(), val + rhs.zz() };
}

template<typename T>
constexpr Matrix3x3<T> operator-(const Matrix3x3<T> &lhs, T val)
{
    return Matrix3x3<T>{ lhs.xx() - val, lhs.xy() - val, lhs.xz() - val,
        lhs.yx() - val, lhs.yy() - val, lhs.yz() - val,
        lhs.zx() - val, lhs.zy() - val, lhs.zz() - val };
}

template<typename T>
constexpr Matrix3x3<T> operator-(T val, const Matrix3x3<T> &rhs)
{
    return Matrix3x3<T>{ val - rhs.xx(), val - rhs.xy(), val - rhs.xz(),
        val - rhs.yx(), val - rhs.yy(), val - rhs.yz(),
        val - rhs.zx(), val - rhs.zy(), val - rhs.zz() };
}

template<typename T>
constexpr Matrix3x3<T> operator*(const Matrix3x3<T> &lhs, T val)
{
    return Matrix3x3<T>{ lhs.xx() * val, lhs.xy() * val, lhs.xz() * val,
        lhs.yx() * val, lhs.yy() * val, lhs.yz() * val,
        lhs.zx() * val, lhs.zy() * val, lhs.zz() * val };
}

template<typename T>
constexpr Matrix3x3<T> operator*(T val, const Matrix3x3<T> &rhs)
{
    return Matrix3x3<T>{ val * rhs.xx(), val * rhs.xy(), val * rhs.xz(),
        val * rhs.yx(), val * rhs.yy(), val * rhs.yz(),
        val * rhs.zx(), val * rhs.zy(), val * rhs.zz() };
}

template<typename T>
constexpr Matrix3x3<T> operator/(const Matrix3x3<T> &lhs, T val)
{
    return Matrix3x3<T>{ lhs.xx() / val, lhs.xy() / val, lhs.xz() / val,
        lhs.yx() / val, lhs.yy() / val, lhs.yz() / val,
        lhs.zx() / val, lhs.zy() / val, lhs.zz() / val };
}

template<typename T>
inline std::ostream& operator<<(std::ostream &os, const Matrix3x3<T> &m)
{
    os << m.xx() << ' ' << m.xy() << ' ' << m.xz() << '\n'
       << m.yx() << ' ' << m.yy() << ' ' << m.yz() << '\n'
       << m.zx() << ' ' << m.zy() << ' ' << m.zz();
    return os;
}

template<typename T>
inline std::istream &operator>>(std::istream &is, Matrix3x3<T> &m) 
{
    // (reads in col-major format)
    T xx{ 0. }, yx{ 0. }, zx{ 0. }, 
           xy{ 0. }, yy{ 0. }, zy{ 0. }, 
           xz{ 0. }, yz{ 0. }, zz{ 0. };
    is >> xx >> yx >> zx >> xy >> yy >> zy >> xz >> yz >> zz;

    m = Matrix3x3<T>{ xx, xy, xz,
                   yx, yy, yz,
                   zx, zy, zz };
    return is;
}


template<typename T>
constexpr Vector3d<T> operator*(const Matrix3x3<T> &m, const Vector3d<T> &v)
{
    return Vector3d<T>{ m.xx() * v.x() + m.xy() * v.y() + m.xz() * v.z(),
        m.yx() * v.x() + m.yy() * v.y() + m.yz() * v.z(),
        m.zx() * v.x() + m.zy() * v.y() + m.zz() * v.z() };
}

template<typename T>
constexpr Vector3d<T> operator*(const Vector3d<T> &v, const Matrix3x3<T> &m)
{
    return m.t() * v;
}

template<typename T>
inline Vector3d<T> linsolve(const Matrix3x3<T> &m_lhs, const Vector3d<T> &v)
{
    const auto d = m_lhs.det();
    if (fabs(d) <= eps) {
        std::cerr << "tr7::linsolve:: warning, input Matrix3x3<T> may "
            "be singular or badly scaled" << std::endl;
        return Vector3d<T>{ inf, inf, inf };
    }

    return Vector3d<T>{ v.x() * m_lhs.yy() * m_lhs.zz() +
        v.y() * m_lhs.zy() * m_lhs.xz() +
        m_lhs.xy() * m_lhs.yz() * v.z() -
        v.z() * m_lhs.xz() * m_lhs.yy() -
        m_lhs.zy() * m_lhs.yz() * v.x() -
        v.y() * m_lhs.xy() * m_lhs.zz(),

        m_lhs.xx() * v.y() * m_lhs.zz() +
        m_lhs.yx() * v.z() * m_lhs.xz() +
        v.x() * m_lhs.yz() * m_lhs.zx() -
        m_lhs.zx() * m_lhs.xz() * v.y() -
        v.z() * m_lhs.yz() * m_lhs.xx() -
        m_lhs.yx() * v.x() * m_lhs.zz(),

        m_lhs.xx() * m_lhs.yy() * v.z() +
        m_lhs.yx() * m_lhs.zy() * v.x() +
        m_lhs.xy() * v.y() * m_lhs.zx() -
        m_lhs.zx() * v.x() * m_lhs.yy() -
        m_lhs.zy() * v.y() * m_lhs.xx() -
        m_lhs.yx() * m_lhs.xy() * v.z() } / d;

}

template<typename T>
inline Vector3d<T> linsolve(const Vector3d<T> &v, const Matrix3x3<T> &m_rhs)
{
    return linsolve(m_rhs.t(), v);
}

template<typename T>
constexpr Matrix3x3<T> outer(const Vector3d<T> &v_lhs, const Vector3d<T> &v_rhs)
{
    return Matrix3x3<T>{ v_lhs.x() * v_rhs.x(), v_lhs.x() * v_rhs.y(), v_lhs.x() * v_rhs.z(),
        v_lhs.y() * v_rhs.x(), v_lhs.y() * v_rhs.y(), v_lhs.y() * v_rhs.z(),
        v_lhs.z() * v_rhs.x(), v_lhs.z() * v_rhs.y(), v_lhs.z() * v_rhs.z() };
}

template<typename T>
constexpr Matrix3x3<T> diag(const Vector3d<T> &v)
{
    return Matrix3x3<T>{ v.x(), 0., 0.,
                      0., v.y(), 0.,
                      0., 0., v.z() };
}

template<typename T>
constexpr Matrix3x3<T> crossmat(const Vector3d<T> &v)
{
    return Matrix3x3<T>{ 0., -v.z(), v.y(),
                      v.z(), 0., -v.x(),
                      -v.y(), v.x(), 0. };
}

