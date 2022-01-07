//
// tr7::Vector3d<T> implementation header, #included by "tr7/base/vector3d.h"
//
template<typename T>
constexpr Vector3d<T>::Vector3d(T x, T y, T z) :
    base_type{ x, y, z } {}

template<typename T>
constexpr Vector3d<T>::Vector3d(T val) :
    base_type{ val, val, val } {}

template<typename T>
constexpr Vector3d<T>::Vector3d(const base_type &arr) :
    base_type(arr) {}

template<typename T>
inline Vector3d<T>::Vector3d(const std::vector<T> &v) noexcept(false)
{
    const bool size_ok{ v.size() == Vector3d<T>::size };
    if (size_ok) {
        *this = Vector3d<T>{ v[0], v[1], v[2] };
    }
    else {
        std::runtime_error excp("tr7::Vector3d<T>::Vector3d<T>: invalid argument, the input std::vector is not of size 3.");
        std::cerr << excp.what() << std::endl;
        throw excp;
    }
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::normalize()
{
    return *this /= norm();
}

template<typename T>
constexpr T Vector3d<T>::norm() const
{
    return std::sqrt(x() * x() + y() * y() + z() * z());
}

template<typename T>
constexpr T Vector3d<T>::safe_norm(T tol) const
{
    return std::max(norm(), std::abs(tol));
}

template<typename T>
constexpr Vector3d<T> Vector3d<T>::zeros()
{
    return Vector3d<T>{};
}

template<typename T>
constexpr Vector3d<T> Vector3d<T>::ones()
{
    return Vector3d<T>{ 1. };
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator+=(const Vector3d<T> &other)
{
    x() += other.x();
    y() += other.y();
    z() += other.z();
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator-=(const Vector3d<T> &other)
{
    x() -= other.x();
    y() -= other.y();
    z() -= other.z();
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator*=(const Vector3d<T> &other)
{
    x() *= other.x();
    y() *= other.y();
    z() *= other.z();
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator/=(const Vector3d<T> &other)
{
    x() /= other.x();
    y() /= other.y();
    z() /= other.z();
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator=(T val)
{
    x() = val;
    y() = val;
    z() = val;
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator+=(T val)
{
    x() += val;
    y() += val;
    z() += val;
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator-=(T val)
{
    x() -= val;
    y() -= val;
    z() -= val;
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator*=(T val)
{
    x() *= val;
    y() *= val;
    z() *= val;
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator/=(T val)
{
    const auto inv_val{ 1. / val };
    x() *= inv_val;
    y() *= inv_val;
    z() *= inv_val;
    return *this;
}

template<typename T>
constexpr Vector3d<T>& Vector3d<T>::operator=(const base_type &arr)
{
    x() = arr[0];
    y() = arr[1];
    z() = arr[2];
    return *this;
}


template<typename T>
constexpr Vector3d<T> normalize(Vector3d<T> arg)
{
    return arg.normalize();
}

template<typename T>
constexpr T norm(const Vector3d<T> &arg)
{
    return arg.norm();
}

template<typename T>
constexpr T angle(const Vector3d<T> &lhs, const Vector3d<T> &rhs)
{
    return wrap_to_pi<T>(2.0 *
        std::atan2(cross(lhs, rhs).norm(), dot(lhs, rhs) + lhs.safe_norm() * rhs.safe_norm()));
}

template<typename T>
constexpr T dot(const Vector3d<T> &lhs, const Vector3d<T> &rhs)
{
    return lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z();
}

template<typename U, typename V, typename W>
constexpr Vector3d<W> cross(const Vector3d<U> &lhs, const Vector3d<V> &rhs)
{
    return Vector3d<W>{ lhs.y() * rhs.z() - lhs.z() * rhs.y(),
        lhs.z() * rhs.x() - lhs.x() * rhs.z(),
        lhs.x() * rhs.y() - lhs.y() * rhs.x() };
}

template<typename T>
constexpr T distance(const Vector3d<T> &lhs, const Vector3d<T> &rhs)
{
    return norm(lhs - rhs);
}

template<typename T>
constexpr Vector3d<T> sqrt(const Vector3d<T> &arg)
{
    return Vector3d<T>{ std::sqrt(arg.x()),
        std::sqrt(arg.y()),
        std::sqrt(arg.z()) };
}

template<typename T>
constexpr Vector3d<T> abs(const Vector3d<T> &arg)
{
    return Vector3d<T>{ std::abs(arg.x()),
        std::abs(arg.y()),
        std::abs(arg.z()) };
}

template<typename T>
constexpr T sum(const Vector3d<T> &arg)
{
    return arg.x() + arg.y() + arg.z();
}

template<typename T>
constexpr Vector3d<T> operator+(const Vector3d<T> &lhs, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ lhs.x() + rhs.x(),
        lhs.y() + rhs.y(),
        lhs.z() + rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator-(const Vector3d<T> &rhs)
{
    return Vector3d<T>{ -rhs.x(),
        -rhs.y(),
        -rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator-(const Vector3d<T> &lhs, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ lhs.x() - rhs.x(),
        lhs.y() - rhs.y(),
        lhs.z() - rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator*(const Vector3d<T> &lhs, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ lhs.x() * rhs.x(),
        lhs.y() * rhs.y(),
        lhs.z() * rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator/(const Vector3d<T> &lhs, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ lhs.x() / rhs.x(),
        lhs.y() / rhs.y(),
        lhs.z() / rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator+(const Vector3d<T> &lhs, T val)
{
    return Vector3d<T>{ lhs.x() + val,
        lhs.y() + val,
        lhs.z() + val };
}

template<typename T>
constexpr Vector3d<T> operator+(T val, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ val + rhs.x(),
        val + rhs.y(),
        val + rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator-(const Vector3d<T> &lhs, T val)
{
    return Vector3d<T>{ lhs.x() - val,
        lhs.y() - val,
        lhs.z() - val };
}

template<typename T>
constexpr Vector3d<T> operator-(T val, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ val - rhs.x(),
        val - rhs.y(),
        val - rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator*(const Vector3d<T> &lhs, T val)
{
    return Vector3d<T>{ lhs.x() * val,
        lhs.y() * val,
        lhs.z() * val };
}

template<typename T>
constexpr Vector3d<T> operator*(T val, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ val * rhs.x(),
        val * rhs.y(),
        val * rhs.z() };
}

template<typename T>
constexpr Vector3d<T> operator/(const Vector3d<T> &lhs, T val)
{
    const auto inv_val{ 1. / val };
    return Vector3d<T>{ lhs.x() * inv_val,
        lhs.y() * inv_val,
        lhs.z() * inv_val };
}

template<typename T>
constexpr Vector3d<T> operator/(T val, const Vector3d<T> &rhs)
{
    return Vector3d<T>{ val / rhs.x(),
        val / rhs.y(),
        val / rhs.z() };
}

template<typename T>
inline std::ostream& operator<<(std::ostream &os, const Vector3d<T> &v)
{
    os << v.x() << ' ' << v.y() << ' ' << v.z();
    return os;
}

template<typename T>
inline std::istream& operator>>(std::istream &is, Vector3d<T> &v)
{
    is >> v.x() >> v.y() >> v.z();
    return is;
}


template<typename T>
typename std::enable_if<!std::is_same<T,scalar>::value,Vector3d<T>>::type operator+(Vector3d<T> lhs, const sVector3d& rhs)
{
    lhs[0] += rhs[0];
    lhs[1] += rhs[1];
    lhs[2] += rhs[2];

    return lhs;
}

