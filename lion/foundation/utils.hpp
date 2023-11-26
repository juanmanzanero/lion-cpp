#ifndef LION_FOUNDATION_UTILS_HPP
#define LION_FOUNDATION_UTILS_HPP
#pragma once


#include <regex>
#include <iostream>

#include "lion/foundation/types.h"
#include "lion/foundation/constants.h"
#include "lion/foundation/utils.h"
#include "lion/math/vector3d.hpp"


constexpr double& Value(double& val) { return val; }
constexpr const double& Value(const double& val) { return val; }


template <typename T0, typename T1, typename T2>
constexpr bool in_or_on_interval(const T0 &x, const T1 &lo, const T2 &hi)
{
    //
    // Returns true if "x" is inside or on the extremes of
    // the closed interval "[lo, hi]", false otherwise.
    //

    return x >= lo && x <= hi;
}


template <typename T> 
constexpr T sign(const T &x)
{
    //
    // Strict sign function:
    // if x == 0, sign(x) = 0,
    // elseif x > 0, sign(x) = 1,
    // elseif x < 0 , sign(x) = -1.
    //

    return static_cast<T>((T{ 0 } < x) - (x < T{ 0 }));
}

template<typename T>
constexpr bool samesign(T x, T y)
{
    //
    // Returns true if "x" and "y" have the same sign (0 considered positive).
    //

    return (x >= T{ 0 }) == (y >= T{ 0 });
}

template<typename T>
constexpr bool samesign(T x, T y, T z)
{
    //
    // Returns true if "x", "y" and "z" have the same sign (0 considered positive).
    //

    const auto sx = x >= T{ 0 };

    return (sx == (y >= T{ 0 })) && (sx == (z >= T{ 0 }));
}


template<typename T, typename TolType>
constexpr bool eq_fp(T x, T y, TolType tol)
{
    //
    // Kopriva's AlmostEqual routine,
    // tests equality of two floating-point
    // numbers. Input tol is a multiple of
    // "std::numeric_limits<T>::epsilon()".
    //

    using std::abs;

    const auto teps = tol * std::numeric_limits<TolType>::epsilon();
    const auto diff = abs(x - y);
    const auto ax = abs(x);
    const auto ay = abs(y);

    if (ax <= teps || ay <= teps) {
        return diff <= TolType{ 2 } * teps;
    }
    else {
        return diff <= teps * ax && diff <= teps * ay;
    }
}


template<typename T>
constexpr T mod_mat(T x, T y)
{
    // Matlab's "mod" function, should only be called for floating-point variables. Cannot be
    // constexpr since it uses std::fmod

    // NEW VERSION FROM MATLAB CODE GENERATION
    auto ret{ x };
    if (std::isfinite(x) && std::isfinite(y)) {
        if (eq_fp(x, static_cast<T>(0.))) {
            ret = static_cast<T>(0.);
        }
        else {
            if (!eq_fp(y, static_cast<T>(0.))) {
                ret = std::fmod(x, y);
                bool rEQ0{ eq_fp(ret, static_cast<T>(0.)) }; 
                if ((!rEQ0) && (y > std::floor(y))) {
                    auto q{ fabs(x / y) }; 
                    rEQ0 = (fabs(q - std::floor(q + static_cast<T>(0.5))) <= static_cast<T>(2.2204460492503131E-16) * q); // 2.220446049250313E-16 == eps, but we'll leave it like this since Matlab hardcodes it
                }

                if (rEQ0) {
                    ret = static_cast<T>(0.);
                }
                else {
                    if ((x < static_cast<T>(0.)) != (y < static_cast<T>(0.))) {
                        ret += y;
                    }
                }
            }
        }
    }
    else {
        if (!eq_fp(y, static_cast<T>(0.))) {
            ret = nan_T<T>;
        }
    }

    return ret;
}


template<typename T>
constexpr T wrap_to_pi(T angle_rad)
{
    //
    // Wraps the input angle, in radians, to the
    // open interval [-pi, pi).
    //

    // return mod_mat(angle_rad + pi_T<T>, T{ 2 } * pi_T<T>) - pi_T<T>;
    return angle_rad - T{ 2 } * pi_T<T> *
        std::floor(T{ 0.5 } * angle_rad / pi_T<T> + T{ 0.5 });
}

template<typename T>
constexpr T wrap_to_2pi(T angle_rad)
{
    //
    // Wraps the input angle, in radians, to the
    // open interval [0, 2 * pi).
    //

    // return mod_mat(ang, T{ 2 } * pi_T<T>);
    return angle_rad - T{ 2 } * pi_T<T> *
        std::floor(T{ 0.5 } * angle_rad / pi_T<T>);
}

template<typename T>
constexpr T wrap_to_180(T angle_deg)
{
    //
    // Wraps the input angle, in degrees, to the
    // open interval [-180, 180).
    //

    // return mod_mat(angle_deg + T{ 180 }, T{ 360 }) - T{ 180 };
    return angle_deg - T{ 360 } *
        std::floor(angle_deg / T{ 360 } + T{ 0.5 });
}

template<typename T>
constexpr T wrap_to_360(T angle_deg)
{
    //
    // Wraps the input angle, in degrees, to the
    // open interval [0, 360).
    //

    // return mod_mat(angle_deg, T{ 360 });
    return angle_deg - T{ 360 } *
        std::floor(angle_deg / T{ 360 });
}


template<bool ActuallySmooth, typename T, typename S>
constexpr T smooth_pos(const T &x, S eps2)
{
    //
    // Returns the positive part of input "x", smoothed
    // out near 0 with parameter "eps2" (unless template
    // parameter "ActuallySmooth" is false, in that case
    // the function returns the strict positive part
    // of "x").
    //

    if constexpr (ActuallySmooth) {
        using std::sqrt;

        return S{ 0.5 } * (x + sqrt(x * x + eps2));
    }
    else {
        (void)eps2;

        return x > T{ 0 } ? x : T{ 0 };
    }
}


template<bool ActuallySmooth, typename T, typename S>
constexpr T smooth_neg(const T &x, S eps2)
{
    //
    // Returns the negative part of input "x", smoothed
    // out near 0 with parameter "eps2" (unless template
    // parameter "ActuallySmooth" is false, in that case
    // the function returns the strict negative part
    // of "x").
    //

    return -smooth_pos<ActuallySmooth>(-x, eps2);
}


template<bool ActuallySmooth, bool UseSinAtanFormula, typename T, typename S>
constexpr T smooth_sign(const T &x, S eps)
{
    //
    // Returns the sign of input "x" (where sign(0) = 0),
    // smoothed out near 0 with parameter "eps" (unless
    // template parameter "ActuallySmooth" is false, in that
    // case the function is the strict sign of "x").
    //

    if constexpr (ActuallySmooth) {
        if constexpr (UseSinAtanFormula) {

            return x / smooth_abs(x, eps * eps);
        }
        else {
            using std::tanh;

            return tanh(x / eps);
        }
    }
    else {
        (void)eps;

        return sign(x);
    }
}


template<bool ActuallySmooth, bool UseSinAtanFormula, typename T, typename S>
constexpr T smooth_sign_derivative(const T &x, S eps)
{
    //
    // Returns the sign derivative of input "x" (where sign(0) = 0),
    // smoothed out near 0 with parameter "eps" (unless
    // template parameter "ActuallySmooth" is false, in that
    // case the function is the strict sign of "x").
    //

    if constexpr (ActuallySmooth) {
        if constexpr (UseSinAtanFormula) {
            using std::sqrt;

            return eps * eps / ((eps * eps + x * x) * sqrt(eps * eps + x * x));
        }
        else {
            using std::tanh;

            const auto tanh_x_div_eps = tanh(x / eps);
            return -(tanh_x_div_eps * tanh_x_div_eps - 1.0) / eps;
        }
    }
    else {
        (void)eps;

        return (x >= 0.0 ? 1.0 : -1.0);
    }
}



template<bool ActuallySmooth, typename T, typename S>
constexpr T smooth_abs(const T &x, S eps2)
{
    //
    // Returns the absolute value of input "x", smoothed
    // out near 0 with parameter "eps2" (unless template
    // parameter "ActuallySmooth" is false, in that case
    // the function returns the strict absolute value
    // of "x").
    //

    if constexpr (ActuallySmooth) {
        using std::sqrt;

        return sqrt(x * x + eps2);
    }
    else {
        using std::abs;

        (void)eps2;

        return abs(x);
    }
}


template<bool ActuallySmooth, typename T, typename S>
constexpr T smooth_hypot(const T &x, const T &y, S eps2)
{
    //
    // Returns "sqrt(x^2 + y^2)", smoothed out near [0, 0]
    // with parameter "eps2" (unless template parameter
    // "ActuallySmooth" is false, in that case the function
    // returns the strict hypot of "[x, y]").
    //

    if constexpr (ActuallySmooth) {
        using std::sqrt;

        return sqrt(x * x + y * y + eps2);
    }
    else {
        //using std::hypot; -> not implemented in CppAD
        using std::sqrt;

        (void)eps2;

        //return hypot(x, y);
        return sqrt(x * x + y * y);
    }
}

template<bool ActuallySmooth, typename T, typename T1, typename S>
constexpr T smooth_max(const T& x, const T1 &lo, S eps2)
{
    //
    // Clamps "x" to "[lo, inf)", using a squared root smooth modulation,
    // unless template parameter "ActuallySmooth" is false, in
    // that case we strictly clamp the value.
    //

    if constexpr (ActuallySmooth) {
        using std::sqrt;

        const auto Dxlo = x - lo;
        return S{ 0.5 } * (Dxlo + sqrt(Dxlo * Dxlo + eps2)) + lo;
    }
    else {
        return x > lo ? x : static_cast<T>(lo);
    }
}

template<bool ActuallySmooth, typename T, typename T1, typename S>
constexpr T smooth_min(const T& x, const T1 &hi, S eps2)
{
    //
    // Clamps "x" to "(-inf, hi]", using a squared root smooth modulation,
    // unless template parameter "ActuallySmooth" is false, in
    // that case we strictly clamp the value.
    //

    if constexpr (ActuallySmooth) {
        using std::sqrt;

        const auto Dxhi = x - hi;
        return S{ 0.5 } * (Dxhi - sqrt(Dxhi * Dxhi + eps2)) + hi;
    }
    else {
        return x < hi ? x : static_cast<T>(hi);
    }
}

template<bool ActuallySmooth, typename T, typename T1, typename T2, typename S>
constexpr T smooth_clamp(const T &x, const T1 &lo, const T2 &hi, S eps2)
{
    //
    // Clamps "x" to "[x_lower, x_upper]", using a squared root smooth modulation,
    // unless template parameter "ActuallySmooth" is false, in
    // that case we strictly clamp the value.
    //

    if constexpr (ActuallySmooth) {
        using std::sqrt;

        const auto Dxlo = x - lo;
        const auto Dxhi = x - hi;
        return S{ 0.5 } * (sqrt(Dxlo * Dxlo + eps2) - sqrt(Dxhi * Dxhi + eps2) + lo + hi);
    }
    else {
        return x < hi ? (x > lo ? x : static_cast<T>(lo)) : static_cast<T>(hi);
    }
}

template<bool ActuallySmooth, bool UseSinAtanFormula, typename T, typename S>
constexpr T smooth_step(const T &x, S eps)
{
    //
    // Implements the "smooth step" function (i.e., a smoothed
    // version of the Heavyside step function, with parameter "eps"),
    // unless template parameter "ActuallySmooth" is false, in
    // that case it just evaluates the Heavyside step function.
    //

    if constexpr (ActuallySmooth) {
        return S{ 0.5 } * (S{ 1 } + smooth_sign<ActuallySmooth, UseSinAtanFormula>(x, eps));
    }
    else {
        return x >= S{ 0 } ? T{ 1 } : T{ 0 };
    }
}


template <typename T>
constexpr std::vector<T> linspace(T lo, T hi, std::size_t num_points)
{
    //
    // Returns a vector of size "num_points", holding
    // linearly equally spaced points in the (closed)
    // interval "[lo, hi]".
    //

    static_assert(!std::is_integral_v<T>, "linspace:"
        " the input interval limits should not be"
        " of an integral type.");

    std::vector<T> xs(num_points);
    if (num_points != 0u) {
        // use Matlab's formula rigorously, so that
        // we get the exact same numbers
        const auto nm1 = num_points - 1u;
        const auto extent = hi - lo;
        for (auto count = 1u; count < nm1; ++count) {
            xs[count] = lo + (count * extent) / nm1;
        }

        xs.front() = lo;
        xs.back() = hi;
    }

    return xs;
}


template<typename T>
constexpr std::vector<T> iota(T lo, std::size_t num_points, T increment)
{
    //
    // Returns a vector of size "num_points" whose first
    // element is equal to "lo" and its rest of elements
    // are sequentially increased by "increment".
    //

    std::vector<T> ret(num_points);
    for (auto &&r : ret) {
        r = lo;
        lo += increment;
    }

    return ret;
}


template<typename T>
inline T trapz(const std::vector<T>& x, const std::vector<T>& y)
{
    assert(x.size() == y.size());

    T val{ 0.0 };

    for (size_t i = 1; i < x.size(); ++i)
    {
        val += T{ 0.5 } * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }

    return val;
}


template<typename ArrayType, typename ValueType>
constexpr ValueType sumabs(const ArrayType &x)
{
    //
    // Returns the sum of the absolute values of the
    // elements in the input array.
    //

    return std::accumulate(x.cbegin(), x.cend(),
        ValueType{ 0 },
        [](auto accum, auto xi) { return accum + std::abs(xi); });
}


template<typename ArrayType, typename ValueType>
constexpr ValueType sumsqr(const ArrayType &x)
{
    //
    // Returns the sum of the squares of the
    // elements in the input array.
    //

    return std::accumulate(x.cbegin(), x.cend(),
        ValueType{ 0 },
        [](auto accum, auto xi) { return accum + xi * xi; });
}


template<typename T>
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_closest_point(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, bool closed, const size_t i_start, const scalar maximum_distance)
{
    // Get number of points
    const auto N = xy_polygon.size();

    // (1) do a search along the poly vertices
    size_t i_closest_poly    = i_start;
    sVector3d x_closest_poly = xy_polygon[i_start];
    T d2_closest_poly        = dot(x_closest_poly-x0,x_closest_poly-x0);
    T d2_current             = d2_closest_poly;

    scalar current_distance = 0.0;

    for (size_t i = i_start+1; i < N; ++i)
    {
        current_distance += norm(xy_polygon[i]-xy_polygon[i-1]);

        if ( (current_distance <= maximum_distance) || (i < i_start + 5) )
        {
            d2_current = dot(xy_polygon[i]-x0,xy_polygon[i]-x0);
            if ( d2_current < d2_closest_poly )
            {
                i_closest_poly = i;
                d2_closest_poly = d2_current;
                x_closest_poly = xy_polygon[i];
            }
        }
    }
    
    T d2_closest = d2_closest_poly;
    Vector3d<T> x_closest = Vector3d<T>(x_closest_poly);
    std::array<size_t,2> i_closest = {i_closest_poly, i_closest_poly};
    
    // (2) compute the minimum distance to the next forward face
    bool skip_calculation = false;
    size_t i_next_fwd;
    if ( i_closest_poly == N - 1)
    {
        if ( closed )
            i_next_fwd = 0;
        else
            skip_calculation = true;
    }
    else
        i_next_fwd = i_closest_poly + 1;

    if ( !skip_calculation )
    {
        const sVector3d x_next_fwd = xy_polygon[i_next_fwd];
        sVector3d p = x_next_fwd - x_closest_poly;
    
        T s_next_fwd_closest = dot(p,x0-x_closest_poly)/dot(p,p);
        
        if ( (s_next_fwd_closest > 0.0) && (s_next_fwd_closest < 1.0) )
        {
            const auto x_next_fwd_closest = x_closest_poly + p*s_next_fwd_closest;
            const auto d2_next_fwd_closest = dot(x_next_fwd_closest - x0,x_next_fwd_closest - x0);
            if ( d2_next_fwd_closest < d2_closest_poly )
            {
                x_closest = x_next_fwd_closest;
                d2_closest = d2_next_fwd_closest;
                i_closest = {i_closest_poly, i_next_fwd};
            }
        }
    }
        
    // (3) compute the minimum distance to the next backward face
    size_t    i_next_bwd;
    skip_calculation = false;
    if ( i_closest_poly == 0 )
    {
        if ( closed ) 
            i_next_bwd = xy_polygon.size()-1;
        else
            skip_calculation = true;
    }
    else
        i_next_bwd = i_closest_poly - 1;
    
    if ( !skip_calculation )
    {
        const sVector3d x_next_bwd = xy_polygon[i_next_bwd];
        const sVector3d p = x_next_bwd - x_closest_poly;
        T s_next_bwd_closest = dot(p,x0-x_closest_poly)/dot(p,p);
    
        if ( (s_next_bwd_closest > 0.0) && (s_next_bwd_closest < 1.0) )
        {
            const auto x_next_bwd_closest = x_closest_poly + p*s_next_bwd_closest;
            const auto d2_next_bwd_closest = dot(x_next_bwd_closest - x0, x_next_bwd_closest - x0);
            if ( d2_next_bwd_closest < d2_closest )
            {
                x_closest = x_next_bwd_closest;
                d2_closest = d2_next_bwd_closest;
                i_closest = {i_next_bwd, i_closest_poly};
            }
        }
    }

    return {x_closest,d2_closest,i_closest};
}


template<typename T>
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_intersection(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, const Vector3d<T>& p0, bool closed)
{
    Vector3d<T> v_intersection = {0.0, 0.0, 0.0};
    T            d_intersection2 = 1.0e18;
    std::array<size_t,2> i_intersection = {0, 0};

    for (size_t i = 0; i < xy_polygon.size(); ++i)
    {
        const sVector3d& r1 = xy_polygon[i];
        const sVector3d p1 = xy_polygon[i+1]-xy_polygon[i];

        auto [s0,s1] = find_intersection(x0, p0, r1, p1);

        // We only accept solutions with s > 0, and within [xy_polygon[i], xy_polygon[i+1]]
        if ( s0 < 0.0 || s1 < 0.0 || s1 > 1.0)
            continue;
        else if ( s0*s0 < d_intersection2 )
        {
            v_intersection = r1 + p1*s1;
            d_intersection2 = s0*s0;
            i_intersection = {i,i+1}; 
        }
    } 

    if ( closed )
    {
        const sVector3d& r1 = xy_polygon.back();
        const sVector3d p1 = xy_polygon.front()-xy_polygon.back();

        auto [s0,s1] = find_intersection(x0, p0, r1, p1);

        // We only accept solutions with s > 0, and within [xy_polygon[i], xy_polygon[i+1]]
        if ( s0 < 0.0 || s1 < 0.0 || s1 > 1.0)
        {
        }
        else if ( s0*s0 < d_intersection2 )
        {
            v_intersection = r1 + p1*s1;
            d_intersection2 = s0*s0;
            i_intersection = {xy_polygon.size()-1,0}; 
        }
    }

#ifndef NDEBUG
    if ( (i_intersection.front() == 0) && (i_intersection.back() == 0) )
        throw lion_exception("The straight given does not intersect the polygon");
#endif

    return {v_intersection, d_intersection2, i_intersection};
}

template<typename T>
inline std::tuple<T,T> find_intersection(const Vector3d<T>& r0, const Vector3d<T>& p0, const sVector3d& r1, const sVector3d& p1)
{
    // (1) Check if the straights are parallel 
    const T detA = -p0[0]*p1[1] + p0[1]*p1[0];
    if ( std::abs(detA) < 1.0e-10 )
        return {1.0e18, 1.0e18};

    // (2) Compute the intersection otherwise
    const auto b = r1-r0;
    T s0 = (-b[0]*p1[1]+b[1]*p1[0])/detA;
    T s1 = (p0[0]*b[1] - p0[1]*b[0])/detA;

    return {s0, s1};
}


template<typename SizeType>
constexpr SizeType nchoosek(SizeType n, SizeType k)
{
    //
    // Calculates the binomial coefficient "C(n, k) = n! / ((n - k)! * k!)".
    // The value is calculated in such a way as to avoid overflow and
    // roundoff, using integer arithmetic. The code is taken from function
    // "i4_choose" in file "polynomial.cpp" at
    // "https://people.math.sc.edu/Burkardt/cpp_src/polynomial/polynomial.html"
    //

    const auto mn = std::min(k, n - k);
    if (mn < SizeType{ 0 }) {
        return SizeType{ 0 };
    }
    else if (mn == SizeType{ 0 }) {
        return SizeType{ 1 };
    }
    else {
        const auto mx = std::max(k, n - k);
        auto ret = mx + SizeType{ 1 };
        for (auto i = SizeType{ 2 }; i <= mn; ++i) {
            ret = (ret * (mx + i)) / i;
        }

        return ret;
    }
}


template<typename T, typename Array2Type>
constexpr std::pair<Array2Type, bool> sin_cos_solve(T lhs_s, T lhs_c, T rhs,
    T tolzero)
{
    //
    // Solves "lhs_s * sin(x) + lhs_c * cos(x) = rhs" for "x" (2 solutions),
    // only if the equation has real solutions. Returns an array holding both
    // solutions and an extra flag that is true if they're valid (AND real),
    // false otherwise (in that case, the output array will be equal to
    // "[nan, nan]").
    //
    // NOTE: the output solutions of this function will be within the CLOSED
    // interval "[-pi, pi]". For the solutions to be in "[-pi, pi)", one can
    // for example call function "wrap_to_pi" on them (we don't do it here so
    // that the function can be used with CppAD).
    //

    using std::abs;
    using std::sqrt;
    using std::atan2;
    using std::asin;
    using std::acos;

    if (abs(lhs_s) > tolzero && abs(lhs_c) > tolzero) {
        // lhs_s * sin(x) + lhs_c * cos(x) = d * sin(x + p) = rhs
        // d = sqrt(lhs_s^2 + lhs_c^2)
        // p = atan2(lhs_c, lhs_s)
        const auto arg = rhs / sqrt(lhs_s * lhs_s + lhs_c * lhs_c);
        if (abs(arg) <= T{ 1 }) {
            const auto p = atan2(lhs_c, lhs_s);
            const auto x0 = asin(arg);
            return { Array2Type{ x0 - p, T{ pi } - x0 - p },
                true };
        }
    }
    else if (abs(lhs_s) < tolzero && abs(lhs_c) > tolzero) {
        // lhs_c * cos(x) = rhs
        const auto arg = rhs / lhs_c;
        if (abs(arg) <= T{ 1 }) {
            const auto x0 = acos(arg);
            return { Array2Type{ x0, -x0 },
                true };
        }
    }
    else if (abs(lhs_s) > tolzero && abs(lhs_c) < tolzero) {
        // lhs_s * sin(x) = rhs
        const auto arg = rhs / lhs_s;
        if (abs(arg) <= T{ 1 }) {
            const auto x0 = asin(arg);
            return { Array2Type{ x0, T{ pi } - x0 },
                true };
        }
    }

    // if we've arrived here, the equation was either "0 = rhs"
    // or had complex roots... invalid solutions
    return { Array2Type{ std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN() },
        false };
}


template<typename ContainerOfGridVectorsType,
         typename ScalarType>
std::vector<ScalarType> grid_vectors2points_rowmaj(const ContainerOfGridVectorsType &grid_vectors)
{
    //
    // Expands a collection of "grid vectors" into a
    // matrix of "points" (returned as an std::vector
    // in row-major order), e.g., if
    //
    //   "grid_vectors = { { 0, 1 },
    //                     { 3, 4 },
    //                     { 5, 6 } }",
    //
    // then
    //
    //   "points = { 0, 3, 5,
    //               1, 3, 5,
    //               0, 4, 5,
    //               1, 4, 5,
    //               0, 3, 6,
    //               1, 3, 6,
    //               0, 4, 6,
    //               1, 4, 6 }".
    //

    const auto num_grid_vectors = grid_vectors.size();

    std::vector<std::size_t> grid_vector_sizes(num_grid_vectors);
    std::transform(grid_vectors.cbegin(), grid_vectors.cend(), grid_vector_sizes.begin(),
                   [](auto &gv) { return gv.size(); });

    std::vector<std::size_t> grid_vectors_accumulated_sizes(num_grid_vectors);
    std::exclusive_scan(grid_vector_sizes.cbegin(), grid_vector_sizes.cend(),
                        grid_vectors_accumulated_sizes.begin(), std::size_t{ 1u }, std::multiplies<>{});

    const auto total_num_points = (*grid_vectors_accumulated_sizes.rbegin()) *
                                  (*grid_vector_sizes.rbegin());

    std::vector<ScalarType> points_rowmaj(total_num_points * num_grid_vectors);
    for (auto dim = 0u; dim < num_grid_vectors; ++dim) {
        for (auto p = std::size_t{ 0u }; p < total_num_points; ) {
            for (const auto &grid_value : grid_vectors[dim]) {
                for (auto rep = 0u; rep < grid_vectors_accumulated_sizes[dim]; ++rep) {
                    points_rowmaj[num_grid_vectors * (p + rep) + dim] = grid_value;
                }
                p += grid_vectors_accumulated_sizes[dim];
            }
        }
    }
    return points_rowmaj;
}


template<typename Container, typename ValueType>
constexpr ValueType median(Container cont)
{
    //
    // Returns the median of the values stored in a container
    // (NOTE: the container is copied in order to avoid
    // altering its order).
    //

    const auto size = cont.size();
    if (size == 0u) {
      return std::numeric_limits<ValueType>::max();
    }

    const auto middle = std::next(cont.begin(), size >> 1u);
    std::nth_element(cont.begin(), middle, cont.end());

    if (size % 2u) {
      return *middle;
    }
    else {
      return (*middle + *std::max_element(cont.begin(), middle)) / ValueType{ 2 };
    }
}


template<typename It, typename ValueType>
constexpr It nearest_in_sorted_range(It first, It last, ValueType val)
{
    //
    // Returns an iterator to the element in the SORTED
    // range "[first, last)" that is closest to the input
    // value "val", calculated with a binary search.
    // Distance ties are rounded up. The result is undefined
    // if the range is empty, and incoprrect if it is not
    // sorted.
    //

    using std::abs;

    const auto nearest = std::lower_bound(first, last, val);
    if (nearest == last ||
        (nearest != first && abs(*std::prev(nearest) - val) < abs(*nearest - val))) {

        return std::prev(nearest);
    }
    else {
        return nearest;
    }
}

#endif
