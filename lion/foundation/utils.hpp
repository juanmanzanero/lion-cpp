#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "utils.h"
#include "constants.h"
#include <regex>
#include<iostream>
#include "types.h"
#include "lion/math/vector3d.hpp"

constexpr double& Value(double& val) { return val; }
constexpr const double& Value(const double& val) { return val; }

template <typename T> 
constexpr int sign(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template<typename T>
constexpr bool eq_fp(T x, T y, T tol)
{
    // Kopriva's AlmostEqual routine
    const T teps = tol * eps_T<T>;
    const T diff = fabs(x - y);

    const T ax = fabs(x);
    const T ay = fabs(y);

    if (ax <= teps || ay <= teps) {
        return diff <= static_cast<T>(2.) * teps;
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
            ret = nan_T<T>::nan;
        }
    }

    return ret;
}

template<typename T>
constexpr T wrap_to_pi(T ang)
{
    // wrap to [-pi, pi)
    return mod_mat(ang + pi_T<T>, static_cast<T>(2.0) * pi_T<T>) - pi_T<T>;
}

template<typename T>
constexpr T wrap_to_2pi(T ang)
{
    // wrap to [0, 2 * pi)
    return mod_mat(ang, static_cast<T>(2.0) * pi_T<T>);
}

template<typename T>
constexpr T wrap_to_180(T ang)
{
    // wrap to [-180, 180)
    return mod_mat(ang + static_cast<T>(180.0), static_cast<T>(360.0)) - static_cast<T>(180.);
}

template<typename T>
constexpr T wrap_to_360(T ang)
{
    // wrap to [0, 360)
    return mod_mat(ang, static_cast<T>(360.0));
}



template<typename T>
constexpr T smooth_pos(T a, scalar eps2)
{
    return 0.5*(a + sqrt(a*a + eps2));
}

template<typename T>
constexpr T smooth_sign(T a, scalar eps2)
{
    return tanh(eps2*a);
}


template <typename T>
inline std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;

    // Make sure the last element numerically fits the endpoint
    xs.back() = b;
    return xs;
}


template<typename T>
inline T trapz(const std::vector<T>& x, const std::vector<T>& y)
{
    assert(x.size() == y.size());

    T val{ 0.0 };

    for (size_t i = 1; i < x.size(); ++i)
    {
        val += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }

    return val;
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

#endif
