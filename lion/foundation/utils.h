#ifndef __UTILS_H__
#define __UTILS_H__

#include "constants.h"
#include <regex>
#include<iostream>
#include "lion/thirdparty/include/logger.hpp"
#include "types.h"
#include "lion/math/vector3d.h"

#define PRINTVARIABLE(MSG,VAR) out(2) << "[" << #MSG << "] " << #VAR << ": " << VAR << std::endl

constexpr double& Value(double& val) { return val; }
constexpr const double& Value(const double& val) { return val; }

template <typename T> 
constexpr int sign(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template<typename T = double>
constexpr bool eq_fp(T x, T y, T tol = static_cast<T>(1.))
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

template<typename T = double>
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

template<typename T = double>
constexpr T wrap_to_pi(T ang)
{
    // wrap to [-pi, pi)
    return mod_mat(ang + pi_T<T>, static_cast<T>(2.0) * pi_T<T>) - pi_T<T>;
}

template<typename T = double>
constexpr T wrap_to_2pi(T ang)
{
    // wrap to [0, 2 * pi)
    return mod_mat(ang, static_cast<T>(2.0) * pi_T<T>);
}

template<typename T = double>
constexpr T wrap_to_180(T ang)
{
    // wrap to [-180, 180)
    return mod_mat(ang + static_cast<T>(180.0), static_cast<T>(360.0)) - static_cast<T>(180.);
}

template<typename T = double>
constexpr T wrap_to_360(T ang)
{
    // wrap to [0, 360)
    return mod_mat(ang, static_cast<T>(360.0));
}


inline std::vector<double> string_to_double_vector(std::string s)
{
    std::vector<double> result;
    s = std::regex_replace(s, std::regex(","), " ");
    s = std::regex_replace(s, std::regex("\n"), " ");
    std::string::size_type sz;
    while ( s.find_first_of("0123456789") != std::string::npos )
    {
        result.push_back(std::stod(s,&sz));
        s = s.substr(sz);
    }

    return result;
}

template<typename T = scalar>
constexpr T smooth_pos(T a, scalar eps2)
{
    return 0.5*(a + sqrt(a*a + eps2));
}

template<typename T = scalar>
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
inline std::pair<Vector3d<T>,T> find_closest_point(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0)
{
    // Get number of points
    const auto N = xy_polygon.size() - 1;

    // (1) do a search along the poly vertices
    size_t i_closest_poly      = 0;
    sVector3d x_closest_poly = xy_polygon.front();
    T d2_closest_poly           = dot(x_closest_poly-x0,x_closest_poly-x0);
    T d2_current                = d2_closest_poly;
    
    for (size_t i = 1; i < N; ++i)
    {
        d2_current = dot(xy_polygon[i]-x0,xy_polygon[i]-x0);
        if ( d2_current < d2_closest_poly )
        {
            i_closest_poly = i;
            d2_closest_poly = d2_current;
            x_closest_poly = xy_polygon[i];
        }
    }
    
    T d2_closest = d2_closest_poly;
    Vector3d<T> x_closest = Vector3d<T>(x_closest_poly);
    
    // (2) compute the minimum distance to the next forward face
    sVector3d x_next_fwd;
    if ( i_closest_poly == N - 1)
        x_next_fwd = xy_polygon.front();
    else
        x_next_fwd = xy_polygon[i_closest_poly+1];

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
        }
    }
    
    // (3) compute the minimum distance to the next backward face
    sVector3d x_next_bwd;
    if ( i_closest_poly == 0 )
        x_next_bwd = xy_polygon.back();
    else
        x_next_bwd = xy_polygon[i_closest_poly-1];
    
    p = x_next_bwd - x_closest_poly;
    T s_next_bwd_closest = dot(p,x0-x_closest_poly)/dot(p,p);
    
    if ( (s_next_bwd_closest > 0.0) && (s_next_bwd_closest < 1.0) )
    {
        const auto x_next_bwd_closest = x_closest_poly + p*s_next_bwd_closest;
        const auto d2_next_bwd_closest = dot(x_next_bwd_closest - x0, x_next_bwd_closest - x0);
        if ( d2_next_bwd_closest < d2_closest )
        {
            x_closest = x_next_bwd_closest;
            d2_closest = d2_next_bwd_closest;
        }
    }

    return {x_closest,d2_closest};
}

#endif
