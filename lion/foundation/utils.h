#ifndef __UTILS_H__
#define __UTILS_H__

#include "constants.h"
#include <regex>
#include<iostream>
#include "lion/thirdparty/include/logger.hpp"

#define PRINTVARIABLE(MSG,VAR) out(2) << "[" << #MSG << "] " << #VAR << ": " << VAR << std::endl

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
    const T diff = std::abs(x - y);

    const T ax = std::abs(x);
    const T ay = std::abs(y);

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
                    auto q{ std::abs(x / y) }; 
                    rEQ0 = (std::abs(q - std::floor(q + static_cast<T>(0.5))) <= static_cast<T>(2.2204460492503131E-16) * q); // 2.220446049250313E-16 == eps, but we'll leave it like this since Matlab hardcodes it
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

template<typename T = double>
constexpr T smooth_pos(T a, T eps2)
{
    return 0.5*(a + std::sqrt(a*a + eps2));
}


template <typename T>
inline std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}


#endif
