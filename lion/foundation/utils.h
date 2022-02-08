#ifndef __UTILS_H__
#define __UTILS_H__

#include "constants.h"
#include <regex>
#include<iostream>
#include "lion/thirdparty/include/logger.hpp"
#include "types.h"
#include "lion/math/vector3d.h"

#define PRINTVARIABLE(MSG,VAR) out(2) << "[" << #MSG << "] " << #VAR << ": " << VAR << std::endl

constexpr double& Value(double& val);
constexpr const double& Value(const double& val);

template <typename T> 
constexpr int sign(T val); 

template<typename T = double>
constexpr bool eq_fp(T x, T y, T tol = static_cast<T>(1.));

template<typename T = double>
constexpr T mod_mat(T x, T y);

template<typename T = double>
constexpr T wrap_to_pi(T ang);

template<typename T = double>
constexpr T wrap_to_2pi(T ang);

template<typename T = double>
constexpr T wrap_to_180(T ang);

template<typename T = double>
constexpr T wrap_to_360(T ang);

inline std::vector<double> string_to_double_vector(std::string s);

template<typename T = scalar>
constexpr T smooth_pos(T a, scalar eps2);

template<typename T = scalar>
constexpr T smooth_sign(T a, scalar eps2);


template <typename T>
inline std::vector<T> linspace(T a, T b, size_t N);  


template<typename T>
inline std::pair<Vector3d<T>,T> find_closest_point(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, bool closed);

#endif
