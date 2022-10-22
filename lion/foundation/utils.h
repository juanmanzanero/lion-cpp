#ifndef __UTILS_H__
#define __UTILS_H__

#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <regex>
#include <iostream>
#include <time.h>
#include "constants.h"
#include "lion/thirdparty/include/logger.hpp"
#include "types.h"
#include "lion/math/vector3d.h"
#include "lion/foundation/lion_exception.h"

#define PRINTVARIABLE(MSG,VAR) out(2) << "[" << #MSG << "] " << #VAR << ": " << VAR << std::endl

constexpr double& Value(double& val);
constexpr const double& Value(const double& val);

inline bool to_bool(std::string str) 
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}


inline std::string get_current_date_and_time()
{   
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%dT%X", &tstruct);
    
    return buf;
}

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
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_closest_point(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, bool closed, const size_t i_start, const scalar maximum_distance);

template<typename T>
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_intersection(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, const Vector3d<T>& p0, bool closed);

template<typename T>
inline std::tuple<T,T> find_intersection(const Vector3d<T>& r0, const Vector3d<T>& p0, const sVector3d& r1, const sVector3d& p1);

#endif
