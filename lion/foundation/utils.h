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
#include "types.h"
#include "lion/math/vector3d.h"
#include "lion/foundation/lion_exception.h"

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

template<typename T>
constexpr bool samesign(T x, T y);

template<typename T>
constexpr bool samesign(T x, T y, T z);

template<typename T, typename TolType = T>
constexpr bool eq_fp(T x, T y, TolType tol = TolType{ 1 });

template<typename T>
constexpr T mod_mat(T x, T y);

template<typename T>
constexpr T wrap_to_pi(T ang);

template<typename T>
constexpr T wrap_to_2pi(T ang);

template<typename T>
constexpr T wrap_to_180(T ang);

template<typename T>
constexpr T wrap_to_360(T ang);

template<typename scalar_type = double>
inline std::vector<scalar_type> string_to_double_vector(std::string s)
{
    std::vector<scalar_type> result;
    s = std::regex_replace(s, std::regex(","), " ");
    s = std::regex_replace(s, std::regex("\n"), " ");
    std::string::size_type sz;
    while ( s.find_first_of("0123456789") != std::string::npos )
    {
        if constexpr (std::is_same_v<scalar_type,double>)
            result.push_back(std::stod(s,&sz));
        else if constexpr (std::is_same_v<scalar_type,float>)
            result.push_back(std::stof(s,&sz));

        s = s.substr(sz);
    }

    return result;
}


template<typename T>
constexpr T smooth_pos(T a, scalar eps2);

template<typename T>
constexpr T smooth_sign(T a, scalar eps2);

template <typename T>
constexpr std::vector<T> linspace(T lo, T hi, std::size_t num_points);

template<typename T>
inline T trapz(const std::vector<T>& x, const std::vector<T>& y);

template<typename T>
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_closest_point(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, bool closed, const size_t i_start, const scalar maximum_distance);

template<typename T>
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_intersection(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, const Vector3d<T>& p0, bool closed);

template<typename T>
inline std::tuple<T,T> find_intersection(const Vector3d<T>& r0, const Vector3d<T>& p0, const sVector3d& r1, const sVector3d& p1);


template<typename T>
constexpr std::array<T, 3u> rotmat2ea(const std::array<T, 9u> &M)
{
    //
    // Returns an std::array holding the Euler angles in 'ZYX' sequence
    // "[yaw, pitch, roll]", in rad, from the input rotation matrix,
    // given as an std::array in column-major format.
    //

    // gymbal lock tolerance, calculated with
    // "gymbal_lock_pitch_tol_deg = 0.001;
    //  gymbal_lock_pitch_zone = cos(deg2rad(gymbal_lock_pitch_tol_deg));"
    constexpr auto gymbal_lock_pitch_zone = T{ 0.999999999847691 };
    const auto sinp = -M[2];


    T pitch_rad;
    if (std::fabs(sinp) <= gymbal_lock_pitch_zone) {
        pitch_rad = std::asin(sinp);
    }
    else if (sinp >= T{ 0 }) {
        pitch_rad = T{ 0.5 } * pi_T<T>;
    }
    else {
        pitch_rad = -T{ 0.5 } * pi_T<T>;
    }

    const auto yaw_rad = std::atan2(M[1], M[0]);
    const auto roll_rad = std::atan2(M[5], M[8]);

    return std::array<T, 3u>{ wrap_to_pi(yaw_rad),
                              wrap_to_pi(pitch_rad),
                              wrap_to_pi(roll_rad) };
}


template<typename SizeType>
constexpr SizeType nchoosek(SizeType n, SizeType k);

#endif
