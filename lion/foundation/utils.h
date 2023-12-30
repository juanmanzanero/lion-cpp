#ifndef LION_FOUNDATION_UTILS_H
#define LION_FOUNDATION_UTILS_H
#pragma once


#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <regex>
#include <iostream>
#include <time.h>
#include <numeric>

#include "lion/foundation/types.h"
#include "lion/foundation/constants.h"
#include "lion/foundation/lion_exception.h"
#include "lion/foundation/string_helpers.h"
#include "lion/math/vector3d.h"


constexpr double& Value(double& val);
constexpr const double& Value(const double& val);


inline bool to_bool(const std::string &str)
{
    if (str == "1") {
        return true;
    }
    else if (str == "0") {
        return false;
    }
    else {
        const auto strlo = strtolower(strtrim(str));
        if (strlo == "true") {
            return true;
        }
        else if (strlo == "false") {
            return false;
        }
        else {
            throw std::runtime_error(
                "to_bool: input boolean string \"" + str + "\" is invalid.");
        }
    }
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


template <typename T0, typename T1, typename T2>
constexpr bool in_or_on_interval(const T0 &x, const T1 &lo, const T2 &hi);

template <typename T> 
constexpr T sign(const T &x);

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


template<typename T = double,
         typename std::enable_if_t<std::is_floating_point_v<T> >* = nullptr>
inline std::vector<T> string_to_double_vector(std::string str)
{
    std::vector<T> ret;
    while (true) {
        const auto str_length = str.length();
        if (str_length == 0u) {
            break;
        }

        std::size_t pos;
        if constexpr (std::is_same_v<T, double>) {
            ret.push_back(std::stod(str, &pos));
        }
        else {
            static_assert(std::is_same_v<T, float>);
            ret.push_back(std::stof(str, &pos));
        }

        while (pos != str_length && (str[pos] == ',' || str[pos] == ';' || str[pos] == '\n')) {
            ++pos;
        }
        str = str.substr(pos);
    }

    return ret;
}


template<bool ActuallySmooth = true, typename T, typename S>
constexpr T smooth_pos(const T &x, S eps2);

template<bool ActuallySmooth = true, typename T, typename S>
constexpr T smooth_neg(const T &x, S eps2);

template<bool ActuallySmooth = true, bool UseSinAtanFormula = true, typename T, typename S>
constexpr T smooth_sign(const T &x, S eps);

template<bool ActuallySmooth = true, bool UseSinAtanFormula = true, typename T, typename S>
constexpr T smooth_sign_derivative(const T &x, S eps);

template<bool ActuallySmooth = true, typename T, typename S>
constexpr T smooth_abs(const T &x, S eps2);

template<bool ActuallySmooth = true, typename T, typename S>
constexpr T smooth_hypot(const T &x, const T &y, S eps2);

template<bool ActuallySmooth = true, typename T, typename T1, typename S>
constexpr T smooth_max(const T &x, const T1 &lo, S eps);

template<bool ActuallySmooth = true, typename T, typename T1, typename S>
constexpr T smooth_min(const T &x, const T1 &hi, S eps);

template<bool ActuallySmooth = true, typename T, typename T1, typename T2, typename S>
constexpr T smooth_clamp(const T &x, const T1 &lo, const T2 &hi, S eps2);

template<bool ActuallySmooth = true, bool UseSinAtanFormula = true, typename T, typename S>
constexpr T smooth_step(const T &x, S eps);


template <typename T>
constexpr std::vector<T> linspace(T lo, T hi, std::size_t num_points);

template<typename T>
constexpr std::vector<T> iota(T lo, std::size_t num_points, T increment = T{ 1 });


template<typename T>
inline T trapz(const std::vector<T>& x, const std::vector<T>& y);

template<typename ArrayType,
    typename ValueType = typename ArrayType::value_type>
constexpr ValueType sumabs(const ArrayType &x);

template<typename ArrayType,
    typename ValueType = typename ArrayType::value_type>
constexpr ValueType sumsqr(const ArrayType &x);

template<typename T>
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_closest_point(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, bool closed, const size_t i_start, const scalar maximum_distance);

template<typename T>
inline std::tuple<Vector3d<T>,T,std::array<size_t,2>> find_intersection(const std::vector<sVector3d>& xy_polygon, const Vector3d<T>& x0, const Vector3d<T>& p0, bool closed);

template<typename T>
inline std::tuple<T,T> find_intersection(const Vector3d<T>& r0, const Vector3d<T>& p0, const sVector3d& r1, const sVector3d& p1);

template<typename SizeType>
constexpr SizeType nchoosek(SizeType n, SizeType k);

template<typename T,
    typename Array2Type = std::array<T, 2> >
constexpr std::pair<Array2Type, bool> sin_cos_solve(T lhs_s, T lhs_c, T rhs,
    T tolzero = 1e3 * std::numeric_limits<T>::epsilon());

template<typename ContainerOfGridVectorsType,
         typename ScalarType = typename ContainerOfGridVectorsType::value_type::value_type>
std::vector<ScalarType> grid_vectors2points_rowmaj(const ContainerOfGridVectorsType &grid_vectors);

template<typename Container, typename ValueType = typename Container::value_type>
constexpr ValueType median(Container cont);

template<typename It, typename ValueType>
constexpr It nearest_in_sorted_range(It first, It last, ValueType val);

#endif
