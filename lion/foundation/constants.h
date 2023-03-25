#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <cmath>

enum Axis{ X, Y, Z };

template<typename T = double>
constexpr T pi_T = T{ 3.1415926535897932385L };

template<typename T = double>
constexpr T eps_T = std::numeric_limits<T>::epsilon();

template<typename T = double>
constexpr T inf_T = std::numeric_limits<T>::infinity();

template<typename T = double>
constexpr T nan_T = std::numeric_limits<T>::quiet_NaN();


// handy names when using double, our usual case
constexpr auto pi{ pi_T<double> };
constexpr auto eps{ eps_T<double> };
constexpr auto inf{ inf_T<double> };
constexpr auto Inf{ inf };
constexpr auto NaN{ nan_T<double> };

constexpr double g0{ 9.81 };
constexpr double KMH{ 1.0 / 3.6 };
constexpr double RPM{ pi / 30.0 };
constexpr double CV{ 735.499 };
constexpr double DEG{ pi / 180.0 };
constexpr double RAD{ 180.0 / pi };

#endif
