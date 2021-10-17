#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <cmath>

enum Axis{ X, Y, Z };

template<typename T = double>
constexpr T pi_T = static_cast<T>(3.1415926535897932385L); // const T pi = static_cast<T>(std::acos(static_cast<T>(-1.)));

template<typename T = double>
constexpr T eps_T = std::numeric_limits<T>::epsilon();

template<typename T = double>
constexpr T inf_T = std::numeric_limits<T>::infinity();

template<typename T = double>
struct nan_T
{
    static inline T nan = static_cast<T>(std::nan(""));
};

template<>
struct nan_T<float>
{
    static inline float nan = std::nanf("");
};

template<>
struct nan_T<long double>
{
    static inline double nan = std::nanl("");
};

// handy names when using double, our usual case
constexpr double pi{ pi_T<double> };
constexpr double eps{ eps_T<double> };
constexpr double inf{ inf_T<double> };
constexpr double Inf{ inf };
inline static double NaN{ nan_T<double>::nan };

constexpr inline double g0{9.81};
constexpr inline double KMH{1.0/3.6};
constexpr inline double RPM{pi/30.0};
constexpr inline double CV{735.499};
constexpr inline double DEG{pi/180.0};
constexpr inline double RAD{180.0/pi};
#endif
