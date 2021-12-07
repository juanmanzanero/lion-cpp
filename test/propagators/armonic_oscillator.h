#ifndef __ARMONIC_OSCILLATOR_H__
#define __ARMONIC_OSCILLATOR_H__

class Armonic_oscillator
{
 public:
    std::vector<scalar> operator()(const std::array<scalar,2>& q, const std::array<scalar,2>& u,
                                       const scalar t) { return { q[1], -w0*w0*q[0] }; }

    inline static scalar w0 = 0.5;
};

inline Armonic_oscillator armonic_oscillator;

class No_control
{
 public:
    std::array<scalar,2> operator() (const std::array<scalar,2>& q, const scalar t) const { return {}; }
};

inline No_control no_control;


#endif
