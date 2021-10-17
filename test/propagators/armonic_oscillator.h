#ifndef __ARMONIC_OSCILLATOR_H__
#define __ARMONIC_OSCILLATOR_H__

class Armonic_oscillator
{
 public:
    std::vector<timeseries> operator()(const std::array<timeseries,2>& q, const std::array<timeseries,2>& u,
                                       const timeseries t) { return { q[1], -w0*w0*q[0] }; }

    inline static scalar w0 = 0.5;
};

inline Armonic_oscillator armonic_oscillator;

class No_control
{
 public:
    std::array<timeseries,2> operator() (const timeseries t) const { return {}; }
};

inline No_control no_control;


#endif
