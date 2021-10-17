#ifndef __RK4_H__
#define __RK4_H__

#include <array>
#include "lion/foundation/types.h"

//! Implementation of the rk4 propagator
//! @param F: type of the ODE functor dqdt = f(q,u,t)
//! @param U: type of the control variables functor u = u(t)
//! @param N: order of the ODE system
template<class F, class U, size_t N>
class RK4
{
 public:
    //! Take a step of the explicit euler
    //! @param[inout] f: ODE functor dqdt = f(q,u,t)
    //! @param[inout] u: control variables functor u = u(q,t)
    //! @param[inout] t: time
    //! @param[inout] dt: time step
    static void take_step(F& f, U& u, std::array<timeseries,N>& q, timeseries t, scalar dt);
};

#include "rk4.hpp"

#endif
