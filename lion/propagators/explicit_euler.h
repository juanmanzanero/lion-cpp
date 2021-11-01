#ifndef __EXPLICIT_EULER_H__
#define __EXPLICIT_EULER_H__

#include "lion/foundation/types.h"

//! A class to perform an explicit Euler step
//! @param F: ODE functor dqdt = f(q,u,t)
//! @param U: controls functor u = u(t)
//! @param N: order of the ODE system
template<class F, class U, size_t N>
class Explicit_euler
{
 public:
    //! Take a step of the explicit euler
    //! q = q + dt * f(q, u(q,t), t);
    //! @param[inout] f: ODE functor dqdt = f(q,u,t)
    //! @param[inout] u: control variables functor u = u(q,t)
    //! @param[inout] t: time
    //! @param[inout] dt: time step
    static void take_step(F& f, U& u, std::array<scalar,N>& q, scalar t, scalar dt);
};

#include "explicit_euler.hpp"

#endif
