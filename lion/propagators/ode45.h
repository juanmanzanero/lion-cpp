#ifndef __ODE45_H__
#define __ODE45_H__

#include <algorithm>
#include <string>
#include <array>
#include <iostream>
#include "lion/foundation/types.h"
#include "lion/foundation/constants.h"
#include "lion/math/matrix_extensions.h"

//! Implementation of the ode45 propagator
//! @param F: type of the ODE functor dqdt = f(q,u,t)
//! @param U: type of the control variables functor u = u(t)
//! @param N: order of the ODE system
template<typename F, typename U, size_t N>
class ODE45
{
 public:
    //! Set a property with given name and value
    //! @param[in] name: name of the property
    //! @param[in] value: value to be given to the property
    static void set(const std::string& name, const double value);

    //! Get a property with given name
    //! @param[in] name: name of the property
    static const double& get(const std::string& name);

    //! Get an estimation of an initial time step
    //! @param[in] f: ODE functor, dqdt = f(q,u,t)
    //! @param[in] u: control variables functor, u = u(t)
    //! @param[in] q0: state vector values
    //! @param[in] t_start: simulation initial time
    //! @param[in] t_end: simulation final end
    static double initial_dt_estimation(F& f, U& u, const std::array<scalar,N>& q0, scalar t_start, scalar t_end);

    //! Perform a Runge-Kutta step
    //! @param[inout] f: ODE functor, dqdt = f(q,u,t)
    //! @param[inout] u: control variables functor, u = u(t)
    //! @param[inout] q: state vector values before and after of the step
    //! @param[inout] t: time before and after of the step
    //! @param[inout] dt: time step before and after of the step
    //! @param[in] t_end: simulation final time
    //! @param[in] t_end_reached: whether the simulation is done or not
    static void take_step(F& f, 
                          U& u, 
                          std::array<scalar,N>& q, 
                          scalar& t, 
                          scalar& dt, 
                          const scalar t_end, 
                          bool& t_end_reached );

 private:
    inline static double _rtol = 1.0e-3;        //! Relative tolerance
    inline static double _atol = 1.0e-6;        //! Absolute tolerance
    inline static double _threshold = 1.0e-3;   //! atol / rtol
    inline static double _hmax = 0.05;          //! Maximum time step size

    // Runge-Kutta stages coefficients
    constexpr static double _pow = double(1)/double(5);

    constexpr static double _a2 = double(1)/double(5);
    constexpr static double _a3 = double(3)/double(10);
    constexpr static double _a4 = double(4)/double(5);
    constexpr static double _a5 = double(8)/double(9);
    
    constexpr static double _b11 = double(1)/double(5); 
    constexpr static double _b21 = double(3)/double(40); 
    constexpr static double _b31 = double(44)/double(45);
    constexpr static double _b41 = double(19372)/double(6561);
    constexpr static double _b51 = double(9017)/double(3168);
    constexpr static double _b61 = double(35)/double(384);
    constexpr static double _b22 = double(9)/double(40);
    constexpr static double _b32 = double(-56)/double(15);
    constexpr static double _b42 = double(-25360)/double(2187);
    constexpr static double _b52 = double(-355)/double(33);
    constexpr static double _b33 = double(32)/double(9);
    constexpr static double _b43 = double(64448)/double(6561);
    constexpr static double _b53 = double(46732)/double(5247);
    constexpr static double _b63 = double(500)/double(1113);
    constexpr static double _b44 = double(-212)/double(729);
    constexpr static double _b54 = double(49)/double(176);
    constexpr static double _b64 = double(125)/double(192);
    constexpr static double _b55 = double(-5103)/double(18656);
    constexpr static double _b65 = double(-2187)/double(6784);
    constexpr static double _b66 = double(11)/double(84);
    
    // Error estimation coefficients
    constexpr static double _e1 = double(71)/double(57600);
    constexpr static double _e3 = double(-71)/double(16695);
    constexpr static double _e4 = double(71)/double(1920);
    constexpr static double _e5 = double(-17253)/double(339200);
    constexpr static double _e6 = double(22)/double(525);
    constexpr static double _e7 = double(-1)/double(40);
};

#include "ode45.hpp"

#endif

