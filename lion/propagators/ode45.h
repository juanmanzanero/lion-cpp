#ifndef __ODE45_H__
#define __ODE45_H__

#include <algorithm>
#include <string>
#include <array>
#include <iostream>
#include <iomanip>

#include "lion/foundation/types.h"
#include "lion/foundation/constants.h"
#include "lion/math/matrix_extensions.h"
#include "lion/foundation/lion_exception.h"

namespace lioncpp{
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

template<typename F, typename U, size_t N>
void ODE45<F,U,N>::set(const std::string& name, const double value)
{
    if ( name == "relative error" )
    {
        _rtol = value;
        _threshold = _atol / _rtol;
    }
    else if ( name == "absolute error" )
    {
        _atol = value;
        _threshold = _atol / _rtol;
    }
    else if ( name == "max h" )
        _hmax = value;

    else
        throw lion_exception("Property \"" + name + "\" is not recognized");
}

template<typename F, typename U, size_t N>
const double& ODE45<F,U,N>::get(const std::string& name) 
{
    if ( name == "relative error" )
        return _rtol;

    else if ( name == "max h" )
        return _hmax;

    else
        throw lion_exception("Property \"" + name + "\" is not recognized");
}



template<typename F, typename U, size_t N>
double ODE45<F,U,N>::initial_dt_estimation(F& f, U& u, const std::array<scalar,N>& q0, scalar t_start, scalar t_end) 
{
    const double hmin = 16*std::abs(t_start)*eps;
    const double htspan = t_end - t_start;
    double absh = std::min(_hmax, htspan);
    const auto f0 = ext::abs(f.ode(q0, u(q0,t_start), t_start) / ext::max(ext::abs(q0),_threshold));
    const double rh = *std::max_element(f0.cbegin(), f0.cend())/ (0.8*std::pow(_rtol,_pow));

    if ( absh * rh > 1.0 )
        absh = 1.0 / rh;
    
    return std::max(absh,hmin);
}


template<typename F, typename U, size_t N>
void ODE45<F,U,N>::take_step(F& f, U& u, std::array<scalar,N>& q, scalar& t, scalar& dt, 
                             const scalar t_end, bool& t_end_reached )
{
    t_end_reached = false;
    // Clip the dt with the minimum and maximum values
    const double hmin = 16.0*std::abs(t)*eps;
    double absh = std::min(_hmax, std::max(hmin, dt)); 
    double h = absh;

    // Stretch the step if within 10% of t_end - t
    if ( 1.1*absh >= (t_end - t) )
    {
        h = t_end - t;
        absh = h;
        t_end_reached = true;
    }

    // loop to advance one step
    bool nofailed = true;
    scalar tnew;
    typename std::remove_reference<decltype(q)>::type qnew;
    scalar err;
    while (true)
    {
        const auto   q1 = q;
        const scalar t1 = t;
        const auto   f1 = f.ode(q1, u(q1,t1), t1);
        
        const auto   q2 = q + h * _b11 * f1;
        const scalar t2 = t + h * _a2;
        const auto   f2 = f.ode(q2, u(q2,t2), t2);
        
        const auto   q3 = q + h * (_b21 * f1 + _b22 * f2 );
        const scalar t3 = t + h * _a3;
        const auto   f3 = f.ode(q3, u(q3,t3), t3);
        
        const auto   q4 = q + h * (_b31 * f1 + _b32 * f2 + _b33 * f3 );
        const scalar t4 = t + h * _a4;
        const auto   f4 = f.ode(q4, u(q4,t4), t4);
        
        const auto   q5 = q + h  * (_b41 * f1 + _b42 * f2 + _b43 * f3 + _b44 * f4 );
        const scalar t5 = t + h * _a5;
        const auto   f5 = f.ode(q5, u(q5,t5), t5);
       
        const auto   q6 = q + h * (_b51 * f1 + _b52 * f2 + _b53 * f3 + _b54 * f4 + _b55 * f5 );
        const scalar t6 = t + h;
        const auto   f6 = f.ode(q6, u(q6,t6), t6);

        tnew = t + h;

        if ( t_end_reached )
            tnew = t_end;

        // Free "h" from roundoff errors
        h = tnew - t;

        qnew = q + h * (_b61 * f1 + _b63 * f3 + _b64 * f4 + _b65 * f5 + _b66 * f6 );
        const auto f7 = f.ode(qnew, u(qnew,tnew), tnew);

        // Error estimation
        const auto fE = ext::abs((f1 * _e1 + f3 * _e3 + f4 * _e4 + f5 * _e5 + f6 * _e6 + f7 * _e7)/ext::max(ext::max(ext::abs(q),ext::abs(qnew)),_threshold));

        err = absh * (*std::max_element(fE.cbegin(), fE.cend()));

        if ( err <= _rtol )
            // Successful step
            break;

        else
        {
            if (absh <= hmin)
                throw lion_exception("The current timestep is lower than the minimum timestep allowed");
            
            if (nofailed)
            {
                nofailed = false;
                absh = std::max(hmin, absh * std::max(0.1, 0.8*std::pow(_rtol/err, _pow)));
            }
            else
            {
                absh = std::max(hmin, 0.5*absh);
            }

            h = absh;
            t_end_reached = false;
        }
    }

    // If there were no failures, try to extend the time step
    if (nofailed)
    {
        if ( !t_end_reached )
        {
            const double temp = 1.25*std::pow(err/_rtol,_pow);
            if ( temp > 0.2 )
                absh *= (1.0 / temp);
            else
                absh *= 5.0;
        }
        else
        {
            absh = dt;
        }
    }
    
    t = tnew;
    q = qnew;
    dt = absh;

    return;
}



}
#endif

