#ifndef __ODE45_HPP__
#define __ODE45_HPP__

#include <iomanip>

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
        throw std::runtime_error("Property \"" + name + "\" is not recognized");
}

template<typename F, typename U, size_t N>
const double& ODE45<F,U,N>::get(const std::string& name) 
{
    if ( name == "relative error" )
        return _rtol;

    else if ( name == "max h" )
        return _hmax;

    else
        throw std::runtime_error("Property \"" + name + "\" is not recognized");
}



template<typename F, typename U, size_t N>
double ODE45<F,U,N>::initial_dt_estimation(F& f, U& u, const std::array<scalar,N>& q0, scalar t_start, scalar t_end) 
{
    const double hmin = 16*std::abs(t_start)*eps;
    const double htspan = t_end - t_start;
    double absh = std::min(_hmax, htspan);
    const auto f0 = ext::abs(f(q0, u(q0,t_start), t_start) / ext::max(ext::abs(q0),_threshold));
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
        const auto   f1 = f(q1, u(q1,t1), t1);
        
        const auto   q2 = q + h * _b11 * f1;
        const scalar t2 = t + h * _a2;
        const auto   f2 = f(q2, u(q2,t2), t2);
        
        const auto   q3 = q + h * (_b21 * f1 + _b22 * f2 );
        const scalar t3 = t + h * _a3;
        const auto   f3 = f(q3, u(q3,t3), t3);
        
        const auto   q4 = q + h * (_b31 * f1 + _b32 * f2 + _b33 * f3 );
        const scalar t4 = t + h * _a4;
        const auto   f4 = f(q4, u(q4,t4), t4);
        
        const auto   q5 = q + h  * (_b41 * f1 + _b42 * f2 + _b43 * f3 + _b44 * f4 );
        const scalar t5 = t + h * _a5;
        const auto   f5 = f(q5, u(q5,t5), t5);
       
        const auto   q6 = q + h * (_b51 * f1 + _b52 * f2 + _b53 * f3 + _b54 * f4 + _b55 * f5 );
        const scalar t6 = t + h;
        const auto   f6 = f(q6, u(q6,t6), t6);

        tnew = t + h;

        if ( t_end_reached )
            tnew = t_end;

        // Free "h" from roundoff errors
        h = tnew - t;

        qnew = q + h * (_b61 * f1 + _b63 * f3 + _b64 * f4 + _b65 * f5 + _b66 * f6 );
        const auto f7 = f(qnew, u(qnew,tnew), tnew);

        // Error estimation
        const auto fE = ext::abs((f1 * _e1 + f3 * _e3 + f4 * _e4 + f5 * _e5 + f6 * _e6 + f7 * _e7)/ext::max(ext::max(ext::abs(q),ext::abs(qnew)),_threshold));

        err = absh * (*std::max_element(fE.cbegin(), fE.cend()));

        if ( err <= _rtol )
            // Successful step
            break;

        else
        {
            if (absh <= hmin)
                throw std::runtime_error("The current timestep is lower than the minimum timestep allowed");
            
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





#endif
