#ifndef __RK4_HPP__
#define __RK4_HPP__

template<class F, class U, size_t N>
inline void RK4<F,U,N>::take_step(F& f, U& u, std::array<timeseries,N>& q, scalar t, scalar dt)
{
    // 1st stage
    const auto u_1 = u(t);
    const auto k1 = f(q,u_1,t);

    // 2nd stage
    const scalar t2 = t + 0.5*dt;
    const auto u_2 = u(t2);

    auto q2 = q;
    for (size_t i = 0; i < q.size(); ++i)
        q2[i] += 0.5*dt*k1[i];

    const auto k2 = f(q2,u_2,t2);

    // 3rd stage
    const scalar t3 = t + 0.5*dt;
    const auto u_3 = u(t3);

    auto q3 = q;
    for (size_t i = 0; i < q.size(); ++i)
        q3[i] += 0.5*dt*k2[i];

    const auto k3 = f(q3,u_3,t3);

    // 4th stage
    const scalar t4 = t;

    const auto u_4 = u(t4);

    auto q4 = q;
    for (size_t i = 0; i < q.size(); ++i)
        q4[i] += dt*k3[i];

    const auto k4 = f(q4,u_4,t4);

    for (size_t i = 0; i < q.size(); ++i)
        q[i] += (1.0/6.0)*dt*(k1[i]+2.0*(k2[i]+k3[i])+k4[i]);
    
}


#endif
