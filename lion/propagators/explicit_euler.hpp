#ifndef __EXPLICIT_EULER_HPP__
#define __EXPLICIT_EULER_HPP__

template<class F, class U, size_t N>
inline void Explicit_euler<F,U,N>::take_step(F& f, U& u, std::array<scalar,N>& q, scalar t, scalar dt)
{
    const auto u_t = u(q,t);

    auto dqdt = f(q,u_t,t);
    
    for ( auto [iq,idq] = std::tuple{ q.begin(), dqdt.begin()} ; iq != q.end(); iq++, idq++)
        *iq += ( (*idq) *= dt );
}

#endif
