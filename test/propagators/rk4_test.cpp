#include "gtest/gtest.h"
#include "lion/propagators/rk4.h"
#include <cmath>
#include "armonic_oscillator.h"

TEST(RK4_test, Armonic_oscillator)
{
    std::array<scalar,2> q = {1.0, 0.0};
    const scalar dt = 1.0e-2;
    const size_t n_timesteps = 50000;
    const scalar tf = dt*n_timesteps;

    for (size_t i = 0; i < n_timesteps; ++i)
        RK4<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator,no_control,q,i*dt,dt);

    EXPECT_NEAR(q[0], cos(armonic_oscillator.w0*tf), 1.0e-8);
    EXPECT_NEAR(q[1],-armonic_oscillator.w0*sin(armonic_oscillator.w0*tf), 1.0e-8);
}
