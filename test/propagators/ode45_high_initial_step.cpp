#include "lion/propagators/ode45.h"
#include "gtest/gtest.h"
#include "armonic_oscillator.h"

using namespace lioncpp;

class ODE45_high_initial_step_test : public ::testing::Test
{
 protected:
    ODE45_high_initial_step_test() 
    { 
        ODE45<Armonic_oscillator,No_control,2>::set("max h",100.0); 
        ODE45<Armonic_oscillator,No_control,2>::set("relative error",1.0e-3);
        ODE45<Armonic_oscillator,No_control,2>::set("absolute error",1.0e-6);
    }

    const std::array<std::array<scalar,2>,8> ref_solution = 
                    {{
                        {{  1.0000000000000000e+00,   0.0000000000000000e+00}},
                        {{  8.3003661091622583e-01,  -2.7884781925664337e-01}},
                        {{  3.7793635024455147e-01,  -4.6290779771432911e-01}},
                        {{ -4.2901474076792595e-01,  -4.5162632241534983e-01}},
                        {{ -9.2710058393628225e-01,  -1.8732628528380596e-01}},
                        {{ -9.0579972982766399e-01,   2.1175241679784929e-01}},
                        {{ -3.8283971024204255e-01,   4.6185084201058341e-01}},
                        {{  2.8397359100141184e-01,   4.7935428202062896e-01}} 
                    }};

    const std::array<scalar,8> ref_t = {
                          0.0000000000000000e+00,
                          1.1832019643370637e+00,
                          2.3664039286741274e+00,
                          4.0281749120717691e+00,
                          5.5147302063985526e+00,
                          7.1573241876406506e+00,
                          8.6382922272292255e+00,
                          1.0000000000000000e+01 
                     };

    static constexpr double t0 = 0.0;
    static constexpr double tf = 10.0;

    static inline std::array<scalar,2> q = {1.0, 0.0};
    static inline double t = t0;
    static inline double dt;
    static inline bool t_end_reached = false;
};


TEST_F(ODE45_high_initial_step_test, step1)
{
    dt = 10.0;
    constexpr size_t i = 1;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}

TEST_F(ODE45_high_initial_step_test, step2)
{
    constexpr size_t i = 2;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_high_initial_step_test, step3)
{
    constexpr size_t i = 3;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_high_initial_step_test, step4)
{
    constexpr size_t i = 4;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_NEAR(q[1],ref_solution[i][1],3.0e-14);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_high_initial_step_test, step5)
{
    constexpr size_t i = 5;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_NEAR(q[1],ref_solution[i-1][1],3.0e-14);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_NEAR(q[0],ref_solution[i][0],3.0e-14);
    EXPECT_NEAR(q[1],ref_solution[i][1],3.0e-14);
    EXPECT_NEAR(t,ref_t[i],3.0e-14);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_high_initial_step_test, step6)
{
    constexpr size_t i = 6;
    EXPECT_NEAR(q[0],ref_solution[i-1][0],3.0e-14);
    EXPECT_NEAR(q[1],ref_solution[i-1][1],3.0e-14);
    EXPECT_NEAR(t,ref_t[i-1],3.0e-14);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_NEAR(q[0],ref_solution[i][0],3.0e-14);
    EXPECT_NEAR(q[1],ref_solution[i][1],3.0e-14);
    EXPECT_NEAR(t,ref_t[i],3.0e-14);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_high_initial_step_test, step7)
{
    constexpr size_t i = 7;
    EXPECT_NEAR(q[0],ref_solution[i-1][0],3.0e-14);
    EXPECT_NEAR(q[1],ref_solution[i-1][1],3.0e-14);
    EXPECT_NEAR(t,ref_t[i-1],3.0e-14);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_NEAR(q[0],ref_solution[i][0],3.0e-14);
    EXPECT_NEAR(q[1],ref_solution[i][1],3.0e-14);
    EXPECT_NEAR(t,ref_t[i],3.0e-14);
    EXPECT_EQ(t_end_reached,true);
}

