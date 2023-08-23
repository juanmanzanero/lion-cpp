#include "lion/propagators/ode45.h"
#include "gtest/gtest.h"
#include "armonic_oscillator.h"

using namespace lioncpp;

class ODE45_test : public ::testing::Test
{
 protected:
    ODE45_test() 
    { 
        ODE45<Armonic_oscillator,No_control,2>::set("max h",0.5); 
        ODE45<Armonic_oscillator,No_control,2>::set("relative error",1.0e-3); 
    }

    const std::array<std::array<scalar,2>,25> ref_solution = 
                    {{
                        {{  1.0000000000000000e+00,   0.0000000000000000e+00}},
                        {{  9.9999991923746101e-01,  -2.0095090911099770e-04}},
                        {{  9.9999709254996549e-01,  -1.2057043186148851e-03}},
                        {{  9.9992238820292245e-01,  -6.2293171889733438e-03}},
                        {{  9.9803520655934508e-01,  -3.1327808968785967e-02}},
                        {{  9.5150739217346336e-01,  -1.5381283794546333e-01}},
                        {{  8.4581945283238025e-01,  -2.6673451261293701e-01}},
                        {{  6.8754256125780244e-01,  -3.6307191021931351e-01}},
                        {{  4.8651761621965495e-01,  -4.3683523909110661e-01}},
                        {{  2.5524338260102408e-01,  -4.8343825479957736e-01}},
                        {{  8.0993782713614237e-03,  -4.9998341049602707e-01}},
                        {{ -2.3954817507291715e-01,  -4.8544201215739169e-01}},
                        {{ -4.7230174943006609e-01,  -4.4071817757293752e-01}},
                        {{ -6.7568985436886342e-01,  -3.6859262242605262e-01}},
                        {{ -8.3706680424327007e-01,  -2.7354976859603514e-01}},
                        {{ -9.4639896557649461e-01,  -1.6149892426477611e-01}},
                        {{ -9.9688859956653719e-01,  -3.9406871515075909e-02}},
                        {{ -9.8539651222185543e-01,   8.5135294640867415e-02}},
                        {{ -9.1263723380985096e-01,   2.0438414393884485e-01}},
                        {{ -7.8313459233210869e-01,   3.1092535995018156e-01}},
                        {{ -6.0494044328320562e-01,   3.9813472640091829e-01}},
                        {{ -3.8913404374713267e-01,   4.6058998882392072e-01}},
                        {{ -1.4913319735806810e-01,   4.9440798398319380e-01}},
                        {{  1.0014000017378277e-01,   4.9748607600163713e-01}},
                        {{  2.8366236686220403e-01,   4.7946147954704305e-01}}
                    }};

    const std::array<scalar,25> ref_t = {
                      0.0000000000000000e+00,
                      8.0380365808306554e-04,
                      4.8228219484983941e-03,
                      2.4917913400575038e-02,
                      1.2539337066095824e-01,
                      6.2539337066095824e-01,
                      1.1253933706609582e+00,
                      1.6253933706609582e+00,
                      2.1253933706609582e+00,
                      2.6253933706609582e+00,
                      3.1253933706609582e+00,
                      3.6253933706609582e+00,
                      4.1253933706609587e+00,
                      4.6253933706609587e+00,
                      5.1253933706609587e+00,
                      5.6253933706609587e+00,
                      6.1253933706609587e+00,
                      6.6253933706609587e+00,
                      7.1253933706609587e+00,
                      7.6253933706609587e+00,
                      8.1253933706609587e+00,
                      8.6253933706609587e+00,
                      9.1253933706609587e+00,
                      9.6253933706609587e+00,
                      1.0000000000000000e+01
                     };

    static constexpr double t0 = 0.0;
    static constexpr double tf = 10.0;

    static inline std::array<scalar,2> q = {1.0, 0.0};
    static inline double t = t0;
    static inline double dt;
    static inline bool t_end_reached = false;
};


TEST_F(ODE45_test, initial_step_estimation)
{
    dt = ODE45<Armonic_oscillator,No_control,2>::initial_dt_estimation(armonic_oscillator, no_control, q, t0, tf);
    EXPECT_DOUBLE_EQ(dt, 8.038036580830655e-04);
}

TEST_F(ODE45_test, step1)
{
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


TEST_F(ODE45_test, step2)
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


TEST_F(ODE45_test, step3)
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


TEST_F(ODE45_test, step4)
{
    constexpr size_t i = 4;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step5)
{
    constexpr size_t i = 5;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step6)
{
    constexpr size_t i = 6;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step7)
{
    constexpr size_t i = 7;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step8)
{
    constexpr size_t i = 8;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step9)
{
    constexpr size_t i = 9;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step10)
{
    constexpr size_t i = 10;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step11)
{
    constexpr size_t i = 11;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step12)
{
    constexpr size_t i = 12;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step13)
{
    constexpr size_t i = 13;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step14)
{
    constexpr size_t i = 14;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step15)
{
    constexpr size_t i = 15;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step16)
{
    constexpr size_t i = 16;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step17)
{
    constexpr size_t i = 17;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step18)
{
    constexpr size_t i = 18;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step19)
{
    constexpr size_t i = 19;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step20)
{
    constexpr size_t i = 20;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step21)
{
    constexpr size_t i = 21;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step22)
{
    constexpr size_t i = 22;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step23)
{
    constexpr size_t i = 23;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,false);
}


TEST_F(ODE45_test, step24)
{
    constexpr size_t i = 24;
    EXPECT_DOUBLE_EQ(q[0],ref_solution[i-1][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i-1][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i-1]);

    ODE45<Armonic_oscillator,No_control,2>::take_step(armonic_oscillator, no_control, q, t, dt, tf, t_end_reached);

    EXPECT_DOUBLE_EQ(q[0],ref_solution[i][0]);
    EXPECT_DOUBLE_EQ(q[1],ref_solution[i][1]);
    EXPECT_DOUBLE_EQ(t,ref_t[i]);
    EXPECT_EQ(t_end_reached,true);
}
