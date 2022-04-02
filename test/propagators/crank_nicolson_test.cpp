#include "lion/propagators/crank_nicolson.h"
#include "lion/foundation/utils.hpp"
#include "gtest/gtest.h"
#include "lion/io/Xml_document.h"

class Armonic_oscillator_DAE
{
 public:
    template<typename T>
    std::tuple<std::array<T,2>,std::array<T,0>> operator()
        (const std::array<T,2>& q, const std::array<T,0>& qa, const std::array<T,0>& u, const T t) 
    { return { { q[1], -w0*w0*q[0] }, {}}; }

    inline static scalar w0 = 0.5;
};

class Robertson_equation
{
 public:
    template<typename T>
    std::tuple<std::array<T,2>,std::array<T,1>> operator()
        (const std::array<T,2>& q, const std::array<T,1>& qa, const std::array<T,0>& u, const T t)
    {
        std::array<T,2> dqdt = {-0.04*q[0] + 1.0e4*q[1]*qa[0], 0.04*q[0] - 3.0e7*q[1]*q[1] - 1.0e4*q[1]*qa[0]}; 
        std::array<T,1> dqa  = {q[0] + q[1] + qa[0] - 1.0};
        
        return {dqdt, dqa};
    }
};


TEST(Crank_nicolson_test, ode_armonic_oscillator)
{
    std::array<scalar,2> q = {1.0, 0.0};
    std::array<scalar,0> qa = {};
    const scalar dt = 1.0e-1;
    const size_t n_timesteps = 100;
    const scalar tf = dt*n_timesteps;
    scalar t = 0.0;

    Armonic_oscillator_DAE armonic_oscillator;

    for (size_t i = 0; i < n_timesteps; ++i)
        Crank_nicolson<Armonic_oscillator_DAE,2,0,0>::take_step(armonic_oscillator, {}, {}, q, qa, t, dt);

    EXPECT_NEAR(tf, t, 2.0e-14);

    EXPECT_NEAR(q[0], cos(armonic_oscillator.w0*tf), 1.0e-3);
    EXPECT_NEAR(q[1],-armonic_oscillator.w0*sin(armonic_oscillator.w0*tf), 1.0e-3);
}


TEST(Crank_nicolson_test, robertson_equation)
{
    std::array<scalar,2> q = {1.0, 0.0};
    std::array<scalar,1> qa = {0.0};

    std::vector<std::array<scalar,2>> q_values;
    std::vector<std::array<scalar,1>> qa_values;

    const size_t n_timesteps = 100;
    std::vector<scalar> t = linspace(-6.0,5.0,n_timesteps+1);
    std::transform(t.begin(), t.end(), t.begin(), [](const auto& t) -> auto { return pow(10.0,t); });

    Robertson_equation robertson_equation;

    q_values.push_back(q);
    qa_values.push_back(qa);
    for (size_t i = 0; i < n_timesteps; ++i)
    {
        scalar dt = t[i+1] - t[i];
        scalar t_val = t[i];
        Crank_nicolson<Robertson_equation,2,1,0>::take_step(robertson_equation, {}, {}, q, qa, t_val, dt);

        q_values.push_back(q);
        qa_values.push_back(qa);
    }

    Xml_document results_saved("data/crank_nicolson_test.xml",true);

    const auto t_saved = results_saved.get_element("robertson_equation_test/t").get_value(std::vector<scalar>());
    const auto q0_saved = results_saved.get_element("robertson_equation_test/q0").get_value(std::vector<scalar>());
    const auto q1_saved = results_saved.get_element("robertson_equation_test/q1").get_value(std::vector<scalar>());
    const auto qa0_saved = results_saved.get_element("robertson_equation_test/qa0").get_value(std::vector<scalar>());

    for (size_t i = 0; i < n_timesteps; ++i)
    {
        EXPECT_NEAR(t[i], t_saved[i], 1.0e-6);
        EXPECT_NEAR(q_values[i][0], q0_saved[i], 1.0e-6);
        EXPECT_NEAR(q_values[i][1], q1_saved[i], 1.0e-6);
        EXPECT_NEAR(qa_values[i][0], qa0_saved[i], 1.0e-6);
    }
}
