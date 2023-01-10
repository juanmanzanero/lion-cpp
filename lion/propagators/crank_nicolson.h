#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include <array>
#include "lion/foundation/types.h"
#include "lion/math/linear_algebra.h"
#include "lion/math/ipopt_cppad_handler.hpp"
#include "lion/math/matrix_extensions.h"
#include <cppad/cppad.hpp>

template<typename F>
class Crank_nicolson
{
 public:

    struct Options
    {
        scalar sigma = 0.5;
        size_t max_iter = 100000;
        scalar error_tolerance = 1.0e-12;
        scalar relaxation_factor = 0.25;
    };

    //! Perform a Crank-Nicolson step for a differential-algebraic equations system
    //! @param[inout] f: ODE functor, [states,dstates_dt] = f(inputs,controls,time)
    //! @param[inout] controls_ini: initial control vector
    //! @param[inout] controls_fin: final control vector
    //! @param[inout] inputs: state vector values before and after of the step
    //! @param[inout] time: time before and after of the step
    //! @param[inout] delta_time: time step before and after of the step
    static void take_step(F& f, const std::array<scalar,F::number_of_controls>& controls_ini, const std::array<scalar,F::number_of_controls>& controls_fin, std::array<scalar,F::number_of_inputs>& inputs, 
        scalar& time, scalar delta_time, Options opts)
    {
        constexpr const auto& number_of_inputs           = F::number_of_inputs;
        constexpr const auto& number_of_states           = F::number_of_states;
        constexpr const auto& number_of_controls         = F::number_of_controls;
        static_assert(number_of_inputs == number_of_states);

        // (1) Evaluate f at the initial point: CppAD variables are needed since operator() takes CppAD arrays

        // (1.1) Transform inputs to CppAD
        std::array<CppAD::AD<scalar>, number_of_inputs>   inputs_ini_cppad;
        std::array<CppAD::AD<scalar>, number_of_controls> controls_ini_cppad;

        std::copy(inputs.cbegin(), inputs.cend(), inputs_ini_cppad.begin());
        std::copy(controls_ini.cbegin(), controls_ini.cend(), controls_ini_cppad.begin());
        
        // (1.2) Evaluate functor
        const auto [states_ini_cppad, dstates_dt_ini_cppad] 
            = f(inputs_ini_cppad,controls_ini_cppad, time);

        // (1.3) Return CppAD outputs back to scalar
        std::array<scalar, number_of_states> states_ini;
        std::array<scalar, number_of_states> dstates_dt_ini;

        std::transform(states_ini_cppad.cbegin(), states_ini_cppad.cend(), states_ini.begin(), [](const auto& q) -> auto { return Value(q); });
        std::transform(dstates_dt_ini_cppad.cbegin(), dstates_dt_ini_cppad.cend(), dstates_dt_ini.begin(), [](const auto& q) -> auto { return Value(q); });

        // (2) Start from the initial point 
        auto inputs_fin     = inputs;

        // (3) Iterate
        bool success = false;
        for (size_t iter = 0; iter < opts.max_iter; ++iter)
        {
            // (3.1) Evaluate the ODE, and compute its Jacobian at the new point
            const auto f_and_jac = evaluate_f_and_jac(f, controls_fin, inputs_fin, time + delta_time); 

            // (3.2) Aliases to values
            const auto& states_fin = f_and_jac.states;
            const auto& dstates_dt_fin = f_and_jac.dstates_dt;

            // (3.3) Aliases to Jacobians
            const auto& jac_states_wrt_inputs_fin                  = f_and_jac.jac_states_wrt_inputs;
            const auto& jac_dstates_dt_wrt_inputs_fin              = f_and_jac.jac_dstates_dt_wrt_inputs;

            // (3.1.1) Check errors
            scalar max_error = 0.0;
            for (size_t i = 0; i < number_of_states; ++i)
            {
                max_error = max(max_error, std::abs(states_fin[i] - states_ini[i] - delta_time*((1.0-opts.sigma)*dstates_dt_ini[i] + opts.sigma*dstates_dt_fin[i])));
            }

            if ( max_error < opts.error_tolerance )
            {
                success = true;
                break;
            }

            // (3.2) Assemble rhs
            std::array<scalar,number_of_inputs> rhs;
    
            // (3.2.1) For ODEs, rhs = q + 0.5.dt.dqdt0 + 0.5.dt.dqdt_new - 0.5.dt.jac_dqdt_q_new.q_new - 0.5.dt.jac_dqdt_qa_new.qa_new
            for (size_t i = 0; i < number_of_states; ++i)
            {
                rhs[i] = -states_fin[i] + states_ini[i] + delta_time*((1.0-opts.sigma)*dstates_dt_ini[i] + opts.sigma*dstates_dt_fin[i]);
            }

            // (3.3) Assemble A
            //
            //               |  jac_states_wrt_inputs - 0.5.dt.jac_dqdt_q       - 0.5.dt.jac_dqdt_qa |
            //          A =  |                                                                             |
            //               |    - jac_dqa_q                                         - jac_dqa_qa         |
            //               
            std::array<scalar,number_of_inputs*number_of_inputs> A;

            // (3.3.1) Block A_qq
            for (size_t j = 0; j < number_of_inputs; ++j)
                for (size_t i = 0; i < number_of_states; ++i)
                    A[i + (number_of_inputs)*j] = jac_states_wrt_inputs_fin[i+ number_of_states * j] 
                        - opts.sigma * delta_time * jac_dstates_dt_wrt_inputs_fin[i + number_of_states * j];

            // (3.4) Solve the linear system
            lusolve(rhs.data(), A.data(), number_of_inputs, 1);

            // (3.5) Update q_new and qa_new
            std::array<scalar,number_of_inputs> delta_inputs;
            std::copy(rhs.cbegin(), rhs.cbegin() + number_of_inputs, delta_inputs.begin());

            inputs_fin = inputs_fin + opts.relaxation_factor*(delta_inputs);
        }

        // (4) Check status
        if (!success)
            throw lion_exception("[ERROR] Crank-Nicolson::take_step -> maximum number of iterations exceeded");

        // (5) Copy solution
        std::copy(inputs_fin.begin(), inputs_fin.end(), inputs.begin());

        // (6) Update t
        time += delta_time;

        return;
    }

 private:

    struct Evaluate_f_and_jac
    {
        std::array<scalar,F::number_of_states> states;
        std::array<scalar,F::number_of_states> dstates_dt;
        
        std::array<scalar,F::number_of_states*F::number_of_inputs> jac_states_wrt_inputs;
        std::array<scalar,F::number_of_states*F::number_of_inputs> jac_dstates_dt_wrt_inputs;
    };

    static Evaluate_f_and_jac evaluate_f_and_jac(F& f, const std::array<scalar,F::number_of_controls>& controls, const std::array<scalar, F::number_of_inputs>& inputs, 
        scalar time)
    {
        constexpr const auto& number_of_inputs           = F::number_of_inputs;
        constexpr const auto& number_of_states           = F::number_of_states;
        constexpr const auto& number_of_controls         = F::number_of_controls;
        static_assert(number_of_inputs == number_of_states);

        // (1) Put the states into a single vector, which will be declared as independent variables
        std::vector<CppAD::AD<scalar>> x0(number_of_inputs);
        std::copy(inputs.cbegin(), inputs.cend(), x0.begin());
    
        // (2) Declare the contents of x0 as the independent variables
        CppAD::Independent(x0);
    
        // (3) Create new inputs to the operator() of the vehicle 
        std::array<CppAD::AD<scalar>,number_of_inputs>     inputs_cppad;
        std::array<CppAD::AD<scalar>,number_of_controls>    controls_cppad;

        std::copy(x0.cbegin()         , x0.cbegin() + number_of_inputs, inputs_cppad.begin());
        std::copy(controls.cbegin()   , controls.cend()            , controls_cppad.begin());
    
        // (4) Call operator(), transform arrays to vectors
        auto [states_cppad, dstates_dt_cppad] = f(inputs_cppad, controls_cppad, time);
        std::vector<CppAD::AD<scalar>> all_outputs_cppad(states_cppad.cbegin(), states_cppad.cend());
        all_outputs_cppad.insert(all_outputs_cppad.end(), dstates_dt_cppad.cbegin(), dstates_dt_cppad.cend());
    
        // (5) Create the AD functions and stop the recording
        CppAD::ADFun<scalar> f_cppad;
        f_cppad.Dependent(x0,all_outputs_cppad);
    
        // (6) Transform x0 to scalar, to evaluate the functions
        std::vector<scalar> x0_scalar(x0.size());
    
        for (size_t i = 0; i < x0.size(); ++i)
            x0_scalar[i] = Value(x0[i]);
    
        // (7) Evaluate y = f(q0,u0,0)
        auto all_outputs = f_cppad.Forward(0, x0_scalar);
        auto all_jacobians = f_cppad.Jacobian(x0_scalar);
        
        // (8) Fill the solution struct
        Evaluate_f_and_jac solution;

        // (8.1) Solution
        std::copy(all_outputs.cbegin()           , all_outputs.cbegin() + number_of_states  , solution.states.begin());
        std::copy(all_outputs.cbegin() + number_of_states  , all_outputs.cbegin() + 2*number_of_states, solution.dstates_dt.begin());
    
        // (8.2) Jacobians
        //          - The CppAD jacobian is sorted as [dy1/dx1, ..., dy1/dxn, dy2/dx1, ..., dy2/dxn, ...]
        // (row major, we use col major, its what MATLAB uses)
        for (size_t i = 0; i < number_of_states; ++i)
            for (size_t j = 0; j < number_of_inputs; ++j)
                solution.jac_states_wrt_inputs[i + number_of_states*j] = all_jacobians[j + number_of_inputs*i];

        for (size_t i = 0; i < number_of_states; ++i)
            for (size_t j = 0; j < number_of_inputs; ++j)
                solution.jac_dstates_dt_wrt_inputs[i + number_of_states*j] = all_jacobians[j + number_of_inputs*(i+number_of_states)];

        return solution;
    }
};


#endif
