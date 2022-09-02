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
    //! @param[inout] f: ODE functor, [states,dstates_dt,algebraic_equations] = f(input_states,algebraic_states,controls,time)
    //! @param[inout] controls_ini: initial control vector
    //! @param[inout] controls_fin: final control vector
    //! @param[inout] input_states: state vector values before and after of the step
    //! @param[inout] algebraic_states: algebraic state vector values before and after of the step
    //! @param[inout] time: time before and after of the step
    //! @param[inout] delta_time: time step before and after of the step
    static void take_step(F& f, const std::array<scalar,F::NCONTROL>& controls_ini, const std::array<scalar,F::NCONTROL>& controls_fin, std::array<scalar,F::NSTATE>& input_states, 
        std::array<scalar,F::NALGEBRAIC>& algebraic_states, scalar& time, scalar delta_time, Options opts)
    {
        constexpr const auto NSTATE     = F::NSTATE;
        constexpr const auto NALGEBRAIC = F::NALGEBRAIC;
        constexpr const auto NCONTROL   = F::NCONTROL;

        // (1) Evaluate f at the initial point: CppAD variables are needed since operator() takes CppAD arrays

        // (1.1) Transform inputs to CppAD
        std::array<CppAD::AD<scalar>, NSTATE>     input_states_ini_cppad;
        std::array<CppAD::AD<scalar>, NALGEBRAIC> algebraic_states_ini_cppad;
        std::array<CppAD::AD<scalar>, NCONTROL>   controls_ini_cppad;

        std::copy(input_states.cbegin(), input_states.cend(), input_states_ini_cppad.begin());
        std::copy(algebraic_states.cbegin(), algebraic_states.cend(), algebraic_states_ini_cppad.begin());
        std::copy(controls_ini.cbegin(), controls_ini.cend(), controls_ini_cppad.begin());
        
        // (1.2) Evaluate functor
        const auto [states_ini_cppad, dstates_dt_ini_cppad, algebraic_equations_ini_cppad] 
            = f(input_states_ini_cppad,algebraic_states_ini_cppad, controls_ini_cppad, time);

        // (1.3) Return CppAD outputs back to scalar
        std::array<scalar, NSTATE> states_ini;
        std::array<scalar, NSTATE> dstates_dt_ini;
        std::array<scalar, NALGEBRAIC> algebraic_equations_ini;

        std::transform(states_ini_cppad.cbegin(), states_ini_cppad.cend(), states_ini.begin(), [](const auto& q) -> auto { return Value(q); });
        std::transform(dstates_dt_ini_cppad.cbegin(), dstates_dt_ini_cppad.cend(), dstates_dt_ini.begin(), [](const auto& q) -> auto { return Value(q); });
        std::transform(algebraic_equations_ini_cppad.cbegin(), algebraic_equations_ini_cppad.cend(), algebraic_equations_ini.begin(), [](const auto& q) -> auto { return Value(q); });

        // (2) Start from the initial point 
        auto input_states_fin     = input_states;
        auto algebraic_states_fin = algebraic_states;

        // (3) Iterate
        bool success = false;
        for (size_t iter = 0; iter < opts.max_iter; ++iter)
        {
            // (3.1) Evaluate the ODE, and compute its Jacobian at the new point
            const auto f_and_jac = evaluate_f_and_jac(f, controls_fin, input_states_fin, algebraic_states_fin, time + delta_time); 

            // (3.2) Aliases to values
            const auto& states_fin = f_and_jac.states;
            const auto& dstates_dt_fin = f_and_jac.dstates_dt;
            const auto& algebraic_equations_fin = f_and_jac.algebraic_equations;

            // (3.3) Aliases to Jacobians
            const auto& jac_states_wrt_input_states_fin                  = f_and_jac.jac_states_wrt_input_states;
            const auto& jac_dstates_dt_wrt_input_states_fin              = f_and_jac.jac_dstates_dt_wrt_input_states;
            const auto& jac_dstates_dt_wrt_algebraic_states_fin          = f_and_jac.jac_dstates_dt_wrt_algebraic_states;
            const auto& jac_algebraic_equations_wrt_input_states_fin     = f_and_jac.jac_algebraic_equations_wrt_input_states;
            const auto& jac_algebraic_equations_wrt_algebraic_states_fin = f_and_jac.jac_algebraic_equations_wrt_algebraic_states;

            // (3.1.1) Check errors
            scalar max_error = 0.0;
            for (size_t i = 0; i < NSTATE; ++i)
            {
                max_error = max(max_error, std::abs(states_fin[i] - states_ini[i] - delta_time*((1.0-opts.sigma)*dstates_dt_ini[i] + opts.sigma*dstates_dt_fin[i])));
            }

            for (size_t i = 0; i < NALGEBRAIC; ++i)
            {
                max_error = max(max_error, std::abs(algebraic_equations_fin[i]));
            }

            if ( max_error < opts.error_tolerance )
            {
                success = true;
                break;
            }

            // (3.2) Assemble rhs
            std::array<scalar,NSTATE+NALGEBRAIC> rhs;
    
            // (3.2.1) For ODEs, rhs = q + 0.5.dt.dqdt0 + 0.5.dt.dqdt_new - 0.5.dt.jac_dqdt_q_new.q_new - 0.5.dt.jac_dqdt_qa_new.qa_new
            for (size_t i = 0; i < NSTATE; ++i)
            {
                rhs[i] = -states_fin[i] + states_ini[i] + delta_time*((1.0-opts.sigma)*dstates_dt_ini[i] + opts.sigma*dstates_dt_fin[i]);
            }

            // (3.2.2) For AEs, rhs = dqa_new
            if constexpr (NALGEBRAIC > 0)
            {
                for (size_t i = 0; i < NALGEBRAIC; ++i)
                {
                    rhs[i + NSTATE] = algebraic_equations_fin[i];
                }
            }

            // (3.3) Assemble A
            //
            //               |  jac_states_wrt_input_states - 0.5.dt.jac_dqdt_q       - 0.5.dt.jac_dqdt_qa |
            //          A =  |                                                                             |
            //               |    - jac_dqa_q                                         - jac_dqa_qa         |
            //               
            std::array<scalar,(NSTATE+NALGEBRAIC)*(NSTATE+NALGEBRAIC)> A;

            // (3.3.1) Block A_qq
            for (size_t j = 0; j < NSTATE; ++j)
                for (size_t i = 0; i < NSTATE; ++i)
                    A[i + (NSTATE+NALGEBRAIC)*j] = jac_states_wrt_input_states_fin[i+ NSTATE * j] 
                        - opts.sigma * delta_time * jac_dstates_dt_wrt_input_states_fin[i + NSTATE * j];

            // (3.3.2) Block A_qqa
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                for (size_t i = 0; i < NSTATE; ++i)
                    A[i + (NSTATE+NALGEBRAIC)*(j + NSTATE)] = -opts.sigma*delta_time*jac_dstates_dt_wrt_algebraic_states_fin[i + NSTATE*j]; 

            // (3.3.3) Block A_qaq
            for (size_t j = 0; j < NSTATE; ++j)
                for (size_t i = 0; i < NALGEBRAIC; ++i)
                    A[i + NSTATE + (NSTATE+NALGEBRAIC)*j] = - jac_algebraic_equations_wrt_input_states_fin[i + NALGEBRAIC*j];

            // (3.3.3) Block A_qaqa
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                for (size_t i = 0; i < NALGEBRAIC; ++i)
                    A[i + NSTATE + (NSTATE+NALGEBRAIC)*(j + NSTATE)] = - jac_algebraic_equations_wrt_algebraic_states_fin[i + NALGEBRAIC*j];

    
            // (3.4) Solve the linear system
            lusolve(rhs.data(), A.data(), NSTATE+NALGEBRAIC, 1);

            // (3.5) Update q_new and qa_new
            std::array<scalar,NSTATE> delta_input_states;
            std::array<scalar,NALGEBRAIC> delta_algebraic_states;
            std::copy(rhs.cbegin(), rhs.cbegin() + NSTATE, delta_input_states.begin());
            std::copy(rhs.cbegin() + NSTATE, rhs.cend(), delta_algebraic_states.begin());

            input_states_fin = input_states_fin + opts.relaxation_factor*(delta_input_states);
            algebraic_states_fin = algebraic_states_fin + opts.relaxation_factor*(delta_algebraic_states);
        }

        // (4) Check status
        if (!success)
            throw lion_exception("[ERROR] Crank-Nicolson::take_step -> maximum number of iterations exceeded");

        // (5) Copy solution
        std::copy(input_states_fin.begin(), input_states_fin.end(), input_states.begin());
        std::copy(algebraic_states_fin.begin(), algebraic_states_fin.end(), algebraic_states.begin());

        // (6) Update t
        time += delta_time;

        return;
    }


    /*  
    static void take_step_ipopt(F& f, const std::array<scalar,F::NCONTROL> controls_ini, const std::array<scalar,F::NCONTROL>& controls_fin, 
        std::array<scalar,F::NSTATE>& input_states, std::array<scalar,F::NALGEBRAIC>& algebraic_states, scalar& time, scalar delta_time, Options opts)
    {
        const size_t n_total = NSTATE + NALGEBRAIC;

        std::vector<scalar> x_lb(n_total,-1.0e3);
        std::vector<scalar> x_ub(n_total,+1.0e3);
        std::vector<scalar> c_lb(n_total,0.0);
        std::vector<scalar> c_ub(n_total,0.0);

        std::vector<scalar> x0(n_total);
        std::copy(input_states.cbegin(), input_states.cend(), x0.begin());
        std::copy(algebraic_states.cbegin(), algebraic_states.cend(), x0.begin() + NSTATE);

        std::string ipoptoptions;
        // turn off any printing
        ipoptoptions += "Integer print_level  ";
        ipoptoptions += std::to_string(5);
        ipoptoptions += "\n";
        ipoptoptions += "Integer max_iter ";
        ipoptoptions += std::to_string(opts.max_iter);
        ipoptoptions += "\n";
        ipoptoptions += "String  sb           yes\n";
        ipoptoptions += "Sparse true forward\n";
        ipoptoptions += "Numeric tol          1e-10\n";
        ipoptoptions += "Numeric constr_viol_tol  1e-10\n";
        ipoptoptions += "Numeric acceptable_tol  1e-8\n";
    
        // place to return solution
        CppAD::ipopt_cppad_result<std::vector<scalar>> result;

        Crank_nicolson_fitness fg(f, input_states, algebraic_states, controls_ini, time, delta_time, controls_fin, opts.sigma);
    
        // solve the problem
        CppAD::ipopt_cppad_solve<std::vector<scalar>, Crank_nicolson_fitness>(ipoptoptions, x0, x_lb, x_ub, c_lb, c_ub, fg, result);
    
        bool success = result.status == CppAD::ipopt_cppad_result<std::vector<scalar>>::success; 
    
        if ( !success )
        {
            throw lion_exception("Optimization did not succeed");
        }

        // Return the new solution
        std::copy(result.x.begin(), result.x.begin() + NSTATE, input_states.begin());
        std::copy(result.x.begin() + NSTATE, result.x.end(), algebraic_states.begin());
    }
    */

 private:

    struct Evaluate_f_and_jac
    {
        std::array<scalar,F::NSTATE> states;
        std::array<scalar,F::NSTATE> dstates_dt;
        std::array<scalar,F::NALGEBRAIC> algebraic_equations;
        
        std::array<scalar,F::NSTATE*F::NSTATE> jac_states_wrt_input_states;
        std::array<scalar,F::NSTATE*F::NSTATE> jac_dstates_dt_wrt_input_states;
        std::array<scalar,F::NALGEBRAIC*F::NSTATE> jac_algebraic_equations_wrt_input_states;

        std::array<scalar,F::NSTATE*F::NALGEBRAIC> jac_states_wrt_algebraic_states;
        std::array<scalar,F::NSTATE*F::NALGEBRAIC> jac_dstates_dt_wrt_algebraic_states;
        std::array<scalar,F::NALGEBRAIC*F::NALGEBRAIC> jac_algebraic_equations_wrt_algebraic_states;
    };

    static Evaluate_f_and_jac evaluate_f_and_jac(F& f, const std::array<scalar,F::NCONTROL>& controls, const std::array<scalar, F::NSTATE>& input_states, 
        const std::array<scalar, F::NALGEBRAIC>& algebraic_states, scalar time)
    {
        constexpr const auto NSTATE = F::NSTATE;
        constexpr const auto NALGEBRAIC = F::NALGEBRAIC;
        constexpr const auto NCONTROL = F::NCONTROL;

        // (1) Put the states into a single vector, which will be declared as independent variables
        constexpr const size_t n_total = NSTATE + NALGEBRAIC;
        std::vector<CppAD::AD<scalar>> x0(n_total);
        std::copy(input_states.cbegin(), input_states.cend(), x0.begin());
        std::copy(algebraic_states.cbegin(), algebraic_states.cend(), x0.begin() + NSTATE);
    
        // (2) Declare the contents of x0 as the independent variables
        CppAD::Independent(x0);
    
        // (3) Create new inputs to the operator() of the vehicle 
        std::array<CppAD::AD<scalar>,NSTATE>     input_states_cppad;
        std::array<CppAD::AD<scalar>,NALGEBRAIC> algebraic_states_cppad;
        std::array<CppAD::AD<scalar>,NCONTROL>    controls_cppad;

        std::copy(x0.cbegin()         , x0.cbegin() + NSTATE, input_states_cppad.begin());
        std::copy(x0.cbegin() + NSTATE, x0.cend()           , algebraic_states_cppad.begin());
        std::copy(controls.cbegin()   , controls.cend()            , controls_cppad.begin());
    
        // (4) Call operator(), transform arrays to vectors
        auto [states_cppad, dstates_dt_cppad,algebraic_equations_cppad] = f(input_states_cppad, algebraic_states_cppad, controls_cppad, time);
        std::vector<CppAD::AD<scalar>> all_outputs_cppad(states_cppad.cbegin(), states_cppad.cend());
        all_outputs_cppad.insert(all_outputs_cppad.end(), dstates_dt_cppad.cbegin(), dstates_dt_cppad.cend());
        all_outputs_cppad.insert(all_outputs_cppad.end(), algebraic_equations_cppad.cbegin(), algebraic_equations_cppad.cend());
    
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
        std::copy(all_outputs.cbegin()           , all_outputs.cbegin() + NSTATE  , solution.states.begin());
        std::copy(all_outputs.cbegin() + NSTATE  , all_outputs.cbegin() + 2*NSTATE, solution.dstates_dt.begin());
        std::copy(all_outputs.cbegin() + 2*NSTATE, all_outputs.cend()             , solution.algebraic_equations.begin());
    
        // (8.2) Jacobians
        //          - The CppAD jacobian is sorted as [dy1/dx1, ..., dy1/dxn, dy2/dx1, ..., dy2/dxn, ...]
        // (row major, we use col major, its what MATLAB uses)
        for (size_t i = 0; i < NSTATE; ++i)
            for (size_t j = 0; j < NSTATE; ++j)
                solution.jac_states_wrt_input_states[i + NSTATE*j] = all_jacobians[j + n_total*i];

        for (size_t i = 0; i < NSTATE; ++i)
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                solution.jac_states_wrt_algebraic_states[i + NSTATE*j] = all_jacobians[j + NSTATE + n_total*i];

        for (size_t i = 0; i < NSTATE; ++i)
            for (size_t j = 0; j < NSTATE; ++j)
                solution.jac_dstates_dt_wrt_input_states[i + NSTATE*j] = all_jacobians[j + n_total*(i+NSTATE)];

        for (size_t i = 0; i < NSTATE; ++i)
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                solution.jac_dstates_dt_wrt_algebraic_states[i + NSTATE*j] = all_jacobians[j + NSTATE + n_total*(i+NSTATE)];

        for (size_t i = 0; i < NALGEBRAIC; ++i)
            for (size_t j = 0; j < NSTATE; ++j)
                solution.jac_algebraic_equations_wrt_algebraic_states[i + NALGEBRAIC*j] = all_jacobians[j + n_total*(i + 2*NSTATE)];

        for (size_t i = 0; i < NALGEBRAIC; ++i)
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                solution.jac_algebraic_equations_wrt_algebraic_states[i + NALGEBRAIC*j] = all_jacobians[j + NSTATE + n_total*(i + 2*NSTATE)];

        return solution;
    }

    /*
    class Crank_nicolson_fitness
    {
     public:
        using ADvector = std::vector<CppAD::AD<scalar>>;

        Crank_nicolson_fitness(F& f, const std::array<scalar,NSTATE>& q, const std::array<scalar,NALGEBRAIC>& qa, 
            const std::array<scalar,NCONTROL>& u_ini, const scalar t, const scalar dt, const std::array<scalar,NCONTROL>& u_fin, const scalar sigma) : _f(f), _q_ini(q),
            _qa_ini(qa), _u_ini(u_ini), _t(t), _dt(dt), _sigma(sigma) 
        {
            // Transform final controls to cppad vector
            std::copy(u_fin.cbegin(), u_fin.cend(), _u_fin.begin());

            // Evaluate initial dqdt
            std::array<CppAD::AD<scalar>,NSTATE> q_ini_cppad;
            std::array<CppAD::AD<scalar>,NALGEBRAIC> qa_ini_cppad;
            std::array<CppAD::AD<scalar>,NCONTROL> u_ini_cppad;

            std::copy(q.cbegin(), q.cend(), q_ini_cppad.begin());
            std::copy(qa.cbegin(), qa.cend(), qa_ini_cppad.begin());
            std::copy(u_ini.cbegin(), u_ini.cend(), u_ini_cppad.begin());

            auto [dqdt_ini_cppad, dqa_ini] = _f(q_ini_cppad,qa_ini_cppad,u_ini_cppad,t);

            std::transform(dqdt_ini_cppad.begin(), dqdt_ini_cppad.end(), _dqdt_ini.begin(), [](const auto& dqdt_i) -> auto { return Value(dqdt_i); });
        }

        void operator()(ADvector& fg, const ADvector& x)
        {
            // Load q and qa
            std::array<CppAD::AD<scalar>,NSTATE> q_fin;
            std::array<CppAD::AD<scalar>,NALGEBRAIC> qa_fin;

            std::copy(x.begin(), x.begin() + NSTATE, q_fin.begin());
            std::copy(x.begin() + NSTATE, x.end(), qa_fin.begin());

            // Fitness function: minimize distance to the original point
            fg[0] = 0.0;
            for (size_t i = 0; i < NSTATE; ++i)
                fg[0] += (q_fin[i]-_q_ini[i])*(q_fin[i]-_q_ini[i]);

            for (size_t i = 0; i < NALGEBRAIC; ++i)
                fg[0] += (qa_fin[i]-_qa_ini[i])*(qa_fin[i]-_qa_ini[i]);
    

            auto [dqdt, dqa] = _f(q_fin, qa_fin, _u_fin, _t + _dt);

            for (size_t i = 0; i < NSTATE; ++i)
                fg[i] = q_fin[i] - _q_ini[i] - _dt*((1.0-_sigma)*_dqdt_ini[i]+ _sigma*dqdt[i]);

            for (size_t i = 0; i < NALGEBRAIC; ++i)
                fg[i + NSTATE] = dqa[i];
        }

     private:
        F _f;
        std::array<scalar,NSTATE> _q_ini;
        std::array<scalar,NALGEBRAIC> _qa_ini;
        std::array<scalar,NCONTROL> _u_ini;
        scalar _t;
        scalar _dt;
        scalar _sigma;
        std::array<scalar,NSTATE> _dqdt_ini;

        std::array<CppAD::AD<scalar>,NCONTROL> _u_fin;

    };
    */

 private:
};


#endif
