#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include <array>
#include "lion/foundation/types.h"
#include "lion/math/linear_algebra.h"
#include "lion/math/ipopt_cppad_handler.hpp"
#include "lion/math/matrix_extensions.h"
#include <cppad/cppad.hpp>

template<typename F, size_t NSTATE, size_t NALGEBRAIC, size_t NCONTROL>
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
    //! @param[inout] f: ODE functor, [dqdt,dqa] = f(q,qa,u,t)
    //! @param[inout] u_ini: initial control vector
    //! @param[inout] u_fin: final control vector
    //! @param[inout] q: state vector values before and after of the step
    //! @param[inout] qa: algebraic state vector values before and after of the step
    //! @param[inout] t: time before and after of the step
    //! @param[inout] dt: time step before and after of the step
    static void take_step(F& f, const std::array<scalar,NCONTROL> u_ini, const std::array<scalar,NCONTROL>& u_fin, std::array<scalar,NSTATE>& q, std::array<scalar,NALGEBRAIC>& qa, scalar& t, scalar dt, Options opts)
    {
        // (1) Evaluate f at the initial point: CppAD variables are needed
        std::array<CppAD::AD<scalar>,NSTATE> q0_cppad;          std::copy(q.cbegin(), q.cend(), q0_cppad.begin());
        std::array<CppAD::AD<scalar>,NALGEBRAIC> qa0_cppad;     std::copy(qa.cbegin(), qa.cend(), qa0_cppad.begin());
        std::array<CppAD::AD<scalar>,NCONTROL> u0_cppad;        std::copy(u_ini.cbegin(), u_ini.cend(), u0_cppad.begin());
        
        const auto [dqdt0_cppad,dqa0_cppad] = f(q0_cppad,qa0_cppad,u0_cppad,t);

        std::array<scalar,NSTATE> dqdt0;    std::transform(dqdt0_cppad.cbegin(), dqdt0_cppad.cend(), dqdt0.begin(), [](const auto& q) -> auto { return Value(q); });
        std::array<scalar,NALGEBRAIC> dqa0;  std::transform(dqa0_cppad.cbegin(), dqa0_cppad.cend(), dqa0.begin(), [](const auto& q) -> auto { return Value(q); });

        // (2) Start from the initial point 
        std::array<scalar,NSTATE> q_new;
        std::array<scalar,NALGEBRAIC> qa_new;

        std::copy(q.cbegin(), q.cend(), q_new.begin());
        std::copy(qa.cbegin(), qa.cend(), qa_new.begin());

        // (3) Iterate
        bool success = false;
        for (size_t iter = 0; iter < opts.max_iter; ++iter)
        {
            // (3.1) Evaluate the ODE, and compute its Jacobian at the new point
            const auto f_and_jac = evaluate_f_and_jac(f, u_fin, q_new, qa_new, t+dt); 

            const auto& dqdt_new        = f_and_jac.dqdt;
            const auto& dqa_new         = f_and_jac.dqa;
            const auto& jac_dqdt_q_new  = f_and_jac.jac_dqdt_q;
            const auto& jac_dqdt_qa_new = f_and_jac.jac_dqdt_qa;
            const auto& jac_dqa_q_new   = f_and_jac.jac_dqa_q;
            const auto& jac_dqa_qa_new  = f_and_jac.jac_dqa_qa;

            // (3.1.1) Check errors
            scalar max_error = 0.0;
            for (size_t i = 0; i < NSTATE; ++i)
            {
                max_error = max(max_error, std::abs(q_new[i] - q[i] - dt*((1.0-opts.sigma)*dqdt0[i] + opts.sigma*dqdt_new[i])));
            }

            for (size_t i = 0; i < NALGEBRAIC; ++i)
            {
                max_error = max(max_error, std::abs(dqa_new[i]));
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
                rhs[i] = -q_new[i] + q[i] + dt*((1.0-opts.sigma)*dqdt0[i] + opts.sigma*dqdt_new[i]);
            }

            // (3.2.2) For AEs, rhs = dqa_new
            if constexpr (NALGEBRAIC > 0)
            {
                for (size_t i = 0; i < NALGEBRAIC; ++i)
                {
                    rhs[i + NSTATE] = dqa_new[i];
                }
            }

            // (3.3) Assemble A
            //
            //               |  I - 0.5.dt.jac_dqdt_q       - 0.5.dt.jac_dqdt_qa |
            //          A =  |                                                   |
            //               |    - jac_dqa_q               - jac_dqa_qa         |
            //               
            std::array<scalar,(NSTATE+NALGEBRAIC)*(NSTATE+NALGEBRAIC)> A;

            // (3.3.1) Block A_qq
            for (size_t j = 0; j < NSTATE; ++j)
                for (size_t i = 0; i < NSTATE; ++i)
                    A[i + (NSTATE+NALGEBRAIC)*j] = (i == j ? 1.0 : 0.0) - opts.sigma*dt*jac_dqdt_q_new[i + NSTATE*j]; 

            // (3.3.2) Block A_qqa
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                for (size_t i = 0; i < NSTATE; ++i)
                    A[i + (NSTATE+NALGEBRAIC)*(j + NSTATE)] = -opts.sigma*dt*jac_dqdt_qa_new[i + NSTATE*j]; 

            // (3.3.3) Block A_qaq
            for (size_t j = 0; j < NSTATE; ++j)
                for (size_t i = 0; i < NALGEBRAIC; ++i)
                    A[i + NSTATE + (NSTATE+NALGEBRAIC)*j] = - jac_dqa_q_new[i + NALGEBRAIC*j];

            // (3.3.3) Block A_qaqa
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                for (size_t i = 0; i < NALGEBRAIC; ++i)
                    A[i + NSTATE + (NSTATE+NALGEBRAIC)*(j + NSTATE)] = - jac_dqa_qa_new[i + NALGEBRAIC*j];

    
            // (3.4) Solve the linear system
            lusolve(rhs.data(), A.data(), NSTATE+NALGEBRAIC, 1);

            // (3.5) Update q_new and qa_new
            std::array<scalar,NSTATE> dq;
            std::array<scalar,NALGEBRAIC> dqa;
            std::copy(rhs.cbegin(), rhs.cbegin() + NSTATE, dq.begin());
            std::copy(rhs.cbegin() + NSTATE, rhs.cend(), dqa.begin());

            q_new = q_new + opts.relaxation_factor*(dq);
            qa_new = qa_new + opts.relaxation_factor*(dqa);
        }

        // (4) Check status
        if (!success)
            throw lion_exception("[ERROR] Crank-Nicolson::take_step -> maximum number of iterations exceeded");

        // (5) Copy solution
        std::copy(q_new.begin(), q_new.end(), q.begin());
        std::copy(qa_new.begin(), qa_new.end(), qa.begin());

        // (6) Update t
        t += dt;

        return;
    }


    static void take_step_ipopt(F& f, const std::array<scalar,NCONTROL> u_ini, const std::array<scalar,NCONTROL>& u_fin, std::array<scalar,NSTATE>& q, std::array<scalar,NALGEBRAIC>& qa, scalar& t, scalar dt, Options opts)
    {

        const size_t n_total = NSTATE + NALGEBRAIC;

        std::vector<scalar> x_lb(n_total,-1.0e3);
        std::vector<scalar> x_ub(n_total,+1.0e3);
        std::vector<scalar> c_lb(n_total,0.0);
        std::vector<scalar> c_ub(n_total,0.0);

        std::vector<scalar> x0(n_total);
        std::copy(q.cbegin(), q.cend(), x0.begin());
        std::copy(qa.cbegin(), qa.cend(), x0.begin() + NSTATE);

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

        Crank_nicolson_fitness fg(f, q, qa, u_ini, t, dt, u_fin, opts.sigma);
    
        // solve the problem
        CppAD::ipopt_cppad_solve<std::vector<scalar>, Crank_nicolson_fitness>(ipoptoptions, x0, x_lb, x_ub, c_lb, c_ub, fg, result);
    
        bool success = result.status == CppAD::ipopt_cppad_result<std::vector<scalar>>::success; 
    
        if ( !success )
        {
            throw lion_exception("Optimization did not succeed");
        }

        // Return the new solution
        std::copy(result.x.begin(), result.x.begin() + NSTATE, q.begin());
        std::copy(result.x.begin() + NSTATE, result.x.end(), qa.begin());
    }

 private:

    struct Evaluate_f_and_jac
    {
        std::array<scalar,NSTATE> dqdt;
        std::array<scalar,NALGEBRAIC> dqa;
        
        std::array<scalar,NSTATE*NSTATE> jac_dqdt_q;
        std::array<scalar,NSTATE*NALGEBRAIC> jac_dqdt_qa;
    
        std::array<scalar,NALGEBRAIC*NSTATE> jac_dqa_q;
        std::array<scalar,NALGEBRAIC*NALGEBRAIC> jac_dqa_qa;
    };

    static Evaluate_f_and_jac evaluate_f_and_jac(F& f, const std::array<scalar,NCONTROL>& u, const std::array<scalar,NSTATE>& q, const std::array<scalar,NALGEBRAIC>& qa, scalar t)
    {
        // (1) Put the states into a single vector, which will be declared as independent variables
        constexpr const size_t n_total = NSTATE + NALGEBRAIC;
        std::vector<CppAD::AD<scalar>> x0(n_total);
        std::copy(q.cbegin(), q.cend(), x0.begin());
        std::copy(qa.cbegin(), qa.cend(), x0.begin() + NSTATE);
    
        // (2) Declare the contents of x0 as the independent variables
        CppAD::Independent(x0);
    
        // (3) Create new inputs to the operator() of the vehicle 
        std::array<CppAD::AD<scalar>,NSTATE>     q0;
        std::array<CppAD::AD<scalar>,NALGEBRAIC> qa0;
        std::array<CppAD::AD<scalar>,NCONTROL>    u0;

        std::copy(x0.cbegin()         , x0.cbegin() + NSTATE, q0.begin());
        std::copy(x0.cbegin() + NSTATE, x0.cend()           , qa0.begin());
        std::copy(u.cbegin()          , u.cend()            , u0.begin());
    
        // (4) Call operator(), transform arrays to vectors
        auto [dqdt_out,dqa_out] = f(q0,qa0,u0,t);
        std::vector<CppAD::AD<scalar>> out_vector(dqdt_out.cbegin(), dqdt_out.cend());
        out_vector.insert(out_vector.end(), dqa_out.cbegin(), dqa_out.cend());
    
        // (5) Create the AD functions and stop the recording
        CppAD::ADFun<scalar> f_cppad;
        f_cppad.Dependent(x0,out_vector);
    
        // (6) Transform x0 to scalar, to evaluate the functions
        std::vector<scalar> x0_scalar(x0.size());
    
        for (size_t i = 0; i < x0.size(); ++i)
            x0_scalar[i] = Value(x0[i]);
    
        // (7) Evaluate y = f(q0,u0,0)
        auto out0 = f_cppad.Forward(0, x0_scalar);
        auto out0_jacobian = f_cppad.Jacobian(x0_scalar);
        
        // (8) Fill the solution struct
        Evaluate_f_and_jac solution;

        // (8.1) Solution
        std::copy(out0.cbegin(), out0.cbegin() + NSTATE, solution.dqdt.begin());
        std::copy(out0.cbegin() + NSTATE, out0.cend(),   solution.dqa.begin());
    
        // (8.2) Jacobians
        //          - The CppAD jacobian is sorted as [dy1/dx1, ..., dy1/dxn, dy2/dx1, ..., dy2/dxn, ...]
        for (size_t i = 0; i < NSTATE; ++i)
            for (size_t j = 0; j < NSTATE; ++j)
                solution.jac_dqdt_q[i + NSTATE*j] = out0_jacobian[j + n_total*i];

        for (size_t i = 0; i < NSTATE; ++i)
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                solution.jac_dqdt_qa[i + NSTATE*j] = out0_jacobian[j + NSTATE + n_total*i];

        for (size_t i = 0; i < NALGEBRAIC; ++i)
            for (size_t j = 0; j < NSTATE; ++j)
                solution.jac_dqa_q[i + NALGEBRAIC*j] = out0_jacobian[j + n_total*(i + NSTATE)];

        for (size_t i = 0; i < NALGEBRAIC; ++i)
            for (size_t j = 0; j < NALGEBRAIC; ++j)
                solution.jac_dqa_qa[i + NALGEBRAIC*j] = out0_jacobian[j + NSTATE + n_total*(i + NSTATE)];

        return solution;
    }

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

 private:
};


#endif
