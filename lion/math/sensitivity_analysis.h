#ifndef SENSITIVITY_ANALYSIS_H
#define SENSITIVITY_ANALYSIS_H


#include "lion/foundation/types.h"
#include "lion/thirdparty/include/cppad/cppad.hpp"
#include <vector>
#include "lion/math/matrix_extensions.h"
#include "lion/thirdparty/include/logger.hpp"

template<typename FG>
class Sensitivity_analysis
{
 public:

    struct Options
    {
        scalar max_error_dual_problem = 1.0e-9;
        bool skip_optimality_check = false;
    };


    //! Constructor that only checks optimality of a solution
    Sensitivity_analysis(const FG& fg, 
                     const std::vector<scalar>& x, 
                     const std::vector<scalar>& s,
                     const std::vector<scalar>& lambda, 
                     const std::vector<scalar>& zl, 
                     const std::vector<scalar>& zu, 
                     const std::vector<scalar>& vl, 
                     const std::vector<scalar>& vu, 
                     const std::vector<scalar>& x_lb, 
                     const std::vector<scalar>& x_ub,
                     const std::vector<scalar>& c_lb,  
                     const std::vector<scalar>& c_ub,
                     const Options& opts
                    )
    : optimality_check({.grad_f{},.id_not_ok{},.success{false}}), _fg(fg), _np(0), _n(x.size()), _nc(lambda.size()), 
      _p(), _x(x), _s(s), _lambda(lambda), 
      _zl(zl), _zu(zu), _vl(vl), _vu(vu), 
      _x_lb(x_lb), _x_ub(x_ub), _c_lb(c_lb), _c_ub(c_ub), _equality_constraints(_nc, false), _inequality_constraints_lb(_nc, false),
      _inequality_constraints_ub(_nc,false),
      _opts(opts)
    {
        check_inputs();

        classify_constraints();

        check_optimality();
    }

    //! Constructor that checks optimality of the solution and computes its gradient w.r.t. given parameters
    Sensitivity_analysis(const FG& fg, 
                     const std::vector<scalar>& p,
                     const std::vector<scalar>& x, 
                     const std::vector<scalar>& s,
                     const std::vector<scalar>& lambda, 
                     const std::vector<scalar>& zl, 
                     const std::vector<scalar>& zu, 
                     const std::vector<scalar>& vl, 
                     const std::vector<scalar>& vu, 
                     const std::vector<scalar>& x_lb, 
                     const std::vector<scalar>& x_ub,
                     const std::vector<scalar>& c_lb,  
                     const std::vector<scalar>& c_ub,
                     const Options& opts
                    )
    : optimality_check({.grad_f{},.id_not_ok{},.success{false}}), _fg(fg), _np(p.size()), _n(x.size()), _nc(lambda.size()), 
      _p(p), _x(x), _s(s), _lambda(lambda), 
      _zl(zl), _zu(zu), _vl(vl), _vu(vu), 
      _x_lb(x_lb), _x_ub(x_ub), _c_lb(c_lb), _c_ub(c_ub), _equality_constraints(_nc, false), _inequality_constraints_lb(_nc, false),
      _inequality_constraints_ub(_nc,false),
      _opts(opts)
    {
        check_inputs();

        classify_constraints();
    
        check_optimality();

        compute_sensitivity();
    }

    //! Helper class to group the results of the optimality check
    struct
    {
        std::vector<scalar> grad_f;
        std::vector<size_t> id_not_ok;
        bool success;
    } optimality_check;

    std::vector<std::vector<scalar>> dxdp;

 private:

    void check_inputs() const;

    void check_optimality();

    void classify_constraints();
    
    void compute_sensitivity();

    struct Sparsity_pattern
    {
        CppAD::vector<size_t> row_jac;        
        CppAD::vector<size_t> col_jac;        
        CppAD::vectorBool     pattern_jac;
        CppAD::vector<size_t> row_hes;        
        CppAD::vector<size_t> col_hes;        
        CppAD::vectorBool     pattern_hes;
    };

    static Sparsity_pattern compute_sparsity_pattern(const size_t n, const size_t nc, CppAD::ADFun<scalar>& fg_ad, const bool force_diagonal); 

    FG                   _fg;
    size_t               _np;
    size_t               _n;
    size_t               _nc;
    std::vector<scalar>  _p;
    std::vector<scalar>  _x;
    std::vector<scalar>  _s;
    std::vector<scalar>  _lambda;
    std::vector<scalar>  _zl;
    std::vector<scalar>  _zu;
    std::vector<scalar>  _vl;
    std::vector<scalar>  _vu;
    std::vector<scalar>  _x_lb;
    std::vector<scalar>  _x_ub;
    std::vector<scalar>  _c_lb;
    std::vector<scalar>  _c_ub;
    std::vector<bool>    _equality_constraints;
    std::vector<bool>    _inequality_constraints_lb;
    std::vector<bool>    _inequality_constraints_ub;
    std::vector<size_t>  _lb_inequality_positions;
    std::vector<size_t>  _ub_inequality_positions;
    size_t               _n_equality;
    size_t               _n_inequalities_lb;
    size_t               _n_inequalities_ub;
    size_t               _n_inequality;
    Options              _opts;

    std::vector<scalar>  _x_full;
    CppAD::ADFun<scalar> _fg_ad;
    CppAD::ADFun<scalar> _fg_full_adfun;

    //! Helper classes -----------------------------------------------------------

    class FG_full
    {
     public:
        using ADvector = typename FG::ADvector;

        FG_full(FG& fg_original, const size_t np, const size_t n, const size_t nc,
                const std::vector<scalar>& c_lb, 
                const std::vector<bool>& equality_constraints, const std::vector<bool>& inequality_constraints_lb,
                const std::vector<bool>& inequality_constraints_ub) 
        : _fg_original(fg_original), _np(np), _n(n), _nc(nc), _c_lb(c_lb),
          _equality_constraints(equality_constraints), _inequality_constraints_lb(inequality_constraints_lb),
          _inequality_constraints_ub(inequality_constraints_ub),
          _n_equality(std::count(_equality_constraints.cbegin(), _equality_constraints.cend(), true)),
          _n_inequalities_lb(std::count(_inequality_constraints_lb.cbegin(), _inequality_constraints_lb.cend(), true)),
          _n_inequalities_ub(std::count(_inequality_constraints_ub.cbegin(), _inequality_constraints_ub.cend(), true)),
          _n_inequality(nc - _n_equality)
        {}

        //! Functor evaluation
        //! the content of x_full is: [x, s, lambda, params]
        void operator()(ADvector& f_full, const ADvector& x_full)
        {
            assert(f_full.size() == 1);
            assert(x_full.size() == _np + _n + _n_equality + 2*_n_inequality);

            // (1) Unpack the full x into the original x, the slack variables, lambdas, and the parameters
            ADvector x_original(_n), s(_n_inequality), lambda(_n_equality + _n_inequality), p(_np);

            std::copy(x_full.begin()                     , x_full.begin() + _n                , x_original.begin());
            std::copy(x_full.begin() + _n                , x_full.begin() + _n + _n_inequality, s.begin());
            std::copy(x_full.begin() + _n + _n_inequality, x_full.end() - _np                 , lambda.begin());
            std::copy(x_full.end() - _np                 , x_full.end()                       , p.begin());

            ADvector fg_original(_nc+1);

            if ( p.size() == 0 )
            {
                _fg_original(fg_original, x_original);
            }
            else
            {
                _fg_original(fg_original, x_original, p);
            }

            // (2) Initialize the full f to the original f
            f_full[0] = fg_original[0];

            // (3) Add the constraints
            size_t slack_var_counter = 0;
            for (size_t i = 0; i < _nc; ++i)
            {
                if ( _equality_constraints[i] )
                {
                    f_full[0] += lambda[i]*(fg_original[i+1] - _c_lb[i]);
                }
                else
                {
                    f_full[0] += lambda[i]*(fg_original[i+1] - s[slack_var_counter]);
                    ++slack_var_counter;
                }
            }
        }

     private:
        FG                   _fg_original;
        size_t               _np;
        size_t               _n;
        size_t               _nc;
        std::vector<scalar>  _c_lb;
        std::vector<bool>    _equality_constraints;
        std::vector<bool>    _inequality_constraints_lb;
        std::vector<bool>    _inequality_constraints_ub;
        size_t               _n_equality;
        size_t               _n_inequalities_lb;
        size_t               _n_inequalities_ub;
        size_t               _n_inequality;
    };


};


#include "sensitivity_analysis.hpp"


#endif
