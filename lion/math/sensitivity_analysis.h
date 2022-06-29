#ifndef SENSITIVITY_ANALYSIS_H
#define SENSITIVITY_ANALYSIS_H


#include "lion/foundation/types.h"
#include "lion/thirdparty/include/cppad/cppad.hpp"
#include <vector>
#include "lion/math/matrix_extensions.h"
#include "lion/thirdparty/include/logger.hpp"
#include "lion/foundation/lion_exception.h"

template<typename FG>
class Sensitivity_analysis
{
 public:

    struct Options
    {
        scalar max_error_dual_problem = 1.0e-9;
        scalar ipopt_bound_relax_factor = 1.0e-8;
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
    : optimality_check({.nlp_error{},.id_not_ok{},.success{false}}), _fg(fg), _np(0), _n(x.size()), _nc(lambda.size()), 
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
    : optimality_check({.nlp_error{},.id_not_ok{},.success{false}}), _fg(fg), _np(p.size()), _n(x.size()), _nc(lambda.size()), 
      _p(p), _x(x), _s(s), _lambda(lambda), 
      _zl(zl), _zu(zu), _vl(vl), _vu(vu), 
      _x_lb(x_lb), _x_ub(x_ub), _c_lb(c_lb), _c_ub(c_ub), _equality_constraints(_nc, false), _inequality_constraints_lb(_nc, false),
      _inequality_constraints_ub(_nc,false),
      _opts(opts)
    {
        check_inputs();

        classify_constraints();
    
        check_optimality();

        if ( p.size() > 0 )
            compute_sensitivity();
    }


    //! Helper class to group the results of the optimality check
    struct
    {
        std::vector<scalar> nlp_error;
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
        size_t                nnz_jac;
        CppAD::vector<size_t> row_jac;        
        CppAD::vector<size_t> col_jac;        
        CppAD::vectorBool     pattern_jac;
        size_t                nnz_hes;
        CppAD::vector<size_t> row_hes;        
        CppAD::vector<size_t> col_hes;        
        CppAD::vectorBool     pattern_hes;
    };

    static Sparsity_pattern compute_sparsity_pattern(const size_t n, const size_t nc, CppAD::ADFun<scalar>& fg_ad); 

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
    std::vector<size_t>  _inequality_positions;
    std::vector<size_t>  _lb_inequality_positions;
    std::vector<size_t>  _ub_inequality_positions;
    size_t               _n_equality;
    size_t               _n_inequalities_lb;
    size_t               _n_inequalities_ub;
    size_t               _n_inequality;
    Options              _opts;

    CppAD::ADFun<scalar> _fg_adfun;
    Sparsity_pattern     _sparsity_pattern;

    std::vector<scalar> _x_aug;
    std::vector<scalar> _dfdx;
    std::vector<scalar> _dfdp;
    std::vector<scalar> _dcdx;
    std::vector<size_t> _dcdx_rows;
    std::vector<size_t> _dcdx_cols;
    std::vector<scalar> _dcdp;
    std::vector<size_t> _dcdp_rows;
    std::vector<size_t> _dcdp_cols;
};


#include "sensitivity_analysis.hpp"


#endif
