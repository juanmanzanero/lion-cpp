#ifndef CHECK_OPTIMALITY_H
#define CHECK_OPTIMALITY_H


#include "lion/foundation/types.h"
#include "lion/thirdparty/include/cppad/cppad.hpp"
#include <vector>
#include "lion/math/matrix_extensions.h"
#include "lion/thirdparty/include/logger.hpp"

template<typename FG>
class Check_optimality
{
 public:

    struct Options
    {
        scalar constraint_viol_tolerance = 1.0e-9;
    };


    Check_optimality(const FG& fg, 
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
    : id_not_ok(), success(false), _fg(fg), _n(x.size()), _nc(lambda.size()), _x(x), _s(s), _lambda(lambda), _zl(zl), _zu(zu), _vl(vl), _vu(vu), 
      _x_lb(x_lb), _x_ub(x_ub), _c_lb(c_lb), _c_ub(c_ub), _equality_constraints(_nc, false), _inequality_constraints_lb(_nc, false),
      _inequality_constraints_ub(_nc,false),
      _opts(opts)
    {
        check_inputs();

        check_optimality();
    }

    std::vector<scalar> grad_f;
    std::vector<size_t> id_not_ok;
    bool success;

 private:

    class FG_full
    {
     public:
        using ADvector = typename FG::ADvector;

        FG_full(FG& fg_original, const size_t n, const size_t nc, const std::vector<scalar>& x_lb, const std::vector<scalar>& x_ub,
                const std::vector<scalar>& c_lb, const std::vector<scalar>& c_ub, 
                const std::vector<bool>& equality_constraints, const std::vector<bool>& inequality_constraints_lb,
                const std::vector<bool>& inequality_constraints_ub) 
        : _fg_original(fg_original), _n(n), _nc(nc), _x_lb(x_lb), _x_ub(x_ub), _c_lb(c_lb), _c_ub(c_ub),
          _equality_constraints(equality_constraints), _inequality_constraints_lb(inequality_constraints_lb),
          _inequality_constraints_ub(inequality_constraints_ub),
          _n_equality(std::count(_equality_constraints.cbegin(), _equality_constraints.cend(), true)),
          _n_inequalities_lb(std::count(_inequality_constraints_lb.cbegin(), _inequality_constraints_lb.cend(), true)),
          _n_inequalities_ub(std::count(_inequality_constraints_ub.cbegin(), _inequality_constraints_ub.cend(), true)),
          _n_inequality(nc - _n_equality)
        {}

        //! Functor evaluation
        //! the content of x_full is: [x, s, lambda]
        void operator()(ADvector& f_full, const ADvector& x_full)
        {
            assert(f_full.size() == 1);
            assert(x_full.size() == _n + _n_equality + 2*_n_inequality);

            // (1) Unpack the full x
            ADvector x_original(_n), s(_n_inequality), lambda(_n_equality + _n_inequality);

            std::copy(x_full.begin(), x_full.begin() + _n, x_original.begin());
            std::copy(x_full.begin() + _n, x_full.begin() + _n + _n_inequality, s.begin());
            std::copy(x_full.begin() + _n + _n_inequality, x_full.end(), lambda.begin());

            ADvector fg_original(_nc+1);

            _fg_original(fg_original, x_original);

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
        size_t               _n;
        size_t               _nc;
        std::vector<scalar>  _x_lb;
        std::vector<scalar>  _x_ub;
        std::vector<scalar>  _c_lb;
        std::vector<scalar>  _c_ub;
        std::vector<bool>    _equality_constraints;
        std::vector<bool>    _inequality_constraints_lb;
        std::vector<bool>    _inequality_constraints_ub;
        size_t               _n_equality;
        size_t               _n_inequalities_lb;
        size_t               _n_inequalities_ub;
        size_t               _n_inequality;
    };

    void check_inputs() const;

    void check_optimality();

    void classify_constraints();

    struct Sparsity_pattern
    {
        CppAD::vector<size_t> row_jac;        
        CppAD::vector<size_t> col_jac;        
        CppAD::vectorBool     pattern_jac;
    };

    static Sparsity_pattern compute_sparsity_pattern(const size_t n, const size_t nc, CppAD::ADFun<scalar>& fg_ad); 

    FG                   _fg;
    CppAD::ADFun<scalar> _fg_ad;
    size_t               _n;
    size_t               _nc;
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
    size_t               _n_equality;
    size_t               _n_inequalities_lb;
    size_t               _n_inequalities_ub;
    size_t               _n_inequality;
    Options              _opts;
};


#include "check_optimality.hpp"


#endif
