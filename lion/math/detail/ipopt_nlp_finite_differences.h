#ifndef LIONCPP_MATH_DETAIL_IPOPT_NLP_FINITE_DIFFERENCES_H
#define LIONCPP_MATH_DETAIL_IPOPT_NLP_FINITE_DIFFERENCES_H


#include "lion/foundation/types.h"
#include "lion/thirdparty/include/coin-or/IpIpoptApplication.hpp"
#include "lion/math/matrix_extensions.h"
#include <math.h>


namespace lioncpp::detail {

template<typename FitnessType, typename ConstraintsType>
class Ipopt_NLP_finite_differences : public Ipopt::TNLP
{
    //
    // An Ipopt NLP that optimizes a fitness function of type FitnessType with constraints of type "ConstraintsType" that uses finite differences
    // to approximate the gradient of the fitness function and the Jacobian of the constraints
    //
    // All NLP matrices are considered dense in this very generic class
    //

public:

    using fitness_argument_type = typename std::decay_t<FitnessType>::argument_type;
    using constraints_argument_type = typename std::decay_t<ConstraintsType>::argument_type;


    //! Constructor
    Ipopt_NLP_finite_differences(const size_t n,
        const size_t nc,
        const std::vector<scalar>& x0,
        FitnessType &&f,
        ConstraintsType &&c,
        const std::vector<scalar>& x_lb,
        const std::vector<scalar>& x_ub,
        const std::vector<scalar>& c_lb,
        const std::vector<scalar>& c_ub,
        scalar eps_jac,
        scalar eps_hess,
        bool central_finite_differences)
        : _n(n), _nc(nc), _x0(x0), _x(_n, 0.0), _f{ f }, _c{ c }, _x_lb(x_lb), _x_ub(x_ub), _c_lb(c_lb), _c_ub(c_ub),
        _eps_jac{ eps_jac }, _eps_hess{ eps_hess }, _central_finite_differences{ central_finite_differences }
    {}

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    );

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(
        Ipopt::Index   n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index   m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    );

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(
        Ipopt::Index   n,
        bool           init_x,
        Ipopt::Number* x,
        bool           init_z,
        Ipopt::Number* z_L,
        Ipopt::Number* z_U,
        Ipopt::Index   m,
        bool           init_lambda,
        Ipopt::Number* lambda
    );

    /** Method to return the objective value */
    virtual bool eval_f(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool                 new_x,
        Ipopt::Number& obj_value
    );

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool                 new_x,
        Ipopt::Number* grad_f
    );

    /** Method to return the constraint residuals */
    virtual bool eval_g(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool                 new_x,
        Ipopt::Index         m,
        Ipopt::Number* g
    );

    /** Method to return:
     *   1) The structure of the Jacobian (if "values" is NULL)
     *   2) The values of the Jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool          new_x,
        Ipopt::Index         m,
        Ipopt::Index         nele_jac,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values
    );


    virtual bool eval_h(
        Ipopt::Index         n,
        const Ipopt::Number* x,
        bool                 new_x,
        Ipopt::Number        obj_factor,
        Ipopt::Index         m,
        const Ipopt::Number* lambda,
        bool                 new_lambda,
        Ipopt::Index         nele_hess,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values
    );


    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
        Ipopt::SolverReturn               status,
        Ipopt::Index                      n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_L,
        const Ipopt::Number* z_U,
        Ipopt::Index                      m,
        const Ipopt::Number* g,
        const Ipopt::Number* lambda,
        Ipopt::Number                     obj_value,
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq
    );

    const auto& x() const { return _x; }

private:

    size_t _n;
    size_t _nc;
    std::vector<double> _x0;
    std::vector<double> _x;
    FitnessType &_f;
    ConstraintsType &_c;
    std::vector<double> _x_lb;
    std::vector<double> _x_ub;
    std::vector<double> _c_lb;
    std::vector<double> _c_ub;

    scalar _eps_jac;
    scalar _eps_hess;
    bool _central_finite_differences;
};


template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
    Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    // Number of problem variables
    n = _n;

    // Number of constraints
    m = _nc;

    // Number of nonzeros in the jacobian
    nnz_jac_g = _n * _nc;

    nnz_h_lag = _n * (_n + 1) / 2;

    index_style = TNLP::C_STYLE;

    return true;
}

template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
    Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    assert(n == int(_n));
    assert(m == int(_nc));

    for (size_t i = 0; i < _n; ++i)
    {
        x_l[i] = _x_lb[i];
        x_u[i] = _x_ub[i];
    }

    for (size_t i = 0; i < _nc; ++i)
    {
        g_l[i] = _c_lb[i];
        g_u[i] = _c_ub[i];
    }

    return true;
}


template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
    Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    for (size_t i = 0; i < _n; ++i) {
        x[i] = _x0[i];
    }

    return true;
}


template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    fitness_argument_type x_v;
    if constexpr (is_specialization_v<fitness_argument_type, std::vector>) {
        x_v.assign(x, x + _n);
    }
    else {
        std::copy(x, x + _n, x_v.begin());
    }

    obj_value = _f(x_v);

    return true;
}


template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    fitness_argument_type x_v;
    if constexpr (is_specialization_v<fitness_argument_type, std::vector>) {
        x_v.assign(x, x + _n);
    }
    else {
        std::copy(x, x + _n, x_v.begin());
    }

    // Compute the gradient using finite differences
    if (_central_finite_differences) {
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < _n; ++i) {
            const auto h = _eps_jac * (1.0 + std::abs(x_v[i]));
            auto xp = x_v;

            xp[i] += h;
            const auto yi_p = _f(xp);

            xp[i] = x_v[i] - h;
            const auto yi_m = _f(xp);

            grad_f[i] = (yi_p - yi_m) * (scalar { 0.5 } / h);
        }
    }
    else {
        auto yi = _f(x_v);
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < _n; ++i) {
            const double h = _eps_jac * (1.0 + std::abs(x_v[i]));
            auto xp = x_v;

            xp[i] += h;
            auto yi_p = _f(xp);

            grad_f[i] = (yi_p - yi) / h;
        }
    }

    return true;
}


template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
    Ipopt::Number* g)
{
    constraints_argument_type x_v;
    if constexpr (is_specialization_v<constraints_argument_type, std::vector>) {
        x_v.assign(x, x + _n);
    }
    else {
        std::copy(x, x + _n, x_v.begin());
    }

    const auto c = _c(x_v);
    for (size_t i = 0; i < _nc; ++i) {
        g[i] = c[i];
    }

    return true;
}


template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
    Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    if (values == NULL)
    {
        for (size_t i = 0; i < _n; ++i)
            for (size_t c = 0; c < _nc; ++c)
            {
                iRow[i * _nc + c] = c;
                jCol[i * _nc + c] = i;
            }
    }
    else
    {
        constraints_argument_type x_v;
        if constexpr (is_specialization_v<constraints_argument_type, std::vector>) {
            x_v.assign(x, x + _n);
        }
        else {
            std::copy(x, x + _n, x_v.begin());
        }

        if (_central_finite_differences) {
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < _n; ++i) {
                auto xp = x_v;

                xp[i] += _eps_jac;
                auto yi_p = _c(xp);

                xp[i] = x_v[i] - _eps_jac;
                auto yi_m = _c(xp);

                auto Jcol = (yi_p - yi_m) * (scalar{ 0.5 } / _eps_jac);

                for (size_t c = 0; c < _nc; ++c)
                    values[i * _nc + c] = Jcol[c];
            }
        }
        else {
            const auto yi = _c(x_v);

#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < _n; ++i) {
                auto xp = x_v;

                xp[i] += _eps_jac;
                auto yi_p = _c(xp);

                auto Jcol = (yi_p - yi) * (scalar{ 1 } / _eps_jac);

                for (size_t c = 0; c < _nc; ++c)
                    values[i * _nc + c] = Jcol[c];
            }
        }
    }

    return true;
}

template<typename FitnessType, typename ConstraintsType>
inline bool Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::eval_h(
    Ipopt::Index         n,
    const Ipopt::Number* x,
    bool                 new_x,
    Ipopt::Number        obj_factor,
    Ipopt::Index         m,
    const Ipopt::Number* lambda,
    bool                 new_lambda,
    Ipopt::Index         nele_hess,
    Ipopt::Index* iRow,
    Ipopt::Index* jCol,
    Ipopt::Number* values
)
{
    assert(n == ((int)_n));
    assert(m == ((int)_nc));
    assert(nele_hess == ((int)(_n * (_n + 1) / 2)));

    if (values == NULL)
    {
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.

        // the hessian for this problem is actually dense
        Ipopt::Index idx = 0;
        for (Ipopt::Index row = 0; row < ((int)_n); ++row)
        {
            for (Ipopt::Index col = 0; col <= row; ++col)
            {
                iRow[idx] = row;
                jCol[idx] = col;
                ++idx;
            }
        }

        assert(idx == nele_hess);
    }
    else
    {
        // return the values. This is a symmetric matrix, fill the lower left
        // triangle only
        fitness_argument_type x_v;
        if constexpr (is_specialization_v<fitness_argument_type, std::vector>) {
            x_v.assign(x, x + _n);
        }
        else {
            std::copy(x, x + _n, x_v.begin());
        }

        // evaluate fitness & constraints at "x_v"
        const auto two_f_v = scalar{ 2 } * _f(x_v);
        const auto two_c_v = scalar{ 2 } * _c(x_v);

        // calculate the Hessian of the Lagrangian
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
        for (Ipopt::Index row = 0; row < ((int)_n); ++row)
        {
            const double h_row = _eps_hess * (1.0 + std::abs(x_v[row]));

            // central cross-derivatives (off-diagonal)
            for (Ipopt::Index col = 0; col < row; ++col)
            {
                const double h_col = _eps_hess * (1.0 + std::abs(x_v[col]));

                auto xi_pp = x_v;
                xi_pp[row] += h_row;
                xi_pp[col] += h_col;

                auto xi_pm = x_v;
                xi_pm[row] += h_row;
                xi_pm[col] -= h_col;

                auto xi_mp = x_v;
                xi_mp[row] -= h_row;
                xi_mp[col] += h_col;

                auto xi_mm = x_v;
                xi_mm[row] -= h_row;
                xi_mm[col] -= h_col;

                const auto hess_f = (_f(xi_pp) - _f(xi_pm) - _f(xi_mp) + _f(xi_mm)) * (1.0 / (4.0 * h_row * h_col));
                const auto hess_c = (_c(xi_pp) - _c(xi_pm) - _c(xi_mp) + _c(xi_mm)) * (1.0 / (4.0 * h_row * h_col));

                const auto idx = col + n * row;
                values[idx] = obj_factor * hess_f;

                for (size_t i = 0; i < _nc; ++i)
                    values[idx] += lambda[i] * hess_c[i];
            }

            // central second derivatives on-diagonal (row == col)
            auto xi_p = x_v;
            xi_p[row] += h_row;

            auto xi_m = x_v;
            xi_m[row] -= h_row;

            const auto inv_h_row2 = 1. / (h_row * h_row);
            const auto hess_f = (_f(xi_p) - two_f_v + _f(xi_m)) * inv_h_row2;
            const auto hess_c = (_c(xi_p) - two_c_v + _c(xi_m)) * inv_h_row2;

            const auto idx = row * (n + 1);
            values[idx] = obj_factor * hess_f;

            for (size_t i = 0; i < _nc; ++i)
                values[idx] += lambda[i] * hess_c[i];
        }
    }

    return true;
}


/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
template<typename FitnessType, typename ConstraintsType>
inline void Ipopt_NLP_finite_differences<FitnessType, ConstraintsType>::finalize_solution(
    Ipopt::SolverReturn               status,
    Ipopt::Index                      n,
    const Ipopt::Number* x,
    const Ipopt::Number* z_L,
    const Ipopt::Number* z_U,
    Ipopt::Index                      m,
    const Ipopt::Number* g,
    const Ipopt::Number* lambda,
    Ipopt::Number                     obj_value,
    const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq
)
{
    std::copy(x, x + _n, _x.begin());
}

} // end namespace lioncpp::detail

#endif
