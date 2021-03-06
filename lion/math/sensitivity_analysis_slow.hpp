#ifndef SENSITIVITY_ANALYSIS_SLOW_HPP
#define SENSITIVITY_ANALYSIS_SLOW_HPP

#include "lion/math/mumps_interface.h" 

template<typename FG>
void Sensitivity_analysis_slow<FG>::check_inputs() const
{
    if ( _zl.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis_slow<FG>::check_inputs() -> z_lb.size() != n");

    if ( _zu.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis_slow<FG>::check_inputs() -> z_ub.size() != n");

    if ( _x_lb.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis_slow<FG>::check_inputs() -> x_lb.size() != n");

    if ( _x_ub.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis_slow<FG>::check_inputs() -> x_ub.size() != n");

    if ( _c_lb.size() != _nc )
        throw lion_exception("[ERROR] void Sensitivity_analysis_slow<FG>::check_inputs() -> c_lb.size() != nc");

    if ( _c_ub.size() != _nc )
        throw lion_exception("[ERROR] void Sensitivity_analysis_slow<FG>::check_inputs() -> c_ub.size() != nc");
}



template<typename FG>
void Sensitivity_analysis_slow<FG>::check_optimality()
{
    // (1) Construct the full problem
    // (1.1) Construct the full fitness function
    FG_full fg_full(_fg, _np, _n, _nc, _c_lb, _equality_constraints,
                    _inequality_constraints_lb, _inequality_constraints_ub);


  
    // (1.2) Construct the full independent variables: [x, s, lambda, parameters]
    _x_full = std::vector<scalar>(_n + _n_equality + 2*_n_inequality + _np);

    // (1.1.1) Original variables, x
    std::copy(_x.cbegin(), _x.cend(), _x_full.begin());

    // (1.1.2) Slack variables, s
    size_t slack_var_counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] )
        {
            _x_full[_n + slack_var_counter] = _s[slack_var_counter];
            ++slack_var_counter;
        }
    }

    // (1.1.3) Lagrange multipliers
    for (size_t i = 0; i < _nc; ++i)
        _x_full[_n + _n_inequality + i] = _lambda[i];

    // (1.1.4) Parameters
    for (size_t i = 0; i < _np; ++i)
        _x_full[_n + _n_inequality + _nc + i] = _p[i];

    // (2) Construct the ADFun
    typename FG::ADvector _x_full_0(_x_full.size()), fg_full_0(1);
    std::copy(_x_full.cbegin(), _x_full.cend(), _x_full_0.begin());

    CppAD::Independent(_x_full_0);

    fg_full(fg_full_0, _x_full_0);

    _fg_full_adfun.Dependent(_x_full_0, fg_full_0);
    _fg_full_adfun.optimize();

    // (3) Evaluate the Jacobian
    _fg_full_adfun.Forward(0, _x_full);
    auto& grad_f = optimality_check.grad_f;
    grad_f = _fg_full_adfun.Reverse(1, std::vector<scalar>{1.0});

    // (3.1) Add the bound multipliers jacobian
    for (size_t i = 0; i < _n; ++i)
    {
        grad_f[i] -= _zl[i];
        grad_f[i] += _zu[i];
    }

    size_t counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( _inequality_constraints_lb[i] )
        {
            grad_f[_n+counter] -= _vl[counter];
            counter++;
        }
    }

    counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( _inequality_constraints_ub[i] )
        {
            grad_f[_n+counter] += _vu[counter];
            counter++;
        }
    }

    // (4) Look for errors and export the solution
    for (size_t i = 0; i < grad_f.size() - _np; ++i)
    {
        const auto& grad_f_i = grad_f[i];

        if ( std::abs(grad_f_i) > _opts.max_error_dual_problem )
            optimality_check.id_not_ok.push_back(i);
    }

    optimality_check.success = (optimality_check.id_not_ok.size() == 0);
}


template<typename FG>
void Sensitivity_analysis_slow<FG>::compute_sensitivity()
{
    // (1) Compute the Hessian of the extended problem: x_extended = [params, x, s, lambda]
    auto sparsity_patterns = compute_sparsity_pattern(_n + 2*_n_inequality + _n_equality + _np, 0, _fg_full_adfun, true);

    size_t n_hes = sparsity_patterns.row_hes.size();
    std::vector<scalar> w = {1.0};
    std::vector<scalar> hes(n_hes);
    CppAD::sparse_hessian_work      work_hes;
    _fg_full_adfun.SparseHessian(_x_full, w, sparsity_patterns.pattern_hes, sparsity_patterns.row_hes, sparsity_patterns.col_hes, hes, work_hes);

    // (2) Get the location of the main diagonal
    std::vector<size_t> diag_loc(_n + 2*_n_inequality + _n_equality + _np);
    std::vector<size_t> diag_set(_n + 2*_n_inequality + _n_equality + _np, false);

    for (size_t ij = 0; ij < n_hes; ++ij)
    {
        if ( sparsity_patterns.row_hes[ij] == sparsity_patterns.col_hes[ij] )
        {
            diag_loc[sparsity_patterns.row_hes[ij]] = ij;
            diag_set[sparsity_patterns.row_hes[ij]] = true;
        }
    }

    if ( std::count(diag_set.cbegin(), diag_set.cend(), false) > 0 )
        throw lion_exception("Some diagonal positions were not found");

    // (2) Add the bound multipliers Hessian

    // (2.1) Variables
    for (size_t ij = 0; ij < _n; ++ij)
    {
        hes[diag_loc[ij]] += _zl[ij]/(_x[ij] - _x_lb[ij]) + _zu[ij]/(_x_ub[ij] - _x[ij]); 
    }

    // (2.2) Slack variables
    for (size_t ij = 0; ij < _n_inequality; ++ij)
    {
        hes[diag_loc[_n+ij]] += _vl[ij]/(_opts.ipopt_bound_relax_factor + _s[ij] - _c_lb[_lb_inequality_positions[ij]]) + _vu[ij]/(_opts.ipopt_bound_relax_factor + _c_ub[_ub_inequality_positions[ij]] - _s[ij]);
    }

    // (3) Prepare the system to be solved with mumps
    size_t n_total = _n + 2*_n_inequality + _n_equality;
    std::vector<size_t> rows_lhs;
    std::vector<size_t> cols_lhs;
    std::vector<double> lhs;
    std::vector<std::vector<double>> rhs(_np,std::vector<double>(n_total,0.0));

    for (size_t ij = 0; ij < n_hes; ++ij)
    {
        if ( sparsity_patterns.row_hes[ij] < sparsity_patterns.col_hes[ij] )
            throw lion_exception("[ERROR] Matrix was expected to be lower triangular (i>=j)");

        if ( (sparsity_patterns.row_hes[ij] < n_total) )
        {
            // If indexes are less than n_total, it's the lhs
            rows_lhs.push_back(sparsity_patterns.row_hes[ij]);
            cols_lhs.push_back(sparsity_patterns.col_hes[ij]);
            lhs.push_back(hes[ij]);

        }
        else
        {
            size_t i_p = sparsity_patterns.row_hes[ij] - n_total;
            if ( sparsity_patterns.col_hes[ij] < n_total ) 
                rhs[i_p][sparsity_patterns.col_hes[ij]] = -hes[ij];
        }
    }

    // (4) Solve the system
    dxdp = mumps_solve_linear_system(n_total, lhs.size(), rows_lhs, cols_lhs, lhs, rhs, true);
}

template<typename FG>
void Sensitivity_analysis_slow<FG>::classify_constraints() 
{
    // (1) Check equality constraints
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( std::abs(_c_lb[i] - _c_ub[i]) < 1.0e-8 )
            _equality_constraints[i] = true;
    }

    // (2) Check lower bound inequalities
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] && _c_lb[i] > -1.0e18 ) 
        {
            _inequality_constraints_lb[i] = true;
            _lb_inequality_positions.push_back(i);
        }
    }

    // (3) Check upper bound inequalities
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] && _c_ub[i] <  1.0e18 ) 
        {
            _inequality_constraints_ub[i] = true;
            _ub_inequality_positions.push_back(i);
        }
    }

    // (4) Count them
    _n_equality = std::count(_equality_constraints.cbegin(), _equality_constraints.cend(), true);
    _n_inequalities_lb = std::count(_inequality_constraints_lb.cbegin(), _inequality_constraints_lb.cend(), true);
    _n_inequalities_ub = std::count(_inequality_constraints_ub.cbegin(), _inequality_constraints_ub.cend(), true);
    _n_inequality = _nc - _n_equality;

    assert(_lb_inequality_positions.size() == _n_inequalities_lb);
    assert(_ub_inequality_positions.size() == _n_inequalities_ub);
}


template<typename FG>
typename Sensitivity_analysis_slow<FG>::Sparsity_pattern Sensitivity_analysis_slow<FG>::compute_sparsity_pattern(const size_t n, const size_t nc, CppAD::ADFun<scalar>& fg_ad, const bool force_diagonal)
{
    using ADvector = typename FG::ADvector;

    size_t i, j;
    const size_t m = 1 + nc;
    //
    // -----------------------------------------------------------
    // Jacobian
    CppAD::vectorBool pattern_jac(m*n);
    if( n <= m )
    {   
        // (1.a) Compute sparsity pattern with the forward method 

        // number of bits that are packed into one unit in vectorBool
        const size_t n_column = CppAD::vectorBool::bit_per_unit();

        // sparsity patterns for current columns
        CppAD::vectorBool r(n * n_column), s(m * n_column);

        // compute the sparsity pattern n_column columns at a time
        size_t n_loop = (n - 1) / n_column + 1;
        for(size_t i_loop = 0; i_loop < n_loop; i_loop++)
        {   // starting column index for this iteration
            size_t i_column = i_loop * n_column;

            // pattern that picks out the appropriate columns
            for(i = 0; i < n; i++)
            {   for(j = 0; j < n_column; j++)
                    r[i * n_column + j] = (i == i_column + j);
            }
            s = fg_ad.ForSparseJac(n_column, r);

            // fill in the corresponding columns of total_sparsity
            for(i = 0; i < m; i++)
            {   for(j = 0; j < n_column; j++)
                {   if( i_column + j < n )
                        pattern_jac[i * n + i_column + j] =
                            s[i * n_column + j];
                }
            }
        }
    }
    else
    {   
        // (1.b) Compute sparsity pattern with the reverse mode

        // number of bits that are packed into one unit in vectorBool
        size_t n_row = CppAD::vectorBool::bit_per_unit();

        // sparsity patterns for current rows
        CppAD::vectorBool r(n_row * m), s(n_row * n);

        // compute the sparsity pattern n_row row at a time
        size_t n_loop = (m - 1) / n_row + 1;
        for(size_t i_loop = 0; i_loop < n_loop; i_loop++)
        {   // starting row index for this iteration
            size_t i_row = i_loop * n_row;

            // pattern that picks out the appropriate rows
            for(i = 0; i < n_row; i++)
            {   for(j = 0; j < m; j++)
                    r[i * m + j] = (i_row + i ==  j);
            }
            s = fg_ad.RevSparseJac(n_row, r);

            // fill in correspoding rows of total sparsity
            for(i = 0; i < n_row; i++)
            {   for(j = 0; j < n; j++)
                    if( i_row + i < m )
                        pattern_jac[ (i_row + i) * n + j ] =
                            s[ i  * n + j];
            }
        }
    }

    // Set row and column indices in Jacoian of [f(x), g(x)]
    // for Jacobian of g(x). These indices are in row major order.
    CppAD::vector<size_t> row_jac, col_jac;
    for(i = 1; i < 1+nc; i++)
    {   for(j = 0; j < n; j++)
        {   if( pattern_jac[ i * n + j ] )
            {   row_jac.push_back(i);
                col_jac.push_back(j);
            }
        }
    }

    // Set row and column indices in Jacoian of [f(x), g(x)]
    // for Jacobian of g(x). These indices are in row major order.
    // -----------------------------------------------------------
    // Hessian
    CppAD::vectorBool pattern_hes(n*n);

    // number of bits that are packed into one unit in vectorBool
    const size_t n_column = CppAD::vectorBool::bit_per_unit();

    // sparsity patterns for current columns
    CppAD::vectorBool r(n * n_column), h(n * n_column);

    // sparsity pattern for range space of function
    CppAD::vectorBool s(m);
    for(i = 0; i < m; i++)
        s[i] = true;

    // compute the sparsity pattern n_column columns at a time
    size_t n_loop = (n - 1) / n_column + 1;
    for(size_t i_loop = 0; i_loop < n_loop; i_loop++)
    {   // starting column index for this iteration
        size_t i_column = i_loop * n_column;

        // pattern that picks out the appropriate columns
        for(i = 0; i < n; i++)
        {   for(j = 0; j < n_column; j++)
                r[i * n_column + j] = (i == i_column + j);
        }
        fg_ad.ForSparseJac(n_column, r);

        // sparsity pattern corresponding to paritls w.r.t. (theta, u)
        // of partial w.r.t. the selected columns
        bool transpose = true;
        h = fg_ad.RevSparseHes(n_column, s, transpose);

        // fill in the corresponding columns of total_sparsity
        for(i = 0; i < n; i++)
        {   for(j = 0; j < n_column; j++)
            {   if( i_column + j < n )
                    pattern_hes[i * n + i_column + j] =
                        h[i * n_column + j];
            }
        }
    }

    // Force the diagonal
    for (i = 0; i < n; ++i)
        pattern_hes[i*n + i] = true;

    // Set row and column indices for Lower triangle of Hessian
    // of Lagragian.  These indices are in row major order.
    CppAD::vector<size_t> row_hes, col_hes;
    for(i = 0; i < n; i++)
    {   for(j = 0; j < n; j++)
        {   if( pattern_hes[ i * n + j ] )
            if( j <= i )
            {   row_hes.push_back(i);
                col_hes.push_back(j);
            }
        }
    }

    // Column order indirect sort of the Jacobian indices
    CppAD::vector<size_t>  col_order_jac(col_jac.size());
    index_sort( col_jac, col_order_jac );

    return (Sparsity_pattern)
    {
        .row_jac = row_jac,        
        .col_jac = col_jac,
        .pattern_jac = pattern_jac,
        .row_hes = row_hes,
        .col_hes = col_hes,        
        .pattern_hes = pattern_hes
    };
}



#endif
