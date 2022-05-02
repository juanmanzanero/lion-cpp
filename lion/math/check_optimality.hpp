#ifndef CHECK_OPTIMALITY_HPP
#define CHECK_OPTIMALITY_HPP

template<typename FG>
void Check_optimality<FG>::check_inputs() const
{
    if ( _zl.size() != _n )
        throw std::runtime_error("[ERROR] void Check_optimality<FG>::check_inputs() -> z_lb.size() != n");

    if ( _zu.size() != _n )
        throw std::runtime_error("[ERROR] void Check_optimality<FG>::check_inputs() -> z_ub.size() != n");

    if ( _x_lb.size() != _n )
        throw std::runtime_error("[ERROR] void Check_optimality<FG>::check_inputs() -> x_lb.size() != n");

    if ( _x_ub.size() != _n )
        throw std::runtime_error("[ERROR] void Check_optimality<FG>::check_inputs() -> x_ub.size() != n");

    if ( _c_lb.size() != _nc )
        throw std::runtime_error("[ERROR] void Check_optimality<FG>::check_inputs() -> c_lb.size() != nc");

    if ( _c_ub.size() != _nc )
        throw std::runtime_error("[ERROR] void Check_optimality<FG>::check_inputs() -> c_ub.size() != nc");
}



template<typename FG>
void Check_optimality<FG>::check_optimality()
{
    // (1) Classify the constraints in equalities and inequalities
    classify_constraints();

    // (4) Construct the full problem
    // (4.1) Construct the full fitness function
    FG_full fg_full(_fg, _n, _nc, _x_lb, _x_ub, _c_lb, _c_ub, _equality_constraints,
                    _inequality_constraints_lb, _inequality_constraints_ub);

  
    // (4.2) Construct the full independent variables: [x, s, lambda]
    std::vector<scalar> x_full(_n + _n_equality + 2*_n_inequality);

    std::copy(_x.cbegin(), _x.cend(), x_full.begin());

    size_t slack_var_counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] )
        {
            x_full[_n + slack_var_counter] = _s[slack_var_counter];
            ++slack_var_counter;
        }
    }

    for (size_t i = 0; i < _nc; ++i)
        x_full[_n + _n_inequality + i] = _lambda[i];

    // (4) Construct the ADFun
    typename FG::ADvector x_full_0(x_full.size()), fg_full_0(1);
    std::copy(x_full.cbegin(), x_full.cend(), x_full_0.begin());

    CppAD::Independent(x_full_0);

    fg_full(fg_full_0, x_full_0);

    CppAD::ADFun<scalar> fg_full_adfun;

    fg_full_adfun.Dependent(x_full_0, fg_full_0);

    // (5) Evaluate the Jacobian
    fg_full_adfun.Forward(0, x_full);

    grad_f = fg_full_adfun.Reverse(1, std::vector<scalar>{1.0});

    // (5.1) Add the bound multipliers jacobian
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

    for (size_t i = 0; i < grad_f.size(); ++i)
    {
        const auto& grad_f_i = grad_f[i];

        if ( std::abs(grad_f_i) > _opts.constraint_viol_tolerance )
            id_not_ok.push_back(i);
    }

    success = (id_not_ok.size() == 0);
}

template<typename FG>
void Check_optimality<FG>::classify_constraints() 
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
            _inequality_constraints_lb[i] = true;
    }

    // (3) Check upper bound inequalities
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] && _c_ub[i] <  1.0e18 ) 
            _inequality_constraints_ub[i] = true;
    }

    // (4) Count them
    _n_equality = std::count(_equality_constraints.cbegin(), _equality_constraints.cend(), true);
    _n_inequalities_lb = std::count(_inequality_constraints_lb.cbegin(), _inequality_constraints_lb.cend(), true);
    _n_inequalities_ub = std::count(_inequality_constraints_ub.cbegin(), _inequality_constraints_ub.cend(), true);
    _n_inequality = _nc - _n_equality;
}


template<typename FG>
typename Check_optimality<FG>::Sparsity_pattern Check_optimality<FG>::compute_sparsity_pattern(const size_t n, const size_t nc, CppAD::ADFun<scalar>& fg_ad)
{
    using ADvector = typename FG::ADvector;

    size_t i, j;
    const size_t m = n + nc;
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
        .col_jac = col_jac,
        .row_jac = row_jac,
        .pattern_jac = pattern_jac
    };
}



#endif
