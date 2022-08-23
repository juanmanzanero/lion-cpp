#ifndef SENSITIVITY_ANALYSIS_HPP
#define SENSITIVITY_ANALYSIS_HPP

#include "lion/math/mumps_interface.h" 

template<typename FG>
void Sensitivity_analysis<FG>::check_inputs() const
{
    if ( _zl.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis<FG>::check_inputs() -> z_lb.size() != n");

    if ( _zu.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis<FG>::check_inputs() -> z_ub.size() != n");

    if ( _x_lb.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis<FG>::check_inputs() -> x_lb.size() != n");

    if ( _x_ub.size() != _n )
        throw lion_exception("[ERROR] void Sensitivity_analysis<FG>::check_inputs() -> x_ub.size() != n");

    if ( _c_lb.size() != _nc )
        throw lion_exception("[ERROR] void Sensitivity_analysis<FG>::check_inputs() -> c_lb.size() != nc");

    if ( _c_ub.size() != _nc )
        throw lion_exception("[ERROR] void Sensitivity_analysis<FG>::check_inputs() -> c_ub.size() != nc");
}



template<typename FG>
void Sensitivity_analysis<FG>::check_optimality()
{
    // (1) Construct the augmented unknowns vector with x_aug = [x,p]
    _x_aug.resize(_n + _np);
    std::copy(_x.cbegin(), _x.cend(), _x_aug.begin());
    std::copy(_p.cbegin(), _p.cend(), _x_aug.begin() + _n);

    // (2) Construct the CppAD function and record the operations
    typename FG::ADvector x_aug_cppad(_n + _np), fg_0(_nc+1);
    std::copy(_x_aug.cbegin(), _x_aug.cend(), x_aug_cppad.begin());
    CppAD::Independent(x_aug_cppad);

    // (2.1) Divide again the cppad vector into [x_cppad, p_cppad]
    typename FG::ADvector x_cppad(_n);
    typename FG::ADvector p_cppad(_np);
    std::copy(x_aug_cppad.begin(), x_aug_cppad.begin() + _n, x_cppad.begin());
    std::copy(x_aug_cppad.begin() + _n, x_aug_cppad.end(), p_cppad.begin());

    
    _fg(fg_0, x_cppad, p_cppad);
    _fg_adfun.Dependent(x_aug_cppad, fg_0);


    // (2) Compute the sparsity pattern
    _sparsity_pattern = compute_sparsity_pattern(_n+_np, _nc, _fg_adfun);

    // (3) Evaluate the gradient of f(x;p)
    _fg_adfun.Forward(0, _x_aug);
    std::vector<scalar> eqs_to_compute(_nc+1,0.0);
    eqs_to_compute.front() = 1.0;

    _fg_adfun.Forward(0, _x_aug);
    auto grad_f = _fg_adfun.Reverse(1, eqs_to_compute);
    _dfdx = std::vector<scalar>(grad_f.cbegin(), grad_f.cbegin() + _n);
    _dfdp = std::vector<scalar>(grad_f.cbegin() + _n, grad_f.cend());

    // (4) Evaluate the (sparse) gradient of c(x;p)
    std::vector<scalar> grad_c(_sparsity_pattern.nnz_jac);
    CppAD::sparse_jacobian_work work_jac;
    _fg_adfun.SparseJacobianForward(_x_aug , _sparsity_pattern.pattern_jac, _sparsity_pattern.row_jac, _sparsity_pattern.col_jac, grad_c, work_jac);

    // (4.1) Split the full Jacobian into derivatives w.r.t. x and w.r.t. p
    for (size_t i = 0; i < _sparsity_pattern.nnz_jac; ++i)
    {
        if ( _sparsity_pattern.col_jac[i] < _n )
        {
            _dcdx_rows.push_back(_sparsity_pattern.row_jac[i] - 1); 
            _dcdx_cols.push_back(_sparsity_pattern.col_jac[i]);
            _dcdx.push_back(grad_c[i]);
        }
        else
        {
            _dcdp_rows.push_back(_sparsity_pattern.row_jac[i] - 1);
            _dcdp_cols.push_back(_sparsity_pattern.col_jac[i] - _n);
            _dcdp.push_back(grad_c[i]);
        }
    }

    // (5) Compute the NLP error: [dfdx + sum(lambda.dcdx) - sum(zlb) + sum(zub), -lambda_ineq - sum(vlb) + sum(vub), eq ? c - c_lb : c - s]
    optimality_check.nlp_error.resize(_n + _n_equality + 2*_n_inequality);

    // (5.1) Compute gradient w.r.t. x
    
    // (5.1.1) Add gradient of fitness function
    std::copy(_dfdx.cbegin(), _dfdx.cend(), optimality_check.nlp_error.begin());

    // (5.1.2) Add gradient of the constraints 
    for (size_t i = 0; i < _dcdx.size(); ++i)
        optimality_check.nlp_error[_dcdx_cols[i]] += _lambda[_dcdx_rows[i]]*_dcdx[i];

    // (5.1.3) Add the bound multipliers
    for (size_t i = 0; i < _n; ++i)
        optimality_check.nlp_error[i] += _zu[i] - _zl[i];

    // (5.2) Compute gradient w.r.t. s
    
    // (5.2.1) Add the lagrange multipliers corresponding to inequalities
    for (size_t i = 0; i < _n_inequality; ++i)
        optimality_check.nlp_error[_n + i] = -_lambda[_inequality_positions[i]];

    // (5.2.2) Add the lower bound multipliers
    size_t ineq_counter = 0;
    size_t lb_ineq_counter = 0;
    size_t ub_ineq_counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] )
        {
            if ( _inequality_constraints_lb[i] )
            { 
                optimality_check.nlp_error[_n + ineq_counter] -= _vl[lb_ineq_counter];
                ++lb_ineq_counter;
            }

            if ( _inequality_constraints_ub[i] )
            {
                optimality_check.nlp_error[_n + ineq_counter] += _vu[ub_ineq_counter];
                ++ub_ineq_counter;
            }

            ++ineq_counter;
        }
    }

    assert(ineq_counter == _n_inequality);
    assert(lb_ineq_counter == _n_inequalities_lb);
    assert(ub_ineq_counter == _n_inequalities_ub);

    // (5.3) Compute gradient w.r.t. lambda
    ineq_counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( _equality_constraints[i] )
        {
            optimality_check.nlp_error[_n + _n_inequality + i] = Value(fg_0[i+1]) - _c_lb[i];
        }
        else
        {
            optimality_check.nlp_error[_n + _n_inequality + i] = Value(fg_0[i+1]) - _s[ineq_counter++]; 
        }
    }

    assert(ineq_counter == _n_inequality);

    // (6) Look for errors and export the solution
    for (auto it_error = optimality_check.nlp_error.cbegin(); it_error != optimality_check.nlp_error.cend(); ++it_error)
    {
        if ( std::abs(*it_error) > _opts.max_error_dual_problem )
            optimality_check.id_not_ok.push_back(std::distance(optimality_check.nlp_error.cbegin(), it_error));
    }

    optimality_check.success = (optimality_check.id_not_ok.size() == 0);
}


template<typename FG>
void Sensitivity_analysis<FG>::compute_sensitivity()
{
    // (1) Compute the Hessian matrix of f(x) + sum(lambda.c) 
    std::vector<scalar> w(1+_nc);
    w.front() = 1.0;
    std::copy(_lambda.cbegin(), _lambda.cend(), w.begin() + 1);
    std::vector<scalar> hes(_sparsity_pattern.nnz_hes);
    CppAD::sparse_hessian_work work_hes;
    _fg_adfun.SparseHessian(_x_aug, w, _sparsity_pattern.pattern_hes, _sparsity_pattern.row_hes, _sparsity_pattern.col_hes, hes, work_hes);

    // (1.1) Split the Hessian between Hxx and Hxp. Recall that CppAD computes a lower triangular matrix (i >= j)
    std::vector<size_t> hxx_rows;
    std::vector<size_t> hxx_cols;
    std::vector<scalar> hxx;

    std::vector<size_t> hxp_rows;
    std::vector<size_t> hxp_cols;
    std::vector<scalar> hxp;
    for (size_t i = 0; i < _sparsity_pattern.nnz_hes; ++i)
    {
        if ( _sparsity_pattern.row_hes[i] < _n )
        {
            hxx_rows.push_back(_sparsity_pattern.row_hes[i]);
            hxx_cols.push_back(_sparsity_pattern.col_hes[i]);
            hxx.push_back(hes[i]);
        }
        else if ( _sparsity_pattern.col_hes[i] < _n )
        {
            hxp_rows.push_back(_sparsity_pattern.row_hes[i] - _n);
            hxp_cols.push_back(_sparsity_pattern.col_hes[i]);
            hxp.push_back(hes[i]);
        }
    }

    // (2) Write the lhs matrix (Lower triangular)
    std::vector<size_t> lhs_rows;
    std::vector<size_t> lhs_cols;
    std::vector<scalar> lhs;

    // (2.1) Add fitness and constraints hessian
    lhs_rows = hxx_rows;
    lhs_cols = hxx_cols;
    lhs      = hxx;

    // (2.2) Add the lower and upper bound multipliers
    for (size_t i = 0; i < _n; ++i)
    {
        lhs_rows.push_back(i);
        lhs_cols.push_back(i);
        lhs.push_back(_zl[i]/(_x[i] - _x_lb[i]) + _zu[i]/(_x_ub[i] - _x[i]));
    }

    // (2.3) Add the Jacobian of the constraints
    for (size_t i = 0; i < _dcdx.size(); ++i)
    {
        lhs_rows.push_back(_n + _n_inequality + _dcdx_rows[i]);
        lhs_cols.push_back(_dcdx_cols[i]);
        lhs.push_back(_dcdx[i]);
    }

    // (2.4) add the slack vars lower bound multipliers
    size_t lb_ineq_counter = 0;
    size_t ineq_counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] )
        {
            if ( _inequality_constraints_lb[i] )
            {
                lhs_rows.push_back(_n + ineq_counter);
                lhs_cols.push_back(_n + ineq_counter);
                lhs.push_back(_vl[lb_ineq_counter]/(_opts.ipopt_bound_relax_factor + _s[ineq_counter] - _c_lb[i]));
                ++lb_ineq_counter;
            }
       
            ++ineq_counter;
        }
    }

    // (2.5) add the slack vars upper bound multipliers
    size_t ub_ineq_counter = 0;
    ineq_counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] )
        {
            if ( _inequality_constraints_lb[i] )
            {
                lhs_rows.push_back(_n + ineq_counter);
                lhs_cols.push_back(_n + ineq_counter);
                lhs.push_back(_vu[ub_ineq_counter]/(_opts.ipopt_bound_relax_factor + _c_ub[i] - _s[ineq_counter]));
                ++ub_ineq_counter;
            }
            ++ineq_counter;
        }
    }

    // (2.6) Add the identity matrix associated to -\partial(lambda_i.s_j)/\partial lambda_i \partial s_j
    ineq_counter = 0;
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( !_equality_constraints[i] )
        {
            lhs_rows.push_back(_n + _n_inequality + i); 
            lhs_cols.push_back(_n + ineq_counter);
            lhs.push_back(-1.0);
            ++ineq_counter;
        }
    }

    // (3) Write the vector of rhs (dense)
    size_t n_total = _n + 2*_n_inequality + _n_equality;
    std::vector<std::vector<scalar>> rhs(_np, std::vector<scalar>(n_total,0.0));

    for (size_t i = 0; i < hxp.size(); ++i)
        rhs[hxp_rows[i]][hxp_cols[i]] = -hxp[i];

    for (size_t i = 0; i < _dcdp.size(); ++i)
        rhs[_dcdp_cols[i]][_n + _n_inequality + _dcdp_rows[i]] = -_dcdp[i];

    // (4) Solve the system
    dxdp = mumps_solve_linear_system(n_total, lhs.size(), lhs_rows, lhs_cols, lhs, rhs, true);
}

template<typename FG>
void Sensitivity_analysis<FG>::classify_constraints() 
{
    // (1) Check equality constraints
    for (size_t i = 0; i < _nc; ++i)
    {
        if ( std::abs(_c_lb[i] - _c_ub[i]) < 1.0e-8 )
        {
            _equality_constraints[i] = true;
        }
        else
        {
            _inequality_positions.push_back(i);
        }
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
typename Sensitivity_analysis<FG>::Sparsity_pattern Sensitivity_analysis<FG>::compute_sparsity_pattern(const size_t n, const size_t nc, CppAD::ADFun<scalar>& fg_ad)
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

    return Sparsity_pattern
    {
        .nnz_jac = row_jac.size(),
        .row_jac = row_jac,        
        .col_jac = col_jac,
        .pattern_jac = pattern_jac,
        .nnz_hes = row_hes.size(),
        .row_hes = row_hes,
        .col_hes = col_hes,        
        .pattern_hes = pattern_hes
    };
}



#endif
