#include "lion/math/matrix_extensions.h"

template<typename F, typename C>
inline bool Ipopt_handler<F,C>::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index&  nnz_jac_g,
    Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    // Number of problem variables
    n = _n;

    // Number of constraints
    m = _nc;

    // Number of nonzeros in the jacobian
    nnz_jac_g = _n*_nc;

    nnz_h_lag = _n*(_n+1)/2;

    index_style = TNLP::C_STYLE;

    return true;
} 

template<typename F, typename C>
inline bool Ipopt_handler<F,C>::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, 
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


template<typename F, typename C>
inline bool Ipopt_handler<F,C>::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
    Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    for (size_t i = 0; i < _n; ++i)
        x[i] = _x0[i]; 

    return true;
}     


template<typename F, typename C>
inline bool Ipopt_handler<F,C>::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    typename F::argument_type x_v;
    if constexpr (std::is_same<typename F::argument_type, std::vector<double>>::value)
        x_v = std::vector<double>(x,x+_n);    

    else
        std::copy(x, x+_n, x_v.begin());

    obj_value = (*_f)(x_v);

    return true;
}


template<typename F, typename C>
inline bool Ipopt_handler<F,C>::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    typename F::argument_type x_v;
    if constexpr (std::is_same<typename F::argument_type, std::vector<double>>::value)
        x_v = std::vector<double>(x,x+_n);    

    else
        std::copy(x, x+_n, x_v.begin());
    
    // Compute the gradient using finite differences
    auto yi = (*_f)(x_v);
    for (size_t i = 0; i < _n; ++i)
    {
        const double h = _eps_jac*(1.0 + std::abs(x_v[i])); 
        auto xi_eps = x_v;
        xi_eps[i] += h;

        auto yi_peps = (*_f)(xi_eps);
        
        grad_f[i] = (yi_peps - yi)/h;
    }

    return true;
}


template<typename F, typename C>
inline bool Ipopt_handler<F,C>::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
                                                 Ipopt::Number* g)
{
    typename C::argument_type x_v;
    if constexpr (std::is_same<typename C::argument_type, std::vector<double>>::value)
        x_v = std::vector<double>(x,x+_n);    

    else
        std::copy(x, x+_n, x_v.begin());

    auto c = (*_c)(x_v);

    for (size_t i = 0; i < _nc; ++i)
        g[i] = c[i];

    return true;
}

template<typename F, typename C>
inline bool Ipopt_handler<F,C>::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
    Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    if ( values == NULL )
    {
        for (size_t i = 0; i < _n; ++i)
            for (size_t c = 0; c < _nc; ++c)
            {   
                iRow[i*_nc+c] = c;
                jCol[i*_nc+c] = i;
            }
    }
    else
    {
        typename C::argument_type x_v;
        if constexpr (std::is_same<typename C::argument_type, std::vector<double>>::value)
            x_v = std::vector<double>(x,x+_n);    

        else
            std::copy(x, x+_n, x_v.begin());

        for (size_t i = 0; i < _n; ++i)
        {
            auto yi = (*_c)(x_v);
            auto xi_eps = x_v;
            xi_eps[i] += _eps_jac;
    
            auto yi_peps = (*_c)(xi_eps);
            
            //xi_eps[i] -= 2.0*_eps_jac;
            //auto yi_meps = (*_c)(xi_eps);
    
            auto Jcol = (yi_peps - yi)*(1.0/_eps_jac);

            for (size_t c = 0; c < _nc; ++c)
                values[i*_nc+c] = Jcol[c];
        }
    }

    return true;
}

template<typename F, typename C>
inline bool Ipopt_handler<F,C>::eval_h(
           Ipopt::Index         n,
           const Ipopt::Number* x,
           bool                 new_x,
           Ipopt::Number        obj_factor,
           Ipopt::Index         m,
           const Ipopt::Number* lambda,
           bool                 new_lambda,
           Ipopt::Index         nele_hess,
           Ipopt::Index*        iRow,
           Ipopt::Index*        jCol,
           Ipopt::Number*       values
)
{
   assert(n == ((int)_n));
   assert(m == ((int)_nc));
   assert(nele_hess == ((int)(_n*(_n+1)/2)));

   if( values == NULL )
   {
      // return the structure. This is a symmetric matrix, fill the lower left
      // triangle only.

      // the hessian for this problem is actually dense
      Ipopt::Index idx = 0;
      for( Ipopt::Index row = 0; row < ((int)_n); row++ )
      {
         for( Ipopt::Index col = 0; col <= row; col++ )
         {
            iRow[idx] = row;
            jCol[idx] = col;
            idx++;
         }
      }

      assert(idx == nele_hess);
   }
   else
   {
      // return the values. This is a symmetric matrix, fill the lower left
      // triangle only
      typename F::argument_type x_v;
      if constexpr (std::is_same<typename F::argument_type, std::vector<double>>::value)
          x_v = std::vector<double>(x,x+_n);    

      else
          std::copy(x, x+_n, x_v.begin());

      // fill the objective portion
      Ipopt::Index idx = 0;
      for( Ipopt::Index row = 0; row < ((int)_n); row++ )
      {
         const double h_row = _eps_hess*(1.0 + std::abs(x_v[row]));
         for( Ipopt::Index col = 0; col <= row; col++ )
         {
            const double h_col = _eps_hess*(1.0 + std::abs(x_v[col]));

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

            // Evaluate f
            double hess_f = ((*_f)(xi_pp) - (*_f)(xi_pm) - (*_f)(xi_mp) + (*_f)(xi_mm))*(1.0/(4.0*h_row*h_col));

            // Evaluate c
            auto hess_c = ((*_c)(xi_pp) - (*_c)(xi_pm) - (*_c)(xi_mp) + (*_c)(xi_mm))*(1.0/(4.0*h_row*h_col));

            values[idx] = obj_factor*hess_f;  

            for (size_t i = 0; i < hess_c.size(); ++i)
                values[idx] += lambda[i]*hess_c[i];

            idx++;
         }
      }
      assert(idx == nele_hess);
   }

   return true;
}


/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
template<typename F, typename C>
inline void Ipopt_handler<F,C>::finalize_solution(
           Ipopt::SolverReturn               status,
           Ipopt::Index                      n,
           const Ipopt::Number*              x,
           const Ipopt::Number*              z_L,
           const Ipopt::Number*              z_U,
           Ipopt::Index                      m,
           const Ipopt::Number*              g,
           const Ipopt::Number*              lambda,
           Ipopt::Number                     obj_value,
           const Ipopt::IpoptData*           ip_data,
           Ipopt::IpoptCalculatedQuantities* ip_cq
        ) 
{
    std::copy(x, x+_n, _x.begin());
}


