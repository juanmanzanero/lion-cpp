#ifndef __IPOPT_HANDLER_H__
#define __IPOPT_HANDLER_H__

#include "lion/foundation/types.h"
#include "lion/thirdparty/include/coin-or/IpIpoptApplication.hpp"
#include <math.h>

template<typename F, typename C>
class Ipopt_handler : public Ipopt::TNLP
{
 public:
    //! Constructor
    Ipopt_handler(const size_t n, const size_t nc, const std::vector<scalar>& x0, F& f, C& c, const std::vector<double>& x_lb, 
             const std::vector<double>& x_ub, const std::vector<double>& c_lb, const std::vector<double>& c_ub): 
        _n(n), _nc(nc), _x0(x0), _x(_n,0.0), _f(&f), _c(&c), _x_lb(x_lb), _x_ub(x_ub), _c_lb(c_lb), _c_ub(c_ub) {};

    //! Virtual destructor
    virtual ~Ipopt_handler() {};

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(
       Ipopt::Index&          n,
       Ipopt::Index&          m,
       Ipopt::Index&          nnz_jac_g,
       Ipopt::Index&          nnz_h_lag,
       IndexStyleEnum&        index_style
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
       Ipopt::Number&       obj_value
    );
 
    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
       Ipopt::Index         n,
       const Ipopt::Number* x,
       bool                 new_x,
       Ipopt::Number*       grad_f
    );
 
    /** Method to return the constraint residuals */
    virtual bool eval_g(
       Ipopt::Index         n,
       const Ipopt::Number* x,
       bool                 new_x,
       Ipopt::Index         m,
       Ipopt::Number*       g
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
       Ipopt::Index*        iRow,
       Ipopt::Index*        jCol,
       Ipopt::Number*       values
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
       Ipopt::Index*        iRow,
       Ipopt::Index*        jCol,
       Ipopt::Number*       values
    );

 
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
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
    ); 

    const std::vector<double>& x() const { return _x; }

 private:
    size_t _n;
    size_t _nc;
    std::vector<double> _x0;
    std::vector<double> _x;
    F* _f;
    C* _c;
    std::vector<double> _x_lb;
    std::vector<double> _x_ub;
    std::vector<double> _c_lb;
    std::vector<double> _c_ub;

    scalar _eps_jac = sqrt(2.0e-16);
    scalar _eps_hess = cbrt(2.0e-16);
};

#include "ipopt_handler.hpp"

#endif
