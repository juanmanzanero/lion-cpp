#ifndef LIONCPP_MATH_OPTIMIZATION_RESULT_H
#define LIONCPP_MATH_OPTIMIZATION_RESULT_H

#include <coin-or/IpIpoptData.hpp>
#include <coin-or/IpDenseVector.hpp>

#include "lion/foundation/lion_exception.h"

namespace lioncpp {

    template <class double_vector_type>
    struct Optimization_result
    {
        //! possible values for the result status
        enum class status_type {
            not_defined,
            success,
            maxiter_exceeded,
            stop_at_tiny_step,
            stop_at_acceptable_point,
            local_infeasibility,
            user_requested_stop,
            feasible_point_found,
            diverging_iterates,
            restoration_failure,
            error_in_step_computation,
            invalid_number_detected,
            too_few_degrees_of_freedom,
            internal_error,
            unknown
        };

        //! A function to check if the received "x" matches the one Ipopt has stored internally
        static void check_solution_consistency(const Ipopt::IpoptData* ip_data, const double_vector_type& x, const double_vector_type& xl, const double_vector_type& xu);

        //! A function to retrieve Ipopt's slack variables and their bounds Lagrange multipliers
        struct Slack_and_bound_multipliers
        {
            double_vector_type s;
            double_vector_type vl;
            double_vector_type vu;
        };

        static auto compute_slack_variables_and_their_lagrange_multipliers(
            const Ipopt::IpoptData* ip_data,
            const double_vector_type& constraint_lower_bounds,
            const double_vector_type& constraint_upper_bounds)->Slack_and_bound_multipliers;

        status_type status = status_type::not_defined; //! solution status
        Ipopt::Number iter_count;                      //! The number of iterations
        double_vector_type x;                          //! the approximation solution
        double_vector_type zl;                         //! Lagrange multipliers corresponding to lower bounds on x
        double_vector_type zu;                         //! Lagrange multipliers corresponding to upper bounds on x
        double_vector_type g;                          //! value of g(x)
        double_vector_type lambda;                     //! Lagrange multipliers correspondiing constraints on g(x)
        double obj_value;                              //! value of f(x)
        double_vector_type s;                          //! slack variables
        double_vector_type vl;                         //! Lagrange multipliers corresponding to lower bounds on x
        double_vector_type vu;                         //! Lagrange multipliers corresponding to upper bounds on x
    };

    template<typename double_vector_type>
    inline void Optimization_result<double_vector_type>::
        check_solution_consistency(const Ipopt::IpoptData* ip_data, const double_vector_type& x, const double_vector_type& xl, const double_vector_type& xu)
    {
        size_t i_nz_ipopt_internal{ 0u };

        for (size_t j = 0; j < x.size(); ++j)
        {
            if (std::abs(xl[j] - xu[j]) > 2.0e-16)
            {
                if (std::abs(x[j] - (dynamic_cast<const Ipopt::DenseVector*>(GetRawPtr(ip_data->curr()->x()))->Values())[i_nz_ipopt_internal++]) > 1.0e-10)
                {
                    std::cerr << "[ERROR] x is not ip_data->curr()->x()" << std::endl;
                    return;
                    //throw lion_exception("x is not ip_data->curr()->x()");
                }
            }
            else
            {
                if (std::abs(x[j] - xl[j]) > 2.0e-16)
                {
                    std::cerr << "[ERROR] x is expected to be xl = xu" << std::endl;
                    return;
                    //throw lion_exception("x is expected to be xl = xu");
                }
            }
        }
    }


    template<typename double_vector_type>
    inline auto Optimization_result<double_vector_type>::
        compute_slack_variables_and_their_lagrange_multipliers(const Ipopt::IpoptData* ip_data, const double_vector_type& constraint_lower_bounds, const double_vector_type& constraint_upper_bounds) -> Slack_and_bound_multipliers
    {
        assert(constraint_lower_bounds.size() == constraint_upper_bounds.size());
        const auto m = constraint_lower_bounds.size();

        // Return slack variables. We must fetch them from Ipopt's internal storage.
        const Ipopt::Number* s_values = dynamic_cast<const Ipopt::DenseVector*>(GetRawPtr(ip_data->curr()->s()))->Values();
        const auto s = double_vector_type(s_values, s_values + ip_data->curr()->s()->Dim());

        // Return Lagrange multipliers for slack variable bounds
        // Internally, Ipopt only stores multipliers for active bounds (i.e. |bound| < 1.0e18),
        // so we must skip unilateral constraints
        const auto num_inequalities = s.size();
        auto vl = double_vector_type(num_inequalities, 0.0);
        auto vu = double_vector_type(num_inequalities, 0.0);

        const Ipopt::Number* vl_values = dynamic_cast<const Ipopt::DenseVector*>(GetRawPtr(ip_data->curr()->v_L()))->Values();
        const auto num_lower_constraint_bound = ip_data->curr()->v_L()->Dim();

        if (vl_values == nullptr && num_lower_constraint_bound > 0)
        {
            std::cerr << "[ERROR] Inconsistent lagrange multipliers for lower slack variable bounds. Pointer is nullptr but size is greater than zero (" << num_lower_constraint_bound << ")" << std::endl;;
            return {s, vl, vu};
        }

        const Ipopt::Number* vu_values = dynamic_cast<const Ipopt::DenseVector*>(GetRawPtr(ip_data->curr()->v_U()))->Values();
        const auto num_upper_constraint_bound = ip_data->curr()->v_U()->Dim();

        if (vu_values == nullptr && num_upper_constraint_bound > 0)
        {
            std::cerr << "[ERROR] Inconsistent lagrange multipliers for upper slack variable bounds. Pointer is nullptr but size is greater than zero (" << num_upper_constraint_bound << ")" << std::endl;;
            return {s, vl, vu};
        }

        size_t i_lower_inequality{ 0u };
        size_t i_upper_inequality{ 0u };
        size_t i_inequality{ 0u };
        for (size_t i_constraint = 0; i_constraint < m; ++i_constraint)
        {
            if (std::abs(constraint_lower_bounds[i_constraint] - constraint_upper_bounds[i_constraint]) > 2.0e-16)
            {
                // Upper and lower bounds differ, we have an inequality

                if (constraint_lower_bounds[i_constraint] > -1.0e18)
                {
                    // Lower bound is greater than -1e18, the bound is active
                    vl[i_inequality] = vl_values[i_lower_inequality];
                    ++i_lower_inequality;
                }

                if (constraint_upper_bounds[i_constraint] < 1.0e18)
                {
                    // Upper bound is lower than 1e18, the bound is active
                    vu[i_inequality] = vu_values[i_upper_inequality];
                    ++i_upper_inequality;
                }

                ++i_inequality;
            }
        }

        // Check that all constraints were considered
        if (i_lower_inequality != static_cast<std::size_t>(num_lower_constraint_bound))
        {
            std::cerr << "[ERROR] ipopt_cppad_handler::finalize_solution() -> inconsistent number of inequalities with lower bounds" << std::endl;
        }

        if (i_upper_inequality != static_cast<std::size_t>(num_upper_constraint_bound))
        {
            std::cerr << "[ERROR] ipopt_cppad_handler::finalize_solution() -> inconsistent number of inequalities with upper bounds" << std::endl;
        }

        return { s, vl, vu };
    }

}

#endif
