#ifndef LION_MATH_FSOLVE_TRUST_REGION_DODLEG_H
#define LION_MATH_FSOLVE_TRUST_REGION_DODLEG_H
#pragma once


#include <numeric>

#include "cppad/cppad.hpp"

#include "lion/math/linear_algebra.h"


//
// Defines the "fsolve_trust_region_dodleg" nonlinear
// solver for square systems, "inspired" in Matlab's
// "fsolve" function, to be employed with CppAD.
//


namespace CppAD {

template<typename Dvector>
struct fsolve_trust_region_dodleg_result
{
    //
    // Holds the results from calling the "fsolve_trust_region_dodleg"
    // nonlinear solver for square systems "g(x) = 0". This struct of
    // results attempts to mimick that of "CppAD::ipopt::solve()".
    //

    enum status_type
    {
        not_defined = -100,             // initial status, before calling the function

        success = 1,                    // equation solved, first-order optimality is small
        stop_at_tiny_step = 2,          // change in x smaller than the specified tolerance, or Jacobian at x is undefined
        stop_at_acceptable_point = 3,   // change in residual less than the specified tolerance

        maxiter_exceeded = 0,           // max. number of iterations exceeded
        invalid_number_detected = -1,   // NaN or Inf detected
        local_infeasibility = -2,       // equation not solved
        error_in_step_computation = -3, // trust region radius became too small

    };

    // possible values for solution status
    status_type status{ not_defined };

    // the approximation solution
    Dvector x;

    // value of "g(x)"
    Dvector g;

    // jacobian of "g(x)" at the solution point (stored in
    // column-major order)
    Dvector dg_dx_colmaj;

    // number of iterations of the method (here,
    // num_iters == num_fun_evals)
    std::size_t iter_count;
};


template<typename G, typename Dvector, typename... Args>
inline fsolve_trust_region_dodleg_result<Dvector>
    fsolve_trust_region_dodleg(G &&g_eval, bool retape,
        const Dvector &xi, Args&&... args)
{
    //
    // "trust-region-dodleg" method to solve the dense & unbounded
    // SQUARE system of nonlinear constraints "g(x) = 0" (in which
    // size(g) == size(x)), with seed vector "xi". The Jacobian of
    // the "g" equations is calculated via CppAD. These equations
    // MUST define a member typedef "ADvector" and a member function
    // "void operator()(ADvector &g, const ADvector &x)", i.e., the
    // standard ipopt interface for cost functions.
    //

    using scalar_type = typename Dvector::value_type;

    // declare a helper ADfun type that can retape each time
    // we call its member function "Forward"
    struct adfun_with_retape_capability : CppAD::ADFun<scalar_type>
    {
        using base_type = CppAD::ADFun<scalar_type>;
        using ADvector = typename std::decay_t<G>::ADvector;

        adfun_with_retape_capability(bool retape, G &&g_eval, const Dvector &xi)
        :   _retape{ retape },
            _g_eval{ g_eval },
            _a_x(xi.cbegin(), xi.cend()),
            _a_g(xi.size())
        {
            CppAD::Independent(_a_x);
            _g_eval(_a_g, _a_x);
            base_type::Dependent(_a_x, _a_g);
            if (!_retape) {
                base_type::optimize();
            }
        }

        auto Forward(std::size_t q, const Dvector &x)
        {
            if (_retape && q == 0u) {
                std::copy(x.cbegin(), x.cend(), _a_x.begin());
                CppAD::Independent(_a_x);
                _g_eval(_a_g, _a_x);
                base_type::Dependent(_a_x, _a_g);
            }

            return base_type::Forward(q, x);
        }

        const bool _retape;
        G &_g_eval;
        ADvector _a_x;
        ADvector _a_g;
    };

    // call the solver with the ADFun
    return fsolve_trust_region_dodleg(
        adfun_with_retape_capability(retape, std::forward<G>(g_eval), xi),
        xi, std::forward<Args>(args)...);
}


template<typename ADFunType, typename Dvector,
    typename ScalarType = typename Dvector::value_type>
inline fsolve_trust_region_dodleg_result<Dvector>
    fsolve_trust_region_dodleg(ADFunType &&adfun, const Dvector &xi,
        std::size_t max_iters = 1000u,
        ScalarType tolx = ScalarType{ 1e-10 },
        ScalarType tolf = ScalarType{ 1e-10 },
        Dvector typical_x = Dvector{},
        ScalarType DeltaMax = ScalarType{ 1e10 },
        ScalarType eta1 = ScalarType{ 0.05 },
        ScalarType eta2 = ScalarType{ 0.9 },
        ScalarType alpha1 = ScalarType{ 2.5 },
        ScalarType alpha2 = ScalarType{ 0.25 })
{
    //
    // "trust-region-dodleg" method to solve a dense & unbounded
    // SQUARE system of nonlinear constraints "g(x) = 0" (in which
    // size(g) == size(x)), with seed vector "xi". The equations
    // are represented by input "adfun", which MUST define a member
    // function "Dvector Forward(std::size_t q, const Dvector &x)"
    // to evaluate the equations and their first derivatives
    // w.r.t. "x".
    //

    using scalar_type = typename Dvector::value_type;

    // define helper "subroutines" to evaluate the function "g",
    // its dense-jacobian (which is stored in col-major order),
    // the method's gradient "transpose(jac) * g" and its infinity-norm
    const auto n = xi.size();
    Dvector xwork(n);
    Dvector fwork(n);

    const auto update_g_and_dense_jac = [&adfun, &xwork,
        &fwork](auto &g, auto &jac, const auto &x, auto n)
    {
        g = adfun.Forward(0, x);

        for (auto j = 0u; j < n; ++j) {
            xwork[j] = scalar_type{ 1 };
            fwork = adfun.Forward(1, xwork);
            for (auto i = 0u; i < n; ++i) {
                jac[i + n * j] = fwork[i];
            }
            xwork[j] = scalar_type{ 0 };
        }
    };

    const auto update_grad_and_normgradinf = [](auto &grad, auto &normgradinf,
        const auto &jac, const auto &g, auto n)
    {
       // grad = transpose(jac) * g
       for (auto i = 0u; i < n; ++i) {
           grad[i] = scalar_type{ 0 };
           for (auto j = 0u; j < n; ++j) {
               grad[i] += jac[j + n * i] * g[j];
           }
       }

       // normgradinf = max(abs(grad))
       normgradinf = std::abs(*std::max_element(grad.cbegin(), grad.cend(),
           [](auto a, auto b) { return std::abs(a) < std::abs(b); }));
    };

    // initialize & test convergence at initial point
    fsolve_trust_region_dodleg_result<Dvector> result;
    result.iter_count = 1u;
    result.x = xi;
    result.g.resize(n);
    result.dg_dx_colmaj.resize(n * n);
    update_g_and_dense_jac(result.g, result.dg_dx_colmaj, result.x, n);

    auto Delta = scalar_type{ 1 };
    auto step_accept = true;
    Dvector d(n, scalar_type{ 0 });
    auto normd = scalar_type{ 0 };
    Dvector grad(n);
    scalar_type normgradinf;
    update_grad_and_normgradinf(grad, normgradinf, result.dg_dx_colmaj, result.g, n);

    auto objold = scalar_type{ 1 };
    auto obj = scalar_type{ 0.5 } *
        std::inner_product(result.g.cbegin(), result.g.cend(), result.g.cbegin(),
            scalar_type{ 0 });

    const auto test_stop = [](auto &status,
        auto normgradinf, auto tolf, auto tolx,
        auto step_accept, auto iter_count, auto max_iters, auto Delta, auto normd,
        auto obj, auto objold, const auto &d, const auto &x, const auto &jac, auto n)
    {
        using status_type = std::decay_t<decltype(status)>;

        if (!std::all_of(x.cbegin(), x.cend(), [](auto x_i) { return std::isfinite(x_i); }) ||
            !std::all_of(d.cbegin(), d.cend(), [](auto d_i) { return std::isfinite(d_i); })) {

            status = status_type::invalid_number_detected;
            return true;
        }
        else if (!std::all_of(jac.cbegin(), jac.cend(), [](auto jac_ij) { return std::isfinite(jac_ij); })) {
            status = status_type::stop_at_tiny_step;
            return true;
        }
        else if (step_accept && normgradinf < tolf) {
            status = status_type::success;
            return true;
        }
        else {
            auto s = std::numeric_limits<scalar_type>::lowest();
            for (auto i = 0u; i < n; ++i) {
                s = std::max(s,
                    std::abs(d[i]) / (std::abs(x[i]) + scalar_type{ 1 }));
            }

            if (iter_count > 2u &&
                s < std::max(tolx * tolx, std::numeric_limits<scalar_type>::epsilon())) {
                if (scalar_type{ 2 } * obj < std::sqrt(tolf)) {
                    status = status_type::stop_at_tiny_step;
                }
                else {
                    status = status_type::local_infeasibility;
                }
                return true;
            }
            else if (iter_count > 2u &&
                step_accept &&
                normd < scalar_type{ 0.9 } * Delta &&
                std::abs(objold - obj) < std::max(tolf * tolf,
                    std::numeric_limits<scalar_type>::epsilon()) *
                    (scalar_type{ 1 } + std::abs(objold))) {

                if (scalar_type{ 2 } * obj < std::sqrt(tolf)) {
                    status = status_type::stop_at_acceptable_point;
                }
                else {
                    status = status_type::local_infeasibility;
                }
                return true;
            }
            else if (Delta < scalar_type{ 2 } * std::numeric_limits<scalar_type>::epsilon()) {
                status = status_type::error_in_step_computation;
                return true;
            }
            else if (iter_count >= max_iters) {
                status = status_type::maxiter_exceeded;
                return true;
            }
        }

        return false;
    };

    if (!test_stop(result.status,
        normgradinf, tolf, tolx,
        step_accept, result.iter_count, max_iters, Delta, normd,
        obj, objold, d, result.x, result.dg_dx_colmaj, n)) {

        // main iteration loop, set up the method's subroutines
        // that calculate the step
        const auto dogleg = [](auto &step, auto &quadobj, auto &normstep, auto &normstepscal,
            auto &grad_on_entry, auto &stepwork, auto &jacwork,
            auto n, const auto &g, const auto &jac, auto Delta,
            const auto &scales)
        {
            //
            // Approximately solves the trust region subproblem via a dogleg approach,
            // i.e., finding an approximate solution "d" to problem:
            //
            //     min_d      f + g' * d + 0.5 * d' * Bd
            //     subject to || D * d|| <= Delta
            //
            // where "g" is the gradient of "f", "B" is a Hessian approximation of "f",
            // "D" is a diagonal scaling matrix and "Delta" is a given trust region radius.
            //

            auto &gradscal = grad_on_entry;
            auto &gradscal2 = step;
            auto normgradscal = scalar_type{ 0 };
            for (auto i = 0u; i < n; ++i) {
                const auto invscale_i = scalar_type{ 1 } / scales[i];
                gradscal[i] = invscale_i * grad_on_entry[i];
                normgradscal += gradscal[i] * gradscal[i];
                gradscal2[i] = invscale_i * gradscal[i];
            }
            normgradscal = std::sqrt(normgradscal);

            auto &dCauchy = step;
            auto normdCauchy = scalar_type{ 0 };
            auto objCauchy = scalar_type{ 0 };
            if (normgradscal >= std::numeric_limits<scalar_type>::epsilon()) {
                // calculate the Cauchy step (in scaled space)
                auto normsqr_JACvec = scalar_type{ 0 };
                for (auto i = 0u; i < n; ++i) {
                    auto row_JACvec = scalar_type{ 0 };
                    for (auto j = 0u; j < n; ++j) {
                        row_JACvec += jac[i + n * j] * gradscal2[j];
                    }
                    normsqr_JACvec += row_JACvec * row_JACvec;
                }

                const auto stauC = -std::min(scalar_type{ 1 },
                    normgradscal * normgradscal * normgradscal /
                    (Delta * normsqr_JACvec)) * Delta / normgradscal;
                for (auto i = 0u; i < n; ++i) {
                    dCauchy[i] = stauC * gradscal[i];
                }

                // calculate quadratic objective at Cauchy point
                normsqr_JACvec = scalar_type{ 0 };
                auto dot_gradscal_dCauchy = scalar_type{ 0 };
                for (auto i = 0u; i < n; ++i) {
                    auto row_JACvec = scalar_type{ 0 };
                    for (auto j = 0u; j < n; ++j) {
                        row_JACvec += jac[i + n * j] * dCauchy[j] / scales[j];
                    }
                    normsqr_JACvec += row_JACvec * row_JACvec;
                    objCauchy += gradscal[i] * dCauchy[i];
                    normdCauchy += dCauchy[i] * dCauchy[i];
                }

                objCauchy += scalar_type{ 0.5 } * normsqr_JACvec;
                normdCauchy = std::min(std::sqrt(normdCauchy), Delta);

            }
            else {
                std::fill(dCauchy.begin(), dCauchy.end(), scalar_type{ 0 });
            }

            if (Delta - normdCauchy < std::numeric_limits<scalar_type>::epsilon()) {
                // take the Cauchy step if it is at the boundary of the trust region
                quadobj = objCauchy;
            }
            else {
                // calculate the Gauss-Newton step (in scaled space)
                auto &dNewton = stepwork;
                std::copy(g.cbegin(), g.cend(), dNewton.begin());
                std::copy(jac.cbegin(), jac.cend(), jacwork.begin());
                const auto lusolve_failed =
                    lusolve(dNewton.data(), jacwork.data(), static_cast<int>(n), 1) != 0;

                if (lusolve_failed ||
                    !std::all_of(dNewton.cbegin(), dNewton.cend(),
                        [](auto d) { return std::isfinite(d); })) {

                    // take the Cauchy step if the Gauss-Newton step gives bad values
                    quadobj = objCauchy;
                }
                else {
                    for (auto i = 0u; i < n; ++i) {
                        // scale the Gauss-Newton step
                        dNewton[i] *= -scales[i];
                    }

                    const auto normdNewt2 =
                        std::inner_product(dNewton.cbegin(), dNewton.cend(),
                            dNewton.cbegin(), scalar_type{ 0 });

                    const auto Delta2 = Delta * Delta;
                    auto &stored_dCauchy = dNewton;
                    if (normdNewt2 <= Delta2) {
                        // use the Newton direction as the trial step
                        std::swap(step, dNewton);
                    }
                    else {
                        // find the intersect point along dogleg path
                        const auto normdCauchy2 = std::min(normdCauchy * normdCauchy, Delta2);
                        const auto dCdN = std::inner_product(dCauchy.cbegin(), dCauchy.cend(),
                            dNewton.cbegin(), scalar_type{ 0 });
                        const auto dCdNdist2 = normdCauchy2 + normdNewt2 -
                            scalar_type{ 2 } * dCdN;

                        auto tauI = scalar_type{ 0 };
                        if (dCdNdist2 > scalar_type{ 0 }) {
                            // stable method for solving 1-D quadratic
                            const auto b = dCdN - normdCauchy2;
                            if (b != scalar_type{ 0 }) {
                                const auto c = normdCauchy2 - Delta2;
                                const auto disc = std::sqrt(b * b - dCdNdist2 * c);
                                if (b > scalar_type{ 0 }) {
                                    tauI = -c / (b + disc);
                                }
                                else if (b < scalar_type{ 0 }) {
                                    tauI = (disc - b) / dCdNdist2;
                                }

                                // make sure we take a finite step,
                                // (check for poorly scaled/infinite directions).
                                if (!std::isfinite(tauI)) {
                                    // take Cauchy step
                                    tauI = scalar_type{ 0 };
                                }
                            }
                        }

                        const auto onemtauI = scalar_type{ 1 } - tauI;
                        for (auto i = 0u; i < n; ++i) {
                            const auto newton_i = tauI * dNewton[i];
                            dNewton[i] = dCauchy[i];
                            step[i] = newton_i + onemtauI * dCauchy[i];
                        }
                    }

                    // quadratic objective at the trial point
                    quadobj = scalar_type{ 0 };
                    auto normsqr_JACvec = scalar_type{ 0 };
                    for (auto i = 0u; i < n; ++i) {
                        auto row_JACvec = scalar_type{ 0 };
                        for (auto j = 0u; j < n; ++j) {
                            row_JACvec += jac[i + n * j] * step[j] / scales[j];
                        }
                        normsqr_JACvec += row_JACvec * row_JACvec;
                        quadobj += gradscal[i] * step[i];
                    }

                    quadobj += scalar_type{ 0.5 } * normsqr_JACvec;

                    // compare Cauchy step and trial step (Newton or intersection)
                    if (objCauchy < quadobj) {
                        std::swap(step, stored_dCauchy);
                        quadobj = objCauchy;
                    }
                }
            }

            // the calculated step was the scaled one, unscale it
            normstepscal = std::sqrt(std::inner_product(
                step.cbegin(), step.cend(), step.cbegin(), scalar_type{ 0 }));
            for (auto i = 0u; i < n; ++i) {
                step[i] /= scales[i];
            }
            normstep = std::sqrt(std::inner_product(
                step.cbegin(), step.cend(), step.cbegin(), scalar_type{ 0 }));
        };

        const auto update_Delta = [](auto &Delta,
            auto ratio, auto normdscal, auto eta1, auto eta2,
            auto alpha1, auto alpha2, auto DeltaMax, auto current_trial_well_defined)
        {
            // update the trust region's radius
            if (!current_trial_well_defined) {
                // shrink the trust region if any element of the function vector
                // at the new step is not defined(i.e.it's inf/NaN/complex). The
                // update for Delta in this case follows that in snls
                Delta = std::min(std::min(normdscal / scalar_type{ 20 },
                    Delta / scalar_type{ 20 }), DeltaMax);
            }
            else if (ratio < eta1) {
                Delta = std::min(alpha2 * normdscal, DeltaMax);
            }
            else if (ratio >= eta2) {
                Delta = std::min(std::max(Delta, alpha1 * normdscal), DeltaMax);
            }
        };

        // create the vector of scales, and two work vectors to
        // hold results of the trial steps
        if (typical_x.size() < n) {
            typical_x.resize(n, scalar_type{ 1 });
        }

        auto &scales = typical_x;
        std::for_each(scales.begin(), scales.end(),
            [](auto &scale_i)
            {
                if (std::abs(scale_i) < std::numeric_limits<scalar_type>::epsilon()) {
                    scale_i = scalar_type{ 1 };
                }
                else {
                    scale_i = scalar_type{ 1 } / std::abs(scale_i);
                }
            });

        auto x_trial = result.x;
        auto g_trial = result.g;
        Dvector jac_trial(n * n);
        scalar_type pred;
        scalar_type normdscal;
        do {
            // calculate the step "d" using the dogleg approach
            dogleg(d, pred, normd, normdscal,
                grad, g_trial, jac_trial,
                n, result.g, result.dg_dx_colmaj, Delta, scales);

            // calculate the trial point
            for (auto j = 0u; j < n; ++j) {
                x_trial[j] = result.x[j] + d[j];
            }

            // evaluate the equations and the objective at the trial point
            update_g_and_dense_jac(g_trial, jac_trial, x_trial, n);
            ++result.iter_count;

            const auto obj_trial = scalar_type{ 0.5 } *
                std::inner_product(g_trial.cbegin(), g_trial.cend(), g_trial.cbegin(),
                    scalar_type{ 0 });

            // calculate the ratio between the actual reduction given
            // by x_trial and pred
            const auto ratio = pred >= scalar_type{ 0 } ?
                scalar_type{ 0 } : (obj_trial - obj) / pred;

            // update the fault tolerance
            const auto current_trial_well_defined = std::isfinite(obj_trial);

            // accept or reject the current step
            step_accept = ratio > eta1 && current_trial_well_defined;
            if (step_accept) {
                std::swap(result.x, x_trial);
                std::swap(result.g, g_trial);
                std::swap(result.dg_dx_colmaj, jac_trial);
                objold = obj;
                obj = obj_trial;
                update_grad_and_normgradinf(grad, normgradinf, result.dg_dx_colmaj, result.g, n);

            }

            // update the trust region's radius
            update_Delta(Delta,
                ratio, normdscal, eta1, eta2,
                alpha1, alpha2, DeltaMax, current_trial_well_defined);

        } while (!test_stop(result.status,
            normgradinf, tolf, tolx,
            step_accept, result.iter_count, max_iters, Delta, normd,
            obj, objold, d, result.x, result.dg_dx_colmaj, n));
    }

    return result;
}

} // end namespace CppAD

#endif
