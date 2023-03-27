#ifndef LION_MATH_BASIC_ROOT_SOLVERS_H
#define LION_MATH_BASIC_ROOT_SOLVERS_H
#pragma once


#include <tuple>

#include "lion/foundation/utils.h"


//
// Contains useful numerical algorithms to solve for
// roots of different types of (simple) equations.
//


enum class basic_root_solvers_exitflag : int
{
    //
    // A value that we can use to check the result of
    // the algorithms in this header.
    //
    // NOTE: BY DEFINITION, ANY EXIT FLAG REPRESENTING
    // AN ERROR WILL BE <= 0.
    //

    success = 1,                      // equation solved
    stop_at_tiny_step = 2,            // change in x smaller than the specified tolerance

    maxiter_exceeded = 0,             // max. number of iterations (here, == function evaluations) exceeded
    local_infeasibility = -2,         // equation not solved, possibly converged to a singular point
    invalid_number_detected = -3,     // NaN or Inf detected
    no_sign_change_found = -6,        // did not encounter a sign change of the objective function, roots cannot be found (fzero only)
    same_sign_at_domain_limits = -200 // objective function has the same sign at both ends of the domain (fzero only, the routine may find a zero if lucky!)
};


template<typename Fun,
         typename T,
         typename TolType = T>
inline std::tuple<T, T,
    basic_root_solvers_exitflag> fzero(Fun fun, T x0,
        std::size_t max_fun_evals = 10000u,
        TolType tolx = TolType{ 1e-12 },
        TolType tolf = TolType{ 1e-12 })
{
    //
    // Finds the zeros of a scalar function "fun", given a seed "x0", tolerances
    // "tolx" and "tolf", and a maximum number of function evaluations.
    // Returns an std::tuple containing "[x, fval, exitflag]",
    // where "fval = fun(x)" which should be approx 0 when successful.
    //
    // The "exitflag" output can take the following values from enum class
    // "basic_root_solvers_exitflag":
    //
    //     * success                 -> function converged to a solution x.
    //     * invalid_number_detected -> nan or inf function value was encountered.
    //     * local_infeasibility     -> the algorithm might have converged to a singular point.
    //     * no_sign_change_found    -> fzero did not detect a sign change.
    //     * maxiter_exceeded        -> limit of max. numer of function evaluations was reached
    //
    // NOTE: this function comes from code generation of Matlab's "fzero".
    //

#define FZERO_CHECK_FUN_EVALS { ++fun_evals; if (fun_evals > max_fun_evals) \
    { exitflag = basic_root_solvers_exitflag::maxiter_exceeded;  goto ret; } }

    // convert tolf and tolx to multiples of eps
    const auto tolf_e = tolf / std::numeric_limits<TolType>::epsilon();
    const auto tolx_e = tolx / std::numeric_limits<TolType>::epsilon();

    std::size_t fun_evals = 0u;
    auto exitflag = basic_root_solvers_exitflag::success;

    auto x = x0;
    auto fval = fun(x); FZERO_CHECK_FUN_EVALS;

    if (!std::isfinite(fval)) {
        exitflag = basic_root_solvers_exitflag::invalid_number_detected;
        x = std::numeric_limits<T>::quiet_NaN();
        fval = std::numeric_limits<T>::quiet_NaN();
        goto ret;
    }
    else if (!std::isfinite(x)) {
        exitflag = basic_root_solvers_exitflag::no_sign_change_found;
        x = std::numeric_limits<T>::quiet_NaN();
        fval = std::numeric_limits<T>::quiet_NaN();
        goto ret;
    }

    if (eq_fp(fval, T{ 0 }, tolf_e)) {
        goto ret;
    }
    else {
        auto dx = !eq_fp(x, T{ 0 }, tolx_e) ? x / T{ 50 } : T{ 0.02 };
        auto fx = fval;
        auto a = x;
        auto guard = false;
        while (true) {
            if ((fx > T{ 0 }) == (fval > T{ 0 })) { // look for sign change
                dx *= T{ 1.4142135623730951 }; // standard search, dx is expanded by sqrt(2)
                a = x0 - dx;
                fx = fun(a); FZERO_CHECK_FUN_EVALS
                if (!std::isfinite(fx)) {
                    exitflag = basic_root_solvers_exitflag::invalid_number_detected;
                    x = std::numeric_limits<T>::quiet_NaN();
                    fval = std::numeric_limits<T>::quiet_NaN();
                    goto ret;
                }
                else if (!std::isfinite(a)) {
                    exitflag = basic_root_solvers_exitflag::no_sign_change_found;
                    x = std::numeric_limits<T>::quiet_NaN();
                    fval = std::numeric_limits<T>::quiet_NaN();
                    goto ret;
                }
                else if ((fx > T{ 0 }) != (fval > T{ 0 })) { // we found a sign change in the last evaluation
                    guard = true;
                    break;
                }
                else {
                    x = x0 + dx;
                    fval = fun(x); FZERO_CHECK_FUN_EVALS
                    if (!std::isfinite(fval)) {
                        exitflag = basic_root_solvers_exitflag::invalid_number_detected;
                        x = std::numeric_limits<T>::quiet_NaN();
                        fval = std::numeric_limits<T>::quiet_NaN();
                        goto ret;
                    }
                    else if (!std::isfinite(x)) {
                        exitflag = basic_root_solvers_exitflag::no_sign_change_found;
                        x = std::numeric_limits<T>::quiet_NaN();
                        fval = std::numeric_limits<T>::quiet_NaN();
                        goto ret;
                    }
                    else {
                        guard = false;
                    }
                }
            }
            else {
                guard = true;
                break;
            }
        }

        if (guard) {
            const auto savefa = fx;
            const auto savefb = fval;
            auto fc = fval;
            auto c = x;
            auto e = T{ 0 };
            auto d = T{ 0 };
            auto exit = false;
            while (!exit && (!eq_fp(fval, T{ 0 }, tolf_e) && !eq_fp(a, x, tolx_e))) {
                if ((fval > T{ 0 }) == (fc > T{ 0 })) {
                    c = a;
                    fc = fx;
                    d = x - a;
                    e = d;
                }

                if (std::abs(fc) < std::abs(fval)) {
                    a = x;
                    x = c;
                    c = a;
                    fx = fval;
                    fval = fc;
                    fc = fx;
                }

                const auto m = T{ 0.5 } * (c - x);
                auto q  = std::max(T{ 1 }, std::abs(x));
                const auto toler{ T{ 2 } * tolx * q };
                if ((std::abs(m) <= toler) || eq_fp(fval, T{ 0 }, tolf_e)) {
                    exit = true;
                }
                else {
                    if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fval))) {
                        d = m;
                        e = m;
                    }
                    else {
                        const auto s = fval / fx;
                        if (eq_fp(a, c, tolx_e)) {
                            dx = T{ 2 } * m * s;
                            q = T{ 1 } - s;
                        }
                        else {
                            q = fx / fc;
                            fx = fval / fc;
                            dx = s * (T{ 2 } * m * q * (q - fx) - (x - a) * (fx - T{ 1 }));
                            q = (q - T{ 1 }) * (fx - T{ 1 }) * (s - T{ 1 });
                        }

                        if (dx > T{ 0 }) {
                            q = -q;
                        }
                        else {
                            dx = -dx;
                        }

                        if ((T{ 2 } * dx < T{ 3 } * m * q - std::abs(toler * q)) &&
                            (dx < std::abs(T{ 0.5 } * e * q))) {
                            e = d;
                            d = dx / q;
                        }
                        else {
                            d = m;
                            e = m;
                        }
                    }

                    a = x;
                    fx = fval;
                    if (std::abs(d) > toler) {
                        x += d;
                    }
                    else if (x > c) {
                        x -= toler;
                    }
                    else {
                        x += toler;
                    }

                   fval = fun(x); FZERO_CHECK_FUN_EVALS
                   if (!std::isfinite(fval)) {
                       exitflag = basic_root_solvers_exitflag::invalid_number_detected;
                       x = std::numeric_limits<T>::quiet_NaN();
                       fval = std::numeric_limits<T>::quiet_NaN();
                       goto ret;
                   }
                   else if (!std::isfinite(x)) {
                       exitflag = basic_root_solvers_exitflag::no_sign_change_found;
                       x = std::numeric_limits<T>::quiet_NaN();
                       fval = std::numeric_limits<T>::quiet_NaN();
                       goto ret;
                   }
                }
            }

            const auto q = std::abs(savefa);
            dx = std::abs(savefb);
            const auto b_q = (q > dx) || std::isnan(dx) ? q : dx;
            if (std::abs(fval) > b_q) {
                exitflag = basic_root_solvers_exitflag::local_infeasibility;
            }
        }
    }

    ret:
    return { x, fval, exitflag };

    // macro cleanup
#undef FZERO_CHECK_FUN_EVALS
}


template<typename Fun,
         typename T,
         typename TolType = T>
inline std::tuple<T, T,
    basic_root_solvers_exitflag> fzero(Fun fun, T a, T b,
        std::size_t max_fun_evals = 10000u,
        TolType tolx = TolType{ 1e-12 },
        TolType tolf = TolType{ 1e-12 })
{
    //
    // Finds the zeros of a scalar function "fun" within the interval "[a, b]",
    // given tolerances "tolx" and "tolf", and a maximum number of function
    // evaluations.
    //
    // WARNING: fun(a) and fun(b) SHOULD DIFFER IN SIGN.
    //
    // The function returns an std::tuple containing "[x, fval, exitflag]",
    // where "fval = fun(x)", which should be approx 0 when successful.
    //
    // The "exitflag" output can take the following values from enum class
    // "basic_root_solvers_exitflag":
    //
    //     * success                    -> function converged to a solution x.
    //     * invalid_number_detected    -> nan or inf function value was encountered.
    //     * local_infeasibility        -> the algorithm might have converged to a singular point.
    //     * no_sign_change_found       -> fzero did not detect a sign change.
    //     * maxiter_exceeded           -> limit of max. numer of function evaluations was reached
    //     * same_sign_at_domain_limits -> sign(fun(a)) == sign(fun(b)), but the routine may find a zero if lucky!
    //
    // NOTE: this function comes from code generation of Matlab's "fzero".
    //

#define FZERO_CHECK_FUN_EVALS { ++fun_evals; if (fun_evals > max_fun_evals) \
    { exitflag = basic_root_solvers_exitflag::maxiter_exceeded;  goto ret; } }


    // convert tolf and tolx to multiples of eps
    const auto tolf_e = tolf / std::numeric_limits<TolType>::epsilon();
    const auto tolx_e = tolx / std::numeric_limits<TolType>::epsilon();

    std::size_t fun_evals = 0u;
    auto exitflag = basic_root_solvers_exitflag::success;

    T x = b;
    T fval = fun(x);
    T b_a = a;
    T fa = fun(b_a);

    if (eq_fp(fval, T{ 0 }, tolf_e)) {
        goto ret;
    }
    else {
        FZERO_CHECK_FUN_EVALS
        FZERO_CHECK_FUN_EVALS
        if (!std::isfinite(fval)) {
            exitflag = basic_root_solvers_exitflag::invalid_number_detected;
            x = std::numeric_limits<T>::quiet_NaN();
            fval = std::numeric_limits<T>::quiet_NaN();
            goto ret;
        }
        else if (!std::isfinite(x)) {
            exitflag = basic_root_solvers_exitflag::no_sign_change_found;
            x = std::numeric_limits<T>::quiet_NaN();
            fval = std::numeric_limits<T>::quiet_NaN();
            goto ret;
        }
    }

    if (eq_fp(fa, T{ 0 }, tolf_e)) {
        x = a;
        fval = fa;
        goto ret;
    }
    else {
        FZERO_CHECK_FUN_EVALS;
        if (!std::isfinite(fa)) {
            exitflag = basic_root_solvers_exitflag::invalid_number_detected;
            x = std::numeric_limits<T>::quiet_NaN();
            fval = std::numeric_limits<T>::quiet_NaN();
            goto ret;
        }
        else if (!std::isfinite(b_a)) {
            exitflag = basic_root_solvers_exitflag::no_sign_change_found;
            x = std::numeric_limits<T>::quiet_NaN();
            fval = std::numeric_limits<T>::quiet_NaN();
            goto ret;
        }
    }


    // f(a) and f(b) should formally differ in sign, but we'll let the routine continue
    // (there are cases in which it can find a zero nevertheless)
    if (samesign(fval, fa)) {
        exitflag = basic_root_solvers_exitflag::same_sign_at_domain_limits;
    }


    // zero-finding routine
    {
        const auto savefa = fa;
        const auto savefb = fval;
        auto fc = fval;
        auto c = b;
        auto e = T{ 0 };
        auto d = T{ 0 };

        auto exit{ false };
        while (!exit && (!eq_fp(fval, T{ 0 }, tolf_e) && !eq_fp(b_a, x, tolx_e))) {
            if ((fval > T{ 0 }) == (fc > T{ 0 })) {
                c = b_a;
                fc = fa;
                d = x - b_a;
                e = d;
            }

            if (std::abs(fc) < std::abs(fval)) {
                b_a = x;
                x = c;
                c = b_a;
                fa = fval;
                fval = fc;
                fc = fa;
            }

            const auto m = T{ 0.5 } * (c - x);
            auto p = std::max(T{ 1 }, std::abs(x));
            const auto toler = T{ 2 } * tolx * p;
            if ((std::abs(m) <= toler) || eq_fp(fval, T{ 0 }, tolf_e)) {
                exit = true;
            }
            else {
                if ((std::abs(e) < toler) || (std::abs(fa) <= std::abs(fval))) {
                    d = m;
                    e = m;
                }
                else {
                    const auto s = fval / fa;
                    if (eq_fp(b_a, c, tolx_e)) {
                        p = T{ 2 } * m * s;
                        fa = T{ 1 } - s;
                    }
                    else {
                        fa /= fc;
                        const auto r = fval / fc;
                        p = s * (T{ 2 } * m * fa * (fa - r) - (x - b_a) * (r - T{ 1 }));
                        fa = (fa - T{ 1 }) * (r - T{ 1 }) * (s - T{ 1 });
                    }

                    if (p > T{ 0 }) {
                        fa = -fa;
                    }
                    else {
                        p = -p;
                    }

                    if ((T{ 2 } * p < T{ 3 } * m * fa - std::abs(toler * fa)) &&
                        (p < std::abs(T{ 0.5 } * e * fa))) {
                        e = d;
                        d = p / fa;
                    }
                    else {
                        d = m;
                        e = m;
                    }
                }

                b_a = x;
                fa = fval;
                if (std::abs(d) > toler) {
                    x += d;
                }
                else if (x > c) {
                    x -= toler;
                }
                else {
                    x += toler;
                }

                fval = fun(x);
                if (!std::isfinite(fval)) {
                    exitflag = basic_root_solvers_exitflag::invalid_number_detected;
                    x = std::numeric_limits<T>::quiet_NaN();
                    fval = std::numeric_limits<T>::quiet_NaN();
                    goto ret;
                }
                else if (!std::isfinite(x)) {
                    exitflag = basic_root_solvers_exitflag::no_sign_change_found;
                    x = std::numeric_limits<T>::quiet_NaN();
                    fval = std::numeric_limits<T>::quiet_NaN();
                    goto ret;
                }
            }
        }

        const auto p = std::abs(savefa);
        fa = std::abs(savefb);
        const auto b_p = ((p > fa) || std::isnan(fa)) ? p : fa;
        if (std::abs(fval) > b_p) {
            exitflag = basic_root_solvers_exitflag::local_infeasibility;
        }
    }

    ret:
    return { x, fval, exitflag };

    // macro cleanup
#undef FZERO_CHECK_FUN_EVALS
}


template<typename Fun, typename dFun_dx,
         typename T, typename TolType = T>
constexpr std::tuple<T, T,
    basic_root_solvers_exitflag> Newton_method1d(Fun fun, dFun_dx dfun_dx, T x0,
        T alpha = T{ 1 },
        std::size_t max_fun_evals = 10000u,
        TolType tolx = TolType{ 1e-12 },
        TolType tolf = TolType{ 1e-12 })
{
    //
    // Applies the damped Newton method to solve "fun(x) = 0", with
    // damping factor "alpha" (multiplies the step) and seed "x0".
    // Returns an std::tuple containing "[x, fval, exitflag]".
    //
    // The "exitflag" output can take the following values from enum class
    // "basic_root_solvers_exitflag":
    //
    //     * success                 -> function converged to a solution x with abs(fun(x)) <= tolf
    //     * stop_at_tiny_step       -> function converged to a solution x with abs(fun(x)) > tolf (Newton step has become < tolx)
    //     * invalid_number_detected -> NaN or Inf function value was encountered
    //     * maxiter_exceeded        -> limit of max. numer of function evaluations was reached
    //

    const auto tolf_e = tolf / std::numeric_limits<TolType>::epsilon();

    auto exitflag = basic_root_solvers_exitflag::success;

    auto fval = fun(x0);
    if (eq_fp(fval, T{ 0 }, tolf_e)) {
        return { x0, fval, exitflag };
    }
    else if (!std::isfinite(x0) || !std::isfinite(fval)) {
        exitflag = basic_root_solvers_exitflag::invalid_number_detected;
        x0 = std::numeric_limits<T>::quiet_NaN();
        fval = std::numeric_limits<T>::quiet_NaN();
        return { x0, fval, exitflag };
    }

    for (std::size_t fun_evals = 2u; fun_evals <= max_fun_evals; ++fun_evals) {
        const auto dx = fval / dfun_dx(x0);
        x0 -= alpha * dx;
        fval = fun(x0);

        if (eq_fp(fval, T{ 0 }, tolf_e)) {
            return { x0, fval, exitflag };
        }
        else if (std::abs(T{ 0.5 } * dx) <= TolType{ 2 } * tolx * std::max(T{ 1 }, std::abs(x0))) {
            exitflag = basic_root_solvers_exitflag::stop_at_tiny_step;
            return { x0, fval, exitflag };
        }
        else if (!std::isfinite(x0) || !std::isfinite(fval)) {
            exitflag = basic_root_solvers_exitflag::invalid_number_detected;
            x0 = std::numeric_limits<T>::quiet_NaN();
            fval = std::numeric_limits<T>::quiet_NaN();
            return { x0, fval, exitflag };
        }
    }

    exitflag = basic_root_solvers_exitflag::maxiter_exceeded;
    return { x0, fval, exitflag };
}


template<typename Fun,
         typename T, typename TolType = T>
constexpr std::tuple<T,
    basic_root_solvers_exitflag> fixed_point_iteration(Fun fun, T x0,
        std::size_t max_fun_evals = 10000u,
        TolType tolx = TolType{ 1e-12 })
{
    //
    // Applies fixed-point iteration to solve equation "x = fun(x)". Returns an
    // std::tuple containing "[x, exitflag]".
    //
    // The "exitflag" output can take the following values from enum class
    // "basic_root_solvers_exitflag":
    //
    //     * success                 -> Function converged to a solution x.
    //     * invalid_number_detected -> NaN or Inf function value was encountered.
    //     * maxiter_exceeded        -> limit of max. numer of function evaluations was reached
    //

    const auto tolx_e = tolx / std::numeric_limits<TolType>::epsilon();
    auto exitflag = basic_root_solvers_exitflag::success;

    if (!std::isfinite(x0)) {
        exitflag = basic_root_solvers_exitflag::invalid_number_detected;
        x0 = std::numeric_limits<T>::quiet_NaN();
        return { x0, exitflag };
    }

    for (std::size_t fun_evals = 1u; fun_evals <= max_fun_evals; ++fun_evals) {
        const auto x_prev = x0;
        x0 = fun(x_prev);

        if (eq_fp(x0 - x_prev, T{ 0 }, tolx_e)) {
            return { x0, exitflag };

        }
        else if (!std::isfinite(x0)) {
            exitflag = basic_root_solvers_exitflag::invalid_number_detected;
            x0 = std::numeric_limits<T>::quiet_NaN();
            return { x0, exitflag };

        }
    }

    exitflag = basic_root_solvers_exitflag::maxiter_exceeded;
    return { x0, exitflag };
}


template<typename Fun2d,
         typename T, typename TolType = T>
constexpr std::tuple<std::array<T, 2>,
    std::array<T, 2>,
    basic_root_solvers_exitflag> numjac_Newton_method2d(Fun2d fun2d,
        std::array<T, 2> x0,
        std::array<T, 2> deltas_numjac = std::array<T, 2>{ T{ 1e-9 }, T{ 1e-9 } },
        std::array<T, 2> alphas = std::array<T, 2>{ T{ 1 }, T{ 1 } },
        std::size_t max_fun_evals = 10000u,
        TolType tolx = TolType{ 1e-12 },
        TolType tolf = TolType{ 1e-12 })
{
    //
    // Applies the two.dimensional damped Newton method to solve
    // "fun(x) = [0; 0]" ("fun" takes an std::array<T, 2> and returns
    // an std::array<T, 2>), with damping factor "alpha" (multiplies
    // the step) and seed "x0". Returns an std::tuple containing
    // "[x, fval, exitflag]". The jacobian of "fun" is calculated
    // numerically, taking central finite differences of step sizes
    // given by "deltas_numjac".
    //
    // The "exitflag" output can take the following values from enum class
    // "basic_root_solvers_exitflag":
    //
    //     * success                 -> function converged to a solution x with abs(fun(x)) <= tolf
    //     * stop_at_tiny_step       -> function converged to a solution x with abs(fun(x)) > tolf (Newton step has become < tolx)
    //     * invalid_number_detected -> NaN or Inf function value was encountered
    //     * maxiter_exceeded        -> limit of max. numer of function evaluations was reached
    //

    const auto tolf_e = tolf / std::numeric_limits<TolType>::epsilon();

    auto exitflag = basic_root_solvers_exitflag::success;

    auto fval = fun2d(x0);
    if (eq_fp(std::abs(fval[0]) + std::abs(fval[1]), T{ 0 }, tolf_e)) {
        return { x0, fval, exitflag };
    }
    else if (!std::isfinite(x0[0]) || !std::isfinite(fval[0]) ||
             !std::isfinite(x0[1]) || !std::isfinite(fval[1])) {

        exitflag = basic_root_solvers_exitflag::invalid_number_detected;
        std::fill(x0.begin(), x0.end(), std::numeric_limits<T>::quiet_NaN());
        std::fill(fval.begin(), fval.end(), std::numeric_limits<T>::quiet_NaN());
        return { x0, fval, exitflag };
    }


    const auto fun_twovars = [fun2d](auto x, auto y)
    {
        return fun2d(std::array<T, 2u>{ x, y });
    };

    for (std::size_t fun_evals = 2u; fun_evals <= max_fun_evals; ++fun_evals) {
        const auto Fp0 = fun_twovars(x0[0] + deltas_numjac[0], x0[1]);
        const auto Fm0 = fun_twovars(x0[0] - deltas_numjac[0], x0[1]);
        const auto Fp1 = fun_twovars(x0[0], x0[1] + deltas_numjac[1]);
        const auto Fm1 = fun_twovars(x0[0], x0[1] - deltas_numjac[1]);

        const auto numjac_xx = T{ 0.5 } * (Fp0[0] - Fm0[0]) / deltas_numjac[0];
        const auto numjac_yx = T{ 0.5 } * (Fp0[1] - Fm0[1]) / deltas_numjac[0];
        const auto numjac_xy = T{ 0.5 } * (Fp1[0] - Fm1[0]) / deltas_numjac[1];
        const auto numjac_yy = T{ 0.5 } * (Fp1[1] - Fm1[1]) / deltas_numjac[1];

        const auto inv_det_numjac = T{ 1 } / (numjac_xx * numjac_yy - numjac_yx * numjac_xy);

        const auto dx0 = inv_det_numjac * (fval[0] * numjac_yy - fval[1] * numjac_xy);
        const auto dx1 = inv_det_numjac * (fval[1] * numjac_xx - fval[0] * numjac_yx);

        x0[0] -= alphas[0] * dx0;
        x0[1] -= alphas[1] * dx1;

        fval = fun2d(x0);

        if (eq_fp(std::abs(fval[0]) + std::abs(fval[1]), T{ 0 }, tolf_e)) {
            return { x0, fval, exitflag };
        }
        else if (T{ 0.5 } * (std::abs(dx0) + std::abs(dx1)) <=
                 TolType{ 2 } * tolx * std::max(T{ 1 }, std::abs(x0[0]) + std::abs(x0[1]))) {

            exitflag = basic_root_solvers_exitflag::stop_at_tiny_step;
            return { x0, fval, exitflag };
        }
        else if (!std::isfinite(x0[0]) || !std::isfinite(fval[0]) ||
                 !std::isfinite(x0[1]) || !std::isfinite(fval[1])) {

            exitflag = basic_root_solvers_exitflag::invalid_number_detected;
            std::fill(x0.begin(), x0.end(), std::numeric_limits<T>::quiet_NaN());
            std::fill(fval.begin(), fval.end(), std::numeric_limits<T>::quiet_NaN());
            return { x0, fval, exitflag };
        }
    }

    exitflag = basic_root_solvers_exitflag::maxiter_exceeded;
    return { x0, fval, exitflag };
}


template<typename T>
constexpr T real_root_of_cubic(T a, T b, T c, T d)
{
    //
    // Extracts one real root from the cubic equation of real coefficients
    // "a * x ^ 3 + b * x ^ 2 + c * x + d = 0". See E.A. Pritchard, 1995.
    //
    // WARNING: all algorithms based on classic formulas, though fast, have
    // modest stability & precision properties. The most precise way of
    // calculating polynomial roots is using the companion matrix method (see
    // tr7/linear_algebra.h")
    //

    b = b / a;
    c = c / a;
    d = d / a;

    auto K3sq = b * b;
    auto r1 = K3sq / T{ 3 } - c;

    K3sq = K3sq * b / T{ 27 } - d;

    const auto r2 = -b / T{ 3 };
    auto K3 = K3sq + r1 * r2;

    K3sq = K3 * K3;

    auto K4 = r1 * r1 * r1 / T{ 27 };
    K4 = K3sq - 4. * K4;

    if (K4 >= T{ 0 }) {
        K4 = std::sqrt(K4);
        K3sq = (K3 + K4) * T{ 0.5 };
        r1 = std::copysign(std::cbrt(std::abs(K3sq)), K3sq);
        K3sq = (K3 - K4) * T{ 0.5 };
        K3 = std::copysign(std::cbrt(std::abs(K3sq)), K3sq);

        return r1 + K3 + r2;
    }
    else if(K3 > std::numeric_limits<T>::epsilon()) {
        K4 = std::sqrt(-K4);

        return T{ 2 } * std::sqrt(r1 / T{ 3 }) * std::cos((std::atan(K4 / K3)) / T{ 3 }) + r2;
    }
    else if(K3 < -std::numeric_limits<T>::epsilon()) {
        K4 = std::sqrt(-K4);

        return T{ 2 } * std::sqrt(r1 / T{ 3 }) * std::cos((pi_T<T> + (std::atan(K4 / K3))) / T{ 3 }) + r2;
    }
    else {
        return std::sqrt(r1) + r2;
    }
}


template<typename RealType,
         typename ComplexType = std::complex<RealType> >
constexpr std::tuple<std::array<ComplexType, 3>,
    std::size_t> cubic_roots(RealType a, RealType b, RealType c, RealType d)
{
    //
    // Calculates the roots of the cubic equation of real coefficients
    // "a * x^3 + b * x^2 + c * x + d = 0". Returns an std::tuple holding
    // "[roots, num_real_roots]". If there is only one real root
    // ("num_real_roots = 1"), it will be placed at "roots[0]", and "roots[1, 2]"
    // will only take valid (i.e. complex) values when template parameter
    // "ComplexType" is a specialization of std::complex.
    //
    // WARNING: all algorithms based on classic formulas, though fast, have
    // modest stability & precision properties. The most precise way of
    // calculating polynomial roots is using the companion matrix method (see
    // tr7/linear_algebra.h")
    //

    b /= a;
    c /= a;
    d /= a;

    // we start by applying Halley's method to eliminate the real root
    // As seed we take the value of the cubic at the inflection point
    auto Z = -b / RealType{ 3 };
    if (((Z + b) * Z + c) * Z + d < RealType{ 0 }) {
        Z = b + RealType{ 2 };
    }
    else {
        Z = RealType{ 1 } + b;
    }

    constexpr auto max_iters_Halley = std::size_t{ 50u };
    constexpr auto tolx_Halley = std::numeric_limits<RealType>::epsilon();
    constexpr auto tolf_Halley = std::numeric_limits<RealType>::epsilon();
    auto Halley_converged = false;
    for (auto i = 0u; i < max_iters_Halley; ++i) {
        const auto r2 = RealType{ 3 } * Z + b;
        const auto r1 = RealType{ 1 } / (Z * (r2 + b) + c);
        const auto K3 = ((Z + b) * Z + c) * Z + d;

        auto dZ = K3 * r1;
        dZ /= (RealType{ 1 } - dZ * r2 * r1);

        Z -= dZ;

        if (std::abs(dZ) < tolx_Halley) {
            if (std::abs(((Z + b) * Z + c) * Z + d) < tolf_Halley) {
                Halley_converged = true;
                break;
            }
        }
    }

    if (!Halley_converged) {
        // if Halley fails, we extract the real root directly
        Z = real_root_of_cubic(RealType{ 1 }, b, c, d);
    }


    // divide cubic by the real root and look at the discriminant:
    // if it is positive, we'l have 3 real roots in total
    const auto B = Z + b;
    const auto C = Z * B + c;
    const auto D = B * B - RealType{ 4 } * C;
    const auto sqrtD = std::sqrt(ComplexType{ D });

    const auto num_real_roots = D >= RealType{ 0 } ? std::size_t{ 3u } : std::size_t{ 1u };

    // we always place the real root at roots[0], the other two will only hold valid values
    // if T is a specialization of std::complex.
    return { std::array<ComplexType, 3>{ Z,
                RealType{ 0.5 } * (-B + sqrtD),
                RealType{ 0.5 } * (-B - sqrtD) },
             num_real_roots };
}


template<typename RealType,
         typename ComplexType = std::complex<RealType> >
constexpr std::tuple<std::array<ComplexType, 4>,
    std::size_t> quartic_roots(RealType a, RealType b, RealType c, RealType d, RealType e)
{
    //
    // Calculates the roots of the quartic equation of real coefficients
    // "a * x^4 + b * x^3 + c * x^2 + d * x + e = 0" using Ferrari's method.
    // Returns an std::tuple holding "[roots, num_real_roots]". If there are
    // only 2 real roots ("num_real_roots = 2"), they will be placed at
    // "roots[0, 1]", and "roots[2, 3]" will only hold valid (complex) values
    // when template parameter "ComplexType" is a specialization of std::complex.
    //
    // WARNING: all algorithms based on classic formulas, though fast, have
    // modest stability & precision properties. The most precise way of
    // calculating polynomial roots is using the companion matrix method (see
    // tr7/linear_algebra.h")
    //

    b /= a;
    c /= a;
    d /= a;
    e /= a;

    // calculate the depressed cubic
    const auto y1 = real_root_of_cubic(RealType{ 1 },
                                       -c,
                                       b * d - RealType{ 4 } * e,
                                       RealType{ 4 } * c * e - b * b * e - d * d);


    // p
    const auto disc = b * b - RealType{ 4 } * c + RealType{ 4 } * y1;
    const auto Dp = std::sqrt(std::max(RealType{ 0 }, disc));
    auto p_a = b + Dp;
    auto p_b = b - Dp;

    RealType q_a, q_b;
    if (std::abs(Dp) > RealType{ 1e-3 }) {
        // linear system, well conditioned
        const auto k = (c - RealType{ 0.25 } * p_a * p_b);
        q_a = (k * p_a - RealType{ 2 } * d) / Dp;
        q_b = RealType{ 2 } * c - RealType{ 0.5 } * p_a * p_b - q_a;
    }
    else {
        // corner case, use the Ferrari's formula
        const auto Dq = std::sqrt(std::max(RealType{ 0 }, y1 * y1 - RealType{ 4 } * e));
        q_a = y1 + Dq;
        q_b = y1 - Dq;
        p_a = RealType{ 0 };
        p_b = RealType{ 0 };
    }

    const auto Da = p_a * p_a - RealType{ 8 } * q_a;
    const auto Db = p_b * p_b - RealType{ 8 } * q_b;
    const auto sqrtDa = std::sqrt(ComplexType{ Da });
    const auto sqrtDb = std::sqrt(ComplexType{ Db });

    if (Da >= RealType{ 0 } && Db >= RealType{ 0 }) {
        return { std::array<ComplexType, 4>{ (-p_a + sqrtDa) * RealType{ 0.25 },
                    (-p_a - sqrtDa) * RealType{ 0.25 },
                    (-p_b + sqrtDb) * RealType{ 0.25 },
                    (-p_b - sqrtDb) * RealType{ 0.25 } },
                 std::size_t{ 4 } };
    }
    else if (Da >= RealType{ 0 }) {
        return { std::array<ComplexType, 4>{ (-p_a + sqrtDa) * RealType{ 0.25 },
                    (-p_a - sqrtDa) * RealType{ 0.25 },
                    (-p_b + sqrtDb) * RealType{ 0.25 },
                    (-p_b - sqrtDb) * RealType{ 0.25 } },
                 std::size_t{ 2 } };
    }
    else if (Db >= RealType{ 0 }) {
        return { std::array<ComplexType, 4>{ (-p_b + sqrtDb) * RealType{ 0.25 },
                    (-p_b - sqrtDb) * RealType{ 0.25 },
                    (-p_a + sqrtDa) * RealType{ 0.25 },
                    (-p_a - sqrtDa) * RealType{ 0.25 } },
                 std::size_t{ 2 } };
    }
    else {
        return { std::array<ComplexType, 4>{ (-p_a + sqrtDa) * RealType{ 0.25 },
                    (-p_a - sqrtDa) * RealType{ 0.25 },
                    (-p_b + sqrtDb) * RealType{ 0.25 },
                    (-p_b - sqrtDb) * RealType{ 0.25 } },
                 std::size_t{ 0 } };
    }
}

#endif
