#ifndef LION_MATH_LINEAR_ALGEBRA_H
#define LION_MATH_LINEAR_ALGEBRA_H


#include <cmath>
#include <cstring>
#include <limits>


//
// Defines functions that perform linear algebra
// operations with dense matrices. They're written
// in classic C-style, since they come from generating
// Matlab code.
//


template<typename T>
inline void matmul(T *X, const T *A, const T *B, int rows_A, int cols_A_rows_B, int cols_B)
{
    //
    // Fills the result of multiplying two dense matrices in column-major
    // order, i.e., "X = A * B". Solution "X" must be preallocated on
    // entry, to size "rows_A * cols_B", and is also a column-major matrix.
    // The number of columns of "A" must coincide with the number of rows
    // of "B", equal to input "cols_A_rows_B".
    //

    if (cols_A_rows_B != 1) {
        const auto c{ rows_A * (cols_B - 1) };
        auto br{ 0 };
        auto cr{ 0 };

        while (cr <= c) {
            const auto i0{ cr + rows_A };
            for (auto ic = cr; ic < i0; ic++) {
                X[ic] = T{ 0 };
            }
            cr += rows_A;
        }

        cr = 0;
        while (cr <= c) {
            auto ar{ -1 };
            const auto i0{ br + cols_A_rows_B };
            for (auto ib = br; ib < i0; ++ib) {
                if (B[ib] != T{ 0 }) {
                    auto ia{ ar };
                    const auto i1{ cr + rows_A };
                    for (auto ic = cr; ic < i1; ++ic) {
                        ++ia;
                        X[ic] += B[ib] * A[ia];
                    }
                }
                ar += rows_A;
            }
            br += cols_A_rows_B;
            cr += rows_A;
        }
    }
    else {
        for (auto i0 = 0; i0 < rows_A; ++i0) {
            for (auto i1 = 0; i1 < cols_B; ++i1) {
                X[i0 + rows_A * i1] = T{ 0 };
                for (auto cr = 0; cr < cols_A_rows_B; ++cr) {
                    X[i0 + rows_A * i1] +=
                        A[i0 + rows_A * cr] * B[cr + cols_A_rows_B * i1];
                }
            }
        }
    }
}


template<typename T>
inline int lusolve(T *B, T *A, int rows_A, int cols_B)
{
    //
    // Direct solve of the dense square problem "A * X = B", via LU with
    // partial pivoting (i.e., applying a decomposition "A = P * L * U").
    // All matrices are in column-major order. On entry, "B" contains the
    // problem's RHS (which may have multiple colums), and holds solution
    // "X" on return. Matrix "A" must hold the problem's LHS matrix on
    // entry, and exits as the LU factors ("A = P * L * U", where P is a
    // permutation matrix that could be reconstructed using the ipiv local
    // variable, since the i-th row of matrix "A" is interchanged with row
    // "ipiv[i]"). The function returns 0 upon successful exit.
    //

    int ret{ 0 };

    const auto n{ rows_A };
    auto *const ipiv{ new int[n] };
    for (auto i = 0; i < n; ++i) {
        ipiv[i] = i + 1;
    }

    auto jA{ 0 };
    const auto nm1{ n - 1 };
    for (auto j = 0; j < nm1; j++) {
        const auto mmj{ n - j };
        const auto c{ j * (n + 1) };

        if (mmj < 1) {
            jA = -1;
        }
        else {
            jA = 0;
            if (mmj > 1) {
                auto iX{ c };
                auto smax{ std::abs(A[c]) };
                for (auto k = 1; k < mmj; ++k) {
                    ++iX;
                    auto s{ std::abs(A[iX]) };
                    if (s > smax) {
                        jA = k;
                        smax = s;
                    }
                }
            }
        }

        if (A[c + jA] != T{ 0 }) {
            if (jA != 0) {
                ipiv[j] = j + jA + 1;
                auto iX{ j };
                auto iy{ j + jA };
                for (auto k = 1; k <= n; ++k) {
                    auto smax{ A[iX] };
                    A[iX] = A[iy];
                    A[iy] = smax;
                    iX += n;
                    iy += n;
                }
            }
            const auto i0{ c + mmj };
            for (jA = c + 1; jA < i0; ++jA) {
                A[jA] /= A[c];
            }
        }
        else {
            ret = j + 1;
        }

        const auto kAcol{ n - j - 1 };
        jA = c + n;
        auto iy{ c + n };
        for (auto nb = 1; nb <= kAcol; ++nb) {
            const auto smax{ A[iy] };
            if (A[iy] != T{ 0 }) {
                auto iX{ c + 1 };
                const auto i0{ mmj + jA };
                for (auto k = 1 + jA; k < i0; ++k) {
                    A[k] += A[iX] * -smax;
                    ++iX;
                }
            }

            iy += n;
            jA += n;
        }

    }

    for (std::ptrdiff_t iy = 0; iy < nm1; ++iy) { // changed from "auto" to "std::ptrdiff_t"
                                                  // to avoid a useless gcc warning
        if (ipiv[iy] != iy + 1) {
            const auto kAcol{ ipiv[iy] - 1 };
            for (jA = 0; jA < cols_B; ++jA) {
                const auto smax{ B[iy + n * jA] };
                B[iy + n * jA] = B[kAcol + n * jA];
                B[kAcol + n * jA] = smax;
            }
        }
    }

    for (auto j = 1; j <= cols_B; ++j) {
        const auto iy{ n * (j - 1) };
        for (auto k = 0; k < n; ++k) {
            const auto kAcol{ n * k };
            if (B[k + iy] != T{ 0 }) {
                for (jA = k + 1; jA < n; ++jA) {
                    B[jA + iy] -= B[k + iy] * A[jA + kAcol];
                }
            }
        }
    }

    for (auto j = 1; j <= cols_B; ++j) {
        const auto iy{ n * (j - 1) };
        for (auto k = n - 1; k + 1 > 0; --k) {
            const auto kAcol{ n * k };
            if (B[k + iy] != T{ 0 }) {
                B[k + iy] /= A[k + kAcol];
                for (jA = 0; jA < k; ++jA) {
                    B[jA + iy] -= B[k + iy] * A[jA + kAcol];
                }
            }
        }
    }

    if (ret == 0 && !(A[n * n - 1] != T{ 0 })) {
        ret = n;
    }

    delete[] ipiv;

    return ret;
}


template<typename T>
inline T rcond(const T *A, int n)
{
    //
    // Returns an estimate for the reciprocal condition of input
    // square matrix "A" in 1-norm. The matrix has size "n x n" and
    // should be specified in column-major order. If "A" is well
    // conditioned, "rcond(A)" is near 1.0. If "A" is badly
    // conditioned, "rcond(A)" is near 0.
    //

    // define some helper functions
    const auto xgetrf = [](auto *A_, auto n)
    {
        for (auto j = 0; j < n - 1; ++j) {
            const auto mmj = n - j;
            const auto b_tmp = j * (n + 1);
            decltype(n) iy;
            if (mmj < 1) {
                iy = -1;
            }
            else {
                iy = 0;
                if (mmj > 1) {
                    auto ix = b_tmp;
                    auto smax = std::abs(A_[b_tmp]);
                    for (auto jA = 2; jA <= mmj; ++jA) {
                        ++ix;
                        const auto s = std::abs(A_[ix]);
                        if (s > smax) {
                            iy = jA - 1;
                            smax = s;
                        }
                    }
                }
            }

            if (A_[b_tmp + iy] != T{ 0 }) {
                if (iy != 0) {
                    iy += j;
                    auto ix = j;
                    for (auto jA = 0; jA < n; ++jA) {
                        const auto smax = A_[ix];
                        A_[ix] = A_[iy];
                        A_[iy] = smax;
                        ix += n;
                        iy += n;
                    }
                }

                for (iy = b_tmp + 2; iy <= (b_tmp + mmj); ++iy) {
                    A_[iy - 1] /= A_[b_tmp];
                }
            }

            auto ixinc = b_tmp + n;
            auto jA = (b_tmp + n) + 1;
            for (auto jp1j = 0; jp1j <= mmj - 2; ++jp1j) {
                const auto smax = A_[ixinc];
                if (A_[ixinc] != T{ 0 }) {
                    auto ix = b_tmp;
                    for (auto ijA = jA; ijA < (mmj + jA) - 1; ++ijA) {
                        A_[ijA] += A_[ix + 1] * -smax;
                        ++ix;
                    }
                }

                ixinc += n;
                jA += n;
            }
        }
    };

    const auto b_norm = [](const auto *x, auto n)
    {
        auto y = T{ 0 };
        for (auto k = 0; k < n; ++k) {
            y += std::abs(x[k]);
        }
        return y;
    };


    // "rcond1"
    auto result = T{ 0 };
    if (n == 0) {
        result = std::numeric_limits<T>::infinity();
    }
    else if (n == 1) {
        const auto s = std::abs(A[0]);
        if (std::isnan(s)) {
            result = std::numeric_limits<T>::quiet_NaN();
        }
        else if (!std::isinf(s) && s != T{ 0 }) {
            result = T{ 1 };
        }
    }
    else {
        auto normA = T{ 0 };
        auto j{ 0 };
        while (j <= n - 1) {
            auto s = T{ 0 };
            for (auto i = 0; i < n; ++i) {
                s += std::abs(A[i + n * j]);
            }

            if (std::isnan(s)) {
                normA = std::numeric_limits<T>::quiet_NaN();
                break;
            }
            else {
                if (s > normA) {
                    normA = s;
                }

                ++j;
            }
        }

        if (!(normA == T{ 0 })) {
            // copy the input matrix
            auto A_ = new T[n * n];
            memcpy(A_, A, n * n * sizeof(T));

            xgetrf(A_, n);
            auto ix = n;
            while (true) {
                if (ix > 0) {
                    if (A_[(ix + n * (ix - 1)) - 1] == T{ 0 }) {
                        break;
                    }
                    else {
                        ix--;
                    }
                }
                else {
                    auto ainvnm = T{ 0 };
                    auto iter = 2;
                    auto kase = 1;
                    auto jump = 1;
                    auto *x = new T[n];
                    for (auto i1 = 0; i1 < n; ++i1) {
                        x[i1] = T{ 1 } / static_cast<T>(n);
                    }

                    while (true) {
                        if (kase != 0) {
                            if (kase == 1) {
                                for (auto b_j = 0; b_j < n; ++b_j) {
                                    const auto jjA = b_j + b_j * n;
                                    for (auto i = 0; i <= n - b_j - 2; ++i) {
                                        ix = (b_j + i) + 1;
                                        x[ix] -= x[b_j] * A_[(jjA + i) + 1];
                                    }
                                }
                                for (auto b_j = n; b_j >= 1; b_j--) {
                                    const auto jjA = (b_j + (b_j - 1) * n) - 1;
                                    x[b_j - 1] /= A_[jjA];
                                    for (auto i = 0; i <= b_j - 2; ++i) {
                                        ix = (b_j - i) - 2;
                                        x[ix] -= x[b_j - 1] * A_[(jjA - i) - 1];
                                    }
                                }
                            }
                            else {
                                for (auto b_j = 0; b_j < n; ++b_j) {
                                    const auto jA = b_j * n;
                                    auto s = x[b_j];
                                    for (auto i = 0; i < b_j; ++i) {
                                        s -= A_[jA + i] * x[i];
                                    }

                                    x[b_j] = s / A_[jA + b_j];
                                }

                                for (auto b_j = n; b_j >= 1; b_j--) {
                                    const auto jA = (b_j - 1) * n;
                                    auto s = x[b_j - 1];
                                    for (auto i = n; i >= b_j + 1; i--) {
                                        s -= A_[(jA + i) - 1] * x[i - 1];
                                    }

                                    x[b_j - 1] = s;
                                }
                            }

                            T absrexk;
                            if (jump == 1) {
                                ainvnm = b_norm(x, n);
                                if (std::isinf(ainvnm) || std::isnan(ainvnm)) {
                                    result = ainvnm;
                                    break;
                                }
                                else {
                                    for (ix = 0; ix < n; ++ix) {
                                        const auto s = std::abs(x[ix]);
                                        if (s > std::numeric_limits<T>::min()) {
                                            x[ix] /= s;
                                        }
                                        else {
                                            x[ix] = T{ 1 };
                                        }
                                    }

                                    kase = 2;
                                    jump = 2;
                                }
                            }
                            else if (jump == 2) {
                                j = 0;
                                auto s = std::abs(x[0]);
                                for (ix = 0; ix <= n - 2; ++ix) {
                                    absrexk = std::abs(x[ix + 1]);
                                    if (!(absrexk <= s)) {
                                        j = 1 + ix;
                                        s = absrexk;
                                    }
                                }

                                iter = 2;
                                for (auto i1 = 0; i1 < n; ++i1) {
                                    x[i1] = T{ 0 };
                                }

                                x[j] = T{ 1 };
                                kase = 1;
                                jump = 3;
                            }
                            else if (jump == 3) {
                                ainvnm = b_norm(x, n);
                                if (ainvnm <= x[0]) {
                                    auto s = T{ 1 };
                                    for (ix = 0; ix < n; ++ix) {
                                        x[ix] = s * (T{ 1 } +((T{ 1 } +static_cast<T>(ix)) - T{ 1 }) / (static_cast<T>(n) - T{ 1 }));
                                        s = -s;
                                    }

                                    kase = 1;
                                    jump = 5;
                                }
                                else {
                                    for (ix = 0; ix < n; ++ix) {
                                        const auto s = std::abs(x[ix]);
                                        if (s > std::numeric_limits<T>::min()) {
                                            x[ix] /= s;
                                        }
                                        else {
                                            x[ix] = T{ 1 };
                                        }
                                    }

                                    kase = 2;
                                    jump = 4;
                                }
                            }
                            else if (jump == 4) {
                                const auto jA = j;
                                j = 0;
                                auto s = std::abs(x[0]);
                                for (ix = 0; ix <= n - 2; ++ix) {
                                    absrexk = std::abs(x[ix + 1]);
                                    if (!(absrexk <= s)) {
                                        j = 1 + ix;
                                        s = absrexk;
                                    }
                                }

                                if ((std::abs(x[jA]) != std::abs(x[j])) && (iter <= 5)) {
                                    ++iter;
                                    for (auto i1 = 0; i1 < n; ++i1) {
                                        x[i1] = T{ 0 };
                                    }

                                    x[j] = T{ 1 };
                                    kase = 1;
                                    jump = 3;
                                }
                                else {
                                    auto s = T{ 1 };
                                    for (ix = 0; ix < n; ++ix) {
                                        x[ix] = s * (T{ 1 } +((T{ 1 } +static_cast<T>(ix)) - T{ 1 }) / (static_cast<T>(n) - T{ 1 }));
                                        s = -s;
                                    }

                                    kase = 1;
                                    jump = 5;
                                }
                            }
                            else {
                                if (jump == 5) {
                                    const auto s = T{ 2 } *b_norm(x, n) / T{ 3 } / static_cast<T>(n);
                                    if (s > ainvnm) {
                                        ainvnm = s;
                                    }

                                    kase = 0;
                                }
                            }
                        }
                        else {
                            if (ainvnm != T{ 0 }) {
                                result = T{ 1 } / ainvnm / normA;
                            }

                            break;
                        }
                    }

                    delete[] x;
                    break;
                }
            }

            delete[] A_;
        }
    }

    return result;
}

#endif