#ifndef LION_MATH_LINEAR_ALGEBRA_H
#define LION_MATH_LINEAR_ALGEBRA_H


#include <cmath>


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

    for (auto iy = 0; iy < nm1; ++iy) {
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

#endif