#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <cmath>

// matmul
template<typename T>
static void matmul(T *X, const T *A, const T *B, int rows_A, int cols_A_rows_B, int cols_B)
{
    // X = A * B

    if (cols_A_rows_B != 1) {

        auto k{ cols_A_rows_B };
        auto m{ rows_A };
        auto c{ rows_A * (cols_B - 1) };
        auto br{ 0 };
        auto cr{ 0 };

        while (cr <= c) {
            auto i0{ cr + m };
            for (auto ic = cr; ic + 1 <= i0; ic++) {
                X[ic] = static_cast<T>(0.0);
            }
            cr += m;
        }

        cr = 0;
        while (cr <= c) {
            auto ar{ -1 };
            auto i0{ br + k };
            for (auto ib = br; ib + 1 <= i0; ++ib) {
                if (B[ib] != static_cast<T>(0.0)) {
                    auto ia{ ar };
                    auto i1{ cr + m };
                    for (auto ic = cr; ic + 1 <= i1; ++ic) {
                        ++ia;
                        X[ic] += B[ib] * A[ia];
                    }
                }
                ar += m;
            }
            br += k;
            cr += m;
        }
    }
    else {
        auto br{ rows_A };
        auto ar{ cols_B };
        auto c{ cols_A_rows_B };
        for (auto i0 = 0; i0 < br; ++i0) {
            for (auto i1 = 0; i1 < ar; ++i1) {
                X[i0 + rows_A * i1] = static_cast<T>(0.0);
                for (auto cr = 0; cr < c; ++cr) {
                    X[i0 + rows_A * i1] +=
                        A[i0 + rows_A * cr] * B[cr + cols_A_rows_B * i1];
                }
            }
        }
    }
    return;
}

// lusolve
template<typename T>
inline int lusolve(T *B, T *A, int rows_A, int cols_B)
{
    // A * X = B (square problem), solved through LU with partial pivoting (A = P * L * U)
    // ---> in this function, T *B must enter as a copy of B, and exits as X.
    //                        T *A enters as a copy of A, and exits as the LU factors (A = P * L * U, where P 
    //                        is a permutation matrix that could be reconstructed using the ipiv local variable, 
    //                        since the i-th row of matrix A has been interchanged with row ipiv[i]).
    // ---> function returns 0 upon successful exit

    auto info{ 0 }; // ret value

    auto n{ rows_A };
    auto *ipiv{ new int[n] };
    for (auto i = 0; i < n; ++i) {
        ipiv[i] = i + 1;
    }

    auto u1{ n - 1 };
    auto jA{ 0 };
    for (auto j = 0; j + 1 <= u1; j++) {
        auto mmj{ n - j };
        auto c{ j * (n + 1) };

        if (mmj < 1) {
            jA = -1;
        }
        else {
            jA = 0;
            if (mmj > 1) {
                auto iX{ c };
                auto smax{ std::abs(A[c]) };
                for (auto k = 1; k + 1 <= mmj; ++k) {
                    ++iX;
                    auto s{ std::abs(A[iX]) };
                    if (s > smax) {
                        jA = k;
                        smax = s;
                    }
                }
            }
        }

        if (A[c + jA] != static_cast<T>(0.0)) {
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
            auto i0{ c + mmj };
            for (jA = c + 1; jA + 1 <= i0; ++jA) {
                A[jA] /= A[c];
            }
        }
        else {
            info = j + 1;
        }

        auto kAcol{ n - j - 1 };
        jA = c + n;
        auto iy{ c + n };
        for (auto nb = 1; nb <= kAcol; ++nb) {
            auto smax{ A[iy] };
            if (A[iy] != static_cast<T>(0.0)) {
                auto iX{ c + 1 };
                auto i0{ mmj + jA };
                for (auto k = 1 + jA; k + 1 <= i0; ++k) {
                    A[k] += A[iX] * -smax;
                    ++iX;
                }
            }

            iy += n;
            jA += n;
        }

    }

    auto nb{ cols_B };
    for (auto iy = 0; iy + 1 < n; ++iy) {
        if (ipiv[iy] != iy + 1) {
            auto kAcol{ ipiv[iy] - 1 };
            for (jA = 0; jA + 1 <= nb; ++jA) {
                auto smax{ B[iy + n * jA] };
                B[iy + n * jA] = B[kAcol + n * jA];
                B[kAcol + n * jA] = smax;
            }
        }
    }

    for (auto j = 1; j <= nb; ++j) {
        auto iy{ n * (j - 1) };
        for (auto k = 0; k + 1 <= n; ++k) {
            auto kAcol{ n * k };
            if (B[k + iy] != static_cast<T>(0.0)) {
                for (jA = k + 1; jA + 1 <= n; ++jA) {
                    B[jA + iy] -= B[k + iy] * A[jA + kAcol];
                }
            }
        }
    }

    for (auto j = 1; j <= nb; ++j) {
        auto iy{ n * (j - 1) };
        for (auto k = n - 1; k + 1 > 0; --k) {
            auto kAcol{ n * k };
            if (B[k + iy] != static_cast<T>(0.0)) {
                B[k + iy] /= A[k + kAcol];
                for (jA = 0; jA + 1 <= k; ++jA) {
                    B[jA + iy] -= B[k + iy] * A[jA + kAcol];
                }
            }
        }
    }

    if ((info == 0) && (!(A[n * n - 1] != static_cast<T>(0.0)))) {
        info = n;
    }

    delete[] ipiv; ipiv = nullptr;

    return info;
}



#endif
