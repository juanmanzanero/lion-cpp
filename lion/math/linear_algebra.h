#ifndef LION_MATH_LINEAR_ALGEBRA_H
#define LION_MATH_LINEAR_ALGEBRA_H
#pragma once


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
inline void matmul(T *X, const T *A, const T *B, int num_rows_A, int num_cols_A_rows_B, int num_cols_B)
{
    //
    // Fills the result of multiplying two dense matrices in column-major
    // order, i.e., "X = A * B". Solution "X" must be preallocated on
    // entry, to size "num_rows_A x num_cols_B", and is also a column-major
    // matrix. "A" holds the first operand in column-major order and has
    // size "num_rows_A x num_cols_A_rows_B". "B" holds the second operand
    // in column-major order and has size "num_cols_A_rows_B x num_cols_B".
    // Hence, the number of columns of "A" must coincide with the number of
    // rows of "B", equal to input "num_cols_A_rows_B".
    //

    if (num_cols_A_rows_B != 1) {
        const auto c{ num_rows_A * (num_cols_B - 1) };
        auto br{ 0 };
        auto cr{ 0 };

        while (cr <= c) {
            const auto i0{ cr + num_rows_A };
            for (auto ic = cr; ic < i0; ic++) {
                X[ic] = T{ 0 };
            }
            cr += num_rows_A;
        }

        cr = 0;
        while (cr <= c) {
            auto ar{ -1 };
            const auto i0{ br + num_cols_A_rows_B };
            for (auto ib = br; ib < i0; ++ib) {
                if (B[ib] != T{ 0 }) {
                    auto ia{ ar };
                    const auto i1{ cr + num_rows_A };
                    for (auto ic = cr; ic < i1; ++ic) {
                        ++ia;
                        X[ic] += B[ib] * A[ia];
                    }
                }
                ar += num_rows_A;
            }
            br += num_cols_A_rows_B;
            cr += num_rows_A;
        }
    }
    else {
        for (auto i0 = 0; i0 < num_rows_A; ++i0) {
            for (auto i1 = 0; i1 < num_cols_B; ++i1) {
                X[i0 + num_rows_A * i1] = T{ 0 };
                for (auto cr = 0; cr < num_cols_A_rows_B; ++cr) {
                    X[i0 + num_rows_A * i1] +=
                        A[i0 + num_rows_A * cr] * B[cr + num_cols_A_rows_B * i1];
                }
            }
        }
    }
}


template<typename T>
inline int plufactorize(T *A, int *ipiv_A, int n)
{
    //
    // Performs a PLU decomposition of the dense square matrix "A",
    // of size "n x n", such that "A = P * L * U", where "P" is a
    // permutation matrix that can be reconstructed through "ipiv_A",
    // where the i-th row of matrix "A" is interchanged with row
    // "ipiv_A[i]". On entry, "A" should contain the square matrix
    // of size "n x n" in column-major order, and "ipiv_A" should
    // come preallocated to size "n". On exit, "A" will contain the
    // LU factors of the matrix, and "ipiv_A" the necessary index
    // permutations. The function returns 0 if the factorization is
    // successful (i.e., the matrix is non-singular).
    //

    int ret{ 0 };

    for (auto i = 0; i < n; ++i) {
        ipiv_A[i] = i + 1;
    }

    auto jA{ 0 };
    for (auto j = 0; j < n - 1; j++) {
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
                ipiv_A[j] = j + jA + 1;
                auto iX{ j };
                auto iy{ j + jA };
                for (auto k = 0; k < n; ++k) {
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

    if (ret == 0 && !(A[n * n - 1] != T{ 0 })) {
        ret = n;
    }

    return ret;
}


template<typename T>
inline void plusolve(T *B, const T *LU_A, const int *ipiv_A, int num_rows_A_rows_B, int num_cols_B)
{
    //
    // Direct solve of the dense square problem "A * X = B", via LU with
    // partial pivoting (i.e., employing a decomposition "A = P * L * U").
    // Input buffers "LU_A" and "ipiv_A" should hold the PLU factorization
    // of the problem's square matrix "A" (which has size "num_rows_A_rows_B x
    // num_rows_A_rows_B"), as returned by function "plu_factorize". On entry,
    // "B" contains the problem's RHS (which may have multiple colums) in
    // column-major order, and holds solution "X" on return (also in
    // column-major order) both of size "num_rows_A_rows_B x num_cols_B".
    //

    const auto &n = num_rows_A_rows_B;
    if (n < 1) {
        return;
    }

    for (std::size_t iy = 0; iy < static_cast<std::size_t>(n) - 1u; ++iy) {
        if (static_cast<std::size_t>(ipiv_A[iy]) != iy + 1u) {
            const auto kAcol{ ipiv_A[iy] - 1 };
            for (auto jA = 0; jA < num_cols_B; ++jA) {
                const auto smax{ B[iy + n * jA] };
                B[iy + n * jA] = B[kAcol + n * jA];
                B[kAcol + n * jA] = smax;
            }
        }
    }

    for (auto j = 0; j < num_cols_B; ++j) {
        const auto iy{ n * j };
        for (auto k = 0; k < n; ++k) {
            const auto kAcol{ n * k };
            if (B[k + iy] != T{ 0 }) {
                for (auto jA = k + 1; jA < n; ++jA) {
                    B[jA + iy] -= B[k + iy] * LU_A[jA + kAcol];
                }
            }
        }
    }

    for (auto j = 0; j < num_cols_B; ++j) {
        const auto iy{ n * j };
        for (auto k = n - 1; k + 1 > 0; --k) {
            const auto kAcol{ n * k };
            if (B[k + iy] != T{ 0 }) {
                B[k + iy] /= LU_A[k + kAcol];
                for (auto jA = 0; jA < k; ++jA) {
                    B[jA + iy] -= B[k + iy] * LU_A[jA + kAcol];
                }
            }
        }
    }
}


template<typename T>
inline int lusolve(T *B, T *A, int num_rows_A_rows_B, int num_cols_B)
{
    //
    // Direct solve of the dense square problem "A * X = B", via LU with
    // partial pivoting (i.e., applying a decomposition "A = P * L * U").
    // All matrices are in column-major order. On entry, "B" contains the
    // problem's RHS (which may have multiple colums), and holds solution
    // "X" on return, of size "num_rows_A_rows_B x num_cols_B". Matrix "A" must
    // hold the problem's square LHS matrix on entry, of size "num_rows_A_rows_B x
    // num_rows_A_rows_B" and exits as the LU factors ("A = P * L * U", where
    // "P" is a permutation matrix that could be reconstructed using the "ipiv"
    // local variable, since the i-th row of matrix "A" is interchanged with
    // row "ipiv[i]"). The function returns 0 upon successful exit.
    //

    auto *ipiv_A = new int[num_rows_A_rows_B];

    auto ret = plufactorize(A, ipiv_A, num_rows_A_rows_B);
    plusolve(B, A, ipiv_A, num_rows_A_rows_B, num_cols_B);

    delete[] ipiv_A;
    return ret;
}


template<typename T>
inline void tridiagonalsolve(T *B_colmaj, T *A_b, const T *A_a, const T *A_c, int num_rows_A_rows_B, int num_cols_B)
{
    //
    // Solves the tridiagonal system "A * x = B", in which
    // matrix "A" (of size "num_rows_A_rows_B x num_rows_A_rows_B")
    // comes specified via three vectors: "A_a" holding the lower
    // diagonal (of size "num_rows_A_rows_B - 1", starting in the
    // second row and ending in the last one), "A_b" holding the
    // diagonal (of size "num_rows_A_rows_B"), and "A_c" holding the
    // upper diagonal (of size "num_rows_A_rows_B - 1", starting in
    // the first row and ending in the penultimate one):
    //
    //      [b0, c0,  0,  0,  0, 0, ...          0,         0,         0;
    //       a0, b1, c1,  0,  0, 0, ...          0,         0,         0;
    //  A =   0, a1, b2, c2,  0, 0, ...          0,         0,         0;
    //        0,  0, a2, b3, c3, 0, ...          0,         0,         0;
    //                              ...
    //        0,  0,  0,  0,  0, 0, ..., a_{n - 2}, b_{n - 1}, c_{n - 1};
    //        0,  0,  0,  0,  0, 0, ...          0, a_{n - 1},       b_n]
    //
    // NOTE: the "A_b" vector holding the diagonal of the matrix
    // will get overwritten on exit by the function.
    //
    // The system may have multiple right-hand-sides, given as the
    // columns of the column-major matrix "B_colmaj" (of size
    // "num_rows_A_rows_B x num_cols_B"). On exit, the columns of this
    // matrix will hold the solution vectors corresponding to each
    // right-hand-side.
    //

    // forward elimination
    const auto &n = num_rows_A_rows_B;
    A_b[0] = T{ 1 } / A_b[0];
    for (auto i = 1; i < n; ++i) {
        const auto w = A_a[i - 1] * A_b[i - 1];
        A_b[i] = T{ 1 } / (A_b[i] - w * A_c[i - 1]);
        for (auto j = 0; j < num_cols_B; ++j) {
            const auto col = j * n;
            B_colmaj[i + col] -= w * B_colmaj[i - 1 + col];
        }
    }

    // backward substitution
    for (auto j = 0; j < num_cols_B; ++j) {
        const auto col = j * n;
        B_colmaj[n - 1 + col] *= A_b[n - 1];
        for (auto i = n - 2; i >= 0; --i) {
            B_colmaj[i + col] = (B_colmaj[i + col] -
                A_c[i] * B_colmaj[i + 1 + col]) * A_b[i];
        }
    }
}


template<typename T>
inline void qrsolve(T *X, T *A, T *B, int num_rows_A_rows_B, int num_cols_A, int num_cols_B)
{
    //
    // Fills "X" with the result of non-square problem "A * X = B", which
    // may be either under or overdetermined, employing the QR algorithm
    // with column pivoting. All matrices are in column-major order. Input
    // matrix "A" should hold the problem's LHS matrix on entry, of size
    // "num_rows_A_rows_B x num_cols_A", and gets overwritten by this function
    // with its QR fatorization. Input matrix "B" should contain the
    // problem's RHS, of size "num_rows_A_rows_B x num_cols_B" (i.e., it may
    // have multiple columns), and also gets overwritten by the function.
    // The result of this function is matrix "X", which should come
    // prealllocated on entry to size "num_cols_A x num_cols_B".
    //

    // declare some helper constants and functions
    constexpr auto TOL3Z = T{ 2.2204460492503131E-16 };

    const auto xgeqp3 = [](auto *A, auto m, auto n, auto mn, auto *tau, auto *jpvt)
    {
        // performs the column-pivoting QR factorization of the
        // "m x n" matrix A, such that "A * P = Q * R", where "P"
        // is a column permutation matrix. Info about pivoting is
        // stored in the "jpvt" array: on exit, if jpvt[j] = k,
        // then the j-th column of "A * P" was the k-th column
        // of "A"
        constexpr auto sqrt_TOL3Z = T{ 1.4901161193847656e-08 };
        constexpr auto SMALLNUM = T{ 1.0020841800044864E-292 };
        constexpr auto BIGNUM = T{ 9.9792015476736E+291 };

        const auto xnrm2 = [](auto n, const auto *x, auto ix0)
        {
            constexpr auto SMALLSCALE = T{ 3.3121686421112381E-170 };

            auto y = T{ 0 };
            if (n >= 1) {
                if (n == 1) {
                    y = std::abs(x[ix0 - 1]);
                }
                else {
                    auto scale{ SMALLSCALE };
                    auto kend{ ix0 + n };
                    for (auto k = ix0; k < kend; ++k) {
                        auto absxk{ std::abs(x[k - 1]) };
                        if (absxk > scale) {
                            auto t{ scale / absxk };
                            y = T{ 1 } + y * t * t;
                            scale = absxk;
                        }
                        else {
                            auto t{ absxk / scale };
                            y += t * t;
                        }
                    }
                    y = scale * std::sqrt(y);
                }
            }

            return y;
        };

        const auto rt_hypotd_snf = [](auto u0, auto u1)
        {
            auto a{ std::abs(u0) };
            auto b{ std::abs(u1) };
            if (a < b) {
                a /= b;
                return b * std::sqrt(a * a + T{ 1 });
            }
            else if (a > b) {
                b /= a;
                return a * std::sqrt(b * b + T{ 1 });
            }
            else if (std::isnan(b)) {
                return b;
            }
            else {
                return a * std::sqrt(T{ 2 });
            }
        };

        const auto xscal = [](auto n, auto a, auto *x, auto ix0)
        {
           for (auto k = ix0; k < ix0 + n; ++k) {
               x[k - 1] *= a;
           }
        };

        auto *work{ new T[n] };
        auto *vn1{ new T[n] };
        auto *vn2{ new T[n] };
        auto k{ 1 };
        for (auto i = 0; i < n; ++i) {
            jpvt[i] = i + 1;
            work[i] = T{ 0 };

            vn1[i] = xnrm2(m, A, k);
            vn2[i] = vn1[i];
            k += m;
        }

        for (auto i = 0; i < mn; ++i) {
            auto i_i{ i + i * m };
            auto nmi{ n - i };
            auto mmi{ m - i - 1 };
            auto itemp{ 1 };

            if (nmi < 1) {
                itemp = 0;
            }
            else {
                itemp = 1;
                if (nmi > 1) {
                    auto ix{ i };
                    auto smax{ std::abs(vn1[ix]) };
                    for (k = 2; k <= nmi; ++k) {
                        ++ix;
                        auto s{ std::abs(vn1[ix]) };
                        if (s > smax) {
                            itemp = k;
                            smax = s;
                        }
                    }

                }
            }

            auto pvt{ i + itemp - 1 };
            if (pvt != i) {
                auto ix{ m * pvt };
                auto iy{ m * i };
                for (k = 1; k <= m; ++k) {
                    auto smax{ A[ix] };
                    A[ix] = A[iy];
                    A[iy] = smax;
                    ++ix;
                    ++iy;
                }

                itemp = jpvt[pvt];
                jpvt[pvt] = jpvt[i];
                jpvt[i] = itemp;
                vn1[pvt] = vn1[i];
                vn2[pvt] = vn2[i];
            }

            if (i < m - 1) {
                auto atmp{ A[i_i] };
                auto d0 = T{ 0 };
                if (!(1 + mmi <= 0)) {
                    auto smax{ xnrm2(mmi, A, i_i + 2) };
                    if (smax != T{ 0 }) {
                        auto s{ rt_hypotd_snf(A[i_i], smax) };
                        if (A[i_i] >= T{ 0 }) {
                            s = -s;
                        }

                        if (std::abs(s) < SMALLNUM) {
                            itemp = 0;
                            do {
                                itemp++;
                                xscal(mmi, BIGNUM, A, i_i + 2);
                                s *= BIGNUM;
                                atmp *= BIGNUM;
                            } while (std::abs(s) < SMALLNUM);

                            s = rt_hypotd_snf(atmp, xnrm2(mmi, A, i_i + 2));
                            if (atmp >= T{ 0 }) {
                                s = -s;
                            }

                            d0 = (s - atmp) / s;
                            xscal(mmi, T{ 1 } / (atmp - s), A, i_i + 2);
                            for (k = 1; k <= itemp; ++k) {
                                s *= SMALLNUM;
                            }
                            atmp = s;

                        }
                        else {
                            d0 = (s - A[i_i]) / s;
                            smax = T{ 1 } / (A[i_i] - s);
                            xscal(mmi, smax, A, i_i + 2);
                            atmp = s;
                        }
                    }
                }

                tau[i] = d0;
                A[i_i] = atmp;
            }
            else {
                tau[i] = T{ 0 };
            }

            if (i < n - 1) {
                const auto atmp{ A[i_i] };
                auto lastc{ 0 };
                auto lastv{ 0 };
                auto i_ip1{ i + (i + 1) * m + 1 };

                A[i_i] = T{ 1 };

                if (tau[i] != T{ 0 }) {
                    lastc = nmi - 1;
                    lastv = mmi;
                    itemp = i_i + mmi;
                    while ((lastv + 1 > 0) && (A[itemp] == T{ 0 })) {
                        --lastv;
                        --itemp;
                    }

                    auto exitg2{ false };
                    while ((!exitg2) && (lastc > 0)) {
                        itemp = i_ip1 + (lastc - 1) * m;
                        k = itemp;
                        auto exitg1{ 0 };
                        do {
                            if (k <= itemp + lastv) {
                                if (A[k - 1] != T{ 0 }) {
                                    exitg1 = 1;
                                }
                                else {
                                    ++k;
                                }
                            }
                            else {
                                --lastc;
                                exitg1 = 2;
                            }
                        } while (exitg1 == 0);

                        if (exitg1 == 1) {
                            exitg2 = true;
                        }
                    }
                }
                else {
                    lastc = 0;
                    lastv = -1;
                }

                if (lastv >= 0) {
                    if (lastc != 0) {
                        for (auto iy = 1; iy <= lastc; ++iy) {
                            work[iy - 1] = T{ 0 };
                        }

                        auto iy{ 0 };
                        auto i1{ i_ip1 + m * (lastc - 1) };
                        itemp = i_ip1;
                        while (((m > 0) && (itemp <= i1)) || ((m < 0) && (itemp >= i1))) {
                            auto ix{ i_i };
                            auto smax = T{ 0 };
                            pvt = itemp + lastv;
                            for (k = itemp; k <= pvt; ++k) {
                                smax += A[k - 1] * A[ix];
                                ++ix;
                            }

                            work[iy] += smax;
                            ++iy;
                            itemp += m;
                        }
                    }

                    if (-tau[i] != T{ 0.0 }) {
                        itemp = 0;
                        for (nmi = 1; nmi <= lastc; nmi++) {
                            if (work[itemp] != T{ 0 }) {
                                const auto smax{ work[itemp] * -tau[i] };
                                auto ix{ i_i };
                                for (pvt = i_ip1; pvt <= lastv + i_ip1; ++pvt) {
                                    A[pvt - 1] += A[ix] * smax;
                                    ++ix;
                                }
                            }

                            ++itemp;
                            i_ip1 += m;
                        }
                    }
                }

                A[i_i] = atmp;
            }

            for (nmi = i + 1; nmi < n; ++nmi) {
                if (vn1[nmi] != T{ 0 }) {
                    auto smax{ std::abs(A[i + m * nmi]) / vn1[nmi] };
                    smax = T{ 1 } - smax * smax;
                    if (smax < T{ 0 }) {
                        smax = 0.0;
                    }

                    auto s{ vn1[nmi] / vn2[nmi] };
                    s = smax * s * s;
                    if (s <= sqrt_TOL3Z) {
                        if (i + 1 < m) {
                            vn1[nmi] = xnrm2(mmi, A, i + m * nmi + 2);
                            vn2[nmi] = vn1[nmi];
                        }
                        else {
                            vn1[nmi] = T{ 0 };
                            vn2[nmi] = T{ 0 };
                        }
                    }
                    else {
                        vn1[nmi] *= std::sqrt(smax);
                    }
                }
            }
        }

        delete[] work;
        delete[] vn1;
        delete[] vn2;
    };

    // perform the column-pivoting QR factorization of A
    const auto &m{ num_rows_A_rows_B };
    const auto &n{ num_cols_A };
    const auto mn{ std::min(m, n) };
    auto *jpvt{ new int[n] };
    auto *tau{ new T[mn] };
    xgeqp3(A, m, n, mn, tau, jpvt);

    // rank from QR
    auto rankR{ 0 };
    auto tol{ static_cast<T>(std::max(m, n)) * std::abs(A[0]) * TOL3Z };
    while ((rankR < mn) && (!(std::abs(A[rankR + m * rankR]) <= tol))) {
        ++rankR;
    }

    for (auto i = 0; i < n * num_cols_B; ++i) {
        X[i] = T{ 0 };
    }

    for (auto k = 0; k < mn; ++k) {
        if (tau[k] != T{ 0 }) {
            for (auto j = 0; j < num_cols_B; ++j) {
                tol = B[k + m * j];
                for (auto i = k + 1; i < m; ++i) {
                    tol += A[i + m * k] * B[i + m * j];
                }

                tol *= tau[k];
                if (tol != T{ 0 }) {
                    B[k + m * j] -= tol;
                    for (auto i = k + 1; i < m; ++i) {
                        B[i + m * j] -= A[i + m * k] * tol;
                    }
                }
            }
        }
    }

    for (auto k = 0; k < num_cols_B; ++k) {
        for (auto i = 0; i < rankR; ++i) {
            X[jpvt[i] + n * k - 1] = B[i + m * k];
        }

        for (auto j = rankR - 1; j >= 0; --j) {
            X[jpvt[j] + n * k - 1] /= A[j + m * j];
            for (auto i = 0; i < j; ++i) {
                X[jpvt[i] + n * k - 1] -= X[jpvt[j] + n * k - 1] *
                    A[i + m * j];
            }
        }
    }

    delete[] tau;
    delete[] jpvt;
}


template<typename T>
inline T rcond(const T *A, int n)
{
    //
    // Returns an estimate for the reciprocal condition of input
    // square matrix "A" in 1-norm. This input square matrix has
    // size "n x n" and should be specified in column-major order.
    // If "A" is well conditioned, then its "rcond" number will be
    // near 1.0. When "A" is badly conditioned, its "rcond" will
    // be near 0.
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
                                        x[ix] = s * (T{ 1 } + ((T{ 1 } +static_cast<T>(ix)) - T{ 1 }) /
                                            (static_cast<T>(n) - T{ 1 }));
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
                                        x[ix] = s * (T{ 1 } +((T{ 1 } +static_cast<T>(ix)) - T{ 1 }) /
                                            (static_cast<T>(n) - T{ 1 }));
                                        s = -s;
                                    }

                                    kase = 1;
                                    jump = 5;
                                }
                            }
                            else {
                                if (jump == 5) {
                                    const auto s = T{ 2 } * b_norm(x, n) / T{ 3 } /
                                        static_cast<T>(n);
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


template<typename T, typename SizeType>
inline void lufactorize_Doolittle(T *A, SizeType n)
{
    //
    // Does a raw LU factorization of the dense & square column-major
    // matrix "A", of size "n x n". Unsafer than "plufactorize" since
    // it doesn't do pivoting, so use it with caution, preferably with
    // small and trusty matrices... however, this algorithm doesn't have
    // "ifs", which makes it more adequate to tape CppAD::ADFuns. On
    // entry, "A" should contain the square matrix of size "n x n" in
    // column-major order. On exit, "A" will contain the LU factors of
    // the matrix, also in column-major order.
    //

    for (auto i = decltype(n){ 0 }; i < n; ++i) {
        for (auto j = decltype(i){ 0 }; j < i; ++j) {
            const auto jn = j * n;
            const auto ij = i + jn;
            for (auto k = decltype(j){ 0 }; k < j; ++k) {
                A[ij] -= A[i + k * n] * A[k + jn];
            }
            A[ij] /= A[j + jn];
        }

        for (auto j = i; j < n; ++j) {
            const auto jn = j * n;
            const auto ij = i + jn;
            for (auto k = decltype(i){ 0 }; k < i; ++k) {
                A[ij] -= A[i + k * n] * A[k + jn];
            }
        }
    }
}


template<typename T, typename SizeType0, typename SizeType1>
inline void lusolve_Doolittle(T *B, T *LU_A, SizeType0 num_rows_A_rows_B, SizeType1 num_cols_B)
{
    //
    // Direct solve of the dense square problem "A * X = B", with
    // "A = L * U" (i.e., raw LU decomposition without pivoting). Input
    // buffer "LU_A" should hold the column-major LU factorization of the
    // problem's square matrix "A" (which has size "num_rows_A_rows_B x
    // num_rows_A_rows_B"), e.g. as returned by function "lufactorize_Doolittle".
    // On entry, "B" contains the problem's RHS (which may have multiple colums)
    // in column-major order, and holds solution "X" on return (also in
    // column-major order) both of size "num_rows_A_rows_B x num_cols_B".
    //

    for (auto j = decltype(num_cols_B){ 0 }; j < num_cols_B; ++j) {
        const auto jn = j * num_rows_A_rows_B;
        for (auto i = decltype(num_rows_A_rows_B){ 0 }; i < num_rows_A_rows_B; ++i) {
            const auto ij = i + jn;
            for (auto k = decltype(i){ 0 }; k < i; ++k) {
                B[ij] -= LU_A[i + k * num_rows_A_rows_B] * B[k + jn];
            }
        }

        for (auto i1 = num_rows_A_rows_B; i1 > decltype(num_rows_A_rows_B){ 0 }; --i1) {
            const auto i = i1 - decltype(num_rows_A_rows_B){ 1 };
            const auto ij = i + jn;
            for (auto k = i1; k < num_rows_A_rows_B; ++k) {
                B[ij] -= LU_A[i + k * num_rows_A_rows_B] * B[k + jn];
            }
            B[ij] /= LU_A[i + i * num_rows_A_rows_B];
        }
    }
}

#endif
