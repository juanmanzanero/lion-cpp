#ifndef LION_MATH_SMALL_DENSE_INVERSE_MATRICES_H
#define LION_MATH_SMALL_DENSE_INVERSE_MATRICES_H
#pragma once


//
// Defines functions that return the determinants
// and inverses of some small dense matrices, in
// column-major order. I've generated these functions
// using Matlab's symbolic toolbox.
//


template<typename Matrix2x2Type>
constexpr Matrix2x2Type::value_type det2x2(const Matrix2x2Type &A)
{
    //
    // Returns the determinant of the input 2 x 2 matrix,
    // which should be specified in column-major order.
    //

    const auto &a00 = A[0];
    const auto &a10 = A[1];

    const auto &a01 = A[2];
    const auto &a11 = A[3];

    return a00 * a11 - a01 * a10;
}


template<typename Matrix2x2Type>
constexpr Matrix2x2Type inv2x2(const Matrix2x2Type &A)
{
    //
    // Returns the inverse of the input 2 x 2 matrix,
    // both the input matrix and the resulting inverse
    // should be specified in column-major order.
    //

    const auto invdet = typename Matrix2x2Type::value_type{ 1 } /
        det2x2(A);

    const auto &a00 = A[0];
    const auto &a10 = A[1];

    const auto &a01 = A[2];
    const auto &a11 = A[3];

    return Matrix2x2Type{
        a11 * invdet,
        -a10 * invdet,
        -a01 * invdet,
        a00 * invdet };
}


template<typename Matrix3x3Type>
constexpr Matrix3x3Type::value_type det3x3(const Matrix3x3Type &A)
{
    //
    // Returns the determinant of the input 3 x 3 matrix,
    // which should be specified in column-major order.
    //

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];

    const auto &a01 = A[3];
    const auto &a11 = A[4];
    const auto &a21 = A[5];

    const auto &a02 = A[6];
    const auto &a12 = A[7];
    const auto &a22 = A[8];

    return
        a00 * (a11 * a22 - a12 * a21) +
        a01 * (a12 * a20 - a10 * a22) +
        a02 * (a10 * a21 - a11 * a20);
}


template<typename Matrix3x3Type>
constexpr Matrix3x3Type inv3x3(const Matrix3x3Type &A)
{
    //
    // Returns the inverse of the input 3 x 3 matrix,
    // both the input matrix and the resulting inverse
    // should be specified in column-major order.
    //

    const auto invdet = typename Matrix3x3Type::value_type{ 1 } /
        det3x3(A);

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];

    const auto &a01 = A[3];
    const auto &a11 = A[4];
    const auto &a21 = A[5];

    const auto &a02 = A[6];
    const auto &a12 = A[7];
    const auto &a22 = A[8];

    return Matrix3x3Type{
        invdet * (a11 * a22 - a12 * a21),
        invdet * (a12 * a20 - a10 * a22),
        invdet * (a10 * a21 - a11 * a20),
        invdet * (a02 * a21 - a01 * a22),
        invdet * (a00 * a22 - a02 * a20),
        invdet * (a01 * a20 - a00 * a21),
        invdet * (a01 * a12 - a02 * a11),
        invdet * (a02 * a10 - a00 * a12),
        invdet * (a00 * a11 - a01 * a10) };
}


template<typename Matrix4x4Type>
constexpr Matrix4x4Type::value_type det4x4(const Matrix4x4Type &A)
{
    //
    // Returns the determinant of the input 4 x 4 matrix,
    // which should be specified in column-major order.
    //

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];
    const auto &a30 = A[3];

    const auto &a01 = A[4];
    const auto &a11 = A[5];
    const auto &a21 = A[6];
    const auto &a31 = A[7];

    const auto &a02 = A[8];
    const auto &a12 = A[9];
    const auto &a22 = A[10];
    const auto &a32 = A[11];

    const auto &a03 = A[12];
    const auto &a13 = A[13];
    const auto &a23 = A[14];
    const auto &a33 = A[15];

    return
        a00 * (a11 * (a22 * a33 - a23 * a32) +
               a12 * (a23 * a31 - a21 * a33) +
               a13 * (a21 * a32 - a22 * a31)) +
        a01 * (a10 * (a23 * a32 - a22 * a33) +
               a12 * (a20 * a33 - a23 * a30) +
               a13 * (a22 * a30 - a20 * a32)) +
        a02 * (a10 * (a21 * a33 - a23 * a31) +
               a11 * (a23 * a30 - a20 * a33) +
               a13 * (a20 * a31 - a21 * a30)) +
        a03 * (a10 * (a22 * a31 - a21 * a32) +
               a11 * (a20 * a32 - a22 * a30) +
               a12 * (a21 * a30 - a20 * a31));
}


template<typename Matrix4x4Type>
constexpr Matrix4x4Type inv4x4(const Matrix4x4Type &A)
{
    //
    // Returns the inverse of the input 4 x 4 matrix,
    // both the input matrix and the resulting inverse
    // should be specified in column-major order.
    //

    const auto invdet = typename Matrix4x4Type::value_type{ 1 } /
        det4x4(A);

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];
    const auto &a30 = A[3];

    const auto &a01 = A[4];
    const auto &a11 = A[5];
    const auto &a21 = A[6];
    const auto &a31 = A[7];

    const auto &a02 = A[8];
    const auto &a12 = A[9];
    const auto &a22 = A[10];
    const auto &a32 = A[11];

    const auto &a03 = A[12];
    const auto &a13 = A[13];
    const auto &a23 = A[14];
    const auto &a33 = A[15];

    return Matrix4x4Type{
        invdet * (a11 * (a22 * a33 - a23 * a32) +
                  a12 * (a23 * a31 - a21 * a33) +
                  a13 * (a21 * a32 - a22 * a31)),
        invdet * (a10 * (a23 * a32 - a22 * a33) +
                  a12 * (a20 * a33 - a23 * a30) +
                  a13 * (a22 * a30 - a20 * a32)),
        invdet * (a10 * (a21 * a33 - a23 * a31) +
                  a11 * (a23 * a30 - a20 * a33) +
                  a13 * (a20 * a31 - a21 * a30)),
        invdet * (a10 * (a22 * a31 - a21 * a32) +
                  a11 * (a20 * a32 - a22 * a30) +
                  a12 * (a21 * a30 - a20 * a31)),
        invdet * (a01 * (a23 * a32 - a22 * a33) +
                  a02 * (a21 * a33 - a23 * a31) +
                  a03 * (a22 * a31 - a21 * a32)),
        invdet * (a00 * (a22 * a33 - a23 * a32) +
                  a02 * (a23 * a30 - a20 * a33) +
                  a03 * (a20 * a32 - a22 * a30)),
        invdet * (a00 * (a23 * a31 - a21 * a33) +
                  a01 * (a20 * a33 - a23 * a30) +
                  a03 * (a21 * a30 - a20 * a31)),
        invdet * (a00 * (a21 * a32 - a22 * a31) +
                  a01 * (a22 * a30 - a20 * a32) +
                  a02 * (a20 * a31 - a21 * a30)),
        invdet * (a01 * (a12 * a33 - a13 * a32) +
                  a02 * (a13 * a31 - a11 * a33) +
                  a03 * (a11 * a32 - a12 * a31)),
        invdet * (a00 * (a13 * a32 - a12 * a33) +
                  a02 * (a10 * a33 - a13 * a30) +
                  a03 * (a12 * a30 - a10 * a32)),
        invdet * (a00 * (a11 * a33 - a13 * a31) +
                  a01 * (a13 * a30 - a10 * a33) +
                  a03 * (a10 * a31 - a11 * a30)),
        invdet * (a00 * (a12 * a31 - a11 * a32) +
                  a01 * (a10 * a32 - a12 * a30) +
                  a02 * (a11 * a30 - a10 * a31)),
        invdet * (a01 * (a13 * a22 - a12 * a23) +
                  a02 * (a11 * a23 - a13 * a21) +
                  a03 * (a12 * a21 - a11 * a22)),
        invdet * (a00 * (a12 * a23 - a13 * a22) +
                  a02 * (a13 * a20 - a10 * a23) +
                  a03 * (a10 * a22 - a12 * a20)),
        invdet * (a00 * (a13 * a21 - a11 * a23) +
                  a01 * (a10 * a23 - a13 * a20) +
                  a03 * (a11 * a20 - a10 * a21)),
        invdet * (a00 * (a11 * a22 - a12 * a21) +
                  a01 * (a12 * a20 - a10 * a22) +
                  a02 * (a10 * a21 - a11 * a20)) };
}


template<typename Matrix5x5Type>
constexpr Matrix5x5Type::value_type det5x5(const Matrix5x5Type &A)
{
    //
    // Returns the determinant of the input 5 x 5 matrix,
    // which should be specified in column-major order.
    //

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];
    const auto &a30 = A[3];
    const auto &a40 = A[4];

    const auto &a01 = A[5];
    const auto &a11 = A[6];
    const auto &a21 = A[7];
    const auto &a31 = A[8];
    const auto &a41 = A[9];

    const auto &a02 = A[10];
    const auto &a12 = A[11];
    const auto &a22 = A[12];
    const auto &a32 = A[13];
    const auto &a42 = A[14];

    const auto &a03 = A[15];
    const auto &a13 = A[16];
    const auto &a23 = A[17];
    const auto &a33 = A[18];
    const auto &a43 = A[19];

    const auto &a04 = A[20];
    const auto &a14 = A[21];
    const auto &a24 = A[22];
    const auto &a34 = A[23];
    const auto &a44 = A[24];

    return
        a00 * (a11 * (a22 * (a33 * a44 - a34 * a43) +
                      a23 * (a34 * a42 - a32 * a44) +
                      a24 * (a32 * a43 - a33 * a42)) +
               a12 * (a21 * (a34 * a43 - a33 * a44) +
                      a23 * (a31 * a44 - a34 * a41) +
                      a24 * (a33 * a41 - a31 * a43)) +
               a13 * (a21 * (a32 * a44 - a34 * a42) +
                      a22 * (a34 * a41 - a31 * a44) +
                      a24 * (a31 * a42 - a32 * a41)) +
               a14 * (a21 * (a33 * a42 - a32 * a43) +
                      a22 * (a31 * a43 - a33 * a41) +
                      a23 * (a32 * a41 - a31 * a42))) +
        a01 * (a10 * (a22 * (a34 * a43 - a33 * a44) +
                      a23 * (a32 * a44 - a34 * a42) +
                      a24 * (a33 * a42 - a32 * a43)) +
               a12 * (a20 * (a33 * a44 - a34 * a43) +
                      a23 * (a34 * a40 - a30 * a44) +
                      a24 * (a30 * a43 - a33 * a40)) +
               a13 * (a20 * (a34 * a42 - a32 * a44) +
                      a22 * (a30 * a44 - a34 * a40) +
                      a24 * (a32 * a40 - a30 * a42)) +
               a14 * (a20 * (a32 * a43 - a33 * a42) +
                      a22 * (a33 * a40 - a30 * a43) +
                      a23 * (a30 * a42 - a32 * a40))) +
        a02 * (a10 * (a21 * (a33 * a44 - a34 * a43) +
                      a23 * (a34 * a41 - a31 * a44) +
                      a24 * (a31 * a43 - a33 * a41)) +
               a11 * (a20 * (a34 * a43 - a33 * a44) +
                      a23 * (a30 * a44 - a34 * a40) +
                      a24 * (a33 * a40 - a30 * a43)) +
               a13 * (a20 * (a31 * a44 - a34 * a41) +
                      a21 * (a34 * a40 - a30 * a44) +
                      a24 * (a30 * a41 - a31 * a40)) +
               a14 * (a20 * (a33 * a41 - a31 * a43) +
                      a21 * (a30 * a43 - a33 * a40) +
                      a23 * (a31 * a40 - a30 * a41))) +
        a03 * (a10 * (a21 * (a34 * a42 - a32 * a44) +
                      a22 * (a31 * a44 - a34 * a41) +
                      a24 * (a32 * a41 - a31 * a42)) +
               a11 * (a20 * (a32 * a44 - a34 * a42) +
                      a22 * (a34 * a40 - a30 * a44) +
                      a24 * (a30 * a42 - a32 * a40)) +
               a12 * (a20 * (a34 * a41 - a31 * a44) +
                      a21 * (a30 * a44 - a34 * a40) +
                      a24 * (a31 * a40 - a30 * a41)) +
               a14 * (a20 * (a31 * a42 - a32 * a41) +
                      a21 * (a32 * a40 - a30 * a42) +
                      a22 * (a30 * a41 - a31 * a40))) +
        a04 * (a10 * (a21 * (a32 * a43 - a33 * a42) +
                      a22 * (a33 * a41 - a31 * a43) +
                      a23 * (a31 * a42 - a32 * a41)) +
               a11 * (a20 * (a33 * a42 - a32 * a43) +
                      a22 * (a30 * a43 - a33 * a40) +
                      a23 * (a32 * a40 - a30 * a42)) +
               a12 * (a20 * (a31 * a43 - a33 * a41) +
                      a21 * (a33 * a40 - a30 * a43) +
                      a23 * (a30 * a41 - a31 * a40)) +
               a13 * (a20 * (a32 * a41 - a31 * a42) +
                      a21 * (a30 * a42 - a32 * a40) +
                      a22 * (a31 * a40 - a30 * a41)));
}


template<typename Matrix5x5Type>
constexpr Matrix5x5Type inv5x5(const Matrix5x5Type &A)
{
    //
    // Returns the inverse of the input 5 x 5 matrix,
    // both the input matrix and the resulting inverse
    // should be specified in column-major order.
    //

    const auto invdet = typename Matrix5x5Type::value_type{ 1 } /
        det5x5(A);

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];
    const auto &a30 = A[3];
    const auto &a40 = A[4];

    const auto &a01 = A[5];
    const auto &a11 = A[6];
    const auto &a21 = A[7];
    const auto &a31 = A[8];
    const auto &a41 = A[9];

    const auto &a02 = A[10];
    const auto &a12 = A[11];
    const auto &a22 = A[12];
    const auto &a32 = A[13];
    const auto &a42 = A[14];

    const auto &a03 = A[15];
    const auto &a13 = A[16];
    const auto &a23 = A[17];
    const auto &a33 = A[18];
    const auto &a43 = A[19];

    const auto &a04 = A[20];
    const auto &a14 = A[21];
    const auto &a24 = A[22];
    const auto &a34 = A[23];
    const auto &a44 = A[24];

    return Matrix5x5Type{
        invdet * (a11 * (a22 * (a33 * a44 - a34 * a43) +
                         a23 * (a34 * a42 - a32 * a44) +
                         a24 * (a32 * a43 - a33 * a42)) +
                  a12 * (a21 * (a34 * a43 - a33 * a44) +
                         a23 * (a31 * a44 - a34 * a41) +
                         a24 * (a33 * a41 - a31 * a43)) +
                  a13 * (a21 * (a32 * a44 - a34 * a42) +
                         a22 * (a34 * a41 - a31 * a44) +
                         a24 * (a31 * a42 - a32 * a41)) +
                  a14 * (a21 * (a33 * a42 - a32 * a43) +
                         a22 * (a31 * a43 - a33 * a41) +
                         a23 * (a32 * a41 - a31 * a42))),
        invdet * (a10 * (a22 * (a34 * a43 - a33 * a44) +
                         a23 * (a32 * a44 - a34 * a42) +
                         a24 * (a33 * a42 - a32 * a43)) +
                  a12 * (a20 * (a33 * a44 - a34 * a43) +
                         a23 * (a34 * a40 - a30 * a44) +
                         a24 * (a30 * a43 - a33 * a40)) +
                  a13 * (a20 * (a34 * a42 - a32 * a44) +
                         a22 * (a30 * a44 - a34 * a40) +
                         a24 * (a32 * a40 - a30 * a42)) +
                  a14 * (a20 * (a32 * a43 - a33 * a42) +
                         a22 * (a33 * a40 - a30 * a43) +
                         a23 * (a30 * a42 - a32 * a40))),
        invdet * (a10 * (a21 * (a33 * a44 - a34 * a43) +
                         a23 * (a34 * a41 - a31 * a44) +
                         a24 * (a31 * a43 - a33 * a41)) +
                  a11 * (a20 * (a34 * a43 - a33 * a44) +
                         a23 * (a30 * a44 - a34 * a40) +
                         a24 * (a33 * a40 - a30 * a43)) +
                  a13 * (a20 * (a31 * a44 - a34 * a41) +
                         a21 * (a34 * a40 - a30 * a44) +
                         a24 * (a30 * a41 - a31 * a40)) +
                  a14 * (a20 * (a33 * a41 - a31 * a43) +
                         a21 * (a30 * a43 - a33 * a40) +
                         a23 * (a31 * a40 - a30 * a41))),
        invdet * (a10 * (a21 * (a34 * a42 - a32 * a44) +
                         a22 * (a31 * a44 - a34 * a41) +
                         a24 * (a32 * a41 - a31 * a42)) +
                  a11 * (a20 * (a32 * a44 - a34 * a42) +
                         a22 * (a34 * a40 - a30 * a44) +
                         a24 * (a30 * a42 - a32 * a40)) +
                  a12 * (a20 * (a34 * a41 - a31 * a44) +
                         a21 * (a30 * a44 - a34 * a40) +
                         a24 * (a31 * a40 - a30 * a41)) +
                  a14 * (a20 * (a31 * a42 - a32 * a41) +
                         a21 * (a32 * a40 - a30 * a42) +
                         a22 * (a30 * a41 - a31 * a40))),
        invdet * (a10 * (a21 * (a32 * a43 - a33 * a42) +
                         a22 * (a33 * a41 - a31 * a43) +
                         a23 * (a31 * a42 - a32 * a41)) +
                  a11 * (a20 * (a33 * a42 - a32 * a43) +
                         a22 * (a30 * a43 - a33 * a40) +
                         a23 * (a32 * a40 - a30 * a42)) +
                  a12 * (a20 * (a31 * a43 - a33 * a41) +
                         a21 * (a33 * a40 - a30 * a43) +
                         a23 * (a30 * a41 - a31 * a40)) +
                  a13 * (a20 * (a32 * a41 - a31 * a42) +
                         a21 * (a30 * a42 - a32 * a40) +
                         a22 * (a31 * a40 - a30 * a41))),
        invdet * (a01 * (a22 * (a34 * a43 - a33 * a44) +
                         a23 * (a32 * a44 - a34 * a42) +
                         a24 * (a33 * a42 - a32 * a43)) +
                  a02 * (a21 * (a33 * a44 - a34 * a43) +
                         a23 * (a34 * a41 - a31 * a44) +
                         a24 * (a31 * a43 - a33 * a41)) +
                  a03 * (a21 * (a34 * a42 - a32 * a44) +
                         a22 * (a31 * a44 - a34 * a41) +
                         a24 * (a32 * a41 - a31 * a42)) +
                  a04 * (a21 * (a32 * a43 - a33 * a42) +
                         a22 * (a33 * a41 - a31 * a43) +
                         a23 * (a31 * a42 - a32 * a41))),
        invdet * (a00 * (a22 * (a33 * a44 - a34 * a43) +
                         a23 * (a34 * a42 - a32 * a44) +
                         a24 * (a32 * a43 - a33 * a42)) +
                  a02 * (a20 * (a34 * a43 - a33 * a44) +
                         a23 * (a30 * a44 - a34 * a40) +
                         a24 * (a33 * a40 - a30 * a43)) +
                  a03 * (a20 * (a32 * a44 - a34 * a42) +
                         a22 * (a34 * a40 - a30 * a44) +
                         a24 * (a30 * a42 - a32 * a40)) +
                  a04 * (a20 * (a33 * a42 - a32 * a43) +
                         a22 * (a30 * a43 - a33 * a40) +
                         a23 * (a32 * a40 - a30 * a42))),
        invdet * (a00 * (a21 * (a34 * a43 - a33 * a44) +
                         a23 * (a31 * a44 - a34 * a41) +
                         a24 * (a33 * a41 - a31 * a43)) +
                  a01 * (a20 * (a33 * a44 - a34 * a43) +
                         a23 * (a34 * a40 - a30 * a44) +
                         a24 * (a30 * a43 - a33 * a40)) +
                  a03 * (a20 * (a34 * a41 - a31 * a44) +
                         a21 * (a30 * a44 - a34 * a40) +
                         a24 * (a31 * a40 - a30 * a41)) +
                  a04 * (a20 * (a31 * a43 - a33 * a41) +
                         a21 * (a33 * a40 - a30 * a43) +
                         a23 * (a30 * a41 - a31 * a40))),
        invdet * (a00 * (a21 * (a32 * a44 - a34 * a42) +
                         a22 * (a34 * a41 - a31 * a44) +
                         a24 * (a31 * a42 - a32 * a41)) +
                  a01 * (a20 * (a34 * a42 - a32 * a44) +
                         a22 * (a30 * a44 - a34 * a40) +
                         a24 * (a32 * a40 - a30 * a42)) +
                  a02 * (a20 * (a31 * a44 - a34 * a41) +
                         a21 * (a34 * a40 - a30 * a44) +
                         a24 * (a30 * a41 - a31 * a40)) +
                  a04 * (a20 * (a32 * a41 - a31 * a42) +
                         a21 * (a30 * a42 - a32 * a40) +
                         a22 * (a31 * a40 - a30 * a41))),
        invdet * (a00 * (a21 * (a33 * a42 - a32 * a43) +
                         a22 * (a31 * a43 - a33 * a41) +
                         a23 * (a32 * a41 - a31 * a42)) +
                  a01 * (a20 * (a32 * a43 - a33 * a42) +
                         a22 * (a33 * a40 - a30 * a43) +
                         a23 * (a30 * a42 - a32 * a40)) +
                  a02 * (a20 * (a33 * a41 - a31 * a43) +
                         a21 * (a30 * a43 - a33 * a40) +
                         a23 * (a31 * a40 - a30 * a41)) +
                  a03 * (a20 * (a31 * a42 - a32 * a41) +
                         a21 * (a32 * a40 - a30 * a42) +
                         a22 * (a30 * a41 - a31 * a40))),
        invdet * (a01 * (a12 * (a33 * a44 - a34 * a43) +
                         a13 * (a34 * a42 - a32 * a44) +
                         a14 * (a32 * a43 - a33 * a42)) +
                  a02 * (a11 * (a34 * a43 - a33 * a44) +
                         a13 * (a31 * a44 - a34 * a41) +
                         a14 * (a33 * a41 - a31 * a43)) +
                  a03 * (a11 * (a32 * a44 - a34 * a42) +
                         a12 * (a34 * a41 - a31 * a44) +
                         a14 * (a31 * a42 - a32 * a41)) +
                  a04 * (a11 * (a33 * a42 - a32 * a43) +
                         a12 * (a31 * a43 - a33 * a41) +
                         a13 * (a32 * a41 - a31 * a42))),
        invdet * (a00 * (a12 * (a34 * a43 - a33 * a44) +
                         a13 * (a32 * a44 - a34 * a42) +
                         a14 * (a33 * a42 - a32 * a43)) +
                  a02 * (a10 * (a33 * a44 - a34 * a43) +
                         a13 * (a34 * a40 - a30 * a44) +
                         a14 * (a30 * a43 - a33 * a40)) +
                  a03 * (a10 * (a34 * a42 - a32 * a44) +
                         a12 * (a30 * a44 - a34 * a40) +
                         a14 * (a32 * a40 - a30 * a42)) +
                  a04 * (a10 * (a32 * a43 - a33 * a42) +
                         a12 * (a33 * a40 - a30 * a43) +
                         a13 * (a30 * a42 - a32 * a40))),
        invdet * (a00 * (a11 * (a33 * a44 - a34 * a43) +
                         a13 * (a34 * a41 - a31 * a44) +
                         a14 * (a31 * a43 - a33 * a41)) +
                  a01 * (a10 * (a34 * a43 - a33 * a44) +
                         a13 * (a30 * a44 - a34 * a40) +
                         a14 * (a33 * a40 - a30 * a43)) +
                  a03 * (a10 * (a31 * a44 - a34 * a41) +
                         a11 * (a34 * a40 - a30 * a44) +
                         a14 * (a30 * a41 - a31 * a40)) +
                  a04 * (a10 * (a33 * a41 - a31 * a43) +
                         a11 * (a30 * a43 - a33 * a40) +
                         a13 * (a31 * a40 - a30 * a41))),
        invdet * (a00 * (a11 * (a34 * a42 - a32 * a44) +
                         a12 * (a31 * a44 - a34 * a41) +
                         a14 * (a32 * a41 - a31 * a42)) +
                  a01 * (a10 * (a32 * a44 - a34 * a42) +
                         a12 * (a34 * a40 - a30 * a44) +
                         a14 * (a30 * a42 - a32 * a40)) +
                  a02 * (a10 * (a34 * a41 - a31 * a44) +
                         a11 * (a30 * a44 - a34 * a40) +
                         a14 * (a31 * a40 - a30 * a41)) +
                  a04 * (a10 * (a31 * a42 - a32 * a41) +
                         a11 * (a32 * a40 - a30 * a42) +
                         a12 * (a30 * a41 - a31 * a40))),
        invdet * (a00 * (a11 * (a32 * a43 - a33 * a42) +
                         a12 * (a33 * a41 - a31 * a43) +
                         a13 * (a31 * a42 - a32 * a41)) +
                  a01 * (a10 * (a33 * a42 - a32 * a43) +
                         a12 * (a30 * a43 - a33 * a40) +
                         a13 * (a32 * a40 - a30 * a42)) +
                  a02 * (a10 * (a31 * a43 - a33 * a41) +
                         a11 * (a33 * a40 - a30 * a43) +
                         a13 * (a30 * a41 - a31 * a40)) +
                  a03 * (a10 * (a32 * a41 - a31 * a42) +
                         a11 * (a30 * a42 - a32 * a40) +
                         a12 * (a31 * a40 - a30 * a41))),
        invdet * (a01 * (a12 * (a24 * a43 - a23 * a44) +
                         a13 * (a22 * a44 - a24 * a42) +
                         a14 * (a23 * a42 - a22 * a43)) +
                  a02 * (a11 * (a23 * a44 - a24 * a43) +
                         a13 * (a24 * a41 - a21 * a44) +
                         a14 * (a21 * a43 - a23 * a41)) +
                  a03 * (a11 * (a24 * a42 - a22 * a44) +
                         a12 * (a21 * a44 - a24 * a41) +
                         a14 * (a22 * a41 - a21 * a42)) +
                  a04 * (a11 * (a22 * a43 - a23 * a42) +
                         a12 * (a23 * a41 - a21 * a43) +
                         a13 * (a21 * a42 - a22 * a41))),
        invdet * (a00 * (a12 * (a23 * a44 - a24 * a43) +
                         a13 * (a24 * a42 - a22 * a44) +
                         a14 * (a22 * a43 - a23 * a42)) +
                  a02 * (a10 * (a24 * a43 - a23 * a44) +
                         a13 * (a20 * a44 - a24 * a40) +
                         a14 * (a23 * a40 - a20 * a43)) +
                  a03 * (a10 * (a22 * a44 - a24 * a42) +
                         a12 * (a24 * a40 - a20 * a44) +
                         a14 * (a20 * a42 - a22 * a40)) +
                  a04 * (a10 * (a23 * a42 - a22 * a43) +
                         a12 * (a20 * a43 - a23 * a40) +
                         a13 * (a22 * a40 - a20 * a42))),
        invdet * (a00 * (a11 * (a24 * a43 - a23 * a44) +
                         a13 * (a21 * a44 - a24 * a41) +
                         a14 * (a23 * a41 - a21 * a43)) +
                  a01 * (a10 * (a23 * a44 - a24 * a43) +
                         a13 * (a24 * a40 - a20 * a44) +
                         a14 * (a20 * a43 - a23 * a40)) +
                  a03 * (a10 * (a24 * a41 - a21 * a44) +
                         a11 * (a20 * a44 - a24 * a40) +
                         a14 * (a21 * a40 - a20 * a41)) +
                  a04 * (a10 * (a21 * a43 - a23 * a41) +
                         a11 * (a23 * a40 - a20 * a43) +
                         a13 * (a20 * a41 - a21 * a40))),
        invdet * (a00 * (a11 * (a22 * a44 - a24 * a42) +
                         a12 * (a24 * a41 - a21 * a44) +
                         a14 * (a21 * a42 - a22 * a41)) +
                  a01 * (a10 * (a24 * a42 - a22 * a44) +
                         a12 * (a20 * a44 - a24 * a40) +
                         a14 * (a22 * a40 - a20 * a42)) +
                  a02 * (a10 * (a21 * a44 - a24 * a41) +
                         a11 * (a24 * a40 - a20 * a44) +
                         a14 * (a20 * a41 - a21 * a40)) +
                  a04 * (a10 * (a22 * a41 - a21 * a42) +
                         a11 * (a20 * a42 - a22 * a40) +
                         a12 * (a21 * a40 - a20 * a41))),
        invdet * (a00 * (a11 * (a23 * a42 - a22 * a43) +
                         a12 * (a21 * a43 - a23 * a41) +
                         a13 * (a22 * a41 - a21 * a42)) +
                  a01 * (a10 * (a22 * a43 - a23 * a42) +
                         a12 * (a23 * a40 - a20 * a43) +
                         a13 * (a20 * a42 - a22 * a40)) +
                  a02 * (a10 * (a23 * a41 - a21 * a43) +
                         a11 * (a20 * a43 - a23 * a40) +
                         a13 * (a21 * a40 - a20 * a41)) +
                  a03 * (a10 * (a21 * a42 - a22 * a41) +
                         a11 * (a22 * a40 - a20 * a42) +
                         a12 * (a20 * a41 - a21 * a40))),
        invdet * (a01 * (a12 * (a23 * a34 - a24 * a33) +
                         a13 * (a24 * a32 - a22 * a34) +
                         a14 * (a22 * a33 - a23 * a32)) +
                  a02 * (a11 * (a24 * a33 - a23 * a34) +
                         a13 * (a21 * a34 - a24 * a31) +
                         a14 * (a23 * a31 - a21 * a33)) +
                  a03 * (a11 * (a22 * a34 - a24 * a32) +
                         a12 * (a24 * a31 - a21 * a34) +
                         a14 * (a21 * a32 - a22 * a31)) +
                  a04 * (a11 * (a23 * a32 - a22 * a33) +
                         a12 * (a21 * a33 - a23 * a31) +
                         a13 * (a22 * a31 - a21 * a32))),
        invdet * (a00 * (a12 * (a24 * a33 - a23 * a34) +
                         a13 * (a22 * a34 - a24 * a32) +
                         a14 * (a23 * a32 - a22 * a33)) +
                  a02 * (a10 * (a23 * a34 - a24 * a33) +
                         a13 * (a24 * a30 - a20 * a34) +
                         a14 * (a20 * a33 - a23 * a30)) +
                  a03 * (a10 * (a24 * a32 - a22 * a34) +
                         a12 * (a20 * a34 - a24 * a30) +
                         a14 * (a22 * a30 - a20 * a32)) +
                  a04 * (a10 * (a22 * a33 - a23 * a32) +
                         a12 * (a23 * a30 - a20 * a33) +
                         a13 * (a20 * a32 - a22 * a30))),
        invdet * (a00 * (a11 * (a23 * a34 - a24 * a33) +
                         a13 * (a24 * a31 - a21 * a34) +
                         a14 * (a21 * a33 - a23 * a31)) +
                  a01 * (a10 * (a24 * a33 - a23 * a34) +
                         a13 * (a20 * a34 - a24 * a30) +
                         a14 * (a23 * a30 - a20 * a33)) +
                  a03 * (a10 * (a21 * a34 - a24 * a31) +
                         a11 * (a24 * a30 - a20 * a34) +
                         a14 * (a20 * a31 - a21 * a30)) +
                  a04 * (a10 * (a23 * a31 - a21 * a33) +
                         a11 * (a20 * a33 - a23 * a30) +
                         a13 * (a21 * a30 - a20 * a31))),
        invdet * (a00 * (a11 * (a24 * a32 - a22 * a34) +
                         a12 * (a21 * a34 - a24 * a31) +
                         a14 * (a22 * a31 - a21 * a32)) +
                  a01 * (a10 * (a22 * a34 - a24 * a32) +
                         a12 * (a24 * a30 - a20 * a34) +
                         a14 * (a20 * a32 - a22 * a30)) +
                  a02 * (a10 * (a24 * a31 - a21 * a34) +
                         a11 * (a20 * a34 - a24 * a30) +
                         a14 * (a21 * a30 - a20 * a31)) +
                  a04 * (a10 * (a21 * a32 - a22 * a31) +
                         a11 * (a22 * a30 - a20 * a32) +
                         a12 * (a20 * a31 - a21 * a30))),
        invdet * (a00 * (a11 * (a22 * a33 - a23 * a32) +
                         a12 * (a23 * a31 - a21 * a33) +
                         a13 * (a21 * a32 - a22 * a31)) +
                  a01 * (a10 * (a23 * a32 - a22 * a33) +
                         a12 * (a20 * a33 - a23 * a30) +
                         a13 * (a22 * a30 - a20 * a32)) +
                  a02 * (a10 * (a21 * a33 - a23 * a31) +
                         a11 * (a23 * a30 - a20 * a33) +
                         a13 * (a20 * a31 - a21 * a30)) +
                  a03 * (a10 * (a22 * a31 - a21 * a32) +
                         a11 * (a20 * a32 - a22 * a30) +
                         a12 * (a21 * a30 - a20 * a31))) };
}


template<typename Matrix6x6Type>
constexpr Matrix6x6Type::value_type det6x6(const Matrix6x6Type &A)
{
    //
    // Returns the determinant of the input 6 x 6 matrix,
    // which should be specified in column-major order.
    //

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];
    const auto &a30 = A[3];
    const auto &a40 = A[4];
    const auto &a50 = A[5];

    const auto &a01 = A[6];
    const auto &a11 = A[7];
    const auto &a21 = A[8];
    const auto &a31 = A[9];
    const auto &a41 = A[10];
    const auto &a51 = A[11];

    const auto &a02 = A[12];
    const auto &a12 = A[13];
    const auto &a22 = A[14];
    const auto &a32 = A[15];
    const auto &a42 = A[16];
    const auto &a52 = A[17];

    const auto &a03 = A[18];
    const auto &a13 = A[19];
    const auto &a23 = A[20];
    const auto &a33 = A[21];
    const auto &a43 = A[22];
    const auto &a53 = A[23];

    const auto &a04 = A[24];
    const auto &a14 = A[25];
    const auto &a24 = A[26];
    const auto &a34 = A[27];
    const auto &a44 = A[28];
    const auto &a54 = A[29];

    const auto &a05 = A[30];
    const auto &a15 = A[31];
    const auto &a25 = A[32];
    const auto &a35 = A[33];
    const auto &a45 = A[34];
    const auto &a55 = A[35];

    return
        a00 * (a11 * (a22 * (a33 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a53 - a43 * a55) +
                             a35 * (a43 * a54 - a44 * a53)) +
                      a23 * (a32 * (a45 * a54 - a44 * a55) +
                             a34 * (a42 * a55 - a45 * a52) +
                             a35 * (a44 * a52 - a42 * a54)) +
                      a24 * (a32 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a52 - a42 * a55) +
                             a35 * (a42 * a53 - a43 * a52)) +
                      a25 * (a32 * (a44 * a53 - a43 * a54) +
                             a33 * (a42 * a54 - a44 * a52) +
                             a34 * (a43 * a52 - a42 * a53))) +
               a12 * (a21 * (a33 * (a45 * a54 - a44 * a55) +
                             a34 * (a43 * a55 - a45 * a53) +
                             a35 * (a44 * a53 - a43 * a54)) +
                      a23 * (a31 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a54 - a44 * a51)) +
                      a24 * (a31 * (a45 * a53 - a43 * a55) +
                             a33 * (a41 * a55 - a45 * a51) +
                             a35 * (a43 * a51 - a41 * a53)) +
                      a25 * (a31 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a51 - a41 * a54) +
                             a34 * (a41 * a53 - a43 * a51))) +
               a13 * (a21 * (a32 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a52 - a42 * a55) +
                             a35 * (a42 * a54 - a44 * a52)) +
                      a22 * (a31 * (a45 * a54 - a44 * a55) +
                             a34 * (a41 * a55 - a45 * a51) +
                             a35 * (a44 * a51 - a41 * a54)) +
                      a24 * (a31 * (a42 * a55 - a45 * a52) +
                             a32 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a52 - a42 * a51)) +
                      a25 * (a31 * (a44 * a52 - a42 * a54) +
                             a32 * (a41 * a54 - a44 * a51) +
                             a34 * (a42 * a51 - a41 * a52))) +
               a14 * (a21 * (a32 * (a45 * a53 - a43 * a55) +
                             a33 * (a42 * a55 - a45 * a52) +
                             a35 * (a43 * a52 - a42 * a53)) +
                      a22 * (a31 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a53 - a43 * a51)) +
                      a23 * (a31 * (a45 * a52 - a42 * a55) +
                             a32 * (a41 * a55 - a45 * a51) +
                             a35 * (a42 * a51 - a41 * a52)) +
                      a25 * (a31 * (a42 * a53 - a43 * a52) +
                             a32 * (a43 * a51 - a41 * a53) +
                             a33 * (a41 * a52 - a42 * a51))) +
               a15 * (a21 * (a32 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a52 - a42 * a54) +
                             a34 * (a42 * a53 - a43 * a52)) +
                      a22 * (a31 * (a44 * a53 - a43 * a54) +
                             a33 * (a41 * a54 - a44 * a51) +
                             a34 * (a43 * a51 - a41 * a53)) +
                      a23 * (a31 * (a42 * a54 - a44 * a52) +
                             a32 * (a44 * a51 - a41 * a54) +
                             a34 * (a41 * a52 - a42 * a51)) +
                      a24 * (a31 * (a43 * a52 - a42 * a53) +
                             a32 * (a41 * a53 - a43 * a51) +
                             a33 * (a42 * a51 - a41 * a52)))) +
        a01 * (a10 * (a22 * (a33 * (a45 * a54 - a44 * a55) +
                             a34 * (a43 * a55 - a45 * a53) +
                             a35 * (a44 * a53 - a43 * a54)) +
                      a23 * (a32 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a52 - a42 * a55) +
                             a35 * (a42 * a54 - a44 * a52)) +
                      a24 * (a32 * (a45 * a53 - a43 * a55) +
                             a33 * (a42 * a55 - a45 * a52) +
                             a35 * (a43 * a52 - a42 * a53)) +
                      a25 * (a32 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a52 - a42 * a54) +
                             a34 * (a42 * a53 - a43 * a52))) +
               a12 * (a20 * (a33 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a53 - a43 * a55) +
                             a35 * (a43 * a54 - a44 * a53)) +
                      a23 * (a30 * (a45 * a54 - a44 * a55) +
                             a34 * (a40 * a55 - a45 * a50) +
                             a35 * (a44 * a50 - a40 * a54)) +
                      a24 * (a30 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a53 - a43 * a50)) +
                      a25 * (a30 * (a44 * a53 - a43 * a54) +
                             a33 * (a40 * a54 - a44 * a50) +
                             a34 * (a43 * a50 - a40 * a53))) +
               a13 * (a20 * (a32 * (a45 * a54 - a44 * a55) +
                             a34 * (a42 * a55 - a45 * a52) +
                             a35 * (a44 * a52 - a42 * a54)) +
                      a22 * (a30 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a54 - a44 * a50)) +
                      a24 * (a30 * (a45 * a52 - a42 * a55) +
                             a32 * (a40 * a55 - a45 * a50) +
                             a35 * (a42 * a50 - a40 * a52)) +
                      a25 * (a30 * (a42 * a54 - a44 * a52) +
                             a32 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a52 - a42 * a50))) +
               a14 * (a20 * (a32 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a52 - a42 * a55) +
                             a35 * (a42 * a53 - a43 * a52)) +
                      a22 * (a30 * (a45 * a53 - a43 * a55) +
                             a33 * (a40 * a55 - a45 * a50) +
                             a35 * (a43 * a50 - a40 * a53)) +
                      a23 * (a30 * (a42 * a55 - a45 * a52) +
                             a32 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a52 - a42 * a50)) +
                      a25 * (a30 * (a43 * a52 - a42 * a53) +
                             a32 * (a40 * a53 - a43 * a50) +
                             a33 * (a42 * a50 - a40 * a52))) +
               a15 * (a20 * (a32 * (a44 * a53 - a43 * a54) +
                             a33 * (a42 * a54 - a44 * a52) +
                             a34 * (a43 * a52 - a42 * a53)) +
                      a22 * (a30 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a53 - a43 * a50)) +
                      a23 * (a30 * (a44 * a52 - a42 * a54) +
                             a32 * (a40 * a54 - a44 * a50) +
                             a34 * (a42 * a50 - a40 * a52)) +
                      a24 * (a30 * (a42 * a53 - a43 * a52) +
                             a32 * (a43 * a50 - a40 * a53) +
                             a33 * (a40 * a52 - a42 * a50)))) +
        a02 * (a10 * (a21 * (a33 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a53 - a43 * a55) +
                             a35 * (a43 * a54 - a44 * a53)) +
                      a23 * (a31 * (a45 * a54 - a44 * a55) +
                             a34 * (a41 * a55 - a45 * a51) +
                             a35 * (a44 * a51 - a41 * a54)) +
                      a24 * (a31 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a53 - a43 * a51)) +
                      a25 * (a31 * (a44 * a53 - a43 * a54) +
                             a33 * (a41 * a54 - a44 * a51) +
                             a34 * (a43 * a51 - a41 * a53))) +
               a11 * (a20 * (a33 * (a45 * a54 - a44 * a55) +
                             a34 * (a43 * a55 - a45 * a53) +
                             a35 * (a44 * a53 - a43 * a54)) +
                      a23 * (a30 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a54 - a44 * a50)) +
                      a24 * (a30 * (a45 * a53 - a43 * a55) +
                             a33 * (a40 * a55 - a45 * a50) +
                             a35 * (a43 * a50 - a40 * a53)) +
                      a25 * (a30 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a53 - a43 * a50))) +
               a13 * (a20 * (a31 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a54 - a44 * a51)) +
                      a21 * (a30 * (a45 * a54 - a44 * a55) +
                             a34 * (a40 * a55 - a45 * a50) +
                             a35 * (a44 * a50 - a40 * a54)) +
                      a24 * (a30 * (a41 * a55 - a45 * a51) +
                             a31 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a51 - a41 * a50)) +
                      a25 * (a30 * (a44 * a51 - a41 * a54) +
                             a31 * (a40 * a54 - a44 * a50) +
                             a34 * (a41 * a50 - a40 * a51))) +
               a14 * (a20 * (a31 * (a45 * a53 - a43 * a55) +
                             a33 * (a41 * a55 - a45 * a51) +
                             a35 * (a43 * a51 - a41 * a53)) +
                      a21 * (a30 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a53 - a43 * a50)) +
                      a23 * (a30 * (a45 * a51 - a41 * a55) +
                             a31 * (a40 * a55 - a45 * a50) +
                             a35 * (a41 * a50 - a40 * a51)) +
                      a25 * (a30 * (a41 * a53 - a43 * a51) +
                             a31 * (a43 * a50 - a40 * a53) +
                             a33 * (a40 * a51 - a41 * a50))) +
               a15 * (a20 * (a31 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a51 - a41 * a54) +
                             a34 * (a41 * a53 - a43 * a51)) +
                      a21 * (a30 * (a44 * a53 - a43 * a54) +
                             a33 * (a40 * a54 - a44 * a50) +
                             a34 * (a43 * a50 - a40 * a53)) +
                      a23 * (a30 * (a41 * a54 - a44 * a51) +
                             a31 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a51 - a41 * a50)) +
                      a24 * (a30 * (a43 * a51 - a41 * a53) +
                             a31 * (a40 * a53 - a43 * a50) +
                             a33 * (a41 * a50 - a40 * a51)))) +
        a03 * (a10 * (a21 * (a32 * (a45 * a54 - a44 * a55) +
                             a34 * (a42 * a55 - a45 * a52) +
                             a35 * (a44 * a52 - a42 * a54)) +
                      a22 * (a31 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a54 - a44 * a51)) +
                      a24 * (a31 * (a45 * a52 - a42 * a55) +
                             a32 * (a41 * a55 - a45 * a51) +
                             a35 * (a42 * a51 - a41 * a52)) +
                      a25 * (a31 * (a42 * a54 - a44 * a52) +
                             a32 * (a44 * a51 - a41 * a54) +
                             a34 * (a41 * a52 - a42 * a51))) +
               a11 * (a20 * (a32 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a52 - a42 * a55) +
                             a35 * (a42 * a54 - a44 * a52)) +
                      a22 * (a30 * (a45 * a54 - a44 * a55) +
                             a34 * (a40 * a55 - a45 * a50) +
                             a35 * (a44 * a50 - a40 * a54)) +
                      a24 * (a30 * (a42 * a55 - a45 * a52) +
                             a32 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a52 - a42 * a50)) +
                      a25 * (a30 * (a44 * a52 - a42 * a54) +
                             a32 * (a40 * a54 - a44 * a50) +
                             a34 * (a42 * a50 - a40 * a52))) +
               a12 * (a20 * (a31 * (a45 * a54 - a44 * a55) +
                             a34 * (a41 * a55 - a45 * a51) +
                             a35 * (a44 * a51 - a41 * a54)) +
                      a21 * (a30 * (a44 * a55 - a45 * a54) +
                             a34 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a54 - a44 * a50)) +
                      a24 * (a30 * (a45 * a51 - a41 * a55) +
                             a31 * (a40 * a55 - a45 * a50) +
                             a35 * (a41 * a50 - a40 * a51)) +
                      a25 * (a30 * (a41 * a54 - a44 * a51) +
                             a31 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a51 - a41 * a50))) +
               a14 * (a20 * (a31 * (a42 * a55 - a45 * a52) +
                             a32 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a52 - a42 * a51)) +
                      a21 * (a30 * (a45 * a52 - a42 * a55) +
                             a32 * (a40 * a55 - a45 * a50) +
                             a35 * (a42 * a50 - a40 * a52)) +
                      a22 * (a30 * (a41 * a55 - a45 * a51) +
                             a31 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a51 - a41 * a50)) +
                      a25 * (a30 * (a42 * a51 - a41 * a52) +
                             a31 * (a40 * a52 - a42 * a50) +
                             a32 * (a41 * a50 - a40 * a51))) +
               a15 * (a20 * (a31 * (a44 * a52 - a42 * a54) +
                             a32 * (a41 * a54 - a44 * a51) +
                             a34 * (a42 * a51 - a41 * a52)) +
                      a21 * (a30 * (a42 * a54 - a44 * a52) +
                             a32 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a52 - a42 * a50)) +
                      a22 * (a30 * (a44 * a51 - a41 * a54) +
                             a31 * (a40 * a54 - a44 * a50) +
                             a34 * (a41 * a50 - a40 * a51)) +
                      a24 * (a30 * (a41 * a52 - a42 * a51) +
                             a31 * (a42 * a50 - a40 * a52) +
                             a32 * (a40 * a51 - a41 * a50)))) +
        a04 * (a10 * (a21 * (a32 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a52 - a42 * a55) +
                             a35 * (a42 * a53 - a43 * a52)) +
                      a22 * (a31 * (a45 * a53 - a43 * a55) +
                             a33 * (a41 * a55 - a45 * a51) +
                             a35 * (a43 * a51 - a41 * a53)) +
                      a23 * (a31 * (a42 * a55 - a45 * a52) +
                             a32 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a52 - a42 * a51)) +
                      a25 * (a31 * (a43 * a52 - a42 * a53) +
                             a32 * (a41 * a53 - a43 * a51) +
                             a33 * (a42 * a51 - a41 * a52))) +
               a11 * (a20 * (a32 * (a45 * a53 - a43 * a55) +
                             a33 * (a42 * a55 - a45 * a52) +
                             a35 * (a43 * a52 - a42 * a53)) +
                      a22 * (a30 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a53 - a43 * a50)) +
                      a23 * (a30 * (a45 * a52 - a42 * a55) +
                             a32 * (a40 * a55 - a45 * a50) +
                             a35 * (a42 * a50 - a40 * a52)) +
                      a25 * (a30 * (a42 * a53 - a43 * a52) +
                             a32 * (a43 * a50 - a40 * a53) +
                             a33 * (a40 * a52 - a42 * a50))) +
               a12 * (a20 * (a31 * (a43 * a55 - a45 * a53) +
                             a33 * (a45 * a51 - a41 * a55) +
                             a35 * (a41 * a53 - a43 * a51)) +
                      a21 * (a30 * (a45 * a53 - a43 * a55) +
                             a33 * (a40 * a55 - a45 * a50) +
                             a35 * (a43 * a50 - a40 * a53)) +
                      a23 * (a30 * (a41 * a55 - a45 * a51) +
                             a31 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a51 - a41 * a50)) +
                      a25 * (a30 * (a43 * a51 - a41 * a53) +
                             a31 * (a40 * a53 - a43 * a50) +
                             a33 * (a41 * a50 - a40 * a51))) +
               a13 * (a20 * (a31 * (a45 * a52 - a42 * a55) +
                             a32 * (a41 * a55 - a45 * a51) +
                             a35 * (a42 * a51 - a41 * a52)) +
                      a21 * (a30 * (a42 * a55 - a45 * a52) +
                             a32 * (a45 * a50 - a40 * a55) +
                             a35 * (a40 * a52 - a42 * a50)) +
                      a22 * (a30 * (a45 * a51 - a41 * a55) +
                             a31 * (a40 * a55 - a45 * a50) +
                             a35 * (a41 * a50 - a40 * a51)) +
                      a25 * (a30 * (a41 * a52 - a42 * a51) +
                             a31 * (a42 * a50 - a40 * a52) +
                             a32 * (a40 * a51 - a41 * a50))) +
               a15 * (a20 * (a31 * (a42 * a53 - a43 * a52) +
                             a32 * (a43 * a51 - a41 * a53) +
                             a33 * (a41 * a52 - a42 * a51)) +
                      a21 * (a30 * (a43 * a52 - a42 * a53) +
                             a32 * (a40 * a53 - a43 * a50) +
                             a33 * (a42 * a50 - a40 * a52)) +
                      a22 * (a30 * (a41 * a53 - a43 * a51) +
                             a31 * (a43 * a50 - a40 * a53) +
                             a33 * (a40 * a51 - a41 * a50)) +
                      a23 * (a30 * (a42 * a51 - a41 * a52) +
                             a31 * (a40 * a52 - a42 * a50) +
                             a32 * (a41 * a50 - a40 * a51)))) +
        a05 * (a10 * (a21 * (a32 * (a44 * a53 - a43 * a54) +
                             a33 * (a42 * a54 - a44 * a52) +
                             a34 * (a43 * a52 - a42 * a53)) +
                      a22 * (a31 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a51 - a41 * a54) +
                             a34 * (a41 * a53 - a43 * a51)) +
                      a23 * (a31 * (a44 * a52 - a42 * a54) +
                             a32 * (a41 * a54 - a44 * a51) +
                             a34 * (a42 * a51 - a41 * a52)) +
                      a24 * (a31 * (a42 * a53 - a43 * a52) +
                             a32 * (a43 * a51 - a41 * a53) +
                             a33 * (a41 * a52 - a42 * a51))) +
               a11 * (a20 * (a32 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a52 - a42 * a54) +
                             a34 * (a42 * a53 - a43 * a52)) +
                      a22 * (a30 * (a44 * a53 - a43 * a54) +
                             a33 * (a40 * a54 - a44 * a50) +
                             a34 * (a43 * a50 - a40 * a53)) +
                      a23 * (a30 * (a42 * a54 - a44 * a52) +
                             a32 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a52 - a42 * a50)) +
                      a24 * (a30 * (a43 * a52 - a42 * a53) +
                             a32 * (a40 * a53 - a43 * a50) +
                             a33 * (a42 * a50 - a40 * a52))) +
               a12 * (a20 * (a31 * (a44 * a53 - a43 * a54) +
                             a33 * (a41 * a54 - a44 * a51) +
                             a34 * (a43 * a51 - a41 * a53)) +
                      a21 * (a30 * (a43 * a54 - a44 * a53) +
                             a33 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a53 - a43 * a50)) +
                      a23 * (a30 * (a44 * a51 - a41 * a54) +
                             a31 * (a40 * a54 - a44 * a50) +
                             a34 * (a41 * a50 - a40 * a51)) +
                      a24 * (a30 * (a41 * a53 - a43 * a51) +
                             a31 * (a43 * a50 - a40 * a53) +
                             a33 * (a40 * a51 - a41 * a50))) +
               a13 * (a20 * (a31 * (a42 * a54 - a44 * a52) +
                             a32 * (a44 * a51 - a41 * a54) +
                             a34 * (a41 * a52 - a42 * a51)) +
                      a21 * (a30 * (a44 * a52 - a42 * a54) +
                             a32 * (a40 * a54 - a44 * a50) +
                             a34 * (a42 * a50 - a40 * a52)) +
                      a22 * (a30 * (a41 * a54 - a44 * a51) +
                             a31 * (a44 * a50 - a40 * a54) +
                             a34 * (a40 * a51 - a41 * a50)) +
                      a24 * (a30 * (a42 * a51 - a41 * a52) +
                             a31 * (a40 * a52 - a42 * a50) +
                             a32 * (a41 * a50 - a40 * a51))) +
               a14 * (a20 * (a31 * (a43 * a52 - a42 * a53) +
                             a32 * (a41 * a53 - a43 * a51) +
                             a33 * (a42 * a51 - a41 * a52)) +
                      a21 * (a30 * (a42 * a53 - a43 * a52) +
                             a32 * (a43 * a50 - a40 * a53) +
                             a33 * (a40 * a52 - a42 * a50)) +
                      a22 * (a30 * (a43 * a51 - a41 * a53) +
                             a31 * (a40 * a53 - a43 * a50) +
                             a33 * (a41 * a50 - a40 * a51)) +
                      a23 * (a30 * (a41 * a52 - a42 * a51) +
                             a31 * (a42 * a50 - a40 * a52) +
                             a32 * (a40 * a51 - a41 * a50))));
}


template<typename Matrix6x6Type>
constexpr Matrix6x6Type inv6x6(const Matrix6x6Type &A)
{
    //
    // Returns the inverse of the input 6 x 6 matrix,
    // both the input matrix and the resulting inverse
    // should be specified in column-major order.
    //

    const auto invdet = typename Matrix6x6Type::value_type{ 1 } /
        det6x6(A);

    const auto &a00 = A[0];
    const auto &a10 = A[1];
    const auto &a20 = A[2];
    const auto &a30 = A[3];
    const auto &a40 = A[4];
    const auto &a50 = A[5];

    const auto &a01 = A[6];
    const auto &a11 = A[7];
    const auto &a21 = A[8];
    const auto &a31 = A[9];
    const auto &a41 = A[10];
    const auto &a51 = A[11];

    const auto &a02 = A[12];
    const auto &a12 = A[13];
    const auto &a22 = A[14];
    const auto &a32 = A[15];
    const auto &a42 = A[16];
    const auto &a52 = A[17];

    const auto &a03 = A[18];
    const auto &a13 = A[19];
    const auto &a23 = A[20];
    const auto &a33 = A[21];
    const auto &a43 = A[22];
    const auto &a53 = A[23];

    const auto &a04 = A[24];
    const auto &a14 = A[25];
    const auto &a24 = A[26];
    const auto &a34 = A[27];
    const auto &a44 = A[28];
    const auto &a54 = A[29];

    const auto &a05 = A[30];
    const auto &a15 = A[31];
    const auto &a25 = A[32];
    const auto &a35 = A[33];
    const auto &a45 = A[34];
    const auto &a55 = A[35];

    return Matrix6x6Type{
        invdet * (a11 * (a22 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a23 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a24 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a25 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53))) +
                  a12 * (a21 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a23 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a24 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a25 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51))) +
                  a13 * (a21 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a22 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a24 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a25 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52))) +
                  a14 * (a21 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a22 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a23 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a25 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51))) +
                  a15 * (a21 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52)) +
                         a22 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53)) +
                         a23 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51)) +
                         a24 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52)))),
        invdet * (a10 * (a22 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a23 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a24 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a25 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52))) +
                  a12 * (a20 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a23 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a24 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a25 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53))) +
                  a13 * (a20 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a22 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a24 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a25 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50))) +
                  a14 * (a20 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a22 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a25 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52))) +
                  a15 * (a20 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53)) +
                         a22 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52)) +
                         a24 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50)))),
        invdet * (a10 * (a21 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a23 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a24 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a25 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53))) +
                  a11 * (a20 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a23 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a24 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a25 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50))) +
                  a13 * (a20 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a21 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a24 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a25 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51))) +
                  a14 * (a20 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a21 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a25 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50))) +
                  a15 * (a20 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51)) +
                         a21 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50)) +
                         a24 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51)))),
        invdet * (a10 * (a21 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a22 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a24 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a25 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51))) +
                  a11 * (a20 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a22 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a24 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a25 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52))) +
                  a12 * (a20 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a21 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a24 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a25 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50))) +
                  a14 * (a20 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a21 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a22 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a25 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51))) +
                  a15 * (a20 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52)) +
                         a21 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50)) +
                         a22 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51)) +
                         a24 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50)))),
        invdet * (a10 * (a21 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a22 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a23 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a25 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52))) +
                  a11 * (a20 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a22 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a25 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50))) +
                  a12 * (a20 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a21 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a25 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51))) +
                  a13 * (a20 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a21 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a22 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a25 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50))) +
                  a15 * (a20 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51)) +
                         a21 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52)) +
                         a22 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50)) +
                         a23 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51)))),
        invdet * (a10 * (a21 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53)) +
                         a22 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51)) +
                         a23 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52)) +
                         a24 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51))) +
                  a11 * (a20 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52)) +
                         a22 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50)) +
                         a24 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52))) +
                  a12 * (a20 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53)) +
                         a21 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51)) +
                         a24 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50))) +
                  a13 * (a20 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51)) +
                         a21 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52)) +
                         a22 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50)) +
                         a24 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51))) +
                  a14 * (a20 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52)) +
                         a21 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50)) +
                         a22 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51)) +
                         a23 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50)))),
        invdet * (a01 * (a22 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a23 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a24 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a25 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52))) +
                  a02 * (a21 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a23 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a24 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a25 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53))) +
                  a03 * (a21 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a22 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a24 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a25 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51))) +
                  a04 * (a21 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a22 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a23 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a25 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52))) +
                  a05 * (a21 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53)) +
                         a22 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51)) +
                         a23 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52)) +
                         a24 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51)))),
        invdet * (a00 * (a22 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a23 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a24 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a25 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53))) +
                  a02 * (a20 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a23 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a24 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a25 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50))) +
                  a03 * (a20 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a22 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a24 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a25 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52))) +
                  a04 * (a20 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a22 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a25 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50))) +
                  a05 * (a20 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52)) +
                         a22 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50)) +
                         a24 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52)))),
        invdet * (a00 * (a21 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a23 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a24 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a25 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51))) +
                  a01 * (a20 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a23 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a24 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a25 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53))) +
                  a03 * (a20 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a21 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a24 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a25 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50))) +
                  a04 * (a20 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a21 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a25 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51))) +
                  a05 * (a20 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53)) +
                         a21 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51)) +
                         a24 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50)))),
        invdet * (a00 * (a21 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a22 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a24 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a25 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52))) +
                  a01 * (a20 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a22 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a24 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a25 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50))) +
                  a02 * (a20 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a21 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a24 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a25 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51))) +
                  a04 * (a20 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a21 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a22 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a25 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50))) +
                  a05 * (a20 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51)) +
                         a21 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52)) +
                         a22 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50)) +
                         a24 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51)))),
        invdet * (a00 * (a21 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a22 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a23 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a25 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51))) +
                  a01 * (a20 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a22 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a25 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52))) +
                  a02 * (a20 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a21 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a25 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50))) +
                  a03 * (a20 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a21 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a22 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a25 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51))) +
                  a05 * (a20 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52)) +
                         a21 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50)) +
                         a22 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51)) +
                         a23 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50)))),
        invdet * (a00 * (a21 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52)) +
                         a22 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53)) +
                         a23 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51)) +
                         a24 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52))) +
                  a01 * (a20 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53)) +
                         a22 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50)) +
                         a23 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52)) +
                         a24 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50))) +
                  a02 * (a20 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51)) +
                         a21 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53)) +
                         a23 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50)) +
                         a24 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51))) +
                  a03 * (a20 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52)) +
                         a21 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50)) +
                         a22 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51)) +
                         a24 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50))) +
                  a04 * (a20 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51)) +
                         a21 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52)) +
                         a22 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50)) +
                         a23 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51)))),
        invdet * (a01 * (a12 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a13 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a14 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a15 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53))) +
                  a02 * (a11 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a13 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a14 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a15 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51))) +
                  a03 * (a11 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a12 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a14 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a15 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52))) +
                  a04 * (a11 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a12 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a13 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a15 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51))) +
                  a05 * (a11 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52)) +
                         a12 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53)) +
                         a13 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51)) +
                         a14 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52)))),
        invdet * (a00 * (a12 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a13 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a14 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a15 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52))) +
                  a02 * (a10 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a13 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a14 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a15 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53))) +
                  a03 * (a10 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a12 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a14 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a15 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50))) +
                  a04 * (a10 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a12 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a13 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a15 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52))) +
                  a05 * (a10 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53)) +
                         a12 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50)) +
                         a13 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52)) +
                         a14 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50)))),
        invdet * (a00 * (a11 * (a33 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a53 - a43 * a55) +
                                a35 * (a43 * a54 - a44 * a53)) +
                         a13 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a14 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a15 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53))) +
                  a01 * (a10 * (a33 * (a45 * a54 - a44 * a55) +
                                a34 * (a43 * a55 - a45 * a53) +
                                a35 * (a44 * a53 - a43 * a54)) +
                         a13 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a14 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a15 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50))) +
                  a03 * (a10 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a11 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a14 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a15 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51))) +
                  a04 * (a10 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a11 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a13 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a15 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50))) +
                  a05 * (a10 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51)) +
                         a11 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53)) +
                         a13 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50)) +
                         a14 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51)))),
        invdet * (a00 * (a11 * (a32 * (a45 * a54 - a44 * a55) +
                                a34 * (a42 * a55 - a45 * a52) +
                                a35 * (a44 * a52 - a42 * a54)) +
                         a12 * (a31 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a54 - a44 * a51)) +
                         a14 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a15 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51))) +
                  a01 * (a10 * (a32 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a54 - a44 * a52)) +
                         a12 * (a30 * (a45 * a54 - a44 * a55) +
                                a34 * (a40 * a55 - a45 * a50) +
                                a35 * (a44 * a50 - a40 * a54)) +
                         a14 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a15 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52))) +
                  a02 * (a10 * (a31 * (a45 * a54 - a44 * a55) +
                                a34 * (a41 * a55 - a45 * a51) +
                                a35 * (a44 * a51 - a41 * a54)) +
                         a11 * (a30 * (a44 * a55 - a45 * a54) +
                                a34 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a54 - a44 * a50)) +
                         a14 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a15 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50))) +
                  a04 * (a10 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a11 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a12 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a15 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51))) +
                  a05 * (a10 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52)) +
                         a11 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50)) +
                         a12 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51)) +
                         a14 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50)))),
        invdet * (a00 * (a11 * (a32 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a52 - a42 * a55) +
                                a35 * (a42 * a53 - a43 * a52)) +
                         a12 * (a31 * (a45 * a53 - a43 * a55) +
                                a33 * (a41 * a55 - a45 * a51) +
                                a35 * (a43 * a51 - a41 * a53)) +
                         a13 * (a31 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a52 - a42 * a51)) +
                         a15 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52))) +
                  a01 * (a10 * (a32 * (a45 * a53 - a43 * a55) +
                                a33 * (a42 * a55 - a45 * a52) +
                                a35 * (a43 * a52 - a42 * a53)) +
                         a12 * (a30 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a53 - a43 * a50)) +
                         a13 * (a30 * (a45 * a52 - a42 * a55) +
                                a32 * (a40 * a55 - a45 * a50) +
                                a35 * (a42 * a50 - a40 * a52)) +
                         a15 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50))) +
                  a02 * (a10 * (a31 * (a43 * a55 - a45 * a53) +
                                a33 * (a45 * a51 - a41 * a55) +
                                a35 * (a41 * a53 - a43 * a51)) +
                         a11 * (a30 * (a45 * a53 - a43 * a55) +
                                a33 * (a40 * a55 - a45 * a50) +
                                a35 * (a43 * a50 - a40 * a53)) +
                         a13 * (a30 * (a41 * a55 - a45 * a51) +
                                a31 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a51 - a41 * a50)) +
                         a15 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51))) +
                  a03 * (a10 * (a31 * (a45 * a52 - a42 * a55) +
                                a32 * (a41 * a55 - a45 * a51) +
                                a35 * (a42 * a51 - a41 * a52)) +
                         a11 * (a30 * (a42 * a55 - a45 * a52) +
                                a32 * (a45 * a50 - a40 * a55) +
                                a35 * (a40 * a52 - a42 * a50)) +
                         a12 * (a30 * (a45 * a51 - a41 * a55) +
                                a31 * (a40 * a55 - a45 * a50) +
                                a35 * (a41 * a50 - a40 * a51)) +
                         a15 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50))) +
                  a05 * (a10 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51)) +
                         a11 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52)) +
                         a12 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50)) +
                         a13 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51)))),
        invdet * (a00 * (a11 * (a32 * (a44 * a53 - a43 * a54) +
                                a33 * (a42 * a54 - a44 * a52) +
                                a34 * (a43 * a52 - a42 * a53)) +
                         a12 * (a31 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a53 - a43 * a51)) +
                         a13 * (a31 * (a44 * a52 - a42 * a54) +
                                a32 * (a41 * a54 - a44 * a51) +
                                a34 * (a42 * a51 - a41 * a52)) +
                         a14 * (a31 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a51 - a41 * a53) +
                                a33 * (a41 * a52 - a42 * a51))) +
                  a01 * (a10 * (a32 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a52 - a42 * a54) +
                                a34 * (a42 * a53 - a43 * a52)) +
                         a12 * (a30 * (a44 * a53 - a43 * a54) +
                                a33 * (a40 * a54 - a44 * a50) +
                                a34 * (a43 * a50 - a40 * a53)) +
                         a13 * (a30 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a52 - a42 * a50)) +
                         a14 * (a30 * (a43 * a52 - a42 * a53) +
                                a32 * (a40 * a53 - a43 * a50) +
                                a33 * (a42 * a50 - a40 * a52))) +
                  a02 * (a10 * (a31 * (a44 * a53 - a43 * a54) +
                                a33 * (a41 * a54 - a44 * a51) +
                                a34 * (a43 * a51 - a41 * a53)) +
                         a11 * (a30 * (a43 * a54 - a44 * a53) +
                                a33 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a53 - a43 * a50)) +
                         a13 * (a30 * (a44 * a51 - a41 * a54) +
                                a31 * (a40 * a54 - a44 * a50) +
                                a34 * (a41 * a50 - a40 * a51)) +
                         a14 * (a30 * (a41 * a53 - a43 * a51) +
                                a31 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a51 - a41 * a50))) +
                  a03 * (a10 * (a31 * (a42 * a54 - a44 * a52) +
                                a32 * (a44 * a51 - a41 * a54) +
                                a34 * (a41 * a52 - a42 * a51)) +
                         a11 * (a30 * (a44 * a52 - a42 * a54) +
                                a32 * (a40 * a54 - a44 * a50) +
                                a34 * (a42 * a50 - a40 * a52)) +
                         a12 * (a30 * (a41 * a54 - a44 * a51) +
                                a31 * (a44 * a50 - a40 * a54) +
                                a34 * (a40 * a51 - a41 * a50)) +
                         a14 * (a30 * (a42 * a51 - a41 * a52) +
                                a31 * (a40 * a52 - a42 * a50) +
                                a32 * (a41 * a50 - a40 * a51))) +
                  a04 * (a10 * (a31 * (a43 * a52 - a42 * a53) +
                                a32 * (a41 * a53 - a43 * a51) +
                                a33 * (a42 * a51 - a41 * a52)) +
                         a11 * (a30 * (a42 * a53 - a43 * a52) +
                                a32 * (a43 * a50 - a40 * a53) +
                                a33 * (a40 * a52 - a42 * a50)) +
                         a12 * (a30 * (a43 * a51 - a41 * a53) +
                                a31 * (a40 * a53 - a43 * a50) +
                                a33 * (a41 * a50 - a40 * a51)) +
                         a13 * (a30 * (a41 * a52 - a42 * a51) +
                                a31 * (a42 * a50 - a40 * a52) +
                                a32 * (a40 * a51 - a41 * a50)))),
        invdet * (a01 * (a12 * (a23 * (a45 * a54 - a44 * a55) +
                                a24 * (a43 * a55 - a45 * a53) +
                                a25 * (a44 * a53 - a43 * a54)) +
                         a13 * (a22 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a52 - a42 * a55) +
                                a25 * (a42 * a54 - a44 * a52)) +
                         a14 * (a22 * (a45 * a53 - a43 * a55) +
                                a23 * (a42 * a55 - a45 * a52) +
                                a25 * (a43 * a52 - a42 * a53)) +
                         a15 * (a22 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a52 - a42 * a54) +
                                a24 * (a42 * a53 - a43 * a52))) +
                  a02 * (a11 * (a23 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a53 - a43 * a55) +
                                a25 * (a43 * a54 - a44 * a53)) +
                         a13 * (a21 * (a45 * a54 - a44 * a55) +
                                a24 * (a41 * a55 - a45 * a51) +
                                a25 * (a44 * a51 - a41 * a54)) +
                         a14 * (a21 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a53 - a43 * a51)) +
                         a15 * (a21 * (a44 * a53 - a43 * a54) +
                                a23 * (a41 * a54 - a44 * a51) +
                                a24 * (a43 * a51 - a41 * a53))) +
                  a03 * (a11 * (a22 * (a45 * a54 - a44 * a55) +
                                a24 * (a42 * a55 - a45 * a52) +
                                a25 * (a44 * a52 - a42 * a54)) +
                         a12 * (a21 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a54 - a44 * a51)) +
                         a14 * (a21 * (a45 * a52 - a42 * a55) +
                                a22 * (a41 * a55 - a45 * a51) +
                                a25 * (a42 * a51 - a41 * a52)) +
                         a15 * (a21 * (a42 * a54 - a44 * a52) +
                                a22 * (a44 * a51 - a41 * a54) +
                                a24 * (a41 * a52 - a42 * a51))) +
                  a04 * (a11 * (a22 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a52 - a42 * a55) +
                                a25 * (a42 * a53 - a43 * a52)) +
                         a12 * (a21 * (a45 * a53 - a43 * a55) +
                                a23 * (a41 * a55 - a45 * a51) +
                                a25 * (a43 * a51 - a41 * a53)) +
                         a13 * (a21 * (a42 * a55 - a45 * a52) +
                                a22 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a52 - a42 * a51)) +
                         a15 * (a21 * (a43 * a52 - a42 * a53) +
                                a22 * (a41 * a53 - a43 * a51) +
                                a23 * (a42 * a51 - a41 * a52))) +
                  a05 * (a11 * (a22 * (a44 * a53 - a43 * a54) +
                                a23 * (a42 * a54 - a44 * a52) +
                                a24 * (a43 * a52 - a42 * a53)) +
                         a12 * (a21 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a51 - a41 * a54) +
                                a24 * (a41 * a53 - a43 * a51)) +
                         a13 * (a21 * (a44 * a52 - a42 * a54) +
                                a22 * (a41 * a54 - a44 * a51) +
                                a24 * (a42 * a51 - a41 * a52)) +
                         a14 * (a21 * (a42 * a53 - a43 * a52) +
                                a22 * (a43 * a51 - a41 * a53) +
                                a23 * (a41 * a52 - a42 * a51)))),
        invdet * (a00 * (a12 * (a23 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a53 - a43 * a55) +
                                a25 * (a43 * a54 - a44 * a53)) +
                         a13 * (a22 * (a45 * a54 - a44 * a55) +
                                a24 * (a42 * a55 - a45 * a52) +
                                a25 * (a44 * a52 - a42 * a54)) +
                         a14 * (a22 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a52 - a42 * a55) +
                                a25 * (a42 * a53 - a43 * a52)) +
                         a15 * (a22 * (a44 * a53 - a43 * a54) +
                                a23 * (a42 * a54 - a44 * a52) +
                                a24 * (a43 * a52 - a42 * a53))) +
                  a02 * (a10 * (a23 * (a45 * a54 - a44 * a55) +
                                a24 * (a43 * a55 - a45 * a53) +
                                a25 * (a44 * a53 - a43 * a54)) +
                         a13 * (a20 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a54 - a44 * a50)) +
                         a14 * (a20 * (a45 * a53 - a43 * a55) +
                                a23 * (a40 * a55 - a45 * a50) +
                                a25 * (a43 * a50 - a40 * a53)) +
                         a15 * (a20 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a53 - a43 * a50))) +
                  a03 * (a10 * (a22 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a52 - a42 * a55) +
                                a25 * (a42 * a54 - a44 * a52)) +
                         a12 * (a20 * (a45 * a54 - a44 * a55) +
                                a24 * (a40 * a55 - a45 * a50) +
                                a25 * (a44 * a50 - a40 * a54)) +
                         a14 * (a20 * (a42 * a55 - a45 * a52) +
                                a22 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a52 - a42 * a50)) +
                         a15 * (a20 * (a44 * a52 - a42 * a54) +
                                a22 * (a40 * a54 - a44 * a50) +
                                a24 * (a42 * a50 - a40 * a52))) +
                  a04 * (a10 * (a22 * (a45 * a53 - a43 * a55) +
                                a23 * (a42 * a55 - a45 * a52) +
                                a25 * (a43 * a52 - a42 * a53)) +
                         a12 * (a20 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a53 - a43 * a50)) +
                         a13 * (a20 * (a45 * a52 - a42 * a55) +
                                a22 * (a40 * a55 - a45 * a50) +
                                a25 * (a42 * a50 - a40 * a52)) +
                         a15 * (a20 * (a42 * a53 - a43 * a52) +
                                a22 * (a43 * a50 - a40 * a53) +
                                a23 * (a40 * a52 - a42 * a50))) +
                  a05 * (a10 * (a22 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a52 - a42 * a54) +
                                a24 * (a42 * a53 - a43 * a52)) +
                         a12 * (a20 * (a44 * a53 - a43 * a54) +
                                a23 * (a40 * a54 - a44 * a50) +
                                a24 * (a43 * a50 - a40 * a53)) +
                         a13 * (a20 * (a42 * a54 - a44 * a52) +
                                a22 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a52 - a42 * a50)) +
                         a14 * (a20 * (a43 * a52 - a42 * a53) +
                                a22 * (a40 * a53 - a43 * a50) +
                                a23 * (a42 * a50 - a40 * a52)))),
        invdet * (a00 * (a11 * (a23 * (a45 * a54 - a44 * a55) +
                                a24 * (a43 * a55 - a45 * a53) +
                                a25 * (a44 * a53 - a43 * a54)) +
                         a13 * (a21 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a54 - a44 * a51)) +
                         a14 * (a21 * (a45 * a53 - a43 * a55) +
                                a23 * (a41 * a55 - a45 * a51) +
                                a25 * (a43 * a51 - a41 * a53)) +
                         a15 * (a21 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a51 - a41 * a54) +
                                a24 * (a41 * a53 - a43 * a51))) +
                  a01 * (a10 * (a23 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a53 - a43 * a55) +
                                a25 * (a43 * a54 - a44 * a53)) +
                         a13 * (a20 * (a45 * a54 - a44 * a55) +
                                a24 * (a40 * a55 - a45 * a50) +
                                a25 * (a44 * a50 - a40 * a54)) +
                         a14 * (a20 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a53 - a43 * a50)) +
                         a15 * (a20 * (a44 * a53 - a43 * a54) +
                                a23 * (a40 * a54 - a44 * a50) +
                                a24 * (a43 * a50 - a40 * a53))) +
                  a03 * (a10 * (a21 * (a45 * a54 - a44 * a55) +
                                a24 * (a41 * a55 - a45 * a51) +
                                a25 * (a44 * a51 - a41 * a54)) +
                         a11 * (a20 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a54 - a44 * a50)) +
                         a14 * (a20 * (a45 * a51 - a41 * a55) +
                                a21 * (a40 * a55 - a45 * a50) +
                                a25 * (a41 * a50 - a40 * a51)) +
                         a15 * (a20 * (a41 * a54 - a44 * a51) +
                                a21 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a51 - a41 * a50))) +
                  a04 * (a10 * (a21 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a53 - a43 * a51)) +
                         a11 * (a20 * (a45 * a53 - a43 * a55) +
                                a23 * (a40 * a55 - a45 * a50) +
                                a25 * (a43 * a50 - a40 * a53)) +
                         a13 * (a20 * (a41 * a55 - a45 * a51) +
                                a21 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a51 - a41 * a50)) +
                         a15 * (a20 * (a43 * a51 - a41 * a53) +
                                a21 * (a40 * a53 - a43 * a50) +
                                a23 * (a41 * a50 - a40 * a51))) +
                  a05 * (a10 * (a21 * (a44 * a53 - a43 * a54) +
                                a23 * (a41 * a54 - a44 * a51) +
                                a24 * (a43 * a51 - a41 * a53)) +
                         a11 * (a20 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a53 - a43 * a50)) +
                         a13 * (a20 * (a44 * a51 - a41 * a54) +
                                a21 * (a40 * a54 - a44 * a50) +
                                a24 * (a41 * a50 - a40 * a51)) +
                         a14 * (a20 * (a41 * a53 - a43 * a51) +
                                a21 * (a43 * a50 - a40 * a53) +
                                a23 * (a40 * a51 - a41 * a50)))),
        invdet * (a00 * (a11 * (a22 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a52 - a42 * a55) +
                                a25 * (a42 * a54 - a44 * a52)) +
                         a12 * (a21 * (a45 * a54 - a44 * a55) +
                                a24 * (a41 * a55 - a45 * a51) +
                                a25 * (a44 * a51 - a41 * a54)) +
                         a14 * (a21 * (a42 * a55 - a45 * a52) +
                                a22 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a52 - a42 * a51)) +
                         a15 * (a21 * (a44 * a52 - a42 * a54) +
                                a22 * (a41 * a54 - a44 * a51) +
                                a24 * (a42 * a51 - a41 * a52))) +
                  a01 * (a10 * (a22 * (a45 * a54 - a44 * a55) +
                                a24 * (a42 * a55 - a45 * a52) +
                                a25 * (a44 * a52 - a42 * a54)) +
                         a12 * (a20 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a54 - a44 * a50)) +
                         a14 * (a20 * (a45 * a52 - a42 * a55) +
                                a22 * (a40 * a55 - a45 * a50) +
                                a25 * (a42 * a50 - a40 * a52)) +
                         a15 * (a20 * (a42 * a54 - a44 * a52) +
                                a22 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a52 - a42 * a50))) +
                  a02 * (a10 * (a21 * (a44 * a55 - a45 * a54) +
                                a24 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a54 - a44 * a51)) +
                         a11 * (a20 * (a45 * a54 - a44 * a55) +
                                a24 * (a40 * a55 - a45 * a50) +
                                a25 * (a44 * a50 - a40 * a54)) +
                         a14 * (a20 * (a41 * a55 - a45 * a51) +
                                a21 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a51 - a41 * a50)) +
                         a15 * (a20 * (a44 * a51 - a41 * a54) +
                                a21 * (a40 * a54 - a44 * a50) +
                                a24 * (a41 * a50 - a40 * a51))) +
                  a04 * (a10 * (a21 * (a45 * a52 - a42 * a55) +
                                a22 * (a41 * a55 - a45 * a51) +
                                a25 * (a42 * a51 - a41 * a52)) +
                         a11 * (a20 * (a42 * a55 - a45 * a52) +
                                a22 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a52 - a42 * a50)) +
                         a12 * (a20 * (a45 * a51 - a41 * a55) +
                                a21 * (a40 * a55 - a45 * a50) +
                                a25 * (a41 * a50 - a40 * a51)) +
                         a15 * (a20 * (a41 * a52 - a42 * a51) +
                                a21 * (a42 * a50 - a40 * a52) +
                                a22 * (a40 * a51 - a41 * a50))) +
                  a05 * (a10 * (a21 * (a42 * a54 - a44 * a52) +
                                a22 * (a44 * a51 - a41 * a54) +
                                a24 * (a41 * a52 - a42 * a51)) +
                         a11 * (a20 * (a44 * a52 - a42 * a54) +
                                a22 * (a40 * a54 - a44 * a50) +
                                a24 * (a42 * a50 - a40 * a52)) +
                         a12 * (a20 * (a41 * a54 - a44 * a51) +
                                a21 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a51 - a41 * a50)) +
                         a14 * (a20 * (a42 * a51 - a41 * a52) +
                                a21 * (a40 * a52 - a42 * a50) +
                                a22 * (a41 * a50 - a40 * a51)))),
        invdet * (a00 * (a11 * (a22 * (a45 * a53 - a43 * a55) +
                                a23 * (a42 * a55 - a45 * a52) +
                                a25 * (a43 * a52 - a42 * a53)) +
                         a12 * (a21 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a53 - a43 * a51)) +
                         a13 * (a21 * (a45 * a52 - a42 * a55) +
                                a22 * (a41 * a55 - a45 * a51) +
                                a25 * (a42 * a51 - a41 * a52)) +
                         a15 * (a21 * (a42 * a53 - a43 * a52) +
                                a22 * (a43 * a51 - a41 * a53) +
                                a23 * (a41 * a52 - a42 * a51))) +
                  a01 * (a10 * (a22 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a52 - a42 * a55) +
                                a25 * (a42 * a53 - a43 * a52)) +
                         a12 * (a20 * (a45 * a53 - a43 * a55) +
                                a23 * (a40 * a55 - a45 * a50) +
                                a25 * (a43 * a50 - a40 * a53)) +
                         a13 * (a20 * (a42 * a55 - a45 * a52) +
                                a22 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a52 - a42 * a50)) +
                         a15 * (a20 * (a43 * a52 - a42 * a53) +
                                a22 * (a40 * a53 - a43 * a50) +
                                a23 * (a42 * a50 - a40 * a52))) +
                  a02 * (a10 * (a21 * (a45 * a53 - a43 * a55) +
                                a23 * (a41 * a55 - a45 * a51) +
                                a25 * (a43 * a51 - a41 * a53)) +
                         a11 * (a20 * (a43 * a55 - a45 * a53) +
                                a23 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a53 - a43 * a50)) +
                         a13 * (a20 * (a45 * a51 - a41 * a55) +
                                a21 * (a40 * a55 - a45 * a50) +
                                a25 * (a41 * a50 - a40 * a51)) +
                         a15 * (a20 * (a41 * a53 - a43 * a51) +
                                a21 * (a43 * a50 - a40 * a53) +
                                a23 * (a40 * a51 - a41 * a50))) +
                  a03 * (a10 * (a21 * (a42 * a55 - a45 * a52) +
                                a22 * (a45 * a51 - a41 * a55) +
                                a25 * (a41 * a52 - a42 * a51)) +
                         a11 * (a20 * (a45 * a52 - a42 * a55) +
                                a22 * (a40 * a55 - a45 * a50) +
                                a25 * (a42 * a50 - a40 * a52)) +
                         a12 * (a20 * (a41 * a55 - a45 * a51) +
                                a21 * (a45 * a50 - a40 * a55) +
                                a25 * (a40 * a51 - a41 * a50)) +
                         a15 * (a20 * (a42 * a51 - a41 * a52) +
                                a21 * (a40 * a52 - a42 * a50) +
                                a22 * (a41 * a50 - a40 * a51))) +
                  a05 * (a10 * (a21 * (a43 * a52 - a42 * a53) +
                                a22 * (a41 * a53 - a43 * a51) +
                                a23 * (a42 * a51 - a41 * a52)) +
                         a11 * (a20 * (a42 * a53 - a43 * a52) +
                                a22 * (a43 * a50 - a40 * a53) +
                                a23 * (a40 * a52 - a42 * a50)) +
                         a12 * (a20 * (a43 * a51 - a41 * a53) +
                                a21 * (a40 * a53 - a43 * a50) +
                                a23 * (a41 * a50 - a40 * a51)) +
                         a13 * (a20 * (a41 * a52 - a42 * a51) +
                                a21 * (a42 * a50 - a40 * a52) +
                                a22 * (a40 * a51 - a41 * a50)))),
        invdet * (a00 * (a11 * (a22 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a52 - a42 * a54) +
                                a24 * (a42 * a53 - a43 * a52)) +
                         a12 * (a21 * (a44 * a53 - a43 * a54) +
                                a23 * (a41 * a54 - a44 * a51) +
                                a24 * (a43 * a51 - a41 * a53)) +
                         a13 * (a21 * (a42 * a54 - a44 * a52) +
                                a22 * (a44 * a51 - a41 * a54) +
                                a24 * (a41 * a52 - a42 * a51)) +
                         a14 * (a21 * (a43 * a52 - a42 * a53) +
                                a22 * (a41 * a53 - a43 * a51) +
                                a23 * (a42 * a51 - a41 * a52))) +
                  a01 * (a10 * (a22 * (a44 * a53 - a43 * a54) +
                                a23 * (a42 * a54 - a44 * a52) +
                                a24 * (a43 * a52 - a42 * a53)) +
                         a12 * (a20 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a53 - a43 * a50)) +
                         a13 * (a20 * (a44 * a52 - a42 * a54) +
                                a22 * (a40 * a54 - a44 * a50) +
                                a24 * (a42 * a50 - a40 * a52)) +
                         a14 * (a20 * (a42 * a53 - a43 * a52) +
                                a22 * (a43 * a50 - a40 * a53) +
                                a23 * (a40 * a52 - a42 * a50))) +
                  a02 * (a10 * (a21 * (a43 * a54 - a44 * a53) +
                                a23 * (a44 * a51 - a41 * a54) +
                                a24 * (a41 * a53 - a43 * a51)) +
                         a11 * (a20 * (a44 * a53 - a43 * a54) +
                                a23 * (a40 * a54 - a44 * a50) +
                                a24 * (a43 * a50 - a40 * a53)) +
                         a13 * (a20 * (a41 * a54 - a44 * a51) +
                                a21 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a51 - a41 * a50)) +
                         a14 * (a20 * (a43 * a51 - a41 * a53) +
                                a21 * (a40 * a53 - a43 * a50) +
                                a23 * (a41 * a50 - a40 * a51))) +
                  a03 * (a10 * (a21 * (a44 * a52 - a42 * a54) +
                                a22 * (a41 * a54 - a44 * a51) +
                                a24 * (a42 * a51 - a41 * a52)) +
                         a11 * (a20 * (a42 * a54 - a44 * a52) +
                                a22 * (a44 * a50 - a40 * a54) +
                                a24 * (a40 * a52 - a42 * a50)) +
                         a12 * (a20 * (a44 * a51 - a41 * a54) +
                                a21 * (a40 * a54 - a44 * a50) +
                                a24 * (a41 * a50 - a40 * a51)) +
                         a14 * (a20 * (a41 * a52 - a42 * a51) +
                                a21 * (a42 * a50 - a40 * a52) +
                                a22 * (a40 * a51 - a41 * a50))) +
                  a04 * (a10 * (a21 * (a42 * a53 - a43 * a52) +
                                a22 * (a43 * a51 - a41 * a53) +
                                a23 * (a41 * a52 - a42 * a51)) +
                         a11 * (a20 * (a43 * a52 - a42 * a53) +
                                a22 * (a40 * a53 - a43 * a50) +
                                a23 * (a42 * a50 - a40 * a52)) +
                         a12 * (a20 * (a41 * a53 - a43 * a51) +
                                a21 * (a43 * a50 - a40 * a53) +
                                a23 * (a40 * a51 - a41 * a50)) +
                         a13 * (a20 * (a42 * a51 - a41 * a52) +
                                a21 * (a40 * a52 - a42 * a50) +
                                a22 * (a41 * a50 - a40 * a51)))),
        invdet * (a01 * (a12 * (a23 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a53 - a33 * a55) +
                                a25 * (a33 * a54 - a34 * a53)) +
                         a13 * (a22 * (a35 * a54 - a34 * a55) +
                                a24 * (a32 * a55 - a35 * a52) +
                                a25 * (a34 * a52 - a32 * a54)) +
                         a14 * (a22 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a52 - a32 * a55) +
                                a25 * (a32 * a53 - a33 * a52)) +
                         a15 * (a22 * (a34 * a53 - a33 * a54) +
                                a23 * (a32 * a54 - a34 * a52) +
                                a24 * (a33 * a52 - a32 * a53))) +
                  a02 * (a11 * (a23 * (a35 * a54 - a34 * a55) +
                                a24 * (a33 * a55 - a35 * a53) +
                                a25 * (a34 * a53 - a33 * a54)) +
                         a13 * (a21 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a54 - a34 * a51)) +
                         a14 * (a21 * (a35 * a53 - a33 * a55) +
                                a23 * (a31 * a55 - a35 * a51) +
                                a25 * (a33 * a51 - a31 * a53)) +
                         a15 * (a21 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a51 - a31 * a54) +
                                a24 * (a31 * a53 - a33 * a51))) +
                  a03 * (a11 * (a22 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a52 - a32 * a55) +
                                a25 * (a32 * a54 - a34 * a52)) +
                         a12 * (a21 * (a35 * a54 - a34 * a55) +
                                a24 * (a31 * a55 - a35 * a51) +
                                a25 * (a34 * a51 - a31 * a54)) +
                         a14 * (a21 * (a32 * a55 - a35 * a52) +
                                a22 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a52 - a32 * a51)) +
                         a15 * (a21 * (a34 * a52 - a32 * a54) +
                                a22 * (a31 * a54 - a34 * a51) +
                                a24 * (a32 * a51 - a31 * a52))) +
                  a04 * (a11 * (a22 * (a35 * a53 - a33 * a55) +
                                a23 * (a32 * a55 - a35 * a52) +
                                a25 * (a33 * a52 - a32 * a53)) +
                         a12 * (a21 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a53 - a33 * a51)) +
                         a13 * (a21 * (a35 * a52 - a32 * a55) +
                                a22 * (a31 * a55 - a35 * a51) +
                                a25 * (a32 * a51 - a31 * a52)) +
                         a15 * (a21 * (a32 * a53 - a33 * a52) +
                                a22 * (a33 * a51 - a31 * a53) +
                                a23 * (a31 * a52 - a32 * a51))) +
                  a05 * (a11 * (a22 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a52 - a32 * a54) +
                                a24 * (a32 * a53 - a33 * a52)) +
                         a12 * (a21 * (a34 * a53 - a33 * a54) +
                                a23 * (a31 * a54 - a34 * a51) +
                                a24 * (a33 * a51 - a31 * a53)) +
                         a13 * (a21 * (a32 * a54 - a34 * a52) +
                                a22 * (a34 * a51 - a31 * a54) +
                                a24 * (a31 * a52 - a32 * a51)) +
                         a14 * (a21 * (a33 * a52 - a32 * a53) +
                                a22 * (a31 * a53 - a33 * a51) +
                                a23 * (a32 * a51 - a31 * a52)))),
        invdet * (a00 * (a12 * (a23 * (a35 * a54 - a34 * a55) +
                                a24 * (a33 * a55 - a35 * a53) +
                                a25 * (a34 * a53 - a33 * a54)) +
                         a13 * (a22 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a52 - a32 * a55) +
                                a25 * (a32 * a54 - a34 * a52)) +
                         a14 * (a22 * (a35 * a53 - a33 * a55) +
                                a23 * (a32 * a55 - a35 * a52) +
                                a25 * (a33 * a52 - a32 * a53)) +
                         a15 * (a22 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a52 - a32 * a54) +
                                a24 * (a32 * a53 - a33 * a52))) +
                  a02 * (a10 * (a23 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a53 - a33 * a55) +
                                a25 * (a33 * a54 - a34 * a53)) +
                         a13 * (a20 * (a35 * a54 - a34 * a55) +
                                a24 * (a30 * a55 - a35 * a50) +
                                a25 * (a34 * a50 - a30 * a54)) +
                         a14 * (a20 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a53 - a33 * a50)) +
                         a15 * (a20 * (a34 * a53 - a33 * a54) +
                                a23 * (a30 * a54 - a34 * a50) +
                                a24 * (a33 * a50 - a30 * a53))) +
                  a03 * (a10 * (a22 * (a35 * a54 - a34 * a55) +
                                a24 * (a32 * a55 - a35 * a52) +
                                a25 * (a34 * a52 - a32 * a54)) +
                         a12 * (a20 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a54 - a34 * a50)) +
                         a14 * (a20 * (a35 * a52 - a32 * a55) +
                                a22 * (a30 * a55 - a35 * a50) +
                                a25 * (a32 * a50 - a30 * a52)) +
                         a15 * (a20 * (a32 * a54 - a34 * a52) +
                                a22 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a52 - a32 * a50))) +
                  a04 * (a10 * (a22 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a52 - a32 * a55) +
                                a25 * (a32 * a53 - a33 * a52)) +
                         a12 * (a20 * (a35 * a53 - a33 * a55) +
                                a23 * (a30 * a55 - a35 * a50) +
                                a25 * (a33 * a50 - a30 * a53)) +
                         a13 * (a20 * (a32 * a55 - a35 * a52) +
                                a22 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a52 - a32 * a50)) +
                         a15 * (a20 * (a33 * a52 - a32 * a53) +
                                a22 * (a30 * a53 - a33 * a50) +
                                a23 * (a32 * a50 - a30 * a52))) +
                  a05 * (a10 * (a22 * (a34 * a53 - a33 * a54) +
                                a23 * (a32 * a54 - a34 * a52) +
                                a24 * (a33 * a52 - a32 * a53)) +
                         a12 * (a20 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a53 - a33 * a50)) +
                         a13 * (a20 * (a34 * a52 - a32 * a54) +
                                a22 * (a30 * a54 - a34 * a50) +
                                a24 * (a32 * a50 - a30 * a52)) +
                         a14 * (a20 * (a32 * a53 - a33 * a52) +
                                a22 * (a33 * a50 - a30 * a53) +
                                a23 * (a30 * a52 - a32 * a50)))),
        invdet * (a00 * (a11 * (a23 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a53 - a33 * a55) +
                                a25 * (a33 * a54 - a34 * a53)) +
                         a13 * (a21 * (a35 * a54 - a34 * a55) +
                                a24 * (a31 * a55 - a35 * a51) +
                                a25 * (a34 * a51 - a31 * a54)) +
                         a14 * (a21 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a53 - a33 * a51)) +
                         a15 * (a21 * (a34 * a53 - a33 * a54) +
                                a23 * (a31 * a54 - a34 * a51) +
                                a24 * (a33 * a51 - a31 * a53))) +
                  a01 * (a10 * (a23 * (a35 * a54 - a34 * a55) +
                                a24 * (a33 * a55 - a35 * a53) +
                                a25 * (a34 * a53 - a33 * a54)) +
                         a13 * (a20 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a54 - a34 * a50)) +
                         a14 * (a20 * (a35 * a53 - a33 * a55) +
                                a23 * (a30 * a55 - a35 * a50) +
                                a25 * (a33 * a50 - a30 * a53)) +
                         a15 * (a20 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a53 - a33 * a50))) +
                  a03 * (a10 * (a21 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a54 - a34 * a51)) +
                         a11 * (a20 * (a35 * a54 - a34 * a55) +
                                a24 * (a30 * a55 - a35 * a50) +
                                a25 * (a34 * a50 - a30 * a54)) +
                         a14 * (a20 * (a31 * a55 - a35 * a51) +
                                a21 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a51 - a31 * a50)) +
                         a15 * (a20 * (a34 * a51 - a31 * a54) +
                                a21 * (a30 * a54 - a34 * a50) +
                                a24 * (a31 * a50 - a30 * a51))) +
                  a04 * (a10 * (a21 * (a35 * a53 - a33 * a55) +
                                a23 * (a31 * a55 - a35 * a51) +
                                a25 * (a33 * a51 - a31 * a53)) +
                         a11 * (a20 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a53 - a33 * a50)) +
                         a13 * (a20 * (a35 * a51 - a31 * a55) +
                                a21 * (a30 * a55 - a35 * a50) +
                                a25 * (a31 * a50 - a30 * a51)) +
                         a15 * (a20 * (a31 * a53 - a33 * a51) +
                                a21 * (a33 * a50 - a30 * a53) +
                                a23 * (a30 * a51 - a31 * a50))) +
                  a05 * (a10 * (a21 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a51 - a31 * a54) +
                                a24 * (a31 * a53 - a33 * a51)) +
                         a11 * (a20 * (a34 * a53 - a33 * a54) +
                                a23 * (a30 * a54 - a34 * a50) +
                                a24 * (a33 * a50 - a30 * a53)) +
                         a13 * (a20 * (a31 * a54 - a34 * a51) +
                                a21 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a51 - a31 * a50)) +
                         a14 * (a20 * (a33 * a51 - a31 * a53) +
                                a21 * (a30 * a53 - a33 * a50) +
                                a23 * (a31 * a50 - a30 * a51)))),
        invdet * (a00 * (a11 * (a22 * (a35 * a54 - a34 * a55) +
                                a24 * (a32 * a55 - a35 * a52) +
                                a25 * (a34 * a52 - a32 * a54)) +
                         a12 * (a21 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a54 - a34 * a51)) +
                         a14 * (a21 * (a35 * a52 - a32 * a55) +
                                a22 * (a31 * a55 - a35 * a51) +
                                a25 * (a32 * a51 - a31 * a52)) +
                         a15 * (a21 * (a32 * a54 - a34 * a52) +
                                a22 * (a34 * a51 - a31 * a54) +
                                a24 * (a31 * a52 - a32 * a51))) +
                  a01 * (a10 * (a22 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a52 - a32 * a55) +
                                a25 * (a32 * a54 - a34 * a52)) +
                         a12 * (a20 * (a35 * a54 - a34 * a55) +
                                a24 * (a30 * a55 - a35 * a50) +
                                a25 * (a34 * a50 - a30 * a54)) +
                         a14 * (a20 * (a32 * a55 - a35 * a52) +
                                a22 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a52 - a32 * a50)) +
                         a15 * (a20 * (a34 * a52 - a32 * a54) +
                                a22 * (a30 * a54 - a34 * a50) +
                                a24 * (a32 * a50 - a30 * a52))) +
                  a02 * (a10 * (a21 * (a35 * a54 - a34 * a55) +
                                a24 * (a31 * a55 - a35 * a51) +
                                a25 * (a34 * a51 - a31 * a54)) +
                         a11 * (a20 * (a34 * a55 - a35 * a54) +
                                a24 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a54 - a34 * a50)) +
                         a14 * (a20 * (a35 * a51 - a31 * a55) +
                                a21 * (a30 * a55 - a35 * a50) +
                                a25 * (a31 * a50 - a30 * a51)) +
                         a15 * (a20 * (a31 * a54 - a34 * a51) +
                                a21 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a51 - a31 * a50))) +
                  a04 * (a10 * (a21 * (a32 * a55 - a35 * a52) +
                                a22 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a52 - a32 * a51)) +
                         a11 * (a20 * (a35 * a52 - a32 * a55) +
                                a22 * (a30 * a55 - a35 * a50) +
                                a25 * (a32 * a50 - a30 * a52)) +
                         a12 * (a20 * (a31 * a55 - a35 * a51) +
                                a21 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a51 - a31 * a50)) +
                         a15 * (a20 * (a32 * a51 - a31 * a52) +
                                a21 * (a30 * a52 - a32 * a50) +
                                a22 * (a31 * a50 - a30 * a51))) +
                  a05 * (a10 * (a21 * (a34 * a52 - a32 * a54) +
                                a22 * (a31 * a54 - a34 * a51) +
                                a24 * (a32 * a51 - a31 * a52)) +
                         a11 * (a20 * (a32 * a54 - a34 * a52) +
                                a22 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a52 - a32 * a50)) +
                         a12 * (a20 * (a34 * a51 - a31 * a54) +
                                a21 * (a30 * a54 - a34 * a50) +
                                a24 * (a31 * a50 - a30 * a51)) +
                         a14 * (a20 * (a31 * a52 - a32 * a51) +
                                a21 * (a32 * a50 - a30 * a52) +
                                a22 * (a30 * a51 - a31 * a50)))),
        invdet * (a00 * (a11 * (a22 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a52 - a32 * a55) +
                                a25 * (a32 * a53 - a33 * a52)) +
                         a12 * (a21 * (a35 * a53 - a33 * a55) +
                                a23 * (a31 * a55 - a35 * a51) +
                                a25 * (a33 * a51 - a31 * a53)) +
                         a13 * (a21 * (a32 * a55 - a35 * a52) +
                                a22 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a52 - a32 * a51)) +
                         a15 * (a21 * (a33 * a52 - a32 * a53) +
                                a22 * (a31 * a53 - a33 * a51) +
                                a23 * (a32 * a51 - a31 * a52))) +
                  a01 * (a10 * (a22 * (a35 * a53 - a33 * a55) +
                                a23 * (a32 * a55 - a35 * a52) +
                                a25 * (a33 * a52 - a32 * a53)) +
                         a12 * (a20 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a53 - a33 * a50)) +
                         a13 * (a20 * (a35 * a52 - a32 * a55) +
                                a22 * (a30 * a55 - a35 * a50) +
                                a25 * (a32 * a50 - a30 * a52)) +
                         a15 * (a20 * (a32 * a53 - a33 * a52) +
                                a22 * (a33 * a50 - a30 * a53) +
                                a23 * (a30 * a52 - a32 * a50))) +
                  a02 * (a10 * (a21 * (a33 * a55 - a35 * a53) +
                                a23 * (a35 * a51 - a31 * a55) +
                                a25 * (a31 * a53 - a33 * a51)) +
                         a11 * (a20 * (a35 * a53 - a33 * a55) +
                                a23 * (a30 * a55 - a35 * a50) +
                                a25 * (a33 * a50 - a30 * a53)) +
                         a13 * (a20 * (a31 * a55 - a35 * a51) +
                                a21 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a51 - a31 * a50)) +
                         a15 * (a20 * (a33 * a51 - a31 * a53) +
                                a21 * (a30 * a53 - a33 * a50) +
                                a23 * (a31 * a50 - a30 * a51))) +
                  a03 * (a10 * (a21 * (a35 * a52 - a32 * a55) +
                                a22 * (a31 * a55 - a35 * a51) +
                                a25 * (a32 * a51 - a31 * a52)) +
                         a11 * (a20 * (a32 * a55 - a35 * a52) +
                                a22 * (a35 * a50 - a30 * a55) +
                                a25 * (a30 * a52 - a32 * a50)) +
                         a12 * (a20 * (a35 * a51 - a31 * a55) +
                                a21 * (a30 * a55 - a35 * a50) +
                                a25 * (a31 * a50 - a30 * a51)) +
                         a15 * (a20 * (a31 * a52 - a32 * a51) +
                                a21 * (a32 * a50 - a30 * a52) +
                                a22 * (a30 * a51 - a31 * a50))) +
                  a05 * (a10 * (a21 * (a32 * a53 - a33 * a52) +
                                a22 * (a33 * a51 - a31 * a53) +
                                a23 * (a31 * a52 - a32 * a51)) +
                         a11 * (a20 * (a33 * a52 - a32 * a53) +
                                a22 * (a30 * a53 - a33 * a50) +
                                a23 * (a32 * a50 - a30 * a52)) +
                         a12 * (a20 * (a31 * a53 - a33 * a51) +
                                a21 * (a33 * a50 - a30 * a53) +
                                a23 * (a30 * a51 - a31 * a50)) +
                         a13 * (a20 * (a32 * a51 - a31 * a52) +
                                a21 * (a30 * a52 - a32 * a50) +
                                a22 * (a31 * a50 - a30 * a51)))),
        invdet * (a00 * (a11 * (a22 * (a34 * a53 - a33 * a54) +
                                a23 * (a32 * a54 - a34 * a52) +
                                a24 * (a33 * a52 - a32 * a53)) +
                         a12 * (a21 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a51 - a31 * a54) +
                                a24 * (a31 * a53 - a33 * a51)) +
                         a13 * (a21 * (a34 * a52 - a32 * a54) +
                                a22 * (a31 * a54 - a34 * a51) +
                                a24 * (a32 * a51 - a31 * a52)) +
                         a14 * (a21 * (a32 * a53 - a33 * a52) +
                                a22 * (a33 * a51 - a31 * a53) +
                                a23 * (a31 * a52 - a32 * a51))) +
                  a01 * (a10 * (a22 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a52 - a32 * a54) +
                                a24 * (a32 * a53 - a33 * a52)) +
                         a12 * (a20 * (a34 * a53 - a33 * a54) +
                                a23 * (a30 * a54 - a34 * a50) +
                                a24 * (a33 * a50 - a30 * a53)) +
                         a13 * (a20 * (a32 * a54 - a34 * a52) +
                                a22 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a52 - a32 * a50)) +
                         a14 * (a20 * (a33 * a52 - a32 * a53) +
                                a22 * (a30 * a53 - a33 * a50) +
                                a23 * (a32 * a50 - a30 * a52))) +
                  a02 * (a10 * (a21 * (a34 * a53 - a33 * a54) +
                                a23 * (a31 * a54 - a34 * a51) +
                                a24 * (a33 * a51 - a31 * a53)) +
                         a11 * (a20 * (a33 * a54 - a34 * a53) +
                                a23 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a53 - a33 * a50)) +
                         a13 * (a20 * (a34 * a51 - a31 * a54) +
                                a21 * (a30 * a54 - a34 * a50) +
                                a24 * (a31 * a50 - a30 * a51)) +
                         a14 * (a20 * (a31 * a53 - a33 * a51) +
                                a21 * (a33 * a50 - a30 * a53) +
                                a23 * (a30 * a51 - a31 * a50))) +
                  a03 * (a10 * (a21 * (a32 * a54 - a34 * a52) +
                                a22 * (a34 * a51 - a31 * a54) +
                                a24 * (a31 * a52 - a32 * a51)) +
                         a11 * (a20 * (a34 * a52 - a32 * a54) +
                                a22 * (a30 * a54 - a34 * a50) +
                                a24 * (a32 * a50 - a30 * a52)) +
                         a12 * (a20 * (a31 * a54 - a34 * a51) +
                                a21 * (a34 * a50 - a30 * a54) +
                                a24 * (a30 * a51 - a31 * a50)) +
                         a14 * (a20 * (a32 * a51 - a31 * a52) +
                                a21 * (a30 * a52 - a32 * a50) +
                                a22 * (a31 * a50 - a30 * a51))) +
                  a04 * (a10 * (a21 * (a33 * a52 - a32 * a53) +
                                a22 * (a31 * a53 - a33 * a51) +
                                a23 * (a32 * a51 - a31 * a52)) +
                         a11 * (a20 * (a32 * a53 - a33 * a52) +
                                a22 * (a33 * a50 - a30 * a53) +
                                a23 * (a30 * a52 - a32 * a50)) +
                         a12 * (a20 * (a33 * a51 - a31 * a53) +
                                a21 * (a30 * a53 - a33 * a50) +
                                a23 * (a31 * a50 - a30 * a51)) +
                         a13 * (a20 * (a31 * a52 - a32 * a51) +
                                a21 * (a32 * a50 - a30 * a52) +
                                a22 * (a30 * a51 - a31 * a50)))),
        invdet * (a01 * (a12 * (a23 * (a35 * a44 - a34 * a45) +
                                a24 * (a33 * a45 - a35 * a43) +
                                a25 * (a34 * a43 - a33 * a44)) +
                         a13 * (a22 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a42 - a32 * a45) +
                                a25 * (a32 * a44 - a34 * a42)) +
                         a14 * (a22 * (a35 * a43 - a33 * a45) +
                                a23 * (a32 * a45 - a35 * a42) +
                                a25 * (a33 * a42 - a32 * a43)) +
                         a15 * (a22 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a42 - a32 * a44) +
                                a24 * (a32 * a43 - a33 * a42))) +
                  a02 * (a11 * (a23 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a43 - a33 * a45) +
                                a25 * (a33 * a44 - a34 * a43)) +
                         a13 * (a21 * (a35 * a44 - a34 * a45) +
                                a24 * (a31 * a45 - a35 * a41) +
                                a25 * (a34 * a41 - a31 * a44)) +
                         a14 * (a21 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a43 - a33 * a41)) +
                         a15 * (a21 * (a34 * a43 - a33 * a44) +
                                a23 * (a31 * a44 - a34 * a41) +
                                a24 * (a33 * a41 - a31 * a43))) +
                  a03 * (a11 * (a22 * (a35 * a44 - a34 * a45) +
                                a24 * (a32 * a45 - a35 * a42) +
                                a25 * (a34 * a42 - a32 * a44)) +
                         a12 * (a21 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a44 - a34 * a41)) +
                         a14 * (a21 * (a35 * a42 - a32 * a45) +
                                a22 * (a31 * a45 - a35 * a41) +
                                a25 * (a32 * a41 - a31 * a42)) +
                         a15 * (a21 * (a32 * a44 - a34 * a42) +
                                a22 * (a34 * a41 - a31 * a44) +
                                a24 * (a31 * a42 - a32 * a41))) +
                  a04 * (a11 * (a22 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a42 - a32 * a45) +
                                a25 * (a32 * a43 - a33 * a42)) +
                         a12 * (a21 * (a35 * a43 - a33 * a45) +
                                a23 * (a31 * a45 - a35 * a41) +
                                a25 * (a33 * a41 - a31 * a43)) +
                         a13 * (a21 * (a32 * a45 - a35 * a42) +
                                a22 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a42 - a32 * a41)) +
                         a15 * (a21 * (a33 * a42 - a32 * a43) +
                                a22 * (a31 * a43 - a33 * a41) +
                                a23 * (a32 * a41 - a31 * a42))) +
                  a05 * (a11 * (a22 * (a34 * a43 - a33 * a44) +
                                a23 * (a32 * a44 - a34 * a42) +
                                a24 * (a33 * a42 - a32 * a43)) +
                         a12 * (a21 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a41 - a31 * a44) +
                                a24 * (a31 * a43 - a33 * a41)) +
                         a13 * (a21 * (a34 * a42 - a32 * a44) +
                                a22 * (a31 * a44 - a34 * a41) +
                                a24 * (a32 * a41 - a31 * a42)) +
                         a14 * (a21 * (a32 * a43 - a33 * a42) +
                                a22 * (a33 * a41 - a31 * a43) +
                                a23 * (a31 * a42 - a32 * a41)))),
        invdet * (a00 * (a12 * (a23 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a43 - a33 * a45) +
                                a25 * (a33 * a44 - a34 * a43)) +
                         a13 * (a22 * (a35 * a44 - a34 * a45) +
                                a24 * (a32 * a45 - a35 * a42) +
                                a25 * (a34 * a42 - a32 * a44)) +
                         a14 * (a22 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a42 - a32 * a45) +
                                a25 * (a32 * a43 - a33 * a42)) +
                         a15 * (a22 * (a34 * a43 - a33 * a44) +
                                a23 * (a32 * a44 - a34 * a42) +
                                a24 * (a33 * a42 - a32 * a43))) +
                  a02 * (a10 * (a23 * (a35 * a44 - a34 * a45) +
                                a24 * (a33 * a45 - a35 * a43) +
                                a25 * (a34 * a43 - a33 * a44)) +
                         a13 * (a20 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a44 - a34 * a40)) +
                         a14 * (a20 * (a35 * a43 - a33 * a45) +
                                a23 * (a30 * a45 - a35 * a40) +
                                a25 * (a33 * a40 - a30 * a43)) +
                         a15 * (a20 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a43 - a33 * a40))) +
                  a03 * (a10 * (a22 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a42 - a32 * a45) +
                                a25 * (a32 * a44 - a34 * a42)) +
                         a12 * (a20 * (a35 * a44 - a34 * a45) +
                                a24 * (a30 * a45 - a35 * a40) +
                                a25 * (a34 * a40 - a30 * a44)) +
                         a14 * (a20 * (a32 * a45 - a35 * a42) +
                                a22 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a42 - a32 * a40)) +
                         a15 * (a20 * (a34 * a42 - a32 * a44) +
                                a22 * (a30 * a44 - a34 * a40) +
                                a24 * (a32 * a40 - a30 * a42))) +
                  a04 * (a10 * (a22 * (a35 * a43 - a33 * a45) +
                                a23 * (a32 * a45 - a35 * a42) +
                                a25 * (a33 * a42 - a32 * a43)) +
                         a12 * (a20 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a43 - a33 * a40)) +
                         a13 * (a20 * (a35 * a42 - a32 * a45) +
                                a22 * (a30 * a45 - a35 * a40) +
                                a25 * (a32 * a40 - a30 * a42)) +
                         a15 * (a20 * (a32 * a43 - a33 * a42) +
                                a22 * (a33 * a40 - a30 * a43) +
                                a23 * (a30 * a42 - a32 * a40))) +
                  a05 * (a10 * (a22 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a42 - a32 * a44) +
                                a24 * (a32 * a43 - a33 * a42)) +
                         a12 * (a20 * (a34 * a43 - a33 * a44) +
                                a23 * (a30 * a44 - a34 * a40) +
                                a24 * (a33 * a40 - a30 * a43)) +
                         a13 * (a20 * (a32 * a44 - a34 * a42) +
                                a22 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a42 - a32 * a40)) +
                         a14 * (a20 * (a33 * a42 - a32 * a43) +
                                a22 * (a30 * a43 - a33 * a40) +
                                a23 * (a32 * a40 - a30 * a42)))),
        invdet * (a00 * (a11 * (a23 * (a35 * a44 - a34 * a45) +
                                a24 * (a33 * a45 - a35 * a43) +
                                a25 * (a34 * a43 - a33 * a44)) +
                         a13 * (a21 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a44 - a34 * a41)) +
                         a14 * (a21 * (a35 * a43 - a33 * a45) +
                                a23 * (a31 * a45 - a35 * a41) +
                                a25 * (a33 * a41 - a31 * a43)) +
                         a15 * (a21 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a41 - a31 * a44) +
                                a24 * (a31 * a43 - a33 * a41))) +
                  a01 * (a10 * (a23 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a43 - a33 * a45) +
                                a25 * (a33 * a44 - a34 * a43)) +
                         a13 * (a20 * (a35 * a44 - a34 * a45) +
                                a24 * (a30 * a45 - a35 * a40) +
                                a25 * (a34 * a40 - a30 * a44)) +
                         a14 * (a20 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a43 - a33 * a40)) +
                         a15 * (a20 * (a34 * a43 - a33 * a44) +
                                a23 * (a30 * a44 - a34 * a40) +
                                a24 * (a33 * a40 - a30 * a43))) +
                  a03 * (a10 * (a21 * (a35 * a44 - a34 * a45) +
                                a24 * (a31 * a45 - a35 * a41) +
                                a25 * (a34 * a41 - a31 * a44)) +
                         a11 * (a20 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a44 - a34 * a40)) +
                         a14 * (a20 * (a35 * a41 - a31 * a45) +
                                a21 * (a30 * a45 - a35 * a40) +
                                a25 * (a31 * a40 - a30 * a41)) +
                         a15 * (a20 * (a31 * a44 - a34 * a41) +
                                a21 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a41 - a31 * a40))) +
                  a04 * (a10 * (a21 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a43 - a33 * a41)) +
                         a11 * (a20 * (a35 * a43 - a33 * a45) +
                                a23 * (a30 * a45 - a35 * a40) +
                                a25 * (a33 * a40 - a30 * a43)) +
                         a13 * (a20 * (a31 * a45 - a35 * a41) +
                                a21 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a41 - a31 * a40)) +
                         a15 * (a20 * (a33 * a41 - a31 * a43) +
                                a21 * (a30 * a43 - a33 * a40) +
                                a23 * (a31 * a40 - a30 * a41))) +
                  a05 * (a10 * (a21 * (a34 * a43 - a33 * a44) +
                                a23 * (a31 * a44 - a34 * a41) +
                                a24 * (a33 * a41 - a31 * a43)) +
                         a11 * (a20 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a43 - a33 * a40)) +
                         a13 * (a20 * (a34 * a41 - a31 * a44) +
                                a21 * (a30 * a44 - a34 * a40) +
                                a24 * (a31 * a40 - a30 * a41)) +
                         a14 * (a20 * (a31 * a43 - a33 * a41) +
                                a21 * (a33 * a40 - a30 * a43) +
                                a23 * (a30 * a41 - a31 * a40)))),
        invdet * (a00 * (a11 * (a22 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a42 - a32 * a45) +
                                a25 * (a32 * a44 - a34 * a42)) +
                         a12 * (a21 * (a35 * a44 - a34 * a45) +
                                a24 * (a31 * a45 - a35 * a41) +
                                a25 * (a34 * a41 - a31 * a44)) +
                         a14 * (a21 * (a32 * a45 - a35 * a42) +
                                a22 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a42 - a32 * a41)) +
                         a15 * (a21 * (a34 * a42 - a32 * a44) +
                                a22 * (a31 * a44 - a34 * a41) +
                                a24 * (a32 * a41 - a31 * a42))) +
                  a01 * (a10 * (a22 * (a35 * a44 - a34 * a45) +
                                a24 * (a32 * a45 - a35 * a42) +
                                a25 * (a34 * a42 - a32 * a44)) +
                         a12 * (a20 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a44 - a34 * a40)) +
                         a14 * (a20 * (a35 * a42 - a32 * a45) +
                                a22 * (a30 * a45 - a35 * a40) +
                                a25 * (a32 * a40 - a30 * a42)) +
                         a15 * (a20 * (a32 * a44 - a34 * a42) +
                                a22 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a42 - a32 * a40))) +
                  a02 * (a10 * (a21 * (a34 * a45 - a35 * a44) +
                                a24 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a44 - a34 * a41)) +
                         a11 * (a20 * (a35 * a44 - a34 * a45) +
                                a24 * (a30 * a45 - a35 * a40) +
                                a25 * (a34 * a40 - a30 * a44)) +
                         a14 * (a20 * (a31 * a45 - a35 * a41) +
                                a21 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a41 - a31 * a40)) +
                         a15 * (a20 * (a34 * a41 - a31 * a44) +
                                a21 * (a30 * a44 - a34 * a40) +
                                a24 * (a31 * a40 - a30 * a41))) +
                  a04 * (a10 * (a21 * (a35 * a42 - a32 * a45) +
                                a22 * (a31 * a45 - a35 * a41) +
                                a25 * (a32 * a41 - a31 * a42)) +
                         a11 * (a20 * (a32 * a45 - a35 * a42) +
                                a22 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a42 - a32 * a40)) +
                         a12 * (a20 * (a35 * a41 - a31 * a45) +
                                a21 * (a30 * a45 - a35 * a40) +
                                a25 * (a31 * a40 - a30 * a41)) +
                         a15 * (a20 * (a31 * a42 - a32 * a41) +
                                a21 * (a32 * a40 - a30 * a42) +
                                a22 * (a30 * a41 - a31 * a40))) +
                  a05 * (a10 * (a21 * (a32 * a44 - a34 * a42) +
                                a22 * (a34 * a41 - a31 * a44) +
                                a24 * (a31 * a42 - a32 * a41)) +
                         a11 * (a20 * (a34 * a42 - a32 * a44) +
                                a22 * (a30 * a44 - a34 * a40) +
                                a24 * (a32 * a40 - a30 * a42)) +
                         a12 * (a20 * (a31 * a44 - a34 * a41) +
                                a21 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a41 - a31 * a40)) +
                         a14 * (a20 * (a32 * a41 - a31 * a42) +
                                a21 * (a30 * a42 - a32 * a40) +
                                a22 * (a31 * a40 - a30 * a41)))),
        invdet * (a00 * (a11 * (a22 * (a35 * a43 - a33 * a45) +
                                a23 * (a32 * a45 - a35 * a42) +
                                a25 * (a33 * a42 - a32 * a43)) +
                         a12 * (a21 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a43 - a33 * a41)) +
                         a13 * (a21 * (a35 * a42 - a32 * a45) +
                                a22 * (a31 * a45 - a35 * a41) +
                                a25 * (a32 * a41 - a31 * a42)) +
                         a15 * (a21 * (a32 * a43 - a33 * a42) +
                                a22 * (a33 * a41 - a31 * a43) +
                                a23 * (a31 * a42 - a32 * a41))) +
                  a01 * (a10 * (a22 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a42 - a32 * a45) +
                                a25 * (a32 * a43 - a33 * a42)) +
                         a12 * (a20 * (a35 * a43 - a33 * a45) +
                                a23 * (a30 * a45 - a35 * a40) +
                                a25 * (a33 * a40 - a30 * a43)) +
                         a13 * (a20 * (a32 * a45 - a35 * a42) +
                                a22 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a42 - a32 * a40)) +
                         a15 * (a20 * (a33 * a42 - a32 * a43) +
                                a22 * (a30 * a43 - a33 * a40) +
                                a23 * (a32 * a40 - a30 * a42))) +
                  a02 * (a10 * (a21 * (a35 * a43 - a33 * a45) +
                                a23 * (a31 * a45 - a35 * a41) +
                                a25 * (a33 * a41 - a31 * a43)) +
                         a11 * (a20 * (a33 * a45 - a35 * a43) +
                                a23 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a43 - a33 * a40)) +
                         a13 * (a20 * (a35 * a41 - a31 * a45) +
                                a21 * (a30 * a45 - a35 * a40) +
                                a25 * (a31 * a40 - a30 * a41)) +
                         a15 * (a20 * (a31 * a43 - a33 * a41) +
                                a21 * (a33 * a40 - a30 * a43) +
                                a23 * (a30 * a41 - a31 * a40))) +
                  a03 * (a10 * (a21 * (a32 * a45 - a35 * a42) +
                                a22 * (a35 * a41 - a31 * a45) +
                                a25 * (a31 * a42 - a32 * a41)) +
                         a11 * (a20 * (a35 * a42 - a32 * a45) +
                                a22 * (a30 * a45 - a35 * a40) +
                                a25 * (a32 * a40 - a30 * a42)) +
                         a12 * (a20 * (a31 * a45 - a35 * a41) +
                                a21 * (a35 * a40 - a30 * a45) +
                                a25 * (a30 * a41 - a31 * a40)) +
                         a15 * (a20 * (a32 * a41 - a31 * a42) +
                                a21 * (a30 * a42 - a32 * a40) +
                                a22 * (a31 * a40 - a30 * a41))) +
                  a05 * (a10 * (a21 * (a33 * a42 - a32 * a43) +
                                a22 * (a31 * a43 - a33 * a41) +
                                a23 * (a32 * a41 - a31 * a42)) +
                         a11 * (a20 * (a32 * a43 - a33 * a42) +
                                a22 * (a33 * a40 - a30 * a43) +
                                a23 * (a30 * a42 - a32 * a40)) +
                         a12 * (a20 * (a33 * a41 - a31 * a43) +
                                a21 * (a30 * a43 - a33 * a40) +
                                a23 * (a31 * a40 - a30 * a41)) +
                         a13 * (a20 * (a31 * a42 - a32 * a41) +
                                a21 * (a32 * a40 - a30 * a42) +
                                a22 * (a30 * a41 - a31 * a40)))),
        invdet * (a00 * (a11 * (a22 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a42 - a32 * a44) +
                                a24 * (a32 * a43 - a33 * a42)) +
                         a12 * (a21 * (a34 * a43 - a33 * a44) +
                                a23 * (a31 * a44 - a34 * a41) +
                                a24 * (a33 * a41 - a31 * a43)) +
                         a13 * (a21 * (a32 * a44 - a34 * a42) +
                                a22 * (a34 * a41 - a31 * a44) +
                                a24 * (a31 * a42 - a32 * a41)) +
                         a14 * (a21 * (a33 * a42 - a32 * a43) +
                                a22 * (a31 * a43 - a33 * a41) +
                                a23 * (a32 * a41 - a31 * a42))) +
                  a01 * (a10 * (a22 * (a34 * a43 - a33 * a44) +
                                a23 * (a32 * a44 - a34 * a42) +
                                a24 * (a33 * a42 - a32 * a43)) +
                         a12 * (a20 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a43 - a33 * a40)) +
                         a13 * (a20 * (a34 * a42 - a32 * a44) +
                                a22 * (a30 * a44 - a34 * a40) +
                                a24 * (a32 * a40 - a30 * a42)) +
                         a14 * (a20 * (a32 * a43 - a33 * a42) +
                                a22 * (a33 * a40 - a30 * a43) +
                                a23 * (a30 * a42 - a32 * a40))) +
                  a02 * (a10 * (a21 * (a33 * a44 - a34 * a43) +
                                a23 * (a34 * a41 - a31 * a44) +
                                a24 * (a31 * a43 - a33 * a41)) +
                         a11 * (a20 * (a34 * a43 - a33 * a44) +
                                a23 * (a30 * a44 - a34 * a40) +
                                a24 * (a33 * a40 - a30 * a43)) +
                         a13 * (a20 * (a31 * a44 - a34 * a41) +
                                a21 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a41 - a31 * a40)) +
                         a14 * (a20 * (a33 * a41 - a31 * a43) +
                                a21 * (a30 * a43 - a33 * a40) +
                                a23 * (a31 * a40 - a30 * a41))) +
                  a03 * (a10 * (a21 * (a34 * a42 - a32 * a44) +
                                a22 * (a31 * a44 - a34 * a41) +
                                a24 * (a32 * a41 - a31 * a42)) +
                         a11 * (a20 * (a32 * a44 - a34 * a42) +
                                a22 * (a34 * a40 - a30 * a44) +
                                a24 * (a30 * a42 - a32 * a40)) +
                         a12 * (a20 * (a34 * a41 - a31 * a44) +
                                a21 * (a30 * a44 - a34 * a40) +
                                a24 * (a31 * a40 - a30 * a41)) +
                         a14 * (a20 * (a31 * a42 - a32 * a41) +
                                a21 * (a32 * a40 - a30 * a42) +
                                a22 * (a30 * a41 - a31 * a40))) +
                  a04 * (a10 * (a21 * (a32 * a43 - a33 * a42) +
                                a22 * (a33 * a41 - a31 * a43) +
                                a23 * (a31 * a42 - a32 * a41)) +
                         a11 * (a20 * (a33 * a42 - a32 * a43) +
                                a22 * (a30 * a43 - a33 * a40) +
                                a23 * (a32 * a40 - a30 * a42)) +
                         a12 * (a20 * (a31 * a43 - a33 * a41) +
                                a21 * (a33 * a40 - a30 * a43) +
                                a23 * (a30 * a41 - a31 * a40)) +
                         a13 * (a20 * (a32 * a41 - a31 * a42) +
                                a21 * (a30 * a42 - a32 * a40) +
                                a22 * (a31 * a40 - a30 * a41)))) };
}

#endif
