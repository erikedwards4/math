#!/usr/bin/env python3
"""
Sets up ffibuilder and compiles.
"""

from cffi import FFI

# Init builder
ffibuilder = FFI()

# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
ffibuilder.cdef(
    """
    // Generate: Constants Other_Gen Rand

    // Generate: Constants

    int zeros_s (float *Y, const size_t N);
    int zeros_d (double *Y, const size_t N);
    int zeros_c (float *Y, const size_t N);
    int zeros_z (double *Y, const size_t N);

    int ones_s (float *Y, const size_t N);
    int ones_d (double *Y, const size_t N);
    int ones_c (float *Y, const size_t N);
    int ones_z (double *Y, const size_t N);

    int twos_s (float *Y, const size_t N);
    int twos_d (double *Y, const size_t N);
    int twos_c (float *Y, const size_t N);
    int twos_z (double *Y, const size_t N);

    int e_s (float *Y, const size_t N);
    int e_d (double *Y, const size_t N);
    int e_c (float *Y, const size_t N);
    int e_z (double *Y, const size_t N);

    int ln2_s (float *Y, const size_t N);
    int ln2_d (double *Y, const size_t N);
    int ln2_c (float *Y, const size_t N);
    int ln2_z (double *Y, const size_t N);

    int ln10_s (float *Y, const size_t N);
    int ln10_d (double *Y, const size_t N);
    int ln10_c (float *Y, const size_t N);
    int ln10_z (double *Y, const size_t N);

    int log2e_s (float *Y, const size_t N);
    int log2e_d (double *Y, const size_t N);
    int log2e_c (float *Y, const size_t N);
    int log2e_z (double *Y, const size_t N);

    int log10e_s (float *Y, const size_t N);
    int log10e_d (double *Y, const size_t N);
    int log10e_c (float *Y, const size_t N);
    int log10e_z (double *Y, const size_t N);

    int sqrt2_s (float *Y, const size_t N);
    int sqrt2_d (double *Y, const size_t N);
    int sqrt2_c (float *Y, const size_t N);
    int sqrt2_z (double *Y, const size_t N);

    int isqrt2_s (float *Y, const size_t N);
    int isqrt2_d (double *Y, const size_t N);
    int isqrt2_c (float *Y, const size_t N);
    int isqrt2_z (double *Y, const size_t N);

    int pi_s (float *Y, const size_t N);
    int pi_d (double *Y, const size_t N);
    int pi_c (float *Y, const size_t N);
    int pi_z (double *Y, const size_t N);

    int ipi_s (float *Y, const size_t N);
    int ipi_d (double *Y, const size_t N);
    int ipi_c (float *Y, const size_t N);
    int ipi_z (double *Y, const size_t N);

    int pi_2_s (float *Y, const size_t N);
    int pi_2_d (double *Y, const size_t N);
    int pi_2_c (float *Y, const size_t N);
    int pi_2_z (double *Y, const size_t N);

    int pi_4_s (float *Y, const size_t N);
    int pi_4_d (double *Y, const size_t N);
    int pi_4_c (float *Y, const size_t N);
    int pi_4_z (double *Y, const size_t N);

    int eps_s (float *Y, const size_t N);
    int eps_d (double *Y, const size_t N);
    int eps_c (float *Y, const size_t N);
    int eps_z (double *Y, const size_t N);

    int realmin_s (float *Y, const size_t N);
    int realmin_d (double *Y, const size_t N);
    int realmin_c (float *Y, const size_t N);
    int realmin_z (double *Y, const size_t N);

    int realmax_s (float *Y, const size_t N);
    int realmax_d (double *Y, const size_t N);
    int realmax_c (float *Y, const size_t N);
    int realmax_z (double *Y, const size_t N);

    int inf_s (float *Y, const size_t N);
    int inf_d (double *Y, const size_t N);
    int inf_c (float *Y, const size_t N);
    int inf_z (double *Y, const size_t N);

    int nan_s (float *Y, const size_t N);
    int nan_d (double *Y, const size_t N);
    int nan_c (float *Y, const size_t N);
    int nan_z (double *Y, const size_t N);

    int fill_s (float *Y, const size_t N, const float val);
    int fill_d (double *Y, const size_t N, const double val);
    int fill_c (float *Y, const size_t N, const float val);
    int fill_z (double *Y, const size_t N, const double val);

    // Generate: Other_Gen

    int eye_s (float *Y, const size_t R, const size_t C, const int iscolmajor);
    int eye_d (double *Y, const size_t R, const size_t C, const int iscolmajor);
    int eye_c (float *Y, const size_t R, const size_t C, const int iscolmajor);
    int eye_z (double *Y, const size_t R, const size_t C, const int iscolmajor);

    int linspace_s (float *Y, const size_t N, const float a, const float b);
    int linspace_d (double *Y, const size_t N, const double a, const double b);
    int linspace_c (float *Y, const size_t N, const float a, const float b);
    int linspace_z (double *Y, const size_t N, const double a, const double b);

    int logspace_s (float *Y, const size_t N, const float a, const float b);
    int logspace_d (double *Y, const size_t N, const double a, const double b);
    int logspace_c (float *Y, const size_t N, const float a, const float b);
    int logspace_z (double *Y, const size_t N, const double a, const double b);

    int randperm_s (float *Y, const size_t M, const size_t N);
    int randperm_d (double *Y, const size_t M, const size_t N);

    // Generate: Rand

    int randi_s (float *Y, const int a, const int b, const size_t N);
    int randi_d (double *Y, const int a, const int b, const size_t N);
    int randi_c (float *Y, const int a, const int b, const size_t N);
    int randi_z (double *Y, const int a, const int b, const size_t N);

    int randu_s (float *Y, const float a, const float b, const size_t N);
    int randu_d (double *Y, const double a, const double b, const size_t N);
    int randu_c (float *Y, const float a, const float b, const size_t N);
    int randu_z (double *Y, const double a, const double b, const size_t N);

    int randn_s (float *Y, const float mu, const float sig, const size_t N);
    int randn_d (double *Y, const double mu, const double sig, const size_t N);
    int randn_c (float *Y, const float mu, const float sig, const size_t N);
    int randn_z (double *Y, const double mu, const double sig, const size_t N);


    // Matmanip: Construct Matsel Rearrange Split_Join

    // Matmanip: Construct

    int diagmat_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k);
    int diagmat_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k);
    int diagmat_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k);
    int diagmat_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k);

    int tril_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
    int tril_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
    int tril_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
    int tril_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);

    int triu_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
    int triu_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
    int triu_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
    int triu_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);

    int repmat_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);
    int repmat_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);
    int repmat_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);
    int repmat_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);

    // Matmanip: Matsel

    int diag_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k);
    int diag_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k);
    int diag_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k);
    int diag_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k);

    int row_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r);
    int row_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r);
    int row_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r);
    int row_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r);

    int col_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);
    int col_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);
    int col_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);
    int col_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);

    // Matmanip: Rearrange

    int transpose_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor);
    int transpose_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor);
    int transpose_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor);
    int transpose_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor);

    int ctranspose_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor);
    int ctranspose_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor);
    int ctranspose_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor);
    int ctranspose_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor);

    int rot90_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int K);
    int rot90_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int K);
    int rot90_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int K);
    int rot90_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int K);

    // Matmanip: Split_Join

    int split2_s (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int split2_d (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int split2_c (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int split2_z (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int split3_s (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int split3_d (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int split3_c (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int split3_z (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int join2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int join2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int join2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int join2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int join3_s (float *Y, const float *X1, const float *X2, const float *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const int iscolmajor, const size_t dim);
    int join3_d (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const int iscolmajor, const size_t dim);
    int join3_c (float *Y, const float *X1, const float *X2, const float *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const int iscolmajor, const size_t dim);
    int join3_z (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const int iscolmajor, const size_t dim);


    // Elementwise1: Operators Trig Exp_Log Round Special Nonlin

    // Elementwise1: Operators

    int plusplus_s (float *X, const size_t N);
    int plusplus_d (double *X, const size_t N);

    int minusminus_s (float *X, const size_t N);
    int minusminus_d (double *X, const size_t N);

    int neg_s (float *Y, const float *X, const size_t N);
    int neg_d (double *Y, const double *X, const size_t N);
    int neg_c (float *Y, const float *X, const size_t N);
    int neg_z (double *Y, const double *X, const size_t N);

    // Elementwise1: Trig

    int sin_s (float *Y, const float *X, const size_t N);
    int sin_d (double *Y, const double *X, const size_t N);
    int sin_c (float *Y, const float *X, const size_t N);
    int sin_z (double *Y, const double *X, const size_t N);

    int cos_s (float *Y, const float *X, const size_t N);
    int cos_d (double *Y, const double *X, const size_t N);
    int cos_c (float *Y, const float *X, const size_t N);
    int cos_z (double *Y, const double *X, const size_t N);

    int tan_s (float *Y, const float *X, const size_t N);
    int tan_d (double *Y, const double *X, const size_t N);
    int tan_c (float *Y, const float *X, const size_t N);
    int tan_z (double *Y, const double *X, const size_t N);

    int asin_s (float *Y, const float *X, const size_t N);
    int asin_d (double *Y, const double *X, const size_t N);
    int asin_c (float *Y, const float *X, const size_t N);
    int asin_z (double *Y, const double *X, const size_t N);

    int acos_s (float *Y, const float *X, const size_t N);
    int acos_d (double *Y, const double *X, const size_t N);
    int acos_c (float *Y, const float *X, const size_t N);
    int acos_z (double *Y, const double *X, const size_t N);

    int atan_s (float *Y, const float *X, const size_t N);
    int atan_d (double *Y, const double *X, const size_t N);
    int atan_c (float *Y, const float *X, const size_t N);
    int atan_z (double *Y, const double *X, const size_t N);

    int sinh_s (float *Y, const float *X, const size_t N);
    int sinh_d (double *Y, const double *X, const size_t N);
    int sinh_c (float *Y, const float *X, const size_t N);
    int sinh_z (double *Y, const double *X, const size_t N);

    int cosh_s (float *Y, const float *X, const size_t N);
    int cosh_d (double *Y, const double *X, const size_t N);
    int cosh_c (float *Y, const float *X, const size_t N);
    int cosh_z (double *Y, const double *X, const size_t N);

    int tanh_s (float *Y, const float *X, const size_t N);
    int tanh_d (double *Y, const double *X, const size_t N);
    int tanh_c (float *Y, const float *X, const size_t N);
    int tanh_z (double *Y, const double *X, const size_t N);

    int asinh_s (float *Y, const float *X, const size_t N);
    int asinh_d (double *Y, const double *X, const size_t N);
    int asinh_c (float *Y, const float *X, const size_t N);
    int asinh_z (double *Y, const double *X, const size_t N);

    int acosh_s (float *Y, const float *X, const size_t N);
    int acosh_d (double *Y, const double *X, const size_t N);
    int acosh_c (float *Y, const float *X, const size_t N);
    int acosh_z (double *Y, const double *X, const size_t N);

    int atanh_s (float *Y, const float *X, const size_t N);
    int atanh_d (double *Y, const double *X, const size_t N);
    int atanh_c (float *Y, const float *X, const size_t N);
    int atanh_z (double *Y, const double *X, const size_t N);

    int rad2deg_s (float *Y, const float *X, const size_t N);
    int rad2deg_d (double *Y, const double *X, const size_t N);

    int deg2rad_s (float *Y, const float *X, const size_t N);
    int deg2rad_d (double *Y, const double *X, const size_t N);

    // Elementwise1: Exp_Log

    int exp_s (float *Y, const float *X, const size_t N);
    int exp_d (double *Y, const double *X, const size_t N);
    int exp_c (float *Y, const float *X, const size_t N);
    int exp_z (double *Y, const double *X, const size_t N);

    int exp2_s (float *Y, const float *X, const size_t N);
    int exp2_d (double *Y, const double *X, const size_t N);
    int exp2_c (float *Y, const float *X, const size_t N);
    int exp2_z (double *Y, const double *X, const size_t N);

    int exp10_s (float *Y, const float *X, const size_t N);
    int exp10_d (double *Y, const double *X, const size_t N);
    int exp10_c (float *Y, const float *X, const size_t N);
    int exp10_z (double *Y, const double *X, const size_t N);

    int log_s (float *Y, const float *X, const size_t N);
    int log_d (double *Y, const double *X, const size_t N);
    int log_c (float *Y, const float *X, const size_t N);
    int log_z (double *Y, const double *X, const size_t N);

    int log2_s (float *Y, const float *X, const size_t N);
    int log2_d (double *Y, const double *X, const size_t N);
    int log2_c (float *Y, const float *X, const size_t N);
    int log2_z (double *Y, const double *X, const size_t N);

    int log10_s (float *Y, const float *X, const size_t N);
    int log10_d (double *Y, const double *X, const size_t N);
    int log10_c (float *Y, const float *X, const size_t N);
    int log10_z (double *Y, const double *X, const size_t N);

    // Elementwise1: Round

    int floor_s (float *Y, const float *X, const size_t N);
    int floor_d (double *Y, const double *X, const size_t N);
    int floor_c (float *Y, const float *X, const size_t N);
    int floor_z (double *Y, const double *X, const size_t N);

    int ceil_s (float *Y, const float *X, const size_t N);
    int ceil_d (double *Y, const double *X, const size_t N);
    int ceil_c (float *Y, const float *X, const size_t N);
    int ceil_z (double *Y, const double *X, const size_t N);

    int round_s (float *Y, const float *X, const size_t N);
    int round_d (double *Y, const double *X, const size_t N);
    int round_c (float *Y, const float *X, const size_t N);
    int round_z (double *Y, const double *X, const size_t N);

    int trunc_s (float *Y, const float *X, const size_t N);
    int trunc_d (double *Y, const double *X, const size_t N);
    int trunc_c (float *Y, const float *X, const size_t N);
    int trunc_z (double *Y, const double *X, const size_t N);

    // Elementwise1: Special

    int erf_s (float *Y, const float *X, const size_t N);
    int erf_d (double *Y, const double *X, const size_t N);

    int erfc_s (float *Y, const float *X, const size_t N);
    int erfc_d (double *Y, const double *X, const size_t N);

    int tgamma_s (float *Y, const float *X, const size_t N);
    int tgamma_d (double *Y, const double *X, const size_t N);

    int lgamma_s (float *Y, const float *X, const size_t N);
    int lgamma_d (double *Y, const double *X, const size_t N);

    // Elementwise1: Nonlin

    int abs_s (float *Y, const float *X, const size_t N);
    int abs_d (double *Y, const double *X, const size_t N);
    int abs_c (float *Y, const float *X, const size_t N);
    int abs_z (double *Y, const double *X, const size_t N);

    int square_s (float *Y, const float *X, const size_t N);
    int square_d (double *Y, const double *X, const size_t N);
    int square_c (float *Y, const float *X, const size_t N);
    int square_z (double *Y, const double *X, const size_t N);

    int cube_s (float *Y, const float *X, const size_t N);
    int cube_d (double *Y, const double *X, const size_t N);
    int cube_c (float *Y, const float *X, const size_t N);
    int cube_z (double *Y, const double *X, const size_t N);

    int sqrt_s (float *Y, const float *X, const size_t N);
    int sqrt_d (double *Y, const double *X, const size_t N);
    int sqrt_c (float *Y, const float *X, const size_t N);
    int sqrt_z (double *Y, const double *X, const size_t N);

    int cbrt_s (float *Y, const float *X, const size_t N);
    int cbrt_d (double *Y, const double *X, const size_t N);
    int cbrt_c (float *Y, const float *X, const size_t N);
    int cbrt_z (double *Y, const double *X, const size_t N);

    int reciprocal_s (float *Y, const float *X, const size_t N);
    int reciprocal_d (double *Y, const double *X, const size_t N);
    int reciprocal_c (float *Y, const float *X, const size_t N);
    int reciprocal_z (double *Y, const double *X, const size_t N);

    int sign_s (float *Y, const float *X, const size_t N);
    int sign_d (double *Y, const double *X, const size_t N);

    int deadzone_s (float *Y, const float *X, const size_t N, const float delta);
    int deadzone_d (double *Y, const double *X, const size_t N, const double delta);


    // Elementwise2

    int plus_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int plus_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int plus_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int plus_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int minus_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int minus_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int minus_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int minus_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int times_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int times_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int times_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int times_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int rdivide_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int rdivide_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int rdivide_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int rdivide_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int pow_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int pow_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int pow_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int pow_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int hypot_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int hypot_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int hypot_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int hypot_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int atan2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int atan2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int adiff_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int adiff_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int adiff_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int adiff_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);


    // Vec2scalar: Sums Prctiles Iprctiles Ranges Norms Moments Other_Means Other_Spreads Other_Stats

    // Vec2scalar: Sums

    int sum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int sum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int sum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int sum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int asum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int asum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int asum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int asum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int cnt_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cnt_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cnt_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cnt_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    // Vec2scalar: Prctiles

    int prctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int prctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);

    int median_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int median_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int max_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int max_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int max_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int max_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int min_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int min_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int min_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int min_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int amax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int amax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int amax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int amax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    // Vec2scalar: Iprctiles

    int iprctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int iprctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);

    int imed_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int imed_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int imax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int imax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int imax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int imax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int imin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int imin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int imin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int imin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int iamax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int iamax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int iamax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int iamax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    // Vec2scalar: Ranges

    int range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int range_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int iqr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int iqr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int idr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int idr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    // Vec2scalar: Norms

    int norm0_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm0_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm0_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm0_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int norm1_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm1_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm1_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm1_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int norm2_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int norm2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int normp_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int normp_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);
    int normp_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int normp_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);

    // Vec2scalar: Moments

    int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int var_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int var_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);

    int skewness_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int skewness_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int skewness_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int skewness_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);

    int kurtosis_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int kurtosis_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int kurtosis_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int kurtosis_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);

    // Vec2scalar: Other_Means

    int trimmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
    int trimmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);

    int winsormean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
    int winsormean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);

    int geomean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int geomean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int geomean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int geomean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int harmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int harmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int harmean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int harmean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int genmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int genmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);
    int genmean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int genmean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);

    int rms_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int rms_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int rms_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int rms_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    // Vec2scalar: Other_Spreads

    int std_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int std_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int std_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int std_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);

    int geostd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int geostd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int trimstd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q, const int biased);
    int trimstd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q, const int biased);

    int trimvar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q, const int biased);
    int trimvar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q, const int biased);

    int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);

    int mad_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mad_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    // Vec2scalar: Other_Stats

    int prod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int prod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int prod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int prod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


    // Vec2vec: Center Scale Normalize Reorder Other_Vec2vec

    // Vec2vec: Center

    int mean0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mean0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mean0_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mean0_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int med0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int med0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int geomean1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int geomean1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int geomean1_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int geomean1_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    // Vec2vec: Scale

    int zscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int zscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int zscore_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int zscore_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);

    int mscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int mscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int gscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int gscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int range1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);
    int range1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);

    int iqr1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);
    int iqr1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);

    int idr1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);
    int idr1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);

    // Vec2vec: Normalize

    int normalize1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int normalize1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int normalize1_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int normalize1_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int normalize2_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int normalize2_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int normalize2_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int normalize2_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int normalizep_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int normalizep_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);
    int normalizep_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p);
    int normalizep_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p);

    // Vec2vec: Reorder

    int flip_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int flip_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int flip_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int flip_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int shift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
    int shift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
    int shift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
    int shift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);

    int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
    int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
    int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
    int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);

    int sort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int sort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int sort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int sort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);

    // Vec2vec: Other_Vec2vec

    int sorti_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int sorti_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int sorti_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int sorti_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);

    int ranks_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int ranks_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int ranks_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);
    int ranks_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend);

    int prctiles_s (float *Y, const float *X, const float *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int prctiles_d (double *Y, const double *X, const double *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int moments_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int moments_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int moments_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
    int moments_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);

    int winsorize_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
    int winsorize_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);

    int trim_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
    int trim_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);

    int cumsum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cumsum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cumsum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cumsum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int cumprod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cumprod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cumprod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int cumprod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int prepad_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const float val);
    int prepad_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const double val);
    int prepad_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const float val);
    int prepad_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const double val);

    int postpad_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const float val);
    int postpad_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const double val);
    int postpad_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const float val);
    int postpad_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const double val);

    int softmax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
    int softmax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

    int betamax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float base);
    int betamax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double base);


    // Vecs2scalar: Similarity Dist WStats

    // Vecs2scalar: Similarity

    int dot_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dot_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dot_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dot_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int cov_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const int biased);
    int cov_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const int biased);
    int cov_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const int biased);
    int cov_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const int biased);

    int corr_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int corr_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int corr_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int corr_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int corr_opus_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int corr_opus_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int corr_opus_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int corr_opus_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int cos2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int cos2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int cos2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int cos2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int cokurtosis_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int cokurtosis_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int spearman_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int spearman_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int kendall_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int kendall_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    // Vecs2scalar: Dist

    int dist0_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist0_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist0_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist0_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int dist1_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist1_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist1_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist1_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int dist2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);

    int distp_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const float p);
    int distp_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const double p);
    int distp_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const float p);
    int distp_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim, const double p);

    int dist_cos_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
    int dist_cos_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);


    // Complex

    int complex_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int complex_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int polar_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int polar_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);

    int real_c (float *Y, const float *X, const size_t N);
    int real_z (double *Y, const double *X, const size_t N);

    int imag_c (float *Y, const float *X, const size_t N);
    int imag_z (double *Y, const double *X, const size_t N);

    int conj_c (float *Y, const float *X, const size_t N);
    int conj_z (double *Y, const double *X, const size_t N);

    int arg_s (float *Y, const float *X, const size_t N);
    int arg_d (double *Y, const double *X, const size_t N);
    int arg_c (float *Y, const float *X, const size_t N);
    int arg_z (double *Y, const double *X, const size_t N);

    int proj_c (float *Y, const float *X, const size_t N);
    int proj_z (double *Y, const double *X, const size_t N);


    // Linalg: Matmul Transform Sim_Mat Dist_Mat Other_Linalg

    // Linalg: Matmul

    int matmul1_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor);
    int matmul1_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor);
    int matmul1_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor);
    int matmul1_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor);

    int matmul1t_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int tr);
    int matmul1t_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int tr);
    int matmul1t_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int tr);
    int matmul1t_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int tr);

    int matmul2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);
    int matmul2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);
    int matmul2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);
    int matmul2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);

    int matmul2t_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);
    int matmul2t_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);
    int matmul2t_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);
    int matmul2t_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);

    int matmul3_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const int iscolmajor);
    int matmul3_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const int iscolmajor);
    int matmul3_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const int iscolmajor);
    int matmul3_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const int iscolmajor);

    int matmul3t_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);
    int matmul3t_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);
    int matmul3t_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);
    int matmul3t_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);

    int kronecker_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int kronecker_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int kronecker_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
    int kronecker_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);


    // Linalg: Transform

    int linear_s (float *Y, const float *X, const float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);
    int linear_d (double *Y, const double *X, const double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);
    int linear_c (float *Y, const float *X, const float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);
    int linear_z (double *Y, const double *X, const double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);

    int affine_s (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);
    int affine_d (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);
    int affine_c (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);
    int affine_z (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim);

    // Linalg: Sim_Mat

    int dotmat_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);
    int dotmat_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);
    int dotmat_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);
    int dotmat_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);

    int cosmat_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);
    int cosmat_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);
    int cosmat_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);
    int cosmat_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const size_t dim);

    // Linalg: Other_Linalg

    int solve_s (float *X, const float *A, const float *B, const size_t R1, const size_t C1, const size_t C2, const int iscolmajor, const int tr);
    int solve_d (double *X, const double *A, const double *B, const size_t R1, const size_t C1, const size_t C2, const int iscolmajor, const int tr);
    int solve_c (float *X, const float *A, const float *B, const size_t R1, const size_t C1, const size_t C2, const int iscolmajor, const int tr);
    int solve_z (double *X, const double *A, const double *B, const size_t R1, const size_t C1, const size_t C2, const int iscolmajor, const int tr);

    int chol_s (float *Y, const float *X, const size_t R, const int iscolmajor, const int upper);
    int chol_d (double *Y, const double *X, const size_t R, const int iscolmajor, const int upper);
    int chol_c (float *Y, const float *X, const size_t R, const int iscolmajor, const int upper);
    int chol_z (double *Y, const double *X, const size_t R, const int iscolmajor, const int upper);

    int eig_s (float *U, float *V, const float *X, const size_t R, const int iscolmajor, const size_t K);
    int eig_d (double *U, double *V, const double *X, const size_t R, const int iscolmajor, const size_t K);
    int eig_c (float *U, float *V, const float *X, const size_t R, const int iscolmajor, const size_t K);
    int eig_z (double *U, double *V, const double *X, const size_t R, const int iscolmajor, const size_t K);

    int svd_s (float *U, float *S, float *Vt, const float *X, const size_t R, const size_t C, const int iscolmajor, const size_t K);
    int svd_d (double *U, double *S, double *Vt, const double *X, const size_t R, const size_t C, const int iscolmajor, const size_t K);
    int svd_c (float *U, float *S, float *Vt, const float *X, const size_t R, const size_t C, const int iscolmajor, const size_t K);
    int svd_z (double *U, double *S, double *Vt, const double *X, const size_t R, const size_t C, const int iscolmajor, const size_t K);

    """
)


# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
ffibuilder.set_source("_math_cffi",
"""
     #include "codee_math.h"   // the C header of the library
""",
    library_dirs = ["."],
    libraries = ['math'],   # library name, for the linker
)


if __name__ == "__main__":
    """
    Compiles to local .c, .o and .so files.
    """
    ffibuilder.compile(verbose=True)
